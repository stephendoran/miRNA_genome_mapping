"""
Map an RNA (or DNA) sequence to genomic coordinates by aligning it against
the CURRENT Ensembl reference genome for a given species, using a local
BLAST search.

Why local BLAST rather than a remote Ensembl call
----------------------------------------------------
Ensembl does not provide a programmatic/remote BLAST API - this is
confirmed directly by Ensembl's own developers ("we don't provide a
remote service for BLASTing - you have to download the data and run your
own BLASTs"). Their website BLAST/BLAT tool is browser-only. So "mapping
via Ensembl" in practice means: download Ensembl's own reference genome
FASTA (keeping the result on Ensembl's current assembly, consistent with
the rest of this pipeline), then align locally yourself.

Requires (installed once, system-wide):
    NCBI BLAST+ (makeblastdb, blastn)
      Ubuntu/Debian: sudo apt-get install ncbi-blast+
      conda:         conda install -c bioconda blast
      macOS:         brew install blast

For species Ensembl doesn't have at all (e.g. Python bivittatus,
Caenorhabditis briggsae, most Drosophila species other than melanogaster),
this falls back to NCBI automatically - via NCBI's Entrez E-utilities
(plain HTTPS, no separate tool or install required).

Python dependency:
    pip install pandas --break-system-packages

Usage - the mmu-mir-1966 example from your message:
    python map_sequence_to_genome.py \
        --species mus_musculus \
        --name mmu-mir-1966 \
        --sequence AUAGUGUUGGAAGGGAGCUGGCUCAGGAGAGAGUCCUGAGAUUUAGGCUCUUUCUGACUCAACUCUCCCUUAGCAAGUCAAGU

For batch use later (many sequences, one species), call
map_sequences_batch() directly with a {name: sequence} dict - it runs a
single blastn call for all of them rather than one process per sequence,
and reuses the downloaded genome/BLAST database across calls.
"""

import argparse
import gzip
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd

REST_HOST = "https://rest.ensembl.org"
CACHE_DIR = Path("genome_cache")

# Ensembl Genomes divisions (plants/metazoa/fungi/protists) are hosted on a
# separate FTP server from the main vertebrates site, under their own
# "current" alias - confirmed directly against Ensembl's FTP documentation
# and an EBI rsync example (ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/
# fasta/...). The main site's alias (ftp.ensembl.org/pub/current_fasta/) is
# a different host entirely and returns 404 for plant/metazoa/fungi/protist
# species - this is exactly what caused rice, poplar, grape, and maize to
# fail.
FASTA_BASE_URL_BY_DIVISION = {
    "vertebrates": "https://ftp.ensembl.org/pub/current_fasta/{species_key}/dna/",
    "plants": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/fasta/{species_key}/dna/",
    "metazoa": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/metazoa/current/fasta/{species_key}/dna/",
    "fungi": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/fungi/current/fasta/{species_key}/dna/",
    "protists": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/protists/current/fasta/{species_key}/dna/",
}

BLAST_FIELDS = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore sstrand"
)


def _fetch_assembly_name(species_key):
    """species_key: Ensembl production name, e.g. 'mus_musculus'."""
    url = f"{REST_HOST}/info/assembly/{species_key}?content-type=application/json"
    req = Request(url, headers={"Content-Type": "application/json"})
    with urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode())
    return data["assembly_name"]


def _strip_patch_suffix(assembly_name):
    """
    Ensembl's REST API reports the full patch-qualified assembly name
    (e.g. 'GRCh38.p14'), but its FTP filenames use only the base assembly
    name without the patch suffix (e.g. 'Homo_sapiens.GRCh38.dna...') -
    confirmed against multiple working Ensembl download examples. Mouse
    happened to work without this fix only because GRCm39 currently has
    no patch suffix attached; human's GRCh38.p14 does, which is what
    caused the human genome download to 404.
    """
    return re.sub(r"\.p\d+$", "", assembly_name)


def _download_genome_fasta_ensembl(species_key, division="vertebrates"):
    """
    Download (once, cached locally) the primary-assembly genome FASTA from
    the correct FTP server for this species' Ensembl division - see
    FASTA_BASE_URL_BY_DIVISION above. Uses primary_assembly (real
    chromosomes only, no patches/haplotypes/scaffolds), matching the
    karyotype-only filtering used elsewhere in this pipeline, so results
    should be labeled with plain chromosome names directly.

    Raises if the species isn't on Ensembl at all (no assembly info) or
    the FASTA can't be found under either naming variant - callers should
    catch this and fall back to _download_genome_fasta_ncbi.
    """
    assembly_name = _strip_patch_suffix(_fetch_assembly_name(species_key))

    CACHE_DIR.mkdir(exist_ok=True)
    display_name = species_key.capitalize()  # 'mus_musculus' -> 'Mus_musculus'
    base_url = FASTA_BASE_URL_BY_DIVISION.get(division, FASTA_BASE_URL_BY_DIVISION["vertebrates"]).format(
        species_key=species_key
    )

    # Ensembl only publishes a separate "primary_assembly" file when it
    # actually differs from "toplevel" (i.e. when there are extra
    # haplotype/patch/unplaced scaffolds to exclude). Many genomes -
    # plants especially, which rarely have alternate haplotypes - only
    # have a "toplevel" file at all, so primary_assembly 404s for them.
    # Try primary_assembly first (cleaner - chromosomes only) and fall
    # back to toplevel if that 404s.
    last_error = None
    for variant in ("primary_assembly", "toplevel"):
        basename = f"{display_name}.{assembly_name}.dna.{variant}.fa"
        fasta_gz = CACHE_DIR / f"{basename}.gz"
        fasta = CACHE_DIR / basename

        if fasta.exists():
            print(f"Using cached genome FASTA: {fasta}")
            return fasta

        url = base_url + f"{basename}.gz"
        print(f"Downloading genome FASTA (large file - this can take a while):\n  {url}")
        try:
            req = Request(url)
            with urlopen(req, timeout=600) as resp, open(fasta_gz, "wb") as out:
                shutil.copyfileobj(resp, out)
        except HTTPError as exc:
            last_error = exc
            if exc.code == 404 and variant == "primary_assembly":
                print("  primary_assembly not available for this species - trying toplevel instead.")
                continue
            raise

        print("Decompressing...")
        with gzip.open(fasta_gz, "rb") as f_in, open(fasta, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        if variant == "toplevel":
            print(
                "Note: using 'toplevel' (no separate primary_assembly file exists "
                "for this species) - this can include unplaced/unlocalized "
                "scaffolds alongside named chromosomes, so a hit could "
                "occasionally land on a scaffold rather than a chromosome name."
            )
        return fasta

    raise RuntimeError(
        f"Could not download genome FASTA for '{species_key}' from Ensembl - tried both "
        f"primary_assembly and toplevel, last error: {last_error}"
    )


ENTREZ_HOST = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _lookup_ncbi_assembly_ftp_path(scientific_name, search_term_override=None):
    """
    Find the FTP directory for a species' genome assembly using NCBI's
    Entrez E-utilities (esearch + esummary against the 'assembly'
    database) - plain HTTP, no separate tool required.

    By default searches by organism name and prefers a flagged reference/
    representative genome (i.e. "give me the current best genome for this
    species"). Pass search_term_override to target a SPECIFIC assembly
    instead - e.g. "GCA_000001635.9[Assembly Accession]" or
    '"GRCm39"[Assembly Name]' - used when a particular build is required
    rather than whichever one is currently the reference (see
    map_all_hairpins_to_genomes.py's --target-builds support).

    Returns (ftp_path, accession, assembly_name). ftp_path is the base
    directory containing that assembly's files (e.g.
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14) -
    the actual genomic FASTA is at
    {ftp_path}/{last_path_segment}_genomic.fna.gz.
    """
    from urllib.parse import quote_plus

    search_term = search_term_override or (scientific_name + "[Organism]")
    search_url = (
        f"{ENTREZ_HOST}/esearch.fcgi?db=assembly&retmode=json&retmax=20"
        f"&term={quote_plus(search_term)}"
    )
    req = Request(search_url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=60) as resp:
        search_data = json.loads(resp.read().decode())
    uids = search_data.get("esearchresult", {}).get("idlist", [])
    if not uids:
        raise RuntimeError(f"No NCBI assembly found for search term '{search_term}'.")

    summary_url = f"{ENTREZ_HOST}/esummary.fcgi?db=assembly&retmode=json&id={','.join(uids)}"
    req = Request(summary_url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=60) as resp:
        summary_data = json.loads(resp.read().decode())
    records = summary_data.get("result", {})

    candidates = [records[uid] for uid in uids if uid in records]
    if not candidates:
        raise RuntimeError(f"NCBI esummary returned no usable records for '{search_term}'.")

    if search_term_override:
        # Targeting a specific assembly - the search term should already
        # be precise enough to return exactly what we want, so just take
        # the first result rather than re-ranking by reference/
        # representative status (that ranking is for "give me the best
        # current genome", not relevant when a specific one is named).
        chosen = candidates[0]
    else:
        # Prefer a flagged reference/representative genome; otherwise fall
        # back to the first result (esearch generally returns the most
        # relevant/recent first, but this isn't guaranteed for every organism).
        def category_rank(rec):
            category = (rec.get("refseq_category") or "").lower()
            if "reference" in category:
                return 0
            if "representative" in category:
                return 1
            return 2

        candidates.sort(key=category_rank)
        chosen = candidates[0]

    ftp_path = chosen.get("ftppath_refseq") or chosen.get("ftppath_genbank")
    if not ftp_path:
        raise RuntimeError(f"NCBI record for '{search_term}' has no FTP path available.")

    return ftp_path.replace("ftp://", "https://"), chosen.get("assemblyaccession"), chosen.get("assemblyname")


def _fetch_ncbi_assembly_report_chrom_map(ftp_path):
    """
    Fetch and parse an NCBI assembly's *_assembly_report.txt (published
    alongside every genome at the same FTP path) to map each sequence's
    GenBank accession to its assigned chromosome name - needed because
    NCBI genomic FASTA headers use accession-style identifiers (e.g.
    "gb|CP126666.1|") rather than clean chromosome names like Ensembl's
    "19".

    Returns dict {genbank_accession: assigned_molecule_name}, restricted
    to rows with Sequence-Role "assembled-molecule" (i.e. actual named
    chromosomes, not unplaced scaffolds/contigs).
    """
    basename = ftp_path.rstrip("/").rsplit("/", 1)[-1]
    url = f"{ftp_path}/{basename}_assembly_report.txt"

    chrom_map = {}
    try:
        req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urlopen(req, timeout=60) as resp:
            columns = None
            for raw_line in resp:
                line = raw_line.decode("utf-8", errors="replace").rstrip("\n")
                if not line:
                    continue
                if line.startswith("#"):
                    # the last comment line before data rows is the
                    # tab-separated column header
                    columns = line.lstrip("#").strip().split("\t")
                    continue
                if columns is None:
                    continue
                row = dict(zip(columns, line.split("\t")))
                if row.get("Sequence-Role") == "assembled-molecule":
                    accession = row.get("GenBank-Accn")
                    molecule = row.get("Assigned-Molecule")
                    if accession and molecule and accession != "na":
                        chrom_map[accession] = molecule
    except (HTTPError, URLError, TimeoutError) as exc:
        print(
            f"[warn] Could not fetch NCBI assembly report from {url}: {exc} - "
            "chromosome names will remain as raw accessions.",
            file=sys.stderr,
        )

    return chrom_map


def _download_ncbi_fasta_from_ftp_path(ftp_path, fasta_path, chrom_map_path):
    """Shared download logic: fetch the chromosome map and genomic FASTA
    for a resolved NCBI ftp_path, writing both to the given cache paths."""
    chrom_map = _fetch_ncbi_assembly_report_chrom_map(ftp_path)
    chrom_map_path.write_text(json.dumps(chrom_map))

    basename = ftp_path.rstrip("/").rsplit("/", 1)[-1]
    fasta_url = f"{ftp_path}/{basename}_genomic.fna.gz"
    fasta_gz = fasta_path.with_suffix(fasta_path.suffix + ".gz")

    print(f"Downloading genome FASTA (large file - this can take a while):\n  {fasta_url}")
    req = Request(fasta_url)
    with urlopen(req, timeout=600) as resp, open(fasta_gz, "wb") as out:
        shutil.copyfileobj(resp, out)

    print("Decompressing...")
    with gzip.open(fasta_gz, "rb") as f_in, open(fasta_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    print(f"NCBI genome FASTA ready: {fasta_path}")
    return fasta_path, chrom_map


def _download_genome_fasta_ncbi(species_key):
    """
    Fallback for species Ensembl doesn't have at all (e.g. Python
    bivittatus, Caenorhabditis briggsae, most Drosophila other than
    melanogaster) - looks up the reference/representative assembly via
    NCBI's Entrez API and downloads its genomic FASTA directly over HTTPS,
    the same way _download_genome_fasta_ensembl does for Ensembl. No
    separate CLI tool required.

    Also fetches that assembly's chromosome-name mapping (see
    _fetch_ncbi_assembly_report_chrom_map) and caches it alongside the
    FASTA, since NCBI's own sequence naming in the genomic FASTA is often
    an accession (e.g. "gb|CP126666.1|") rather than a plain chromosome
    name - the mapping lets callers relabel hits to real chromosome names.

    Returns (fasta_path, chrom_map).
    """
    CACHE_DIR.mkdir(exist_ok=True)
    fasta_path = CACHE_DIR / f"{species_key}.NCBI.dna.fa"
    chrom_map_path = CACHE_DIR / f"{species_key}.NCBI.chrom_map.json"

    if fasta_path.exists():
        print(f"Using cached NCBI genome FASTA: {fasta_path}")
        if chrom_map_path.exists():
            chrom_map = json.loads(chrom_map_path.read_text())
        else:
            # FASTA was cached before this feature existed - backfill the
            # chromosome map without re-downloading the (large) FASTA.
            taxon_name = species_key.replace("_", " ")
            ftp_path, _accession, _assembly_name = _lookup_ncbi_assembly_ftp_path(taxon_name)
            chrom_map = _fetch_ncbi_assembly_report_chrom_map(ftp_path)
            chrom_map_path.write_text(json.dumps(chrom_map))
        return fasta_path, chrom_map

    taxon_name = species_key.replace("_", " ")
    print(f"Ensembl does not have this species - falling back to NCBI for taxon '{taxon_name}'...")
    ftp_path, accession, assembly_name = _lookup_ncbi_assembly_ftp_path(taxon_name)
    print(f"  Found NCBI assembly {accession} ({assembly_name})")

    return _download_ncbi_fasta_from_ftp_path(ftp_path, fasta_path, chrom_map_path)


def _download_genome_fasta_ncbi_target(species_key, accession=None, assembly_name=None):
    """
    Like _download_genome_fasta_ncbi, but fetches a SPECIFIC assembly
    (identified by accession if available, otherwise by assembly name)
    rather than whichever NCBI currently considers the reference genome -
    used when an exact target build is specified (see
    map_all_hairpins_to_genomes.py's --target-builds support).

    Returns (fasta_path, chrom_map).
    """
    if not accession and not assembly_name:
        raise ValueError("Need an accession or assembly_name to target a specific NCBI assembly.")

    CACHE_DIR.mkdir(exist_ok=True)
    cache_key = re.sub(r"[^A-Za-z0-9_.-]", "_", accession or assembly_name)
    fasta_path = CACHE_DIR / f"{species_key}.NCBI_{cache_key}.dna.fa"
    chrom_map_path = CACHE_DIR / f"{species_key}.NCBI_{cache_key}.chrom_map.json"

    if fasta_path.exists():
        print(f"Using cached NCBI genome FASTA: {fasta_path}")
        chrom_map = json.loads(chrom_map_path.read_text()) if chrom_map_path.exists() else {}
        return fasta_path, chrom_map

    search_term = f"{accession}[Assembly Accession]" if accession else f'"{assembly_name}"[Assembly Name]'
    print(f"Fetching specific target assembly for '{species_key}': {accession or assembly_name}")
    ftp_path, found_accession, found_assembly_name = _lookup_ncbi_assembly_ftp_path(
        species_key.replace("_", " "), search_term_override=search_term
    )
    print(f"  Found NCBI assembly {found_accession} ({found_assembly_name})")

    return _download_ncbi_fasta_from_ftp_path(ftp_path, fasta_path, chrom_map_path)


def _get_genome_fasta(species_key, division="vertebrates", target_accession=None, target_assembly=None):
    """
    Get the genome FASTA for a species.

    Default behavior (no target given): prefer Ensembl's CURRENT assembly
    (keeps coordinates on Ensembl's current build, consistent with the
    rest of this pipeline), falling back to NCBI's current
    reference/representative genome if Ensembl doesn't have this species
    at all or the download fails for any reason.

    Targeted behavior (target_accession and/or target_assembly given):
    fetch that SPECIFIC assembly instead of whatever is currently latest.
    Checks first whether Ensembl's current assembly already happens to
    match the target name (cheap, avoids an unnecessary NCBI round-trip
    if so); otherwise fetches the specific assembly directly from NCBI by
    accession (preferred, unambiguous) or assembly name.

    Returns (fasta_path, source, chrom_map). source is 'ensembl' or 'ncbi'.
    chrom_map is {} for Ensembl (its sequence names are already clean
    chromosome names) or {accession: chromosome_name} for NCBI (see
    _fetch_ncbi_assembly_report_chrom_map) - callers use it to relabel
    BLAST hit chromosome names from raw accessions to real chromosome
    names where possible.
    """
    if target_accession or target_assembly:
        if target_assembly:
            try:
                current_ensembl_name = _strip_patch_suffix(_fetch_assembly_name(species_key))
                if current_ensembl_name == target_assembly:
                    fasta_path = _download_genome_fasta_ensembl(species_key, division=division)
                    return fasta_path, "ensembl", {}
            except Exception:
                pass  # Ensembl doesn't have it, or a different build - fall through to NCBI

        fasta_path, chrom_map = _download_genome_fasta_ncbi_target(
            species_key, accession=target_accession, assembly_name=target_assembly
        )
        return fasta_path, "ncbi", chrom_map

    try:
        fasta_path = _download_genome_fasta_ensembl(species_key, division=division)
        return fasta_path, "ensembl", {}
    except Exception as exc:
        print(f"[warn] Ensembl genome unavailable for '{species_key}' ({exc}) - trying NCBI fallback.", file=sys.stderr)
        fasta_path, chrom_map = _download_genome_fasta_ncbi(species_key)
        return fasta_path, "ncbi", chrom_map


def _build_blast_db(fasta_path):
    """Build (once, cached) a nucleotide BLAST database from the genome FASTA."""
    db_path = fasta_path  # blast appends its own extensions (.nsq etc.) to this basename
    if Path(str(db_path) + ".nsq").exists() or Path(str(db_path) + ".00.nsq").exists():
        print(f"Using cached BLAST database: {db_path}")
        return db_path

    print("Building BLAST database (one-time step for this genome)...")
    subprocess.run(
        ["makeblastdb", "-in", str(fasta_path), "-dbtype", "nucl", "-parse_seqids", "-out", str(db_path)],
        check=True,
    )
    return db_path


DEFAULT_CHUNK_SIZE = 200  # sequences per blastn call - keeps memory bounded for large species


def _run_blast_chunk(sequences_chunk, db_path, min_identity, min_coverage):
    """Run one blastn call for a chunk of sequences and return best-hit rows."""
    query_fasta = CACHE_DIR / "query_chunk.fa"
    with open(query_fasta, "w") as fh:
        for name, seq in sequences_chunk.items():
            dna_seq = seq.upper().replace("U", "T")
            fh.write(f">{name}\n{dna_seq}\n")

    # Choose task/word_size by actual query length. blastn-short + a tiny
    # word_size (e.g. 7) is only appropriate for queries under ~50bp - for
    # longer ones (miRNA hairpins are typically 60-120nt) it generates a huge
    # number of low-specificity seed hits against a whole genome, which is
    # exactly what caused an OOM/SIGKILL on a genome-scale run. Standard
    # blastn's default word_size (11) is both faster and far lighter on
    # memory for this length range, while still sensitive enough for the
    # near-exact match expected from a genomic sequence against its own
    # source genome.
    max_len = max(len(s) for s in sequences_chunk.values())
    if max_len < 50:
        task_args = ["-task", "blastn-short", "-word_size", "7"]
    else:
        task_args = ["-task", "blastn", "-word_size", "11"]

    result = subprocess.run(
        [
            "blastn", *task_args,
            "-db", str(db_path), "-query", str(query_fasta),
            "-evalue", "1e-10",
            "-max_target_seqs", "5",  # only need the best few hits per query, bounds memory
            "-outfmt", f"6 {BLAST_FIELDS}",
        ],
        check=True, capture_output=True, text=True,
    )

    columns = BLAST_FIELDS.split()
    rows = [line.split("\t") for line in result.stdout.strip().splitlines() if line.strip()]
    if not rows:
        return pd.DataFrame(columns=columns)

    hits = pd.DataFrame(rows, columns=columns)
    numeric_cols = ["pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    hits[numeric_cols] = hits[numeric_cols].apply(pd.to_numeric)

    query_lengths = {name: len(seq) for name, seq in sequences_chunk.items()}
    hits["query_length"] = hits["qseqid"].map(query_lengths)
    hits["coverage"] = hits["length"] / hits["query_length"]
    return hits


def map_sequences_batch(
    sequences, species_key, division="vertebrates", min_identity=90.0, min_coverage=0.9,
    chunk_size=DEFAULT_CHUNK_SIZE, return_raw_hits=False,
    target_accession=None, target_assembly=None,
):
    """
    Map multiple sequences to genomic coordinates, chunking large batches
    into several blastn calls (chunk_size sequences each) rather than one
    call for everything - bounds peak memory regardless of how many
    sequences a species has.

    sequences: dict of {name: rna_or_dna_sequence}
    species_key: Ensembl production name, e.g. 'mus_musculus'
    division: one of 'vertebrates', 'plants', 'metazoa', 'fungi', 'protists'
        - determines which FTP server the genome FASTA is downloaded from
        (see FASTA_BASE_URL_BY_DIVISION); defaults to 'vertebrates'
    min_identity: minimum percent identity (0-100) to accept a hit
    min_coverage: minimum fraction of the query sequence that must be
        covered by the alignment to accept a hit
    chunk_size: sequences per blastn call - lower this further if you
        still see out-of-memory kills on a very large-genome species
    return_raw_hits: if True, also return the full hit table (all hits
        clearing BLAST's own e-value cutoff, not just the one best hit per
        query) - needed by callers that want to disambiguate duplicate-
        sequence paralogs, where multiple queries share an identical
        sequence and therefore have identical tied best hits that can't be
        told apart from sequence alone (see
        map_all_hairpins_to_genomes.py's duplicate-resolution step).
    target_accession, target_assembly: if given, map against this SPECIFIC
        genome build rather than whichever is currently latest - see
        _get_genome_fasta for the resolution order (accession preferred,
        assembly name as fallback).

    Returns a pandas DataFrame with one row per successfully-mapped
    sequence: name, chromosome, start, end, strand, percent_identity,
    alignment_length, query_coverage, evalue, bitscore. Sequences with no
    hit clearing the thresholds are omitted - check for them by diffing
    against your input names. If return_raw_hits=True, returns a tuple
    (result_df, raw_hits_df) instead.
    """
    fasta_path, genome_source, chrom_map = _get_genome_fasta(
        species_key, division=division, target_accession=target_accession, target_assembly=target_assembly
    )
    db_path = _build_blast_db(fasta_path)

    items = list(sequences.items())
    all_hits = []
    for chunk_start in range(0, len(items), chunk_size):
        chunk = dict(items[chunk_start:chunk_start + chunk_size])
        print(
            f"  blastn chunk {chunk_start // chunk_size + 1}/"
            f"{-(-len(items) // chunk_size)} ({len(chunk)} sequences)...",
            file=sys.stderr,
        )
        all_hits.append(_run_blast_chunk(chunk, db_path, min_identity, min_coverage))

    hits = pd.concat(all_hits, ignore_index=True) if all_hits else pd.DataFrame()
    if hits.empty:
        print("No BLAST hits at all for any sequence.", file=sys.stderr)
        empty = pd.DataFrame(columns=["name", "chromosome", "start", "end", "strand", "genome_source"])
        return (empty, hits) if return_raw_hits else empty

    if chrom_map:
        # NCBI FASTA headers are often accession-style (e.g. "gb|CP126666.1|")
        # rather than a clean chromosome name - relabel using the assembly
        # report mapping fetched alongside the genome (see
        # _fetch_ncbi_assembly_report_chrom_map). Unmapped sequences (e.g.
        # genuine unplaced scaffolds) are left as-is.
        def _relabel(sseqid):
            token = str(sseqid)
            if "|" in token:
                parts = [p for p in token.split("|") if p]
                token = parts[-1] if parts else token
            return chrom_map.get(token, str(sseqid))

        hits["sseqid"] = hits["sseqid"].map(_relabel)

    hits = hits.sort_values(["qseqid", "bitscore", "pident"], ascending=[True, False, False])

    results = []
    for name, group in hits.groupby("qseqid", sort=False):
        best = group.iloc[0]
        if best["pident"] < min_identity or best["coverage"] < min_coverage:
            print(
                f"[{name}] Best hit below threshold "
                f"(pident={best['pident']:.1f}%, coverage={best['coverage']:.1%}, "
                f"chromosome {best['sseqid']}) - not included in results.",
                file=sys.stderr,
            )
            continue

        tied = group[group["bitscore"] == best["bitscore"]]
        if len(tied) > 1:
            print(
                f"[{name}] WARNING: {len(tied)} equally-good hits found - "
                "sequence may map to more than one genomic location "
                "(repeated element/gene family). Using the first; "
                "inspect manually if this matters for your use case.",
                file=sys.stderr,
            )

        start = int(min(best["sstart"], best["send"]))
        end = int(max(best["sstart"], best["send"]))
        strand = "+" if best["sstrand"] == "plus" else "-"

        results.append(
            {
                "name": name,
                "chromosome": str(best["sseqid"]),
                "start": start,
                "end": end,
                "strand": strand,
                "percent_identity": float(best["pident"]),
                "alignment_length": int(best["length"]),
                "query_coverage": float(best["coverage"]),
                "evalue": float(best["evalue"]),
                "bitscore": float(best["bitscore"]),
                "genome_source": genome_source,
            }
        )

    result_df = pd.DataFrame(results)
    return (result_df, hits) if return_raw_hits else result_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True, help="Ensembl production name, e.g. mus_musculus")
    parser.add_argument(
        "--division", default="vertebrates",
        choices=["vertebrates", "plants", "metazoa", "fungi", "protists"],
        help="Ensembl division this species belongs to - determines which FTP server the genome is downloaded from",
    )
    parser.add_argument("--name", required=True, help="Identifier for this sequence, e.g. mmu-mir-1966")
    parser.add_argument("--sequence", required=True, help="RNA or DNA sequence to map")
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=0.9)
    args = parser.parse_args()

    df = map_sequences_batch(
        {args.name: args.sequence},
        args.species,
        division=args.division,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
    )

    if df.empty:
        print("\nNo confident genomic location found.")
        sys.exit(1)

    print("\nResult:")
    print(df.to_string(index=False))
