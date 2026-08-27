"""
Batch-map ALL miRNA hairpin sequences (grouped by species) to genomic
coordinates on the current Ensembl reference genome for each species,
via local BLAST. Builds on map_sequence_to_genome.py's per-species
alignment logic and mirbase_organisms.txt (miRBase's own
organisms reference table, tab-separated - see
_load_mirbase_organisms_table) - keep all three files in the same
directory.

Input file format: tab-separated, columns identified by NAME from the
header (not by position) - requires 'name' and 'sequence'; optionally
'chromosome', 'start', 'end', 'strand' in any order. Any other columns
(e.g. an 'id' column) are simply ignored, and columns may appear in any
order - this is deliberately tolerant of the file's shape changing over
time. One hairpin per line. Lines are grouped by species prefix (text
before the first '-' in the name) regardless of where they appear in the
file - rows for the same species do NOT need to be consecutive, so
appending new entries (e.g. from MirGeneDB) anywhere in the file is safe.
The coordinate columns, when present, are miRBase's own last-mapped
coordinates for each hairpin (on the build recorded in
mirbase_organisms.txt) - used directly (no network fetch needed) both to
disambiguate duplicate-sequence paralogs and to carry forward coordinates
for sequences that haven't changed (see below). If the file doesn't have
the coordinate columns at all, both of those features are simply
unavailable rather than erroring.

Before mapping each species, this looks up both miRBase's own last-mapped
genome build (from mirbase_organisms.txt's genome_assembly/genome_accession
columns - a local lookup, not a network request, and covers far more
species than the ~31 that actually have a downloadable GFF3) and Ensembl's
current one (via REST /info/assembly/:species). If the build hasn't
changed, each individual sequence is compared against miRBase's own
current hairpin.fa (https://www.mirbase.org/download/hairpin.fa, fetched
once per run for all species): sequences that match exactly are NOT
remapped - their coordinate is carried forward directly from the input
file's own coordinate columns, since nothing about them could have
changed. Only sequences that are new, textually changed (e.g. a
trimmed/extended precursor boundary), or belong to a species whose build
HAS moved, are sent to BLAST. This is the main lever for avoiding
redundant work on a full miRBase-scale run - most species' sequences
won't have changed release-to-release once their build is current.

Outputs (written to disk, not the terminal):
  - mirna_coordinates_mapped.csv  - one row per successfully mapped hairpin
  - sequence_mapping_summary.csv  - one row per species block: how many
    sequences were requested/mapped, and status (ok / skipped / error)
Both are appended to incrementally, one species at a time, so partial
progress is preserved even if a later species fails or the run is
interrupted - safer for a job that may run over many species/hours.

Scale warning: a file spanning all of miRBase can include ~300 species.
Not all of them are necessarily in MIRBASE_PREFIX_TO_SPECIES yet -
unmapped prefixes are skipped and logged (not guessed), since assigning a
sequence to the wrong species/genome would be worse than a visible gap.
Downloading and indexing every covered species' genome is a substantial
one-time cost in disk space and time - consider testing on a handful of
species first with --only before a full run.

Usage:
    python map_all_hairpins_to_genomes.py hairpin_seqs.txt
    python map_all_hairpins_to_genomes.py hairpin_seqs.txt --only hsa,mmu,rno
"""

import argparse
import csv
import itertools
import json
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd

from map_sequence_to_genome import map_sequences_batch
from mirna_coords_batch_biomart import MIRBASE_PREFIX_TO_SPECIES

MIRBASE_HAIRPIN_FASTA_URL = "https://www.mirbase.org/download/hairpin.fa"
MIRBASE_ORGANISMS_TABLE_PATH = "mirbase_organisms.txt"
REST_HOST = "https://rest.ensembl.org"

_mirbase_organisms_table_cache = None  # populated once per run, on first use
_ensembl_assembly_cache = {}
_mirbase_hairpin_sequences_cache = None  # populated once per run, on first use

MIRBASE_PREFIX_TO_SPECIES = {
    # --- vertebrates (main Ensembl) ---
    "hsa": ("Homo sapiens", "vertebrates"),
    "mmu": ("Mus musculus", "vertebrates"),
    "rno": ("Rattus norvegicus", "vertebrates"),
    "dre": ("Danio rerio", "vertebrates"),
    "gga": ("Gallus gallus", "vertebrates"),
    "bta": ("Bos taurus", "vertebrates"),
    "ssc": ("Sus scrofa", "vertebrates"),
    "cfa": ("Canis lupus familiaris", "vertebrates"),
    "fca": ("Felis catus", "vertebrates"),
    "eca": ("Equus caballus", "vertebrates"),
    "oar": ("Ovis aries", "vertebrates"),
    "chi": ("Capra hircus", "vertebrates"),
    "ocu": ("Oryctolagus cuniculus", "vertebrates"),
    "mml": ("Macaca mulatta", "vertebrates"),
    "ptr": ("Pan troglodytes", "vertebrates"),
    "mdo": ("Monodelphis domestica", "vertebrates"),
    "oan": ("Ornithorhynchus anatinus", "vertebrates"),
    "xtr": ("Xenopus tropicalis", "vertebrates"),
    "ola": ("Oryzias latipes", "vertebrates"),
    "tni": ("Tetraodon nigroviridis", "vertebrates"),
    "gac": ("Gasterosteus aculeatus", "vertebrates"),
    "aca": ("Anolis carolinensis", "vertebrates"),
    "cja": ("Callithrix jacchus", "vertebrates"),
    "ssa": ("Salmo salar", "vertebrates"),
    "xla": ("Xenopus laevis", "vertebrates"),
    "age": ("Ateles geoffroyi", "vertebrates"),
    "lla": ("Lagothrix lagotricha", "vertebrates"),
    "sla": ("Saguinus labiatus", "vertebrates"),
    "mne": ("Macaca nemestrina", "vertebrates"),
    "pbi": ("Pygathrix bieti", "vertebrates"),
    "ggo": ("Gorilla gorilla", "vertebrates"),
    "ppa": ("Pan paniscus", "vertebrates"),
    "ppy": ("Pongo pygmaeus", "vertebrates"),
    "ssy": ("Symphalangus syndactylus", "vertebrates"),
    "lca": ("Lemur catta", "vertebrates"),
    "cgr": ("Cricetulus griseus", "vertebrates"),
    "fru": ("Fugu rubripes", "vertebrates"),
    "tgu": ("Taeniopygia guttata", "vertebrates"),
    "pma": ("Petromyzon marinus", "vertebrates"),
    "meu": ("Macropus eugenii", "vertebrates"),
    "sha": ("Sarcophilus harrisii", "vertebrates"),
    "pol": ("Paralichthys olivaceus", "vertebrates"),
    "hhi": ("Hippoglossus hippoglossus", "vertebrates"),
    "ccr": ("Cyprinus carpio", "vertebrates"),
    "aja": ("Artibeus jamaicensis", "vertebrates"),
    "ipu": ("Ictalurus punctatus", "vertebrates"),
    "efu": ("Eptesicus fuscus", "vertebrates"),
    "oha": ("Ophiophagus hannah", "vertebrates"),
    "tch": ("Tupaia chinensis", "vertebrates"),
    "cpi": ("Chrysemys picta", "vertebrates"),
    "ami": ("Alligator mississippiensis", "vertebrates"),
    "cli": ("Columba livia", "vertebrates"),
    "pbv": ("Python bivittatus", "vertebrates"),
    "cmi": ("Callorhinchus milii", "vertebrates"),
    "pal": ("Pteropus alecto", "vertebrates"),
    "cau": ("Carassius auratus", "vertebrates"),
    "omy": ("Oncorhynchus mykiss", "vertebrates"),
    "oni": ("Oreochromis niloticus", "vertebrates"),    
    "abu": ("Astatotilapia burtoni", "vertebrates"),
    "mze": ("Metriaclima zebra", "vertebrates"),
    "nbr": ("Neolamprologus brichardi", "vertebrates"),
    "pny": ("Pundamilia nyererei", "vertebrates"),
    "bbu": ("Bubalus bubalis", "vertebrates"),
    "eel": ("Electrophorus electricus", "vertebrates"),
    "gmo": ("Gadus morhua", "vertebrates"),
    "apl": ("Anas platyrhynchos", "vertebrates"),
    "cpo": ("Cavia porcellus", "vertebrates"),
    "tch": ("Tupaia chinensis", "vertebrates"),
    "dno": ("Dasypus novemcinctus", "vertebrates"),
    "ete": ("Echinops telfairi", "vertebrates"),
    "mmr": ("Microcebus murinus", "vertebrates"),
    "dma": ("Daubentonia madagascariensis", "vertebrates"),
    "nle": ("Nomascus leucogenys", "vertebrates"),
    "sbo": ("Saimiri boliviensis", "vertebrates"),
    "oga": ("Otolemur garnettii", "vertebrates"),
    "pha": ("Papio hamadryas", "vertebrates"),
    "nfu": ("Nothobranchius furzeri", "vertebrates"),
    "laf": ("Loxodonta africana", "vertebrates"),
    "lch": ("Latimeria chalumnae", "vertebrates"),
    "loc": ("Lepisosteus oculatus", "vertebrates"),
    "mal": ("Monopterus albus", "vertebrates"),
    "gja": ("Gekko japonicus", "vertebrates"),
    "ebu": ("Eptatretus burgeri", "vertebrates"),
    "mun": ("Microcaecilia unicolor", "vertebrates"),
    "pae": ("Pongo abelii", "vertebrates"),


    # --- metazoa (invertebrates: Ensembl Metazoa) ---
    "dme": ("Drosophila melanogaster", "metazoa"),
    "cel": ("Caenorhabditis elegans", "metazoa"),
    "ame": ("Apis mellifera", "metazoa"),
    "aga": ("Anopheles gambiae", "metazoa"),
    "aae": ("Aedes aegypti", "metazoa"),
    "bmo": ("Bombyx mori", "metazoa"),
    "tca": ("Tribolium castaneum", "metazoa"),
    "nvi": ("Nasonia vitripennis", "metazoa"),
    "dpu": ("Daphnia pulex", "metazoa"),
    "cqu": ("Culex quinquefasciatus", "metazoa"),
    "isc": ("Ixodes scapularis", "metazoa"),
    "spu": ("Strongylocentrotus purpuratus", "metazoa"),
    "cbr": ("Caenorhabditis briggsae", "metazoa"),
    "bmy": ("Brugia malayi", "metazoa"),
    "aqu": ("Amphimedon queenslandica", "metazoa"),
    "nve": ("Nematostella vectensis", "metazoa"),
    "hma": ("Hydra magnipapillata", "metazoa"),
    "sko": ("Saccoglossus kowalevskii", "metazoa"),
    "cin": ("Ciona intestinalis", "metazoa"),
    "csa": ("Ciona savignyi", "metazoa"),
    "odi": ("Oikopleura dioica", "metazoa"),
    "bfl": ("Branchiostoma floridae", "metazoa"),
    "dan": ("Drosophila ananassae", "metazoa"),
    "der": ("Drosophila erecta", "metazoa"),
    "dgr": ("Drosophila grimshawi", "metazoa"),
    "dmo": ("Drosophila mojavensis", "metazoa"),
    "dpe": ("Drosophila persimilis", "metazoa"),
    "dps": ("Drosophila pseudoobscura", "metazoa"),
    "dse": ("Drosophila sechellia", "metazoa"),
    "dsi": ("Drosophila simulans", "metazoa"),
    "dvi": ("Drosophila virilis", "metazoa"),
    "dwi": ("Drosophila willistoni", "metazoa"),
    "dya": ("Drosophila yakuba", "metazoa"),
    "lmi": ("Locusta migratoria", "metazoa"),
    "sja": ("Schistosoma japonicum", "metazoa"),
    "sma": ("Schistosoma mansoni", "metazoa"),
    "sme": ("Schmidtea mediterranea", "metazoa"),
    "cte": ("Capitella teleta", "metazoa"),
    "cla": ("Cerebratulus lacteus", "metazoa"),
    "hru": ("Haliotis rufescens", "metazoa"),
    "lgi": ("Lottia gigantea", "metazoa"),
    "crm": ("Caenorhabditis remanei", "metazoa"),
    "ppc": ("Pristionchus pacificus", "metazoa"),
    "bma": ("Brugia malayi", "metazoa"),
    "api": ("Acyrthosiphon pisum", "metazoa"),
    "ngi": ("Nasonia giraulti", "metazoa"),
    "nlo": ("Nasonia longicornis", "metazoa"),
    "smr": ("Strigamia maritima", "metazoa"),
    "bpa": ("Brugia pahangi", "metazoa"),
    "xbo": ("Xenoturbella bocki", "metazoa"),
    "egr": ("Echinococcus granulosus", "metazoa"),
    "emu": ("Echinococcus multilocularis", "metazoa"),
    "hme": ("Heliconius melpomene", "metazoa"),
    "gpy": ("Glottidia pyramidata", "metazoa"),
    "tre": ("Terebratulina retusa", "metazoa"),
    "rmi": ("Rhipicephalus microplus", "metazoa"),
    "asu": ("Ascaris suum", "metazoa"),
    "tur": ("Tetranychus urticae", "metazoa"),
    "hco": ("Haemonchus contortus", "metazoa"),
    "mse": ("Manduca sexta", "metazoa"),
    "mja": ("Marsupenaeus japonicus", "metazoa"),
    "cbn": ("Caenorhabditis brenneri", "metazoa"),
    "sci": ("Sycon ciliatum", "metazoa"),
    "lco": ("Leucosolenia complicata", "metazoa"),   
    "prd": ("Panagrellus redivivus", "metazoa"),
    "bbe": ("Branchiostoma belcheri", "metazoa"),
    "pmi": ("Patiria miniata", "metazoa"),
    "lva": ("Lytechinus variegatus", "metazoa"),
    "str": ("Strongyloides ratti", "metazoa"),
    "pxy": ("Plutella xylostella", "metazoa"),
    "gsa": ("Gyrodactylus salaris", "metazoa"),
    "mde": ("Mayetiola destructor", "metazoa"),
    "ago": ("Aphis gossypii", "metazoa"),
    "hpo": ("Heligmosomoides polygyrus", "metazoa"),
    "tcf": ("Triops cancriformis", "metazoa"),
    "sfr": ("Spodoptera frugiperda", "metazoa"),
    "fhe": ("Fasciola hepatica", "metazoa"),
    "bib": ("Biston betularia", "metazoa"),
    "bdo": ("Bactrocera dorsalis", "metazoa"),
    "dqu": ("Dinoponera quadriceps", "metazoa"),
    "pca": ("Polistes canadensis", "metazoa"),
    "pte": ("Parasteatoda tepidariorum", "metazoa"),
    "mco": ("Mesocestoides corti", "metazoa"),
    "mle": ("Melibe leonina", "metazoa"),
    "apa": ("Aiptasia pallida", "metazoa"),
    "phw": ("Parhyale hawaiensis", "metazoa"),
    "bla": ("Branchiostoma lanceolatum", "metazoa"),
    "pfl": ("Ptychodera flava", "metazoa"),
    "agr": ("Acanthopleura granulata", "metazoa"),
    "ofu": ("Owenia fusiformis", "metazoa"),
    "mms": ("Mopalia muscosa", "metazoa"),
    "dlo": ("Diachasmimorpha longicaudata", "metazoa"),
    "pau": ("Phoronis australis", "metazoa"),
    "pve": ("Pictodentalium vernedei", "metazoa"),
    "hmi": ("Hofstenia miamia", "metazoa"),
    "sro": ("Symsagittifera roscoffensis", "metazoa"),
    "mgi": ("Magallana gigas", "metazoa"),
    "dgy": ("Dimorphilus gyrociliatus", "metazoa"),
    "nag": ("Neoechinorhynchus agilis", "metazoa"),
    "hpi": ("Halcurias pilatus", "metazoa"),
    "efe": ("Eisenia fetida", "metazoa"),
    "eba": ("Epimenia babai", "metazoa"),
    "llo": ("Lineus longissimus", "metazoa"),
    "bpl": ("Brachionus plicatilis", "metazoa"),
    "ava": ("Adineta vaga", "metazoa"),
    "sne": ("Seison nebaliae", "metazoa"),
    "war": ("Wirenia argentea", "metazoa"),
    "adi": ("Acropora digitifera", "metazoa"),
    "ovu": ("Octopus vulgaris", "metazoa"),
    "gpa": ("Glossina pallidipes", "metazoa"),
    "ple": ("Pomphorhynchus laevis", "metazoa"),
    "hvl": ("Hydra vulgaris", "metazoa"),
    "csc": ("Centruroides sculpturatus", "metazoa"),
    "lpo": ("Limulus polyphemus", "metazoa"),
    "lhy": ("Laevipilina hyalina", "metazoa"),
    "bko": ("Brachionus koreanus", "metazoa"),
    "pcr": ("Prostheceraeus crozeri", "metazoa"),
    "rph": ("Ruditapes philippinarum", "metazoa"),
    "eml": ("Ephydatia muelleri", "metazoa"),
    "esc": ("Euprymna scolopes", "metazoa"),
    "obi": ("Octopus bimaculoides", "metazoa"),
    "bge": ("Blattella germanica", "metazoa"),
    "npo": ("Nautilus pompilius", "metazoa"),
    "dmg": ("Daphnia magna", "metazoa"),
    "snu": ("Sipunculus nudus", "metazoa"),
    "lan": ("Lingula anatina", "metazoa"),
    "pdu": ("Platynereis dumerilii", "metazoa"),


    # --- plants (Ensembl Plants) ---
    "ath": ("Arabidopsis thaliana", "plants"),
    "osa": ("Oryza sativa", "plants"),
    "zma": ("Zea mays", "plants"),
    "gma": ("Glycine max", "plants"),
    "vvi": ("Vitis vinifera", "plants"),
    "sbi": ("Sorghum bicolor", "plants"),
    "stu": ("Solanum tuberosum", "plants"),
    "sly": ("Solanum lycopersicum", "plants"),
    "mtr": ("Medicago truncatula", "plants"),
    "ptc": ("Populus trichocarpa", "plants"),
    "bdi": ("Brachypodium distachyon", "plants"),
    "bna": ("Brassica napus", "plants"),
    "csi": ("Citrus sinensis", "plants"),
    "hvu": ("Hordeum vulgare", "plants"),
    "pta": ("Pinus taeda", "plants"),
    "ppt": ("Physcomitrella patens", "plants"),
    "smo": ("Selaginella moellendorffii", "plants"),
    "bna": ("Brassica napus", "plants"),
    "bol": ("Brassica oleracea", "plants"),
    "bra": ("Brassica rapa", "plants"),
    "cpa": ("Carica papaya", "plants"),
    "lja": ("Lotus japonicus", "plants"),
    "vun": ("Vigna unguiculata", "plants"),
    "ghb": ("Gossypium herbaceum", "plants"),
    "ghr": ("Gossypium hirsutum", "plants"),
    "gra": ("Gossypium raimondii", "plants"),
    "sof": ("Saccharum officinarum", "plants"),
    "tae": ("Triticum aestivum", "plants"),
    "pvu": ("Phaseolus vulgaris", "plants"),
    "mdm": ("Malus domestica", "plants"),
    "aqc": ("Aquilegia caerulea", "plants"),
    "peu": ("Populus euphratica", "plants"),
    "ccl": ("Citrus clementina", "plants"),
    "crt": ("Citrus reticulata", "plants"),
    "ctr": ("Citrus trifoliata", "plants"),
    "rco": ("Ricinus communis", "plants"),
    "gar": ("Gossypium arboreum", "plants"),
    "aly": ("Arabidopsis lyrata", "plants"),
    "ahy": ("Arachis hypogaea", "plants"),
    "gso": ("Glycine soja", "plants"),
    "pab": ("Picea abies", "plants"),
    "ttu": ("Triticum turgidum", "plants"),    
    "ata": ("Aegilops tauschii", "plants"),
    "far": ("Festuca arundinacea", "plants"),
    "tcc": ("Theobroma cacao", "plants"),
    "rgl": ("Rehmannia glutinosa", "plants"),
    "ssp": ("Saccharum sp.", "plants"),
    "bgy": ("Bruguiera gymnorhiza", "plants"),
    "bcy": ("Bruguiera cylindrica", "plants"),
    "cme": ("Cucumis melo", "plants"),
    "amg": ("Acacia mangium", "plants"),
    "aau": ("Acacia auriculiformis", "plants"),
    "ssl": ("Salvia sclarea", "plants"),
    "dpr": ("Digitalis purpurea", "plants"),
    "nta": ("Nicotiana tabacum", "plants"),
    "egu": ("Elaeis guineensis", "plants"),
    "mes": ("Manihot esculenta", "plants"),
    "cca": ("Cynara cardunculus", "plants"),
    "hbr": ("Hevea brasiliensis", "plants"),
    "pde": ("Pinus densata", "plants"),
    "han": ("Helianthus annuus", "plants"),
    "hci": ("Helianthus ciliaris", "plants"),
    "htu": ("Helianthus tuberosus", "plants"),
    "hex": ("Helianthus exilis", "plants"),
    "har": ("Helianthus argophyllus", "plants"),
    "hpe": ("Helianthus petiolaris", "plants"),
    "hpa": ("Helianthus paradoxus", "plants"),
    "cln": ("Cunninghamia lanceolata", "plants"),
    "lus": ("Linum usitatissimum", "plants"),
    "pgi": ("Panax ginseng", "plants"),
    "ppe": ("Prunus persica", "plants"),
    "ama": ("Avicennia marina", "plants"),
    "atr": ("Amborella trichopoda", "plants"),
    "pso": ("Papaver somniferum", "plants"),
    "pmu": ("Prunus mume", "plants"),
    "sit": ("Setaria italica", "plants"),
    "vca": ("Vriesea carinata", "plants"),
    "eun": ("Eugenia uniflora", "plants"),
    "seu": ("Salicornia europaea", "plants"),
    "fve": ("Fragaria vesca", "plants"),
    "cst": ("Cucumis sativus", "plants"),
    "cas": ("Camelina sativa", "plants"),
    "pla": ("Paeonia lactiflora", "plants"),
    "mpo": ("Marchantia polymorpha", "plants"),
    "smi": ("Salvia miltiorrhiza", "plants"),
    "aof": ("Asparagus officinalis", "plants"),


    # --- protists (Ensembl Protists) ---
    "cre": ("Chlamydomonas reinhardtii", "protists"),
    "ddi": ("Dictyostelium discoideum", "protists"),
    "smc": ("Symbiodinium microadriaticum", "protists"),
    "esi": ("Ectocarpus siliculosus", "protists"),
    "pin": ("Phytophthora infestans", "protists"),
    "psj": ("Phytophthora sojae", "protists"),
    "pra": ("Phytophthora ramorum", "protists"),
    "pti": ("Phaeodactylum tricornutum", "protists"),
}

def _load_mirbase_organisms_table(path=MIRBASE_ORGANISMS_TABLE_PATH):
    """
    Load miRBase's own organisms reference table (tab-separated, columns
    include name_abbr, name, genome_assembly, genome_accession, gff) once
    per run. This is the authoritative source for genome-build comparison -
    a local lookup covering ~285 species, versus the per-species GFF3
    header fetch it replaces, which only worked for the ~31 species that
    actually have a downloadable GFF3 (and required a network round-trip
    per species even for those).

    Returns dict {prefix: {"build_id": str|None, "build_accession": str|None,
    "has_gff": bool, "scientific_name": str}}.
    """
    global _mirbase_organisms_table_cache
    if _mirbase_organisms_table_cache is not None:
        return _mirbase_organisms_table_cache

    table = {}
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                prefix = row["name_abbr"].strip()
                assembly = row["genome_assembly"].strip()
                accession = row["genome_accession"].strip()
                table[prefix] = {
                    "build_id": None if assembly == "NULL" else assembly,
                    # format is "NCBI_Assembly:GCA_000001405.15" - keep just the accession
                    "build_accession": None if accession == "NULL" else accession.split(":")[-1].strip(),
                    "has_gff": row["gff"].strip() == "1",
                    "scientific_name": row["name"].strip(),
                }
    except FileNotFoundError:
        print(
            f"[warn] miRBase organisms table not found at '{path}' - genome-build "
            "comparison will be skipped for every species (everything will be "
            "remapped). Place mirbase_organisms.txt alongside this script to enable it.",
            file=sys.stderr,
        )

    _mirbase_organisms_table_cache = table
    return table


def _fetch_mirbase_genome_build(prefix):
    """
    Look up the genome build miRBase last mapped this species against,
    from the local organisms reference table (see
    _load_mirbase_organisms_table) rather than fetching a GFF3 header over
    HTTP - faster, and covers species that don't have a GFF3 file at all.

    Returns (build_id, build_accession) - either may be None if the prefix
    isn't in the table or that field wasn't recorded for it.
    """
    table = _load_mirbase_organisms_table()
    entry = table.get(prefix)
    if entry is None:
        return None, None
    return entry["build_id"], entry["build_accession"]


def _load_species_build_overrides(path):
    """
    Parse a target-builds override file: tab-separated lines of the form
    'UPDATE  SPECIES  <scientific name>  <field>  <value>', typically two
    lines per species (genome_assembly and genome_accession), though a
    species may only have one of the two given.

    Returns dict {scientific_name_as_written: {"assembly": str|None, "accession": str|None}}
    """
    overrides = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5 or parts[0] != "UPDATE" or parts[1] != "SPECIES":
                continue
            species_name, field, value = parts[2].strip(), parts[3].strip(), parts[4].strip()
            entry = overrides.setdefault(species_name, {"assembly": None, "accession": None})
            if field == "genome_assembly":
                entry["assembly"] = value
            elif field == "genome_accession":
                # format "NCBI_Assembly:GCA_..." - occasionally with a stray
                # space after the colon (seen in real data) - strip either way
                entry["accession"] = value.split(":")[-1].strip()
    return overrides


def _resolve_species_overrides(path):
    """
    Load a target-builds override file and resolve each species name to
    its miRBase prefix, via mirbase_organisms.txt's own scientific-name
    column - NOT the hand-curated MIRBASE_PREFIX_TO_SPECIES dict, since
    naming can differ between sources (e.g. this kind of file may use
    "Canis familiaris" while other references use "Canis lupus familiaris" -
    mirbase_organisms.txt uses the former, matching miRBase's own prefix
    convention, so it's the authoritative match target here). Species that
    can't be matched are reported clearly rather than silently dropped.

    Returns dict {prefix: {"assembly": str|None, "accession": str|None,
    "scientific_name": str}}
    """
    raw_overrides = _load_species_build_overrides(path)

    organisms_table = _load_mirbase_organisms_table()
    name_to_prefix = {
        entry["scientific_name"].strip().lower(): prefix
        for prefix, entry in organisms_table.items()
    }

    resolved = {}
    for species_name, values in raw_overrides.items():
        prefix = name_to_prefix.get(species_name.strip().lower())
        if prefix is None:
            print(
                f"[warn] '{species_name}' from {path} could not be matched to a "
                "miRBase prefix (checked mirbase_organisms.txt's scientific "
                "names) - this species will be SKIPPED. Check spelling/synonyms "
                "against mirbase_organisms.txt's 'name' column.",
                file=sys.stderr,
            )
            continue
        resolved[prefix] = {**values, "scientific_name": species_name}

    return resolved



def _fetch_mirbase_hairpin_sequences():
    """
    Download and parse miRBase's complete hairpin.fa (all species, one
    file) ONCE per run, caching the result in memory. Used to detect which
    individual sequences have actually changed since this miRBase release
    for species where the genome build hasn't changed - so only those need
    remapping, rather than every sequence in the species.

    Returns dict {name: sequence} with sequences normalized (uppercase,
    U->T) for direct comparison against the input file's sequences, which
    go through the same normalization before mapping.
    """
    global _mirbase_hairpin_sequences_cache
    if _mirbase_hairpin_sequences_cache is not None:
        return _mirbase_hairpin_sequences_cache

    print(f"Fetching miRBase's current hairpin.fa (all species, one-time download): {MIRBASE_HAIRPIN_FASTA_URL}", file=sys.stderr)
    sequences = {}
    try:
        req = Request(MIRBASE_HAIRPIN_FASTA_URL, headers={"User-Agent": "Mozilla/5.0"})
        with urlopen(req, timeout=300) as resp:
            current_name = None
            current_seq_parts = []
            for raw_line in resp:
                line = raw_line.decode("utf-8", errors="replace").rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if current_name is not None:
                        sequences[current_name] = "".join(current_seq_parts).upper().replace("U", "T")
                    current_name = line[1:].split()[0]  # first token after '>' is the name
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line.strip())
            if current_name is not None:
                sequences[current_name] = "".join(current_seq_parts).upper().replace("U", "T")
    except (HTTPError, URLError, TimeoutError, OSError) as exc:
        print(f"[warn] Could not fetch miRBase hairpin.fa: {exc} - sequence-change detection disabled, everything will be remapped.", file=sys.stderr)

    print(f"Parsed {len(sequences)} sequences from miRBase hairpin.fa.", file=sys.stderr)
    _mirbase_hairpin_sequences_cache = sequences
    return sequences


def _build_matches(mirbase_build_id, mirbase_accession, ensembl_name, ensembl_accession):
    """
    True if miRBase's last-mapped build matches Ensembl's current one.
    Used to decide whether per-sequence change detection applies at all -
    if the build itself has changed, every sequence needs remapping
    regardless of whether its text matches hairpin.fa, since the whole
    coordinate system moved. Prefers comparing INSDC accessions directly;
    falls back to a name-prefix comparison if an accession isn't
    available on either side (handles cosmetic patch-version suffixes,
    e.g. miRBase says "GRCh38" while Ensembl reports "GRCh38.p14").
    """
    if mirbase_accession and ensembl_accession:
        return mirbase_accession == ensembl_accession
    if mirbase_build_id and ensembl_name:
        return ensembl_name.lower().startswith(mirbase_build_id.lower())
    return False


def _normalize_chrom(chrom):
    """miRBase sometimes uses 'chr19', Ensembl uses '19' - strip common
    prefixes ('chr'/'Chr', 'scaffold_'/'Scaffold') so different naming
    conventions for the same sequence can be compared directly. This does
    NOT help when old and new coordinates are on genuinely different
    assembly structures (e.g. an old scaffold-level draft vs a new
    chromosome-level assembly) - see _resolve_duplicate_group's fallback
    for that case."""
    result = str(chrom).strip()
    for prefix in ("chr", "Chr", "scaffold_", "Scaffold_", "scaffold", "Scaffold"):
        if result.startswith(prefix):
            result = result[len(prefix):]
            break
    return result.strip()


def _resolve_duplicate_group(names, hits_for_shared_sequence, old_coords):
    """
    Assign each of several identical-sequence paralogs to a distinct BLAST
    hit, using miRBase's own (old-assembly) coordinates as an anchor.

    Two tiers are tried, in order:
      1. Chromosome-consistent: each name matched to whichever new hit is
         on the SAME chromosome (after normalization) and closest in
         position to where that name used to be. This is the confident
         case - returned with method='chromosome_anchor'.
      2. Order fallback: used only if tier 1 finds nothing at all, which
         happens when old and new coordinates are on genuinely different
         assembly structures (e.g. the old coordinate is on a fragmented
         scaffold-level draft assembly that doesn't correspond 1:1 to the
         new chromosome-level one - no naming fix can bridge that, since
         it's a real restructuring, not just a naming difference). Falls
         back to matching sort order (earliest old position assigned to
         earliest new hit, and so on) - a reasonable bet, but a much
         weaker guarantee than tier 1, so it's flagged distinctly as
         method='order_fallback' so callers/users can tell the difference
         and review these specifically if it matters for their use case.

    names: list of miRBase names sharing an identical sequence
    hits_for_shared_sequence: DataFrame of all tied best-bitscore BLAST
        hits for that shared sequence
    old_coords: dict {name: (chrom, start, end, strand)} - miRBase's own
        coordinates, read directly from the input file's coordinate
        columns (see _read_species_sequence_blocks)

    Returns (assignment, method, reason). assignment is
    {name: hit_row (a pandas Series)} and method is 'chromosome_anchor' or
    'order_fallback' on success (reason is None). On failure, assignment
    is {} and method is None, with reason explaining specifically why:
    either BLAST found fewer genomic copies than names to assign, or at
    least one name had no coordinate in the input file to anchor against
    - callers should flag these for manual review rather than guessing.
    """
    if len(hits_for_shared_sequence) < len(names):
        return {}, None, (
            f"BLAST found only {len(hits_for_shared_sequence)} distinct genomic "
            f"cop{'y' if len(hits_for_shared_sequence) == 1 else 'ies'}, but "
            f"{len(names)} names need assigning"
        )
    missing_anchor = [n for n in names if n not in old_coords]
    if missing_anchor:
        return {}, None, f"no input-file coordinate for: {missing_anchor}"

    hit_rows = list(hits_for_shared_sequence.itertuples())

    # Tier 1: chromosome-consistent, minimum total distance.
    #
    # For small groups (miRNA paralog counts are normally a handful, 2-4),
    # brute-force over every permutation is cheap and exact. But this is
    # FACTORIAL in group size - a highly-conserved miRNA family with, say,
    # 10 genomic copies would mean 10! (3.6 million) combinations, which
    # would look exactly like a hung process (no progress printed inside
    # this loop at all). Confirmed this is a real risk, not theoretical -
    # capped here: above the threshold, use a greedy nearest-available
    # assignment instead (still requires chromosome match, just doesn't
    # guarantee the GLOBALLY minimal total distance across the whole
    # group - a practical trade for staying fast).
    MAX_EXACT_PERMUTATION_SIZE = 8  # 8! = 40,320 - fast; grows very fast beyond this

    best_assignment = None
    if len(names) <= MAX_EXACT_PERMUTATION_SIZE:
        best_total_distance = None
        for hit_subset in itertools.permutations(hit_rows, len(names)):
            total_distance = 0
            valid = True
            for name, hit in zip(names, hit_subset):
                old_chrom, old_start, _old_end, _old_strand = old_coords[name]
                if _normalize_chrom(old_chrom) != _normalize_chrom(hit.sseqid):
                    valid = False
                    break
                total_distance += abs(old_start - min(hit.sstart, hit.send))
            if not valid:
                continue
            if best_total_distance is None or total_distance < best_total_distance:
                best_total_distance = total_distance
                best_assignment = dict(zip(names, hit_subset))
    else:
        print(
            f"[warn] Duplicate group of {len(names)} names ({names[:3]}...) exceeds the "
            f"exact-search size limit ({MAX_EXACT_PERMUTATION_SIZE}) - using a faster greedy "
            "chromosome-matched assignment instead of the exhaustive search (avoids a "
            "factorial-time hang for large paralog families).",
            file=sys.stderr,
        )
        remaining_hits = list(hit_rows)
        assignment = {}
        # Process names in order of how confidently they can be placed
        # (smallest candidate pool first), greedily taking the closest
        # same-chromosome hit available at each step.
        for name in sorted(names, key=lambda n: old_coords[n][1]):
            old_chrom, old_start, _old_end, _old_strand = old_coords[name]
            same_chrom_hits = [h for h in remaining_hits if _normalize_chrom(h.sseqid) == _normalize_chrom(old_chrom)]
            if not same_chrom_hits:
                assignment = None
                break
            best_hit = min(same_chrom_hits, key=lambda h: abs(old_start - min(h.sstart, h.send)))
            assignment[name] = best_hit
            remaining_hits.remove(best_hit)
        best_assignment = assignment

    if best_assignment:
        return best_assignment, "chromosome_anchor", None

    # Tier 2: order fallback - chromosome naming didn't line up at all
    # (e.g. scaffold vs chromosome restructuring), so fall back to
    # matching sort order instead of giving up outright.
    names_by_old_position = sorted(names, key=lambda n: old_coords[n][1])
    hits_by_new_position = sorted(hit_rows, key=lambda h: min(h.sstart, h.send))
    fallback_assignment = dict(zip(names_by_old_position, hits_by_new_position))
    return fallback_assignment, "order_fallback", None


def _find_duplicate_sequence_groups(sequences):
    """Group names that share an identical sequence (case/U-T normalized).
    Returns a list of name-lists, only for groups with more than one
    member - unique sequences need no disambiguation."""
    by_sequence = {}
    for name, seq in sequences.items():
        key = seq.upper().replace("U", "T")
        by_sequence.setdefault(key, []).append(name)
    return [names for names in by_sequence.values() if len(names) > 1]


def _resolve_duplicates_in_result(result_df, sequences, raw_hits, coords):
    """
    Post-process a species' mapping result: for every group of names that
    share an identical sequence (and therefore currently have identical,
    colliding coordinates from the default best-hit-per-query logic),
    reassign each name to a distinct genomic locus using
    _resolve_duplicate_group, or flag the rows as needing manual review if
    that isn't possible.

    coords: {name: (chromosome, start, end, strand)} - miRBase's own
        coordinates for this species, read directly from the input file
        (see _read_species_sequence_blocks) and used as the anchor for
        disambiguation. No network fetch needed.

    Adds a 'disambiguation_status' column to every row:
      - 'unique_sequence'            - no duplicate, nothing to resolve
      - 'resolved_via_mirbase_anchor' - successfully disambiguated
      - 'AMBIGUOUS_DUPLICATE - ...'   - could not be resolved automatically
    """
    result_df = result_df.copy()
    result_df["disambiguation_status"] = "unique_sequence"

    for names in _find_duplicate_sequence_groups(sequences):
        names_in_result = [n for n in names if n in set(result_df["name"])]
        if len(names_in_result) < 2:
            continue  # only one (or none) of this group actually got mapped

        group_hits = raw_hits[raw_hits["qseqid"].isin(names_in_result)]
        if group_hits.empty:
            continue

        # Take the top N distinct genomic loci by bitscore (N = group
        # size) - NOT just hits tied at the exact top bitscore. Confirmed
        # this was the main cause of unnecessary AMBIGUOUS_DUPLICATE
        # flags: paralogous genomic copies of an identical query sequence
        # very often score slightly differently (a single SNP/indel at
        # one locus vs another is enough), so requiring an exact score
        # tie silently excluded real second/third copies and forced
        # ambiguity even with a perfectly good anchor coordinate on hand.
        distinct_loci = group_hits.drop_duplicates(subset=["sseqid", "sstart", "send"])
        tied = distinct_loci.sort_values("bitscore", ascending=False).head(len(names_in_result))

        old_coords = {n: coords[n] for n in names_in_result if n in coords}
        resolved, method, reason = _resolve_duplicate_group(names_in_result, tied, old_coords)

        status_by_method = {
            "chromosome_anchor": "resolved_via_mirbase_anchor",
            "order_fallback": (
                "resolved_via_order_fallback - LOWER CONFIDENCE (old coordinate's "
                "chromosome/scaffold naming didn't match the new assembly at all, "
                "e.g. old scaffold-level vs new chromosome-level - assigned by "
                "matching sort order instead, worth spot-checking)"
            ),
        }

        for name in names_in_result:
            idx = result_df.index[result_df["name"] == name]
            if name in resolved:
                hit = resolved[name]
                result_df.loc[idx, "chromosome"] = str(hit.sseqid)
                result_df.loc[idx, "start"] = int(min(hit.sstart, hit.send))
                result_df.loc[idx, "end"] = int(max(hit.sstart, hit.send))
                result_df.loc[idx, "strand"] = "+" if hit.sstrand == "plus" else "-"
                result_df.loc[idx, "percent_identity"] = float(hit.pident)
                result_df.loc[idx, "alignment_length"] = int(hit.length)
                result_df.loc[idx, "query_coverage"] = float(hit.coverage)
                result_df.loc[idx, "evalue"] = float(hit.evalue)
                result_df.loc[idx, "bitscore"] = float(hit.bitscore)
                result_df.loc[idx, "disambiguation_status"] = status_by_method[method]
            else:
                other_names = [n for n in names_in_result if n != name]
                result_df.loc[idx, "disambiguation_status"] = (
                    f"AMBIGUOUS_DUPLICATE - needs manual review (identical sequence to "
                    f"{other_names}) - reason: {reason}"
                )

    return result_df


def _fetch_ensembl_assembly_info(species_key):
    """Returns (assembly_name, assembly_accession) for the CURRENT Ensembl
    assembly of this species, or (None, None) if the lookup fails."""
    if species_key in _ensembl_assembly_cache:
        return _ensembl_assembly_cache[species_key]

    url = f"{REST_HOST}/info/assembly/{species_key}?content-type=application/json"
    assembly_name, assembly_accession = None, None
    try:
        req = Request(url, headers={"Content-Type": "application/json"})
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
        assembly_name = data.get("assembly_name")
        assembly_accession = data.get("assembly_accession")
    except (HTTPError, URLError, TimeoutError, ValueError) as exc:
        print(f"[warn] Could not fetch Ensembl assembly info for '{species_key}': {exc}", file=sys.stderr)

    _ensembl_assembly_cache[species_key] = (assembly_name, assembly_accession)
    return assembly_name, assembly_accession


COORDS_OUTPUT = "mirna_coordinates_mapped.csv"
SUMMARY_OUTPUT = "sequence_mapping_summary.csv"


def _read_species_sequence_blocks(filepath):
    """
    Parse the tab-separated input file into (prefix, sequences, coords)
    blocks - one block per distinct species prefix, aggregating every
    matching line regardless of where it appears in the file. This does
    NOT require same-species lines to be consecutive (an earlier version
    did, for streaming efficiency, but that assumption is fragile - e.g.
    newly appended entries from a different source, like MirGeneDB
    additions dropped in at the end of the file, would land far from
    their species' main block and silently create extra fragmented,
    separately-processed blocks for the same species - confirmed this
    actually happens in practice). Aggregating fully by prefix first means
    file ordering never matters, at the cost of reading the whole file
    into memory before yielding anything - fine at this scale (tens of
    MB, not GB).

    Columns are read by NAME from the header, not by position - required
    columns are 'name' and 'sequence'; 'chromosome'/'start'/'end'/'strand'
    are read if present, in whatever order, alongside any other columns
    (e.g. an 'id' column) which are simply ignored. This is deliberately
    robust to the file gaining/losing/reordering columns - a positional
    assumption here previously broke silently and completely when an
    'id' column was added to the front of the file (every row's 'id'
    value was read as if it were the miRNA name, since it was in column
    0 - confirmed this happened in practice, not hypothetical).

    sequences: {name: sequence}
    coords: {name: (chromosome, start, end, strand)} - miRBase's own
        coordinates on the build recorded in mirbase_organisms.txt, read
        directly from the input file's chromosome/start/end/strand
        columns if present. This replaces the earlier approach of
        fetching these from miRBase's GFF3 over HTTP on demand - no
        network call needed, and available for every entry (not just the
        ~31 species with a downloadable GFF3). If the input file doesn't
        have these columns at all, coords will be empty for every name
        and callers fall back to the old behavior (disambiguation/carry-
        forward simply won't be possible without this data). A given
        row's coordinate is simply omitted (not an error) if that row's
        coordinate columns are missing/incomplete.
    """
    sequences_by_prefix = {}
    coords_by_prefix = {}
    prefix_order = []  # preserve first-seen order for deterministic output

    with open(filepath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []

        if "name" not in fieldnames or "sequence" not in fieldnames:
            raise ValueError(
                f"'{filepath}' must have 'name' and 'sequence' columns (found: {fieldnames})"
            )
        has_coords = all(col in fieldnames for col in ("chromosome", "start", "end", "strand"))

        for row in reader:
            name = (row.get("name") or "").strip()
            seq = (row.get("sequence") or "").strip()
            if not name:
                continue
            prefix = name.split("-")[0]

            if prefix not in sequences_by_prefix:
                sequences_by_prefix[prefix] = {}
                coords_by_prefix[prefix] = {}
                prefix_order.append(prefix)

            sequences_by_prefix[prefix][name] = seq

            if has_coords:
                chrom = (row.get("chromosome") or "").strip()
                start = (row.get("start") or "").strip()
                end = (row.get("end") or "").strip()
                strand = (row.get("strand") or "").strip()
                if chrom and start and end and strand:
                    try:
                        coords_by_prefix[prefix][name] = (chrom, int(start), int(end), strand)
                    except ValueError:
                        pass  # leave this row's coordinate out rather than erroring

    for prefix in prefix_order:
        yield prefix, sequences_by_prefix[prefix], coords_by_prefix[prefix]


def _append_csv(df, path):
    """Append to an output CSV, writing the header only the first time."""
    file_exists = Path(path).exists()
    df.to_csv(path, mode="a", header=not file_exists, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="TSV file: name, sequence, grouped by species prefix")
    parser.add_argument(
        "--only",
        help="Comma-separated list of prefixes to process (e.g. hsa,mmu,rno) - useful for testing on a subset first",
    )
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=0.9)
    parser.add_argument(
        "--fresh-start",
        action="store_true",
        help="Delete any existing output files before running, instead of appending to them "
             "(use this if rerunning from scratch rather than resuming/extending a prior run)",
    )
    parser.add_argument(
        "--target-builds",
        help="Path to a species-build override file (tab-separated 'UPDATE SPECIES <name> "
             "<field> <value>' lines - see update_species_statements.txt format). When given, "
             "ONLY the species listed in this file are processed, and each is mapped against "
             "the SPECIFIC genome build/accession given for it, rather than whichever build is "
             "currently latest. Every other species in the input file is skipped entirely.",
    )
    args = parser.parse_args()

    only_prefixes = set(args.only.split(",")) if args.only else None

    species_overrides = None
    if args.target_builds:
        species_overrides = _resolve_species_overrides(args.target_builds)
        print(
            f"Target-builds mode: only processing {len(species_overrides)} species resolved "
            f"from {args.target_builds}; every other species will be skipped.",
        )

    if args.fresh_start:
        for path in (COORDS_OUTPUT, SUMMARY_OUTPUT):
            if Path(path).exists():
                Path(path).unlink()

    already_done = set()
    if Path(SUMMARY_OUTPUT).exists():
        prior_summary = pd.read_csv(SUMMARY_OUTPUT)
        already_done = set(prior_summary["species_prefix"].astype(str))
        if already_done:
            print(
                f"Found existing {SUMMARY_OUTPUT} with {len(already_done)} species already "
                "processed - they'll be skipped. Use --fresh-start to reprocess everything.",
                file=sys.stderr,
            )

    for prefix, sequences, coords in _read_species_sequence_blocks(args.input_file):
        if only_prefixes and prefix not in only_prefixes:
            continue
        if prefix in already_done:
            continue

        n_requested = len(sequences)
        entry = MIRBASE_PREFIX_TO_SPECIES.get(prefix)

        if entry is None:
            print(f"[skip] '{prefix}' ({n_requested} sequences): not in MIRBASE_PREFIX_TO_SPECIES.", file=sys.stderr)
            _append_csv(
                pd.DataFrame([{
                    "species_prefix": prefix, "scientific_name": None, "division": None,
                    "n_requested": n_requested, "n_mapped": 0,
                    "status": "unknown prefix - not in MIRBASE_PREFIX_TO_SPECIES",
                }]),
                SUMMARY_OUTPUT,
            )
            continue

        if species_overrides is not None and prefix not in species_overrides:
            print(f"[skip] '{prefix}' ({n_requested} sequences): not in --target-builds list.", file=sys.stderr)
            _append_csv(
                pd.DataFrame([{
                    "species_prefix": prefix, "scientific_name": entry[0], "division": entry[1],
                    "n_requested": n_requested, "n_mapped": 0,
                    "status": "skipped - not in --target-builds list",
                }]),
                SUMMARY_OUTPUT,
            )
            continue

        try:
            _process_one_species(prefix, sequences, coords, entry, species_overrides, args)
        except Exception as exc:
            print(f"[error] '{prefix}' failed unexpectedly: {exc}", file=sys.stderr)
            _append_csv(
                pd.DataFrame([{
                    "species_prefix": prefix, "scientific_name": entry[0], "division": entry[1],
                    "n_requested": n_requested, "n_mapped": 0,
                    "status": f"unexpected error: {exc}",
                }]),
                SUMMARY_OUTPUT,
            )

    print(f"\nAll done. Coordinates: {COORDS_OUTPUT}  Summary: {SUMMARY_OUTPUT}")


def _process_one_species(prefix, sequences, coords, entry, species_overrides, args):
    """One species' full pipeline: build comparison, unchanged-vs-changed
    split, BLAST mapping, duplicate resolution, and writing results.
    Factored out so the caller can wrap the whole thing in a single
    try/except - any unexpected failure here (not just a BLAST error,
    which has its own handling) still produces a summary row instead of
    crashing the entire run silently."""
    n_requested = len(sequences)
    scientific_name, division = entry
    species_key = scientific_name.strip().lower().replace(" ", "_")

    target_accession, target_assembly = None, None
    mirbase_build_id, mirbase_accession = _fetch_mirbase_genome_build(prefix)
    if species_overrides is not None:
        override = species_overrides[prefix]
        target_assembly, target_accession = override["assembly"], override["accession"]
        # Compare against the TARGET build (from the override file),
        # not Ensembl's live current build - the whole point of
        # --target-builds is to pin to a specific build rather than
        # "whatever's latest".
        ensembl_name, ensembl_accession = target_assembly, target_accession
    else:
        ensembl_name, ensembl_accession = _fetch_ensembl_assembly_info(species_key)

    build_note = (
        f"miRBase: {mirbase_accession or mirbase_build_id or 'unknown'}, "
        f"{'target' if species_overrides is not None else 'Ensembl'}: "
        f"{ensembl_accession or ensembl_name or 'unknown'}"
    )
    build_matches = _build_matches(mirbase_build_id, mirbase_accession, ensembl_name, ensembl_accession)

    # Only worth comparing individual sequences if the build hasn't
    # changed - if it has, every sequence needs remapping regardless
    # of whether its text matches hairpin.fa, since the coordinate
    # system itself moved.
    unchanged_names = set()
    if build_matches:
        hairpin_reference = _fetch_mirbase_hairpin_sequences()
        for name, seq in sequences.items():
            normalized = seq.upper().replace("U", "T")
            if hairpin_reference.get(name) == normalized:
                unchanged_names.add(name)

    to_map = {name: seq for name, seq in sequences.items() if name not in unchanged_names}

    print(
        f"[process] {prefix} -> {scientific_name} ({n_requested} sequences) [{build_note}] - "
        f"{len(unchanged_names)} unchanged (carried forward), {len(to_map)} to map via BLAST"
    )

    carried_forward_rows = []
    if unchanged_names:
        for name in unchanged_names:
            if name not in coords:
                # couldn't confirm the still-valid coordinate - safer to
                # map it than to silently omit it from the output
                to_map[name] = sequences[name]
                continue
            chrom, start, end, strand = coords[name]
            carried_forward_rows.append({
                "name": name,
                "chromosome": _normalize_chrom(chrom),
                "start": start,
                "end": end,
                "strand": strand,
                "percent_identity": None,
                "alignment_length": None,
                "query_coverage": None,
                "evalue": None,
                "bitscore": None,
                "genome_source": "mirbase_input_file",
                "disambiguation_status": "unique_sequence",
                "mapping_source": "carried_forward_unchanged",
            })

    result_df = pd.DataFrame(columns=["name", "chromosome", "start", "end", "strand"])
    raw_hits = pd.DataFrame()
    if to_map:
        try:
            result_df, raw_hits = map_sequences_batch(
                to_map, species_key, division=division,
                min_identity=args.min_identity, min_coverage=args.min_coverage,
                return_raw_hits=True,
                target_accession=target_accession, target_assembly=target_assembly,
            )
        except Exception as exc:
            print(f"[error] '{prefix}' ({scientific_name}) failed: {exc}", file=sys.stderr)
            _append_csv(
                pd.DataFrame([{
                    "species_prefix": prefix, "scientific_name": scientific_name, "division": division,
                    "n_requested": n_requested, "n_mapped": 0,
                    "status": f"error: {exc}",
                }]),
                SUMMARY_OUTPUT,
            )
            return

    n_ambiguous = 0
    if not result_df.empty:
        result_df = _resolve_duplicates_in_result(result_df, to_map, raw_hits, coords)
        n_ambiguous = int(result_df["disambiguation_status"].str.startswith("AMBIGUOUS").sum())
        if n_ambiguous:
            print(f"[warn] {prefix}: {n_ambiguous} row(s) flagged AMBIGUOUS_DUPLICATE - see {COORDS_OUTPUT}.", file=sys.stderr)
        result_df["mapping_source"] = "blast_remapped"

    combined_df = pd.concat(
        [pd.DataFrame(carried_forward_rows), result_df], ignore_index=True
    ) if carried_forward_rows else result_df

    if not combined_df.empty:
        combined_df.insert(0, "species_prefix", prefix)
        _append_csv(combined_df, COORDS_OUTPUT)

    _append_csv(
        pd.DataFrame([{
            "species_prefix": prefix, "scientific_name": scientific_name, "division": division,
            "n_requested": n_requested,
            "n_unchanged_carried_forward": len(carried_forward_rows),
            "n_mapped_via_blast": len(result_df),
            "n_ambiguous_duplicates": n_ambiguous,
            "status": "ok",
        }]),
        SUMMARY_OUTPUT,
    )
    print(f"[done] {prefix}: {len(carried_forward_rows)} carried forward + {len(result_df)} BLAST-mapped / {n_requested} requested\n")


if __name__ == "__main__":
    main()
