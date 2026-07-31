import dotenv
import os

from sqlalchemy import URL

from utils.logger import get_logger

# final id prefix in db and on NeDRex
valid_db_to_id_prefix = {
    'SNOMED CT':'snomedct',
    'MONDO':'mondo',
    'MedDRA':'medra',
    'UMLS':'umls',
    'DOID':'doid',
    'NCIt':'ncit',
    'ORPHA net':'orpha'
}

valid_db_to_hpo_prefix = {
    'SNOMED CT':['SNOMEDCT_US','SNOMED_CT'],
    'MONDO':'MONDO',
    'UMLS':'UMLS',
    'DOID':'DOID',
    'NCIt':'NCIT',
    'ORPHA':'ORPHA',
    'MedDRA':'MEDDRA',
    # Not in NeDRex therefor currently never used
    'Fyler Code':'Fyler',
    'PubMed':'PMID',
    'EFO':'EFO',
    'MPATH':'MPATH',
    'MP':'MP',
    'EPCC':'EPCC',
    'ICD-9': 'ICD-9',
    'ICD9':'ICD9',
    'ICD-10': 'ICD-10',
    'ICD10': 'ICD10',
    'ICD-O': 'ICD-O',
    'COHD': 'COHD',
}

logger = get_logger(__name__)


dotenv.load_dotenv()

# Debug settings
DEBUG = os.getenv("DEBUG", "false").lower() == "true"


# Database settings
DB_USER = os.getenv('DATABASE_USER')
DB_PASSWORD = os.getenv('DATABASE_PASS')
DB_HOST = os.getenv('DB_HOST')
DB_PORT = os.getenv('DB_PORT')
DB_NAME = os.getenv('DATABASE_NAME')

url = url_object = URL.create(
    "postgresql",
    username=DB_USER,
    password=DB_PASSWORD,
    host=DB_HOST,
    port=DB_PORT,
    database=DB_NAME,
)


# Name of the cohort study
OBSERVATIONS = os.getenv("OBSERVATION_SOURCE")

# Data paths
PHENO_PATH = os.getenv("PHENOTYPE_META_PATH")
PROTEIN_PATH = os.getenv("PROTEIN_META_PATH")
METABOLITE_PATH = os.getenv("METABOLITE_META_PATH")
EDGES_PATH = os.getenv("CALCULATED_EDGES_PATH")
DATA_DIR = os.getenv("DATA_DIR")
PARAMETRIC_EDGES_PATH = os.getenv("PARAMETRIC_EDGES_PATH") or None
NONPARAMETRIC_EDGES_PATH = os.getenv("NONPARAMETRIC_EDGES_PATH") or None

# Node sources, mirroring the backend's network.utils.data_manager.combine_data(): one
# entry per data source (e.g. phenotypes, proteins, metabolites), contributing label/type
# (required) to the nodes table. Comma-separated, same convention as combine_data(): all
# three lists must have the same number of entries, and DATA_ROOT (if set) is prepended
# to relative entries in DATA_META_PATHS.
DATA_ROOT = os.getenv("DATA_ROOT")
DATA_META_PATHS = os.getenv("DATA_META_PATHS")
DATA_LABEL_COLUMNS = os.getenv("DATA_LABEL_COLUMNS")
DATA_TYPE_COLUMNS = os.getenv("DATA_TYPE_COLUMNS")

# Optional per-source enrichment columns (display_name/description/xref/group) for the
# nodes table - the backend doesn't need these for scoring, but the DB wants them where
# available. Same length as DATA_META_PATHS if set; use an empty entry to skip a source
# that doesn't have that column. Left entirely unset, nodes just get no enrichment.
DATA_DP_NAME_COLUMNS = os.getenv("DATA_DP_NAME_COLUMNS")
DATA_DESCRIPTION_COLUMNS = os.getenv("DATA_DESCRIPTION_COLUMNS")
DATA_XREF_COLUMNS = os.getenv("DATA_XREF_COLUMNS")
DATA_GROUP_COLUMNS = os.getenv("DATA_GROUP_COLUMNS")
DATA_SUBGROUP_COLUMNS = os.getenv("DATA_SUBGROUP_COLUMNS")
EXTRA_EDGES = os.getenv("EXTRA_EDGES")
VARIANT_META_PATH = os.getenv("VARIANT_META_PATH")

INPUT_ID_DB = os.getenv("INPUT_DB_ID")
ID_PREFIX = valid_db_to_id_prefix.get(INPUT_ID_DB, None)
# If User input nedrex ID prefix instead of DB/ Ontology name
# if INPUT_ID_DB in valid_db_to_id_prefix.values():
#     for key, value in valid_db_to_id_prefix.items():
#         if value == ID_PREFIX:
#             INPUT_ID_DB = key
#             ID_PREFIX = value
HPO_ID_PREFIX = valid_db_to_hpo_prefix.get(INPUT_ID_DB, None)

# Cohort Columns
COHORT_COLUMNS = {
    "phenotype": {'unique_id': os.getenv("PHENOTYPE_LABEL_COLUMN"),
                  'display_name': os.getenv("PHENOTYPE_DP_NAME_COLUMN"),
                  'description': os.getenv("PHENOTYPE_DESCRIPTION_COLUMN"),
                  'xref': os.getenv("PHENOTYPE_XREF_COLUMN")},

    "protein": {'unique_id': os.getenv("PROTEIN_LABEL_COLUMN"),
                'display_name': os.getenv("PROTEIN_DP_NAME_COLUMN"),
                'description': os.getenv("PROTEIN_DESCRIPTION_COLUMN"),
                'xref': os.getenv("PROTEIN_XREF_COLUMN")},

    "metabolite": {'unique_id': os.getenv("METABOLITE_LABEL_COLUMN"),
                   'display_name': os.getenv("METABOLITE_DP_NAME_COLUMN"),
                   'description': os.getenv("METABOLITE_DESCRIPTION_COLUMN"),
                   'xref': os.getenv("METABOLITE_XREF_COLUMN")},

    "variant": {'unique_id': os.getenv("VARIANT_LABEL_COLUMN"),
                'display_name': os.getenv("VARIANT_DP_NAME_COLUMN"),
                'description': os.getenv("VARIANT_DESCRIPTION_COLUMN"),
                'xref': os.getenv("VARIANT_XREF_COLUMN")}
}

# Other settings
CHUNK_SIZE = int(os.getenv("CHUNK_SIZE")) if os.getenv("CHUNK_SIZE") else 10_000_000
VISUALIZE = os.getenv("VISUALIZE", "false").lower() == "true"


logger.debug("---------- Database settings ----------")
logger.debug(f"DB_USER: {DB_USER}")
logger.debug(f"DB_PASSWORD: {DB_PASSWORD}")
logger.debug(f"DB_HOST: {DB_HOST}")
logger.debug(f"DB_PORT: {DB_PORT}")
logger.debug(f"DB_NAME: {DB_NAME}")
logger.debug("---------- Input paths ----------")
logger.debug(f"PHENO_PATH: {PHENO_PATH}")
logger.debug(f"PROTEIN_PATH: {PROTEIN_PATH}")
logger.debug(f"METABOLITE_PATH: {METABOLITE_PATH}")
logger.debug(f"VARIANT_META_PATH: {VARIANT_META_PATH}")
logger.debug(f"EDGES_PATH: {EDGES_PATH}")
logger.debug(f"EXTRA_EDGES: {EXTRA_EDGES}")
logger.debug(f"DATA_DIR: {DATA_DIR}")
logger.debug("---------- Cohort Columns ----------")
logger.debug(f"Protein columns: {COHORT_COLUMNS['protein']}")
logger.debug(f"Phenotype columns: {COHORT_COLUMNS['phenotype']}")
logger.debug(f"Metabolite columns: {COHORT_COLUMNS['metabolite']}")
logger.debug(f"Variant columns: {COHORT_COLUMNS['variant']}")
logger.debug(f"Phenotype ID Type: Database: {INPUT_ID_DB}, Prefix (NeDRex API): "
             f"{ID_PREFIX}, Prefix (HPO API): {HPO_ID_PREFIX}")
logger.debug("---------- Miscellaneous ----------")
logger.debug(f"OBSERVATIONS: {OBSERVATIONS}")
logger.debug(f"CHUNK_SIZE: {CHUNK_SIZE:,}")
logger.debug("--------------------\n")
