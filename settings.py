import dotenv
import os
from utils.logger import get_logger

logger = get_logger(__name__)


dotenv.load_dotenv()

# Debug settings
DEBUG = True if os.getenv("DEBUG") == "True" else False


# Database settings
DB_USER = os.getenv('DATABASE_USER')
DB_PASSWORD = os.getenv('DATABASE_PASS')
DB_HOST = os.getenv('DB_HOST')
DB_PORT = os.getenv('DB_PORT')
DB_NAME = os.getenv('DATABASE_NAME')

# Name of the cohort study
OBSERVATIONS = os.getenv("OBSERVATION_SOURCE")

# Data paths
PHENO_PATH = os.getenv("PHENOTYPE_META_PATH")
PROTEIN_PATH = os.getenv("PROTEIN_META_PATH")
METABOLITE_PATH = os.getenv("METABOLITE_META_PATH")
EDGES_PATH = os.getenv("CALCULATED_EDGES_PATH")
DATA_DIR = os.getenv("DATA_DIR")
GWAS_STATS_PATH = os.getenv("GWAS_STATS_PATH")
GENOMIC_VARIANT_META_PATH = os.getenv("GENOMIC_VARIANT_META_PATH")

# Other settings
CHUNK_SIZE = os.getenv("CHUNK_SIZE") if os.getenv("CHUNK_SIZE") else 10_000_000


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
logger.debug(f"EDGES_PATH: {EDGES_PATH}")
logger.debug(f"DATA_DIR: {DATA_DIR}")
logger.debug("---------- Miscellaneous ----------")
logger.debug(f"OBSERVATIONS: {OBSERVATIONS}")
logger.debug(f"CHUNK_SIZE: {CHUNK_SIZE:,}")
logger.debug("--------------------\n")
