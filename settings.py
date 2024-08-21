import dotenv
import os

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

if DEBUG:
    print("Launching with the following settings:")
    print("---------- Database settings ----------")
    print(f"DB_USER: {DB_USER}")
    print(f"DB_PASSWORD: {DB_PASSWORD}")
    print(f"DB_HOST: {DB_HOST}")
    print(f"DB_PORT: {DB_PORT}")
    print(f"DB_NAME: {DB_NAME}")
    print("---------- Input paths ----------")
    print(f"PHENO_PATH: {PHENO_PATH}")
    print(f"PROTEIN_PATH: {PROTEIN_PATH}")
    print(f"METABOLITE_PATH: {METABOLITE_PATH}")
    print(f"EDGES_PATH: {EDGES_PATH}")
    print(f"DATA_DIR: {DATA_DIR}")
    print("---------- Miscellaneous ----------")
    print(f"OBSERVATIONS: {OBSERVATIONS}")
    print(f"CHUNK_SIZE: {CHUNK_SIZE:,}")
    print("--------------------\n")
