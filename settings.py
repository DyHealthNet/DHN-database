import dotenv
import os

dotenv.load_dotenv()

# Debug settings
DEBUG = os.getenv('DEBUG')


# Database settings
DB_USER = os.getenv('DB_USER')
DB_PASSWORD = os.getenv('DB_PASSWORD')
DB_HOST = os.getenv('DB_HOST')
DB_PORT = os.getenv('DB_PORT')
DB_NAME = os.getenv('DB_NAME')


OBSERVATIONS = os.getenv("OBSERVATION_SOURCE")

# Data paths
PHENO_PATH = os.getenv("PHENOTYPE_PATH")
PROTEIN_PATH = os.getenv("PROTEIN_PATH")
METABOLITE_PATH = os.getenv("METABOLITE_PATH")
EDGES_PATH = os.getenv("CALCULATED_EDGES_PATH")
DATA_DIR = os.getenv("DATA_DIR")

# Other settings
CHUNK_SIZE = os.getenv("CHUNK_SIZE") if os.getenv("CHUNK_SIZE") else 10_000_000

if DEBUG:
    print("Launching with the following settings:")
    print(f"DB_USER: {DB_USER}")
    print(f"DB_PASSWORD: {DB_PASSWORD}")
    print(f"DB_HOST: {DB_HOST}")
    print(f"DB_PORT: {DB_PORT}")
    print(f"DB_NAME: {DB_NAME}")
    print("--------------------")
    print(f"OBSERVATIONS: {OBSERVATIONS}")
    print(f"PHENO_PATH: {PHENO_PATH}")
    print(f"PROTEIN_PATH: {PROTEIN_PATH}")
    print(f"METABOLITE_PATH: {METABOLITE_PATH}")
    print(f"EDGES_PATH: {EDGES_PATH}")
    print(f"DATA_DIR: {DATA_DIR}")
    print("--------------------\n")
