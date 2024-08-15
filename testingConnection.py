from sqlalchemy import create_engine

# Define the connection string
DATABASE_URL = "postgresql://postgres:password@localhost:9852/postgres"

# Create an engine instance
engine = create_engine(DATABASE_URL)

# Try connecting
try:
    with engine.connect() as connection:
        print("Connection successful!")
except Exception as e:
    print(f"Error connecting: {e}")
