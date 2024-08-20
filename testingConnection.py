from sqlalchemy import create_engine
import nedrex
from nedrex.core import api_keys_active, get_api_key

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


from nedrex.core import iter_nodes, iter_edges, get_node_types, get_collection_attributes, get_edge_types


#nedrex.config.set_url_base("https://api.nedrex.net/open/")
nedrex.config.set_url_base("https://apps.cosy.bio/licensed")
if api_keys_active():
    api_key = get_api_key(accept_eula=True)
    nedrex.config.set_api_key(api_key)
