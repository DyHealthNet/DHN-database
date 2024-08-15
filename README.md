# DHN-database

This repository handles the creation of all node types as well as edges between the nodes.
Most data is dynamically pulled from the NeDRex API and subsequently inserted into the database.
However, we also connect to the HPO database as well as the HMDB database for additional information.


## Pre-requisites
Before you can run the script, you need to ensure the following:
- you need to download the data from the HMDB database and place it in the data folder before 
  running the setup as this can not be done automatically
- make sure your .env file is set up correctly, follow the .env.example file
- the files need some columns specifically named  for the script to run:
```
# For each file we have:
usage: column_name

# Proteins
unique_id: protein_id
display_name: UniProt # currently we take it from the external db part
description: long_description
cross-ref: UniProt

# Metabolites
unique_id: analyte_name
display_name: biochemical_name
description: analyte_class
cross-ref: hmdb_id

# Phenotypes
unique_id: label
display_name: snomed_term
description: description
cross-ref: snomed_id
```

## Build the database

```bash
docker compose up --build
```
```bash
docker compose up -d
```
## End Docker image & delete database volume (only if you want to start from scratch)
```bash
docker compose down -v
```
## Run

First, create a conda environment using the provided environment.yml file:

```bash
conda env create -f environment.yml
```

Then, activate the environment:

```bash
conda activate DHN-database
```

Finally, run the main script:

```bash
python setup_db.py
```

Alternatively, you can run the script in the background and log the output to a file:

```bash
nohup python -u setup_db.py > setup_db.log 2>&1 &
```
