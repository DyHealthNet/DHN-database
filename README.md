# DHN-database

This repository handles the creation of all node types as well as edges between the nodes.
Most data is dynamically pulled from the NeDRex API and subsequently inserted into the database.
However, we also connect to the HPO database as well as the HMDB database for additional information.


## Pre-requisites
Before you can run the script, you need to ensure the following:
- (tbd)
- you need to download the data from the HMDB database and place it in the data folder before 
  running the setup as this can not be done automatically
## Build &start Docker image

```bash
docker compose up --build
```
```bash
docker compose up -d
```
## End Docker image
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
