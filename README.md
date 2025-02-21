# DyHealthNet Database Repository

This repository handles the creation of all node types as well as edges between the nodes.
Most data is dynamically pulled from the NeDRex API and subsequently inserted into the database.
However, we also connect to the HPO database as well as the HMDB database for additional information.


## Pre-requisites
Before you can run the script, you need to ensure the following:
- you need to download the data from the HMDB database and place it in the data folder before 
  running the setup as this can not be done automatically
- make sure your .env file is set up correctly, follow the .env.example file

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

## Adding extra edges to the database
If you calculated an extra set of edges between your nodes that you would like to add,
you can specify the path(s) to the edge files in the `EXTRA_EDGES` variable in the `.env` file.
To specify multiple files, they must be comma-separated.

Furthermore, make sure they have the following column headers:
```
label1  label2  pval  effsize effsize_type  test  adj_pval
```
Each row should represent an edge between two nodes, where `label1` and `label2` are the cohort ids of the nodes.
Columns need to be separated by tabs.
