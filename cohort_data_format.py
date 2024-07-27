from models import Phenotype, Disorder, Metabolite, Protein, CohortPhenotype, CohortProtein, CohortMetabolite
import pandas as pd


def cohort_phenotype_data(session, phenotype_path: str = None, obs_source: str = None) -> list:
    print("Adding cohort phenotype data to the database.")
    # read the phenotype data
    raw_phenotypes = pd.read_csv(phenotype_path, sep='\t')
    # from the database, get the snomed ids and associated hpo ids/ mondo ids
    phenotypes = session.query(Phenotype).filter(Phenotype.observation_source == obs_source).all()
    disorders = session.query(Disorder).filter(Disorder.observation_source == obs_source).all()

    # map the snomed ids to the hpo ids
    snomed_map = {[y.split('.')[1] for y in x.xrefs if y.startswith('snomedct.')][0]: x.hpo_id for x in phenotypes}
    snomed_map.update(
        {[y.split('.')[1] for y in x.xrefs if y.startswith('snomedct.')][0]: x.mondo_id for x in disorders})

    print(f"Length of snomed map: {len(snomed_map)}")
    # go through the raw phenotype data and add the phenotypes to the database
    phenotypes_to_add = []
    for index, row in raw_phenotypes.iterrows():
        snomed_id = row['snomed_id']
        mondo_id = hpo_id = None
        if snomed_id in snomed_map:
            if snomed_map[snomed_id].startswith('hpo'):
                hpo_id = snomed_map[snomed_id]
            else:
                mondo_id = snomed_map[snomed_id]

        new_phenotype = CohortPhenotype(cohort_id=row['label'], display_name=row['snomed_term'],
                                        description=row['description'], mondo_id=mondo_id, hpo_id=hpo_id)
        phenotypes_to_add.append(new_phenotype)
    return phenotypes_to_add


def cohort_metabolite_data(session, metabolite_path: str = None, obs_source: str = None) -> list:
    print("Adding cohort metabolite data to the database.")
    metabolite_matches = session.query(Metabolite).filter(Metabolite.observation_source == obs_source).all()
    # This method will assign each metabolite in the cohort study one HMDB id (+ description) even though in the data
    # itself, there may be multiple HMDB ids for the same metabolite.
    metabolite_map = {x.display_name: (x.hmdb_id, x.description) for x in metabolite_matches}
    for x in metabolite_matches:
        for syn in x.synonyms:
            metabolite_map[syn] = x.hmdb_id

    # read the metabolite data
    raw_metabolites = pd.read_csv(metabolite_path, sep='\t')
    metabolites_to_add = []
    for index, row in raw_metabolites.iterrows():
        name = row['analyte_name']
        new_metabolite = CohortMetabolite(cohort_id=name,
                                          display_name=name,
                                          description=metabolite_map.get(name, [None, None])[1],
                                          hmdb_id=metabolite_map.get(name, [None, None])[0])
        metabolites_to_add.append(new_metabolite)
    return metabolites_to_add


def cohort_protein_data(session, protein_path: str = None, obs_source: str = None) -> list:
    print("Adding cohort protein data to the database.")
    protein_matches = session.query(Protein).filter(Protein.observation_source == obs_source).all()
    protein_map = {x.uniprot_id: x.display_name for x in protein_matches}

    raw_proteins = pd.read_csv(protein_path, sep='\t')
    proteins_to_add = []
    for index, row in raw_proteins.iterrows():
        name = row['protein_id']
        uniprot = f"uniprot.{row['UniProt']}" if row['UniProt'] else None
        uniprot = uniprot if uniprot in protein_map else None
        new_protein = CohortProtein(cohort_id=name,
                                    display_name=protein_map.get(name, None),
                                    description=row['long_description'],
                                    uniprot_id=uniprot)
        proteins_to_add.append(new_protein)
    return proteins_to_add
