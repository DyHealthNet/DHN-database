from sqlalchemy.orm import DeclarativeBase
from models import (Phenotype, Disorder, Metabolite, Protein, CohortPhenotype, CohortProtein, CohortMetabolite,
                    CohortReferencesMetabolite, CohortReferencesProtein, CohortReferencesPhenotype,
                    CohortReferencesDisease)
import pandas as pd


def cohort_phenotype_data(session, phenotype_path: str = None, obs_source: str = None) -> tuple[list, list, list]:
    print("Adding cohort phenotype data to the database.")
    # from the database, get the snomed ids and associated hpo ids/ mondo ids
    phenotypes = session.query(Phenotype).filter(Phenotype.observation_source == obs_source).all()
    disorders = session.query(Disorder).filter(Disorder.observation_source == obs_source).all()

    # map the snomed ids to the hpo ids
    snomed_map = {[y.split('.')[1] for y in x.xrefs if y.startswith('snomedct.')][0]: x.hpo_id for x in phenotypes}
    snomed_map.update(
        {[y.split('.')[1] for y in x.xrefs if y.startswith('snomedct.')][0]: x.mondo_id for x in disorders})

    print(f"Length of snomed map: {len(snomed_map)}")
    # read the phenotype data
    raw_phenotypes = pd.read_csv(phenotype_path, sep='\t')
    raw_phenotypes['snomed_id'] = raw_phenotypes['snomed_id'].fillna('')
    # go through the raw phenotype data and add the phenotypes to the database
    phenotypes_to_add = []
    disorder_references_to_add = []
    phenotype_references_to_add = []
    missing = set()
    for index, row in raw_phenotypes.iterrows():

        new_phenotype = CohortPhenotype(cohort_id=row['label'], display_name=row['snomed_term'],
                                        description=row['description'], xrefs=row['snomed_id'])
        phenotypes_to_add.append(new_phenotype)

        mondo_id = hpo_id = None
        for snomed_id in row['snomed_id'].split(";"):
            snomed_id = snomed_id.strip()
            if snomed_id in snomed_map:
                if snomed_map[snomed_id].startswith('hpo'):
                    hpo_id = snomed_map[snomed_id]
                else:
                    mondo_id = snomed_map[snomed_id]

            # add the references to the knowledge graph for the phenotypes
            if hpo_id:
                new_reference = CohortReferencesPhenotype(cohort_id=row['label'], hpo_id=hpo_id)
                disorder_references_to_add.append(new_reference)
            elif mondo_id:
                new_reference = CohortReferencesDisease(cohort_id=row['label'], mondo_id=mondo_id)
                phenotype_references_to_add.append(new_reference)
            else:
                missing.add(snomed_id)

    print(f"{len(missing)} snomed ids could not be mapped: {missing}")
    return phenotypes_to_add, disorder_references_to_add, phenotype_references_to_add


def cohort_metabolite_data(session, metabolite_path: str = None, obs_source: str = None) -> tuple[list, list]:
    print("Adding cohort metabolite data to the database.")
    metabolite_matches = session.query(Metabolite).filter(Metabolite.observation_source == obs_source).all()
    metabolite_map = {x.hmdb_id: x.hmdb_id for x in metabolite_matches}
    for x in metabolite_matches:
        for syn in x.synonyms:
            metabolite_map[syn] = x.hmdb_id

    raw_metabolites = pd.read_csv(metabolite_path, sep='\t')
    raw_metabolites['hmdb_id'] = raw_metabolites['hmdb_id'].fillna('')
    metabolites_to_add = []
    references_to_add = []
    missing = set()
    for index, row in raw_metabolites.iterrows():
        new_metabolite = CohortMetabolite(cohort_id=row['analyte_name'],
                                          display_name=row['biochemical_name'],
                                          description=row['analyte_class'],
                                          xrefs=row['hmdb_id'])

        metabolites_to_add.append(new_metabolite)

        # add the references to the knowledge graph for the metabolites
        for hmdb_id in row['hmdb_id'].split(';'):
            hmdb_id = f"hmdb.{hmdb_id}"
            if hmdb_id in metabolite_map:
                new_reference = CohortReferencesMetabolite(cohort_id=row['analyte_name'],
                                                           hmdb_id=metabolite_map[hmdb_id])
                references_to_add.append(new_reference)
            else:
                missing.add(hmdb_id)

    print(f"{len(missing)} hmdb ids could not be mapped: {missing}")
    return metabolites_to_add, references_to_add


def cohort_protein_data(session, protein_path: str = None, obs_source: str = None) -> tuple[list, list]:
    print("Adding cohort protein data to the database.")
    protein_matches = session.query(Protein).filter(Protein.observation_source == obs_source).all()
    protein_map = {x.uniprot_id: x.display_name for x in protein_matches}

    raw_proteins = pd.read_csv(protein_path, sep='\t')
    raw_proteins['UniProt'] = raw_proteins['UniProt'].fillna('')
    proteins_to_add = []
    references_to_add = []
    missing = set()
    for index, row in raw_proteins.iterrows():
        display_name = ", ".join([protein_map.get(f"uniprot.{x}", x) for x in row['UniProt'].split('|')])
        new_protein = CohortProtein(cohort_id=row['protein_id'],
                                    display_name=display_name,
                                    description=row['long_description'],
                                    xrefs=row['UniProt'])
        proteins_to_add.append(new_protein)

        # add the references to the knowledge graph for the proteins
        for uniprot_id in row['UniProt'].split('|'):
            uniprot_id = f"uniprot.{uniprot_id}"
            if uniprot_id in protein_map:
                new_reference = CohortReferencesProtein(cohort_id=row['protein_id'], uniprot_id=uniprot_id)
                references_to_add.append(new_reference)
            else:
                missing.add(uniprot_id)

    print(f"{len(missing)} uniprot ids could not be mapped: {missing}")
    return proteins_to_add, references_to_add
