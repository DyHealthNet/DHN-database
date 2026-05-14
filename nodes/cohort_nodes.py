import os
import timeit

from sqlalchemy.orm import Session, declarative_base

from utils.models import (Phenotype, Disorder, Metabolite, Protein, CohortPhenotype, CohortProtein, CohortMetabolite,
                          CohortVariant, CohortReferencesVariant, GenomicVariant, CohortReferencesMetabolite,
                          CohortReferencesProtein, CohortReferencesPhenotype, CohortReferencesDisease)
from utils.settings import COHORT_COLUMNS, ID_PREFIX, INPUT_ID_DB
from utils.logger import get_logger
import pandas as pd

logger = get_logger(__name__)


def read_data(data_path):
    file_name, ending = os.path.splitext(data_path)
    if ending not in ['.csv', '.tsv']:
        raise ValueError(f"Unsupported file format for phenotypes meta file: {ending}. "
                         f"Only CSV and TSV files are supported.")
    sep = "," if ending == ".csv" else "\t"
    return pd.read_csv(data_path, sep=sep)


def get_cols(node_type: str) -> tuple[str, str, str, str]:
    cols = COHORT_COLUMNS.get(node_type)
    if not cols:
        raise ValueError(f"Node type {node_type} not found in the cohort columns")
    return cols['unique_id'], cols['display_name'], cols['description'], cols['xref']


def cohort_phenotype_data(session: Session, phenotype_path: str = None, obs_source: str = None) \
        -> tuple[list, list, list]:
    logger.debug("Adding cohort phenotype data to the database.")
    u_id, dp_name, desc, xrefs = get_cols('phenotype')

    # from the database, get the phenotype reference ids and associated hpo ids/ mondo ids
    phenotypes = session.query(Phenotype).filter(Phenotype.observation_source == obs_source).all()
    disorders = session.query(Disorder).filter(Disorder.observation_source == obs_source).all()

    # map the phenotype reference ids to the hpo ids
    pheno_id_map = {}
    for x in phenotypes:
        if [y.split('.')[1] for y in x.xrefs if y.startswith(ID_PREFIX+'.')][0] in pheno_id_map:
            pheno_id_map[[y.split('.')[1] for y in x.xrefs if y.startswith(ID_PREFIX+'.')][0]].append(x.hpo_id)
        else:
            pheno_id_map[[y.split('.')[1] for y in x.xrefs if y.startswith(ID_PREFIX + '.')][0]] = [x.hpo_id]
    #pheno_id_map = {[y.split('.')[1] for y in x.xrefs if y.startswith(ID_PREFIX+'.')][0]: x.hpo_id for x in phenotypes}
    pheno_id_map.update(
        {[y.split('.')[1] for y in x.xrefs if y.startswith(ID_PREFIX+'.')][0]: [x.mondo_id] for x in disorders})

    logger.debug(f"Length of pheno id map: {len(pheno_id_map)}")
    # read the phenotype data
    raw_phenotypes = read_data(phenotype_path)
    has_xrefs = bool(xrefs) and xrefs in raw_phenotypes.columns
    if has_xrefs:
        raw_phenotypes[xrefs] = raw_phenotypes[xrefs].fillna('')
    # go through the raw phenotype data and add the phenotypes to the database
    phenotypes_to_add = []
    disorder_references_to_add = []
    phenotype_references_to_add = []
    missing = set()
    for index, row in raw_phenotypes.iterrows():
        display_name = row[dp_name] if dp_name and row[dp_name] and isinstance(row[dp_name], str) else row[COHORT_COLUMNS['phenotype']['unique_id']]   # TODO change 'label' to settings value
        xref_value = row[xrefs] if has_xrefs else ''
        new_phenotype = CohortPhenotype(cohort_id=row[u_id], display_name=display_name,
                                        description=row[desc],
                                        xrefs="|".join([f"{ID_PREFIX}.{x}" for x in str(xref_value).split(';') if x]))
        phenotypes_to_add.append(new_phenotype)

        if not has_xrefs:
            continue

        mondo_id = hpo_id = None
        for pheno_id in str(row[xrefs]).split(";"):
            pheno_id = pheno_id.strip()
            if not pheno_id:
                continue
            if pheno_id in pheno_id_map:
                for id in pheno_id_map[pheno_id]:
                    # add the references to the knowledge graph for the phenotypes
                    if id.startswith('hpo'):
                        new_reference = CohortReferencesPhenotype(cohort_id=row[u_id], hpo_id=id)
                        disorder_references_to_add.append(new_reference)
                    else:
                        new_reference = CohortReferencesDisease(cohort_id=row[u_id], mondo_id=id)
                        phenotype_references_to_add.append(new_reference)
            else:
                missing.add(pheno_id)

    logger.info(f"Some {INPUT_ID_DB} ids could not be mapped {list(missing)[:min(len(missing) - 1, 5)]} "
                f"and {max(len(missing) - 5, 0)} more")
    return phenotypes_to_add, disorder_references_to_add, phenotype_references_to_add


def cohort_metabolite_data(session: Session, metabolite_path: str = None, obs_source: str = None) -> tuple[list, list]:
    logger.debug("Adding cohort metabolite data to the database.")
    u_id, dp_name, desc, xrefs = get_cols('metabolite')

    metabolite_matches = session.query(Metabolite).filter(Metabolite.observation_source == obs_source).all()
    metabolite_map = {x.hmdb_id: x.hmdb_id for x in metabolite_matches}
    for x in metabolite_matches:
        for syn in x.synonyms:
            metabolite_map[syn] = x.hmdb_id

    raw_metabolites = read_data(metabolite_path)
    has_xrefs = bool(xrefs) and xrefs in raw_metabolites.columns
    if has_xrefs:
        raw_metabolites[xrefs] = raw_metabolites[xrefs].fillna('')
    metabolites_to_add = []
    references_to_add = []
    missing = set()
    for index, row in raw_metabolites.iterrows():
        xref_value = row[xrefs] if has_xrefs else ''
        new_metabolite = CohortMetabolite(cohort_id=row[u_id],
                                          display_name=row[dp_name],
                                          description=row[desc],
                                          xrefs="|".join([f"hmdb.{x}" for x in str(xref_value).split(';') if x]))

        metabolites_to_add.append(new_metabolite)

        if not has_xrefs:
            continue

        # add the references to the knowledge graph for the metabolites
        for hmdb_id in row[xrefs].split(';'):
            hmdb_id = f"hmdb.{hmdb_id}"
            if hmdb_id in metabolite_map:
                new_reference = CohortReferencesMetabolite(cohort_id=row[u_id],
                                                           hmdb_id=metabolite_map[hmdb_id])
                references_to_add.append(new_reference)
            else:
                missing.add(hmdb_id)

    logger.info(f"Some HMDB IDs could not be mapped: {list(missing)[:min(len(missing) - 1, 5)]} "
                f"and {max(len(missing) - 5, 0)} more")
    return metabolites_to_add, references_to_add


def cohort_protein_data(session: Session, protein_path: str = None, obs_source: str = None) -> tuple[list, list]:
    logger.debug("Adding cohort protein data to the database.")
    u_id, dp_name, desc, xrefs = get_cols('protein')

    protein_matches = session.query(Protein).filter(Protein.observation_source == obs_source).all()
    protein_map = {x.uniprot_id: x.display_name for x in protein_matches}

    raw_proteins = read_data(protein_path)
    has_xrefs = bool(xrefs) and xrefs in raw_proteins.columns
    if has_xrefs:
        raw_proteins[xrefs] = raw_proteins[xrefs].fillna('')
    proteins_to_add = []
    references_to_add = []
    missing = set()
    for index, row in raw_proteins.iterrows():
        display_name = ", ".join([protein_map.get(f"uniprot.{x}", x) for x in row[dp_name].split('|')]) if dp_name else row[u_id]
        xref_value = row[xrefs] if has_xrefs else ''
        new_protein = CohortProtein(cohort_id=row[u_id],
                                    display_name=display_name,
                                    description=row[desc],
                                    xrefs="|".join([f"uniprot.{x}" for x in str(xref_value).split('|') if x]))
        proteins_to_add.append(new_protein)

        if not has_xrefs:
            continue

        # add the references to the knowledge graph for the proteins
        for uniprot_id in row[xrefs].split('|'):
            uniprot_id = f"uniprot.{uniprot_id}"
            if uniprot_id in protein_map:
                new_reference = CohortReferencesProtein(cohort_id=row[u_id], uniprot_id=uniprot_id)
                references_to_add.append(new_reference)
            else:
                missing.add(uniprot_id)

    logger.info(f"Some UniProt IDs could not be mapped: {list(missing)[:min(len(missing) - 1, 5)]} "
                f"and {max(len(missing) - 5, 0)} more")
    return proteins_to_add, references_to_add


def cohort_variant_data(session: Session, variants_meta_path: str, obs_source: str) -> tuple[set, set]:
    """
    Retrieves the variant data from the file
    """
    variants_meta_df = read_data(variants_meta_path, dtype=str)
    unique_id, dp_name, desc, xrefs = get_cols('variant')
    start = timeit.default_timer()

    def create_variant(row):
        return CohortVariant(
            cohort_id=row[unique_id],
            description=row[desc],
            display_name=row[dp_name],
            xrefs=f"rsid.{row[xrefs]}"
        )

    variant_set = set(variants_meta_df.apply(create_variant, axis=1))
    logger.debug(f"Time taken to process variants: {timeit.default_timer() - start}")

    variant_refs = get_cohort_references_variant(session, variant_set, obs_source)
    return variant_set, variant_refs


def get_cohort_references_variant(session: Session, variants: declarative_base, obs_source: str) -> set:
    new_cohort_references_set = set()

    existing_cohort_id = {(genomic_variant.description, f"{genomic_variant.cohort_id[-1]}")
                          for genomic_variant in variants}

    desc_map = {f"{genomic_variant.description}{genomic_variant.cohort_id[-1]}": str(genomic_variant.cohort_id)
                for genomic_variant in variants}

    for variant in session.query(GenomicVariant).all():
        variant_domain_ids = variant.xrefs
        dbsnp_id = next((variant_id.replace("dbsnp.", "rs") for variant_id in variant_domain_ids
                         if "dbsnp." in variant_id), None)
        clinvar_id = variant.clinvar_id
        alt_seq = variant.alternative_sequence

        if (dbsnp_id, alt_seq) in existing_cohort_id:
            cohort_id = desc_map[f"{dbsnp_id}{alt_seq}"]
            new_cohort_references_variant = CohortReferencesVariant(
                cohort_id=cohort_id,
                clinvar_id=clinvar_id
            )
            new_cohort_references_set.add(new_cohort_references_variant)

    return new_cohort_references_set
