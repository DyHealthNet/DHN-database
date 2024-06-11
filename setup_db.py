import os

import networkx as nx
from sqlalchemy import URL, create_engine
from sqlalchemy.orm import sessionmaker, DeclarativeBase

from metabolite_mapping import read_metabolite_mapping, read_hmdb_data, download_metabolite_data
from models import *
from query_nedrex import get_needed_snomed_ids, domain_id_to_mondo, get_disorder_data, get_edge_associations, \
    get_harmonizome_data, get_gene_data, get_phenotype_data
from hpo_mapping import download_hpo_ontology, read_hpo_ontology, ontology_data_to_network, snomed_from_hpo
from protein_mapping import get_proteinID_neddrex, read_proteinID_chris, add_proteinSet_data, \
    retrieve_interacting_proteins_neo4j
from nedrex.core import get_collection_attributes

# create postrgres db engine in memory
url = url_object = URL.create(
    "postgresql",
    username="postgres",
    password="password",  # plain (unescaped) text
    host="0.0.0.0",
    port=9000,
    database="postgres",
)
engine = create_engine(url)


def create_tables():
    # does not recreate tables if they already exist
    Base.metadata.create_all(engine)


def example_query(session):
    # Querying the database
    results = (session
               .query(Gene, Phenotype)
               .join(GeneAssocDisorder, Gene.entrez_id == GeneAssocDisorder.entrez_id)
               .filter(Gene.entrez_id == "12345")
               .all())
    for gene, phenotype in results:
        print(
            f"Gene entrez_id: {gene.entrez_id}, Phenotype hpo_id: {phenotype.hpo_id},"
            f" SNOMED ID: {phenotype.snomed_id}, OMIM ID: {phenotype.omim_id}")


def mondo_in_association_graph(mondo_id: str, assoc_graph: nx.Graph) -> tuple[list[str], list[str]] | tuple[None, None]:
    """
    Retrieves the genes associated with a mondo id from the association graph
    :param mondo_id: A mondo id, can be None
    :param assoc_graph: The association graph from NEDRex
    :return: genes associated with the mondo id, source of the data, harmonizome data if not in the association graph
    """
    if mondo_id is None:
        return None, None
    if mondo_id in assoc_graph:
        genes = []
        sources = []
        for edge in assoc_graph.edges(mondo_id, data=True):
            genes.append(edge[1])
            sources.append(edge[2]['source'][0])
    else:
        harm = get_harmonizome_data(mondo_id)
        if harm is None:
            return None, None
        genes = []
        sources = []
        for key, value in harm.items():
            genes.extend(value)
            sources.extend([key] * len(value))
    return genes, sources


def retrieve_disorder_data(needed_snomed: set[str], snomed_to_mondo: dict[str, str], descriptions: dict[str, str],
                           xrefs: dict, gene_info: dict, assoc_graph: nx.Graph) -> tuple[set, set, set, int]:
    """
    Queries the needed snomed ids and retrieves the associated genes and disorders from NEDRex
    :param gene_info: Information about the genes needed for the database (display name, synonyms, etc.)
    :param xrefs: cross references for the mondo ids to other databases
    :param descriptions: descriptions for the mondo ids
    :param needed_snomed: snomed ids in the dataset
    :param snomed_to_mondo: map from snomed to mondo ids from nedrex
    :param assoc_graph: association graph from nedrex of mondo ids to genes
    :return: list of genes to add, list of disorders to add, list of gene associations to add, number of snomed ids found
    """
    gene_associations = set()
    genes_to_add = set()
    disorders = set()
    found = 0
    for snomed in needed_snomed:
        # check if ; in snomed id and if so, do this for all ids
        snomed_ids = str(snomed).split(';')
        for snomed_id in snomed_ids:
            mondo_id = snomed_to_mondo.get(snomed_id)
            xref = xrefs.get(mondo_id, None)
            description = descriptions.get(mondo_id, None)
            genes, sources = mondo_in_association_graph(mondo_id, assoc_graph)
            if genes is None:
                continue
            # add genes to set
            for gene in genes:
                gene_data = gene_info.get(gene, None)
                if gene_data is None:
                    new_gene = Gene(entrez_id=gene)
                    genes_to_add.add(new_gene)
                    continue
                new_gene = Gene(entrez_id=gene_data['primaryDomainId'],
                                display_name=gene_data['displayName'],
                                description=gene_data['description'],
                                synonyms=gene_data['synonyms'],
                                chromosome=gene_data['chromosome'])
                genes_to_add.add(new_gene)

            disorders.add(Disorder(mondo_id=mondo_id, xrefs=xref, description=description))
            # add gene associations to set for each source
            gene_associations.update([GeneAssocDisorder(entrez_id=gene, mondo_id=mondo_id, edge_source=source)
                                      for gene, source in zip(genes, sources)])
            found += 1
        continue
    return genes_to_add, disorders, gene_associations, found


def retrieve_phenotype_data(available_ids: dict, additional_data: dict):
    """
    Retrieve the phenotype data for the needed snomed ids
    # HPO conversion: Pathway
    # HPO data (HPO_ID ----> SNOMED_ID) - look for needed SNOMED IDs
    # -> map to Mondo (SNOMED_ID -- OMIM_ID/ORPHA_ID --> Mondo_ID)
    #

    :param additional_data: dictionary with additional data for the hpo ids, must be a dictionary with hpo ids as keys
    :param available_ids: dictionary with snomed ids as keys and hpo ids as values
    :return: dictionary with the phenotype data
    """
    found = 0
    genes_to_add = set()
    phenotypes = set()
    disorder_associations = set()

    available_snomed_ids = available_ids
    # go through all the nodes in the HPO graph and find the ones that have xrefs to SNOMED

    print(f'Found {len(available_snomed_ids)} snomed ids in the HPO ontology')

    # get edge associations for disorder_has_phenotype
    assoc_graph = get_edge_associations(set(available_snomed_ids.values()), edge_type='disorder_has_phenotype')

    # find the genes that are associated with the mondo ids
    for snomed_id, hpo_id in available_snomed_ids.items():
        snomed_id = f"snomedct.{snomed_id}"
        if additional_data.get(hpo_id, None) is None:
            continue
        phenotype_data = additional_data[hpo_id]
        phenotypes.add(Phenotype(hpo_id=hpo_id,
                                 xrefs=set(phenotype_data['domainIds'] + [snomed_id]),
                                 description=phenotype_data['description'],
                                 synonyms=phenotype_data['synonyms'],
                                 display_name=phenotype_data['displayName']))
        if hpo_id not in assoc_graph:
            continue
        # get the disorder ids associated with the hpo id
        for edge in assoc_graph.edges(hpo_id, data=True):
            disorder = edge[1]
            source = edge[2]['source'][0]
            new_assoc = DisorderAssocPhenotype(mondo_id=disorder, hpo_id=hpo_id, edge_source=source)
            disorder_associations.add(new_assoc)

        found += 1
    return genes_to_add, phenotypes, disorder_associations, found


def add_items(session, items: iter, column: type[DeclarativeBase], filter_args: list):
    """
    Adds items to the database if they do not already exist
    :param session: Session object
    :param items: iterable of items to add
    :param column: the type of item to add, must be a class from models.py
    :param filter_args: the attributes to filter by to check if the item already exists
    :return: None
    """
    for item in items:
        filter_values = {key: getattr(item, key) for key in filter_args}
        exists = session.query(column).filter_by(**filter_values).first()
        if exists is not None:
            continue
        try:
            session.add(item)
        except Exception as e:
            print("Exception: ", e)
    session.commit()


def add_disorder_data(session, snomed_id_path: str):
    """
    Adds disorder data to the database given a path to a file with snomed ids
    :param session: Database session object
    :param snomed_id_path: str, path to file with snomed ids
    :return: None
    """
    needed_snomed = get_needed_snomed_ids(snomed_id_path)
    data = get_disorder_data(needed_snomed)
    snomed_to_mondo = domain_id_to_mondo(data)
    data = {x['primaryDomainId']: x for x in data}
    xrefs = {mondo: data[mondo]['domainIds'] for mondo in snomed_to_mondo.values() if 'domainIds' in data[mondo]}
    assoc_graph = get_edge_associations(set(snomed_to_mondo.values()), edge_type='gene_associated_with_disorder')
    # assoc_graph is filtered for ids that we need, now we can get the data for all genes in the graph since they're
    # all associated with the mondo ids
    gene_info = get_gene_data(set(assoc_graph.nodes))
    gene_dict = {x['primaryDomainId']: x for x in gene_info}

    mondo_description = {mondo: data[mondo]['description'] for mondo in snomed_to_mondo.values()}

    genes_to_add, disorders, gene_associations, found = retrieve_disorder_data(needed_snomed, snomed_to_mondo,
                                                                               mondo_description, xrefs, gene_dict,
                                                                               assoc_graph)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    add_items(session, disorders, Disorder, ['mondo_id'])
    add_items(session, gene_associations, GeneAssocDisorder, ['entrez_id', 'mondo_id'])
    session.commit()
    print(f"Found and successfully added {found} snomed ids with diseases to db")


def add_phenotype_data(session, phenotype_path: str, data_dir: str = '../data'):
    """
    Adds phenotype data to the database given a path to a file with phenotype data
    :param session: Database session object
    :param phenotype_path: str, path to file with phenotype data
    :param data_dir: str, path to the data directory
    :return: None
    """
    # data handling
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    needed_files = [f'{data_dir}/hp.json', f'{data_dir}/phenotype.hpoa']
    if not all([os.path.exists(f) for f in needed_files]):
        download_hpo_ontology(data_dir)
    needed_ids = get_needed_snomed_ids(phenotype_path)

    hpo_data = read_hpo_ontology(needed_files[0])
    hpo_graph = ontology_data_to_network(hpo_data)

    available_snomed_ids = snomed_from_hpo(hpo_graph, needed_ids)

    #         hpo_id = hpo_id.replace(':', '.').replace('HP', 'hpo')
    available_snomed_ids = {k: v.replace(':', '.').replace('HP', 'hpo') for k, v in available_snomed_ids.items()}
    pheno_data = get_phenotype_data(set(available_snomed_ids.values()))

    additional_data = {item['primaryDomainId']: item for item in pheno_data}
    genes_to_add, phenotypes, disorder_associations, _ = retrieve_phenotype_data(available_snomed_ids, additional_data)

    # since some phenotypes are subtypes of disorders, we only add phenotypes that are
    # not already in the disorder database
    removable_phenotypes = []
    removable_associations = []

    for phenotype in phenotypes:
        snomed = [x for x in phenotype.xrefs if 'snomedct' in x][0]
        items_disorder = session.query(Disorder).filter(Disorder.xrefs.any(snomed)).first()
        if not items_disorder:
            continue
        removable_phenotypes.append(phenotype)

    # also remove associations to phenotypes that are not in the disorder table
    for assoc in disorder_associations:
        if assoc.hpo_id in [x.hpo_id for x in removable_phenotypes]:
            removable_associations.append(assoc)
        if session.query(Disorder).filter_by(mondo_id=assoc.mondo_id).first() is None:
            removable_associations.append(assoc)

    print(f"Removing {len(removable_phenotypes)} phenotypes and {len(removable_associations)} associations that are "
          f"already in the disorder database")
    phenotypes = phenotypes - set(removable_phenotypes)
    disorder_associations = disorder_associations - set(removable_associations)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    add_items(session, phenotypes, Phenotype, ['hpo_id'])
    add_items(session, disorder_associations, DisorderAssocPhenotype, ['mondo_id', 'hpo_id'])

    session.commit()
    print(f"Found and successfully added {len(phenotypes)} snomed ids with phenotypes to db")


def add_protein_data(session, proteinData_path):
    proteinData = read_proteinID_chris(proteinData_path)
    neddrexProteins = []
    counter = 0
    for proteinID in proteinData:
        print(proteinID)
        if counter == 5:
            break
        entrez_id = get_proteinID_neddrex(proteinID)[0]['geneName']
        print(entrez_id, proteinID)
        newProtein = Protein(
            uniprot_id=proteinID)  #entrez_id = entrez_id, what to do if the id doesnt exist in gene? should i create a new gene?
        neddrexProteins.append(newProtein)
        #get_proteinID_neddrex(proteinID)[0]['geneName']
        counter += 1
    #print(*neddrexProteins)
    add_items(session, neddrexProteins, Protein, ['uniprot_id'])
    session.commit()


def add_metabolite_data(session, metabolite_path, data_dir: str = '../data'):
    hmdb_data_path = f'{data_dir}/hmdb_metabolites.xml'
    download_metabolite_data(data_dir)
    metabolite_mapping = read_metabolite_mapping(metabolite_path)
    unique_metabolites = set()

    # split metabolites that have ; in them
    for metabolite in metabolite_mapping['hmdb_id'].dropna():
        unique_metabolites.update(metabolite.split(';'))

    print(f"Found {len(unique_metabolites)} unique metabolites in the mapping file.")

    hmdb_mapping = read_hmdb_data(hmdb_data_path, unique_metabolites)
    print(f"Found info for {len(hmdb_mapping)} metabolites in the hmdb data file out of "
          f"{len(unique_metabolites)} metabolites in the mapping file.")

    omim_diseases = set()
    for metabolite in hmdb_mapping:
        omim_diseases.update(hmdb_mapping[metabolite]['diseases'])
    omim_diseases = {f"omim.{x}" for x in omim_diseases}
    disease_data = get_disorder_data(omim_diseases)
    omim_mapping = domain_id_to_mondo(disease_data, 'omim')

    metabolites = []
    metabolite_protein_associations = []
    metabolite_disease_associations = []

    for metabolite in hmdb_mapping:
        metabolite_name = f"hmdb.{metabolite}"
        # check if metabolite is in the db already and if so, skip
        exists = session.query(Metabolite).filter_by(hmdb_id=metabolite_name).first()
        # if exists is not None:
        #    continue
        metabolites.append(Metabolite(hmdb_id=metabolite_name, display_name=hmdb_mapping[metabolite]['display_name'],
                                      description=hmdb_mapping[metabolite]['description'],
                                      synonyms=hmdb_mapping[metabolite]['synonyms'],
                                      xrefs=hmdb_mapping[metabolite]['xrefs'], ))

        for protein in hmdb_mapping[metabolite]['proteins']:
            # check if the protein is in the db, if it is not, skip
            if session.query(Protein).filter_by(uniprot_id=protein).first() is None:
                continue
            metabolite_protein_associations.append(ProteinAssocMetabolite(hmdb_id=metabolite_name, uniprot_id=protein))

        for disease in hmdb_mapping[metabolite]['diseases']:
            if session.query(Disorder).filter_by(mondo_id=omim_mapping.get(f"omim.{disease}", None)).first() is None:
                continue
            metabolite_disease_associations.append(MetaboliteAssocDisorder(hmdb_id=metabolite_name,
                                                                           mondo_id=omim_mapping[f"omim.{disease}"]))

    print(f"A total of {len(metabolites)} metabolites were found in the mapping file, "
          f"as well as {len(metabolite_protein_associations)} protein associations and "
          f"{len(metabolite_disease_associations)} disease associations.")
    add_items(session, metabolites, Metabolite, ['hmdb_id'])
    add_items(session, metabolite_protein_associations, ProteinAssocMetabolite, ['hmdb_id', 'uniprot_id'])
    add_items(session, metabolite_disease_associations, MetaboliteAssocDisorder, ['hmdb_id', 'mondo_id'])
    session.commit()


if __name__ == '__main__':
    # Define a session
    Session = sessionmaker(bind=engine)
    session = Session()
    create_tables()
    pheno_data_path = '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv'
    protein_data_path = '../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt'
    metabo_data_path = '../data/DyHealthNet/chris_summary_data/metabolites/CHRIS_biocristes7500SumStats.txt'

    # add_disorder_data(session, pheno_data_path)
    # add_phenotype_data(session, pheno_data_path)
    # add_protein_data(session, protein_data_path)

    add_metabolite_data(session, metabo_data_path)

