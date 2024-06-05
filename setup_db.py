import os
import networkx as nx
from sqlalchemy import URL, create_engine
from sqlalchemy.orm import sessionmaker, DeclarativeBase
from models import *
from query_nedrex import get_needed_snomed_ids, domain_id_to_mondo, get_disorder_data, get_edge_associations, \
    get_harmonizome_data, get_gene_data, get_phenotype_data
from hpo_mapping import download_hpo_ontology, read_hpo_ontology, ontology_data_to_network, snomed_from_hpo
from protein_mapping import get_proteinID_neddrex, read_proteinID_chris, add_proteinSet_data

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
        if counter == 1:
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


# add_items(session, genes_to_add, Gene, ['entrez_id'])


if __name__ == '__main__':
    # Define a session
    from nedrex.core import iter_nodes, iter_edges

    help(iter_nodes)

    Session = sessionmaker(bind=engine)
    session = Session()
    create_tables()
    pheno_data_path = '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv'
    protein_data_path = '../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt'
    # add_disorder_data(session, pheno_data_path)
    #add_phenotype_data(session, pheno_data_path)
    #add_protein_data(session, protein_data_path)
    add_proteinSet_data(session=session,protein_data_path=protein_data_path)


    # 'primaryDomainId': 'uniprot.P43320', 'domainIds': ['uniprot.P43320'], 'sequence': 'MASDHQTQAGKPQSLNPKIIIFEQENFQGHSHELNGPCPNLKETGVEKAGSVLVQAGPWVGYEQANCKGEQFVFEKGEYPRWDSWTSSRRTDSLSSLRPIKVDSQEHKIILYENPNFTGKKMEIIDDDVPSFHAHGYQEKVSSVRVQSGTWVGYQYPGYRGLQYLLEKGDYKDSSDFGAPHPQVQSVRRIRDMQWHQRGAFHPSN', 'displayName': 'CRBB2_HUMAN', 'synonyms': ['Beta-crystallin B2', 'Beta-B2 crystallin', 'Beta-crystallin Bp'], 'comments': 'FUNCTION: Crystallins are the dominant structural components of the vertebrate eye lens.\nSUBUNIT: Homo/heterodimer, or complexes of higher-order. The structure of beta-crystallin oligomers seems to be stabilized through interactions between the N-terminal arms (By similarity). {ECO:0000250}.\nINTERACTION: Self; NbExp=5; IntAct=EBI-974082, EBI-974082;\nDOMAIN: Has a two-domain beta-structure, folded into four very similar Greek key motifs.\nMASS SPECTROMETRY: Mass=23291; Mass_error=3; Method=Electrospray; Evidence={ECO:0000269|PubMed:8999933};\nMASS SPECTROMETRY: Mass=23289; Method=Electrospray; Evidence={ECO:0000269|PubMed:8175657};\nMASS SPECTROMETRY: Mass=23290; Method=Electrospray; Evidence={ECO:0000269|PubMed:10930324};\nDISEASE: Cataract 3, multiple types (CTRCT3) [MIM:601547]: An opacification of the crystalline lens of the eye that frequently results in visual impairment or blindness. Opacities vary in morphology, are often confined to a portion of the lens, and may be static or progressive. CTRCT3 includes congenital cerulean and sutural cataract with punctate and cerulean opacities, among others. Cerulean cataract is characterized by peripheral bluish and white opacifications organized in concentric layers with occasional central lesions arranged radially. The opacities are observed in the superficial layers of the fetal nucleus as well as the adult nucleus of the lens. Involvement is usually bilateral. Visual acuity is only mildly reduced in childhood. In adulthood, the opacifications may progress, making lens extraction necessary. Histologically the lesions are described as fusiform cavities between lens fibers which contain a deeply staining granular material. Although the lesions may take on various colors, a dull blue is the most common appearance and is responsible for the designation cerulean cataract. Sutural cataract with punctate and cerulean opacities is characterized by white opacification around the anterior and posterior Y sutures, and grayish and bluish, spindle shaped, oval punctate and cerulean opacities of various sizes arranged in lamellar form. The spots are more concentrated towards the peripheral layers and do not delineate the embryonal or fetal nucleus. Phenotypic variation with respect to the size and density of the sutural opacities as well as the number and position of punctate and cerulean spots is observed among affected subjects. {ECO:0000269|PubMed:10634616, ECO:0000269|PubMed:9158139}. Note=The disease is caused by mutations affecting the gene represented in this entry.\nSIMILARITY: Belongs to the beta/gamma-crystallin family. {ECO:0000305}.\nWEB RESOURCE: Name=Eye disease Crystallin, beta-B2 (CRYBB2); Note=Leiden Open Variation Database (LOVD); URL="http://www.lovd.nl/CRYBB2";', 'geneName': 'CRYBB2', 'taxid': 9606, 'type': 'Protein'}]
