import os
from sqlalchemy import URL
from sqlalchemy.orm import sessionmaker
from models import *
from query_nedrex import get_needed_snomed_ids, domain_id_to_mondo, get_disorder_data, get_edge_associations, \
    get_harmonizome_data
from hpo_mapping import download_hpo_ontology, read_hpo_ontology, ontology_data_to_network, snomed_from_hpo
from protein_mapping import get_proteinID_neddrex, read_proteinID_chris

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


def example_add(session):
    # example test
    new_gene = Gene(entrez_id="12345")
    new_disorder = Disorder(mondo_id='MONDO:0000001', snomed_id='SNOMED:0000001')
    new_phenotype = Phenotype(hpo_id='HP:0000001', snomed_id='SNOMED:0000002', omim_id='OMIM:0000001')

    session.add(new_gene)
    session.add(new_disorder)
    session.add(new_phenotype)
    session.commit()

    # Create associations
    gene_assoc_disorder = GeneAssocDisorder(entrez_id="12345", mondo_id='MONDO:0000001', edge_source='source1')
    gene_assoc_phenotype = GeneAssocPhenotype(entrez_id="12345", hpo_id='HP:0000001')

    session.add(gene_assoc_disorder)
    session.add(gene_assoc_phenotype)
    session.commit()


def example_query(session):
    # Querying the database
    results = (session
               .query(Gene, Phenotype)
               .join(GeneAssocPhenotype, Gene.entrez_id == GeneAssocPhenotype.entrez_id)
               .join(Phenotype, GeneAssocPhenotype.hpo_id == Phenotype.hpo_id)
               .filter(Gene.entrez_id == "12345")
               .all())
    for gene, phenotype in results:
        print(
            f"Gene entrez_id: {gene.entrez_id}, Phenotype hpo_id: {phenotype.hpo_id},"
            f" SNOMED ID: {phenotype.snomed_id}, OMIM ID: {phenotype.omim_id}")


def mondo_in_association_graph(mondo_id, assoc_graph):
    """
    Retrieves the genes associated with a mondo id from the association graph
    :param mondo_id: A mondo id, can be None
    :param assoc_graph: The association graph from NEDRex
    :return: genes associated with the mondo id, source of the data, harmonizome data if not in the association graph
    """
    harm = None
    if mondo_id is None:
        return None, None
    if mondo_id in assoc_graph:
        genes = []
        sources = []
        for edge in assoc_graph.edges(mondo_id, data=True):
            genes = edge[1]
            sources = edge[2]['source'][0]
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


def retrieve_disorder_data(needed_snomed, snomed_to_mondo, assoc_graph):
    """
    Queries the needed snomed ids and retrieves the associated genes and disorders from NEDRex
    :param needed_snomed: snomed ids in the dataset
    :param snomed_to_mondo: map from snomed to mondo ids from nedrex
    :param assoc_graph: association graph from nedrex of mondo ids to genes
    :return: list of genes to add, list of disorders to add, list of gene associations to add, number of snomed ids found
    """
    gene_associations = set()
    genes_to_add = set()
    disorders = set()
    found = 0
    for snomed_id in needed_snomed:
        # check if ; in snomed id and if so, do this for all ids
        snomed_ids = str(snomed_id).split(';')
        for snomed_id in snomed_ids:
            mondo_id = snomed_to_mondo.get(snomed_id)
            genes, sources = mondo_in_association_graph(mondo_id, assoc_graph)
            if genes is None:
                continue
            # add genes to set
            genes_to_add.update([Gene(entrez_id=x) for x in genes])
            disorders.add(Disorder(mondo_id=mondo_id, snomed_id=snomed_id))
            # add gene associations to set for each source
            gene_associations.update([GeneAssocDisorder(entrez_id=gene, mondo_id=mondo_id, edge_source=source)
                                      for gene, source in zip(genes, sources)])
            found += 1
        continue
    return genes_to_add, disorders, gene_associations, found


def retrieve_phenotype_data(needed_snomed, hpo_graph, data_dir):
    """
    Retrieve the phenotype data for the needed snomed ids
    # HPO conversion: Pathway
    # HPO data (HPO_ID ----> SNOMED_ID) - look for needed SNOMED IDs
    # -> map to external database (SNOMED_ID -- HPO_ID --> OMIM_ID/ORPHA_ID)
    # -> map to Mondo (SNOMED_ID -- OMIM_ID/ORPHA_ID --> Mondo_ID)
    # -> get associated genes (SNOMED_ID -- Mondo_ID --> Genes)

    :param needed_snomed: set of snomed ids that are needed
    :param hpo_graph: Graph of the HPO ontology
    :return: dictionary with the phenotype data
    """
    found = 0
    genes_to_add = set()
    phenotypes = set()
    gene_associations = set()
    needed_ids = set(needed_snomed)

    # go through all the nodes in the HPO graph and find the ones that have xrefs to SNOMED
    available_snomed_ids = snomed_from_hpo(hpo_graph, needed_ids)

    print(f'Found {len(available_snomed_ids)} snomed ids in the HPO ontology')

    # find the genes that are associated with the mondo ids
    for snomed_id, hpo_id in available_snomed_ids.items():
        hpo_id = hpo_id.replace(':', '.').replace('HP', 'hpo')
        snomed_id = f"snomedct.{snomed_id}"
        phenotypes.add(Phenotype(hpo_id=hpo_id, snomed_id=snomed_id))
        # add associations to disorders

        found += 1
    return genes_to_add, phenotypes, gene_associations, found


def add_items(session, items: iter, column: type[Gene | Phenotype | Disorder | GeneAssocDisorder | GeneAssocPhenotype | Protein],
              filter_args: list):
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
    assoc_graph = get_edge_associations(set(snomed_to_mondo.values()), edge_type='gene_associated_with_disorder')

    genes_to_add, disorders, gene_associations, found = retrieve_disorder_data(needed_snomed, snomed_to_mondo, assoc_graph)

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
    genes_to_add, phenotypes, gene_associations, _ = retrieve_phenotype_data(needed_ids, hpo_graph, data_dir)

    # since some phenotypes are subtypes of disorders, we only add phenotypes that are
    # not already in the disorder database
    removable_phenotypes = []

    for phenotype in phenotypes:
        if session.query(Disorder).filter_by(snomed_id=phenotype.snomed_id).first() is None:
            continue
        # remove phenotypes that are already in the disorder database
        removable_phenotypes.append(phenotype)

    print(f"Removing {len(removable_phenotypes)} phenotypes that are already in the disorder database")
    phenotypes = phenotypes - set(removable_phenotypes)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    add_items(session, phenotypes, Phenotype, ['hpo_id'])
    add_items(session, gene_associations, GeneAssocPhenotype, ['entrez_id', 'hpo_id'])

    session.commit()
    print(f"Found and successfully added {len(phenotypes)} snomed ids with phenotypes to db")


def add_protein_data(session, proteinData_path):
    proteinData =  read_proteinID_chris(proteinData_path)
    neddrexProteins = []
    counter = 0
    for proteinID in proteinData:
        print(proteinID)
        if counter == 1:
            break
        entrez_id = get_proteinID_neddrex(proteinID)[0]['geneName']
        print(entrez_id, proteinID)
        newProtein = Protein(uniprot_id = proteinID)#entrez_id = entrez_id, what to do if the id doesnt exist in gene? should i create a new gene?
        neddrexProteins.append(newProtein)
        #get_proteinID_neddrex(proteinID)[0]['geneName']
        counter += 1
    #print(*neddrexProteins)
    add_items(session, neddrexProteins, Protein, ['uniprot_id'])
    session.commit()
   # add_items(session, genes_to_add, Gene, ['entrez_id'])


if __name__ == '__main__':
    # Define a session
    Session = sessionmaker(bind=engine)
    session = Session()
    create_tables()
    pheno_data_path = '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv'
    protein_data_path = '../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt'
    # add_disorder_data(session, pheno_data_path)
    add_phenotype_data(session, pheno_data_path)
    # add_protein_data(session, protein_data_path)
    #'primaryDomainId': 'uniprot.P43320', 'domainIds': ['uniprot.P43320'], 'sequence': 'MASDHQTQAGKPQSLNPKIIIFEQENFQGHSHELNGPCPNLKETGVEKAGSVLVQAGPWVGYEQANCKGEQFVFEKGEYPRWDSWTSSRRTDSLSSLRPIKVDSQEHKIILYENPNFTGKKMEIIDDDVPSFHAHGYQEKVSSVRVQSGTWVGYQYPGYRGLQYLLEKGDYKDSSDFGAPHPQVQSVRRIRDMQWHQRGAFHPSN', 'displayName': 'CRBB2_HUMAN', 'synonyms': ['Beta-crystallin B2', 'Beta-B2 crystallin', 'Beta-crystallin Bp'], 'comments': 'FUNCTION: Crystallins are the dominant structural components of the vertebrate eye lens.\nSUBUNIT: Homo/heterodimer, or complexes of higher-order. The structure of beta-crystallin oligomers seems to be stabilized through interactions between the N-terminal arms (By similarity). {ECO:0000250}.\nINTERACTION: Self; NbExp=5; IntAct=EBI-974082, EBI-974082;\nDOMAIN: Has a two-domain beta-structure, folded into four very similar Greek key motifs.\nMASS SPECTROMETRY: Mass=23291; Mass_error=3; Method=Electrospray; Evidence={ECO:0000269|PubMed:8999933};\nMASS SPECTROMETRY: Mass=23289; Method=Electrospray; Evidence={ECO:0000269|PubMed:8175657};\nMASS SPECTROMETRY: Mass=23290; Method=Electrospray; Evidence={ECO:0000269|PubMed:10930324};\nDISEASE: Cataract 3, multiple types (CTRCT3) [MIM:601547]: An opacification of the crystalline lens of the eye that frequently results in visual impairment or blindness. Opacities vary in morphology, are often confined to a portion of the lens, and may be static or progressive. CTRCT3 includes congenital cerulean and sutural cataract with punctate and cerulean opacities, among others. Cerulean cataract is characterized by peripheral bluish and white opacifications organized in concentric layers with occasional central lesions arranged radially. The opacities are observed in the superficial layers of the fetal nucleus as well as the adult nucleus of the lens. Involvement is usually bilateral. Visual acuity is only mildly reduced in childhood. In adulthood, the opacifications may progress, making lens extraction necessary. Histologically the lesions are described as fusiform cavities between lens fibers which contain a deeply staining granular material. Although the lesions may take on various colors, a dull blue is the most common appearance and is responsible for the designation cerulean cataract. Sutural cataract with punctate and cerulean opacities is characterized by white opacification around the anterior and posterior Y sutures, and grayish and bluish, spindle shaped, oval punctate and cerulean opacities of various sizes arranged in lamellar form. The spots are more concentrated towards the peripheral layers and do not delineate the embryonal or fetal nucleus. Phenotypic variation with respect to the size and density of the sutural opacities as well as the number and position of punctate and cerulean spots is observed among affected subjects. {ECO:0000269|PubMed:10634616, ECO:0000269|PubMed:9158139}. Note=The disease is caused by mutations affecting the gene represented in this entry.\nSIMILARITY: Belongs to the beta/gamma-crystallin family. {ECO:0000305}.\nWEB RESOURCE: Name=Eye disease Crystallin, beta-B2 (CRYBB2); Note=Leiden Open Variation Database (LOVD); URL="http://www.lovd.nl/CRYBB2";', 'geneName': 'CRYBB2', 'taxid': 9606, 'type': 'Protein'}]
