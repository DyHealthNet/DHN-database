import pandas as pd

from utils.query_nedrex import get_edge_associations
from utils.settings import DEBUG
from utils.models import GenomicVariant, Variant_affects_gene, Gene
from nedrex.core import iter_nodes
from utils.logger import get_logger

logger = get_logger(__name__)


def read_rsid_chris(variant_data_path: str) -> set[tuple[str, str]]:
    """
    reads variant IDs from CHRIS dataset
    """
    df = pd.read_csv(variant_data_path, sep='\t')
    logger.debug(f"Number of nans in variant col: {df['rsid'].isna().sum()}")
    df['rsid'] = df['rsid'].fillna('')
    return set(zip(df['rsid'], df['alt']))


def get_genomic_variant_nodes(rs_id_from_cohort: set, obs_source: str) -> list[GenomicVariant]:
    logger.debug(f"# of Genomic Variant IDs: {len(rs_id_from_cohort)}")
    no_ids = len(rs_id_from_cohort)

    found_genomic_variants = []
    genomic_variant_count = 0
    for node in iter_nodes('genomic_variant'):
        if genomic_variant_count >= no_ids:
            break

        rs_id_variant = next((x.replace('dbsnp.', 'rs') for x in node['domainIds'] if 'dbsnp.' in x), None)

        alternate_sequence = node.get('alternativeSequence')
        if (rs_id_variant, alternate_sequence) in rs_id_from_cohort:
            genomic_variant = GenomicVariant(
                clinvar_id=node.get('primaryDomainId'),
                alternative_sequence=node.get('alternativeSequence'),
                chromosome=node.get('chromosome'),
                data_sources=node.get('dataSources'),
                xrefs=node.get('domainIds'),
                position=node.get('position'),
                reference_sequence=node.get('referenceSequence'),
                type=node.get('type'),
                variant_type=node.get('variantType'),
                observation_source=obs_source
            )
            found_genomic_variants.append(genomic_variant)
            genomic_variant_count += 1
        if DEBUG and genomic_variant_count > 1000:
            break

    return found_genomic_variants


def add_variant_affects_gene(clinvar_ids: set[str], obs_source: str = "external"):
    variant_affects_gene_graph = get_edge_associations(node_ids=clinvar_ids, edge_type='variant_affects_gene',
                                                       direction='directed')
    variant_affects_gene_dict = {}
    id_list_total = []  # to check if is in database
    variant_edge_list = []
    variant_affects_gene_to_add = set()
    for edge in variant_affects_gene_graph.edges(data=True):
        variant_affects_gene_dict[edge[0]] = edge[1]
        variant_edge_list.append(edge)
        source_domain_id = edge[0]
        targed_domain_id = edge[1]
        id_list_total.append(source_domain_id)
        id_list_total.append(
            targed_domain_id)  # because graph is directed either first or second entry contains the entrez id

        if source_domain_id.startswith("entrez."):
            entrez_id = source_domain_id
            variant = targed_domain_id
        elif targed_domain_id.startswith("entrez."):
            entrez_id = targed_domain_id
            variant = source_domain_id
        else:
            continue
        variant_affects_gene_edge = Variant_affects_gene(clinvar_id=variant,
                                                         entrez_id=entrez_id)
        variant_affects_gene_to_add.add(variant_affects_gene_edge)

    genes_to_add = set()
    for node in iter_nodes('gene'):
        if node['primaryDomainId'] in id_list_total:
            new_gene = Gene(entrez_id=node['primaryDomainId'],
                            display_name=node['displayName'],
                            description=node['description'],
                            synonyms=node['synonyms'],
                            chromosome=node['chromosome'],
                            observation_source=obs_source)
            genes_to_add.add(new_gene)

    return genes_to_add, variant_affects_gene_to_add
