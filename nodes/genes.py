from utils.query_nedrex import get_gene_data
from utils.models import Gene, ProteinAssocGene


def get_missing_genes(missing_genes: set) -> set[Gene]:
    gene_data = get_gene_data(missing_genes, search_col='displayName')
    genes = set()
    for fetched_gene in gene_data:
        gene = Gene(
            entrez_id=fetched_gene['primaryDomainId'],
            display_name=fetched_gene['displayName'],
            description=fetched_gene['description'],
            synonyms=fetched_gene['synonyms'],
            chromosome=fetched_gene['chromosome'],
            observation_source='external'
        )
        genes.add(gene)
    return genes


def gene_associations(protein_gene_map: set, existing_genes: dict, missing_genes: set) \
        -> tuple[set[Gene], set[ProteinAssocGene]]:
    genes = set()
    if missing_genes:
        genes = get_missing_genes(missing_genes)
    gene_map = {x.display_name: x.entrez_id for x in genes}
    gene_map.update(existing_genes)
    protein_assoc_genes = set()
    for protein, gene in protein_gene_map:
        entrez_id = gene_map.get(gene)
        if not entrez_id:
            continue

        protein_assoc_gene = ProteinAssocGene(
            uniprot_id=protein,
            entrez_id=entrez_id
        )
        protein_assoc_genes.add(protein_assoc_gene)
    return genes, protein_assoc_genes
