import pandas as pd

from utils.query_nedrex import get_edge_associations
from settings import DEBUG
from utils.models import CohortGenomicVariant, Genomic_variant, EffectVariantProtein, EffectVariantMetabolite, \
    EffectVariantPhenotype, Variant_affects_gene, Gene, CohortReferencesVariant
from nedrex.core import iter_nodes


def read_rsid_chris(variant_data_path: str):
    """
    reads Protein IDs from CHRIS dataset
    """
    df = pd.read_csv(variant_data_path, sep='\t')
    # check how many nans in the uniprot ocl
    print("Number of nans in variant col:", df['rsid'].isna().sum())
    df['rsid'] = df['rsid'].fillna('')

    rs_id_list = df[['rsid', 'pos', 'ref', 'alt']]
    return rs_id_list


def get_genomic_variant_nodes(rs_id_from_cohort: pd.DataFrame, obs_source: str):
    print("Genomic_variant IDs:", len(rs_id_from_cohort))
    no_ids = len(rs_id_from_cohort)
    rs_ids_with_alt_seq = set(zip(rs_id_from_cohort['rsid'], rs_id_from_cohort['alt']))

    found_genomic_variants = []
    genomic_variant_count = 0
    for node in iter_nodes('genomic_variant'):
        if genomic_variant_count >= no_ids:
            break
        rs_id_genomic_variant = node['domainIds']
        if len(rs_id_genomic_variant) > 1:
            rs_id_genomic_variant = node['domainIds'][1].replace('dbsnp.', 'rs')
        else:
            continue
        alternate_sequence = node.get('alternativeSequence')
        if (rs_id_genomic_variant, alternate_sequence) in rs_ids_with_alt_seq:
            genomic_variant = Genomic_variant(
                clinvar_id=str(node.get('primaryDomainId')),
                alternativeSequence=str(node.get('alternativeSequence')),
                chromosome=str(node.get('chromosome')),
                dataSources=str(node.get('dataSources')),
                xrefs=str(node.get('domainIds')),
                position=str(node.get('position')),
                referenceSequence=str(node.get('referenceSequence')),
                type=str(node.get('type')),
                variantType=str(node.get('variantType')),
                observation_source=obs_source
            )
            found_genomic_variants.append(genomic_variant)
            genomic_variant_count += 1
        if DEBUG:
            if genomic_variant_count > 100:
                break

    return found_genomic_variants


def read_variant_meta_file(variants_meta_path: str):
    """
    reads Protein IDs from Chris dataset
    """
    variants_meta_df = pd.read_csv(variants_meta_path, sep='\t', dtype=str)
    print("Iterating through gwas file")
    variant_set = set()
    for index, row in variants_meta_df.iterrows():
        new_variant = CohortGenomicVariant(
            cohort_id=f"{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}",
            description=row['rsid'],
            display_name=f"{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}",
            xrefs=f"rsid.{row['rsid']}",
        )
        variant_set.add(new_variant)
        if DEBUG:
            if len(variant_set) > 1000:
                break
    return variant_set


def read_variant_gwas_file(gwas_stats_path: str):
    """
    reads Protein IDs from Chris dataset
    """
    variants_meta_df = pd.read_csv(gwas_stats_path, sep='\t', dtype=str)
    effect_variant_protein_set = set()
    effect_variant_metabolite_set = set()
    effect_variant_phenotype_set = set()

    variant_effect_type = {'pheno': (EffectVariantPhenotype, effect_variant_phenotype_set, "phenotype_id"),
                           'prot': (EffectVariantProtein, effect_variant_protein_set, "protein_id"),
                           'metab': (EffectVariantMetabolite, effect_variant_metabolite_set, "metabolite_id")}

    for index, row in variants_meta_df.iterrows():
        if DEBUG and (len(effect_variant_phenotype_set) > 100 and
                      len(effect_variant_metabolite_set) > 100 and
                      len(effect_variant_protein_set) > 100):
            break
        effect_type = row['type']
        effect_class, effect_set, id_type = variant_effect_type[effect_type]

        effect_values = {id_type: row['label2'],
                         'variant_id': row['label1'].replace("chr", ""),
                         'p_value': float(row['pval']),
                         'effect_size': float(row['effsize']),
                         'effect_size_type': row['effsize_type'],
                         'test_statistic': row['test']}

        new_effect = effect_class(**effect_values)
        effect_set.add(new_effect)

    return effect_variant_protein_set, effect_variant_metabolite_set, effect_variant_phenotype_set


def get_cohort_references_variant(session, obs_source):
    genomic_variants = session.query(Genomic_variant).all()
    new_cohort_references_set = set()
    query_result = session.query(CohortGenomicVariant).all()

    existing_cohort_id = {(genomic_variant.description, f"{genomic_variant.cohort_id[-1]}")
                          for genomic_variant in query_result}

    desc_map = {f"{genomic_variant.description}{genomic_variant.cohort_id[-1]}": genomic_variant.cohort_id
                for genomic_variant in query_result}

    for variant in genomic_variants:
        variant_domain_ids = variant.xrefs.replace(",", "").replace("[", "").replace("]", "").replace("'", "").split()
        dbsnp_id = next((variant_id.replace("dbsnp.", "rs") for variant_id in variant_domain_ids
                         if "dbsnp." in variant_id), None)
        clinvar_id = variant.clinvar_id
        alt_seq = variant.alternativeSequence

        if (dbsnp_id, alt_seq) in existing_cohort_id:
            cohort_id = desc_map[f"{dbsnp_id}{alt_seq}"]
            new_cohort_references_variant = CohortReferencesVariant(
                cohort_id=cohort_id,
                clinvar_id=clinvar_id
            )
            new_cohort_references_set.add(new_cohort_references_variant)

    return new_cohort_references_set


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
        sourceDomainId = edge[0]
        targedDomainId = edge[1]
        id_list_total.append(sourceDomainId)
        id_list_total.append(
            targedDomainId)  # because graph is directed either first or second entry contains the entrez id

        if sourceDomainId.startswith("entrez."):
            entrez_id = sourceDomainId
            variant = targedDomainId
        if targedDomainId.startswith("entrez."):
            entrez_id = targedDomainId
            variant = sourceDomainId
        variant_affects_gene_edge = Variant_affects_gene(clinvar_id=variant,
                                                         entrez_id=entrez_id)
        variant_affects_gene_to_add.add(variant_affects_gene_edge)
    entrez_ids = {id for id in id_list_total if id.startswith("entrez.")}
    genomic_variant_node_generator = iter_nodes('gene')
    genes_to_add = set()
    for node in genomic_variant_node_generator:
        if node['primaryDomainId'] in id_list_total:
            new_gene = Gene(entrez_id=node['primaryDomainId'],
                            display_name=node['displayName'],
                            description=node['description'],
                            synonyms=node['synonyms'],
                            chromosome=node['chromosome'],
                            observation_source=obs_source)
            genes_to_add.add(new_gene)

    # stmt = select([genomic_variants.c.id]).where(genomic_variants.c.id == id_to_check)
    # result = db_session.execute(stmt).first()
    return genes_to_add, variant_affects_gene_to_add


def add_genomic_variants_linking_from_gene(session, entrez_ids: set[str] = None, observation_source: str = None):
    variant_affects_gene_graph = get_edge_associations(node_ids=entrez_ids, edge_type='variant_affects_gene',
                                                       direction='directed')
    print(f"Got {len(variant_affects_gene_graph)} edge associations for variant_affects_gene.")
    variant_affects_gene_dict = {}
    clinvar_ids = set()
    for edge in variant_affects_gene_graph.edges(data=True):
        clinvar_ids.add(edge[0])
        variant_affects_gene_dict[edge[0]] = edge[1]

    pattern = r'^[^.]*\.'
    variants_to_add = []
    genomic_variant_node_generator = iter_nodes('genomic_variant')
    for node in genomic_variant_node_generator:
        if node['primaryDomainId'] in clinvar_ids:
            newVariant = Genomic_variant(clinvar_id=node['primaryDomainId'],
                                         alternativeSequence=node['alternativeSequence'],
                                         chromosome=node['chromosome'],
                                         dataSources=node['dataSources'],
                                         xrefs=node['domainIds'],
                                         position=node['position'],
                                         referenceSequence=node['referenceSequence'],
                                         type=node['type'],
                                         variantType=node['variantType'])
            variants_to_add.append(newVariant)
    variant_affects_gene_to_add = []
    genes = []
    available_variants = {variant.clinvar_id for variant in variants_to_add}
    for variant, gene in variant_affects_gene_dict.items():
        if (variant in available_variants and gene in entrez_ids):
            variant_affects_gene_to_add.append(Variant_affects_gene(clinvar_id=variant, entrez_id=gene))
    return variants_to_add, variant_affects_gene_to_add