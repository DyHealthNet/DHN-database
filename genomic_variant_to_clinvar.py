import pandas as pd
import nedrex
from nedrex.core import api_keys_active, get_api_key
from settings import DEBUG
from models import Genomic_variant
from nedrex.core import iter_nodes, iter_edges, get_node_types,get_collection_attributes
nedrex.config.set_url_base("https://apps.cosy.bio/licensed")
if api_keys_active():
    api_key = get_api_key(accept_eula=True)
    nedrex.config.set_api_key(api_key)

rsidPath='/home/leo/Documents/Uni/Masterpraktikum/toy_data/variants_meta.csv'
def read_rsid_chris(variantDataPath: str):
    """
    reads Protein IDs from Chris dataset
    """
    df = pd.read_csv(variantDataPath, sep='\t')
    # check how many nans in the uniprot ocl
    print("Number of nans in variant col:", df['rsid'].isna().sum())
    df['rsid'] = df['rsid'].fillna('')
    rsids = df['rsid'].unique().tolist()
    #uniprot_ids = set([uniprot for sublist in uniprot_ids for uniprot in sublist])
    return rsids


def get_genomic_variant_nodes(clinvarIds , observation_source):
    print("UniProt IDs:", len(clinvarIds))
    #clinvarIds = ['clinvar.' + str(item) for item in clinvarIds]
    noIds = len(clinvarIds)
    foundGenomicVariants = []
    found_proteins = 0
    for node in iter_nodes('genomic_variant'):
        if (found_proteins >= noIds):
            break
        primary_domain_id = node['domainIds']
        if len(primary_domain_id) > 1:
            primary_domain_id = node['domainIds'][1]
        else:
            continue
        if(primary_domain_id in clinvarIds):
            genomic_variant = Genomic_variant(
                variant_primaryDomainId= str(node.get('primaryDomainId')),
                # Column(String, primary_key=True)  # clinvar.17735
                alternativeSequence=str(node.get('alternativeSequence')),  # Column(String)  # 'T',
                chromosome=str(node.get('chromosome')),  # Column(String)  # 'NW_009646201.1',
                created=str(node.get('created')),  # Column(String)  # '2024-06-17T12:36:21.275000'
                dataSources=str(node.get('dataSources')),  # Column(String)  # ['clinvar'],
                domainIds=str(node.get('domainIds')),  # Column(String)  # ['clinvar.17735', 'dbsnp.1556058284']
                position=str(node.get('position')),  # Column(String)  # 83614,
                referenceSequence=str(node.get('referenceSequence')),  # Column(String)  # 'TC',
                type=str(node.get('type')),  # Column(String)  # 'GenomicVariant'
                variantType=str(node.get('variantType'))  # Column(String)  # 'Deletion'}
                #observation_source=observation_source

                )
            foundGenomicVariants.append(genomic_variant)
            found_proteins += 1
        if DEBUG:
            if found_proteins > 100:
                break

    return foundGenomicVariants

xml_file_path = '/home/leo/Documents/Uni/Masterpraktikum/data/ClinVarVCVRelease_00-latest.xml.gz'
summaryData = '/home/leo/Documents/Uni/Masterpraktikum/variant_summary.txt'

"""
rsIDset = read_rsid_chris(rsidPath)
print(get_collection_attributes("genomic_variant",True))
updated_ids_set = {id_.replace('rs', '') for id_ in rsIDset}
print("#Ids from Dataset" ,len(updated_ids_set))
#updated_ids_set.add("397704705")
columns_to_use = ['RS# (dbSNP)', 'VariationID']
idMappingDf = pd.read_csv(summaryData, sep='\t', usecols=columns_to_use )
print(idMappingDf.head(5))
print(idMappingDf.dtypes)
subset_df = idMappingDf[idMappingDf['RS# (dbSNP)'].astype(str).isin(updated_ids_set)].drop_duplicates()

subset_df.to_csv('/home/leo/Documents/Uni/Masterpraktikum/data/output_file.tsv', sep='\t', index=False)
print("#Ids that got matched", len(subset_df))
print(subset_df.head(10))
clinvarIds = subset_df["VariationID"].tolist()
clinvarIds = {id_.replace('rs', 'dbsnp.') for id_ in rsIDset}
#get_genomic_variant_nodes(clinvarIds)

print(len(clinvarIds))


"""
def read_variant_files(gwas_stats_path: str,variants_meta_path: str):
    """
    reads Protein IDs from Chris dataset
    """
   # print(head(gwas_stats_df))
    gwas_stats_df = pd.read_csv(variants_meta_path, sep='\t', dtype= str)
    print("GWAS: ", gwas_stats_df.shape)
    variants_meta_df = pd.read_csv(gwas_stats_path, sep='\t', dtype =str)
    print("VARIANTS: ", variants_meta_df.shape)
    if 'Unnamed: 0' in variants_meta_df.columns:
        variants_meta_df = variants_meta_df.drop(columns=['Unnamed: 0'])
    if 'Unnamed: 0' in gwas_stats_df.columns:
        gwas_stats_df = gwas_stats_df.drop(columns=['Unnamed: 0'])

    merged_df = pd.merge(variants_meta_df, gwas_stats_df, left_on=['chr', 'pos', 'ref', 'alt'],
                         right_on=['chrom', 'pos', 'ref', 'alt'])
    print("MERGED: ", merged_df.shape)
    """
    df = pd.read_csv(variantDataPath, sep='\t')
    # check how many nans in the uniprot ocl
    print("Number of nans in variant col:", df['rsid'].isna().sum())
    df['rsid'] = df['rsid'].fillna('')
    rsids = df['rsid'].unique().tolist()
    #uniprot_ids = set([uniprot for sublist in uniprot_ids for uniprot in sublist])
    return rsids
    """

# gwas_stats file contains following columns: | chr,	pos,	ref,	alt,	neglog10_pval_meta,	beta_meta,	label,	type
gwas_stats_path = '/home/leo/Documents/Uni/Masterpraktikum/toy_data/medium/gwas_stats.csv'
# variants_meta file contains following columns| chrom, pos, ref, alt, rsid
variants_meta_path = '/home/leo/Documents/Uni/Masterpraktikum/toy_data/variants_meta.csv'

#read_variant_files(gwas_stats_path, variants_meta_path)
















#rs562993331

#example that is also in clinvar: rs121917870
#https://www.ncbi.nlm.nih.gov/snp/?term=rs121917870%5BReference%20SNP%20ID%5D
#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
#$url = $base . "epost.fcgi?db=$db1&id=$id_list";


#print(rsIDset[:2])





