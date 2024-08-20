import pandas as pd
import nedrex
from nedrex.core import api_keys_active, get_api_key
from sqlalchemy.sql.type_api import Variant

from settings import DEBUG
from models import CohortGenomicVariant, Genomic_variant, EffectVariantProtein, EffectVariantMetabolite, EffectVariantPhenotype
from nedrex.core import iter_nodes, iter_edges, get_node_types,get_collection_attributes
#nedrex.config.set_url_base("https://apps.cosy.bio/licensed")
nedrex.config.set_url_base("https://api.nedrex.net/licensed/")
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

def read_variant_meta_file(variants_meta_path: str):
    """
    reads Protein IDs from Chris dataset
    """
   # print(head(gwas_stats_df))

    variants_meta_df = pd.read_csv(variants_meta_path, sep='\t', dtype= str).drop(columns=['Unnamed: 0'])
    variantSet = set()
    for index, row in variants_meta_df.iterrows():
        newVariant=CohortGenomicVariant(
            cohort_id=row['rsid'],
            chrom = row['chrom'],
            pos = row['pos'],
            ref = row['ref'],
            alt = row['alt']

        )
        variantSet.add(newVariant)
        if(DEBUG):
            if(len(variantSet) > 10000):
                break
    return variantSet

def read_variant_gwas_file(gwas_stats_path: str):
    """
    reads Protein IDs from Chris dataset
    """
   # print(head(gwas_stats_df))

    variants_meta_df = pd.read_csv(gwas_stats_path, sep='\t', dtype= str)
    effectVariantProteinSet = set()
    effectVariantMetaboliteSet = set()
    effectVariantPhenotypeSet = set()

    for index, row in variants_meta_df.iterrows():
        if (DEBUG):
            if (len(effectVariantPhenotypeSet) > 100):
                break
        type = row['type']
        if(type == "pheno"):
            newEffectVariantPhenotype = EffectVariantPhenotype(
                phenotype_id =  row['label'],
                variant_id =  row['variant'],
                pvalue =  float(row['pvalue']),
                effect_size = float(row['eff_size'])
            )

            effectVariantPhenotypeSet.add(newEffectVariantPhenotype)
        if (type == "prot"):
            newEffectVariantProtein = EffectVariantProtein(
                protein_id=row['label'],
                variant_id=row['variant'],
                pvalue=float(row['pvalue']),
                effect_size=float(row['eff_size'])
            )

            effectVariantProteinSet.add(newEffectVariantProtein)


    return effectVariantProteinSet, effectVariantMetaboliteSet, effectVariantPhenotypeSet











# gwas_stats file contains following columns: | chr,	pos,	ref,	alt,	neglog10_pval_meta,	beta_meta,	label,	type
gwas_stats_path = '/home/leo/Documents/Uni/Masterpraktikum/data/DyHealthNet/chris_summary_data/variants/fully_simulated_gwas.tsv'
# variants_meta file contains following columns| chrom, pos, ref, alt, rsid
variants_meta_path = '/home/leo/Documents/Uni/Masterpraktikum/data/DyHealthNet/chris_summary_data/variants/variants_meta.csv'

#read_variant_files(gwas_stats_path, variants_meta_path)
















#rs562993331

#example that is also in clinvar: rs121917870
#https://www.ncbi.nlm.nih.gov/snp/?term=rs121917870%5BReference%20SNP%20ID%5D
#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
#$url = $base . "epost.fcgi?db=$db1&id=$id_list";


#print(rsIDset[:2])





