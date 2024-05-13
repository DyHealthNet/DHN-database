from query_nedrex import get_disorder_data, domain_id_to_mondo, needed_snomed_ids


if __name__ == '__main__':
    needed_snomed_ids = needed_snomed_ids('../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv')
    disorder_data = get_disorder_data()