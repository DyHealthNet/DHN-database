#%%
from sqlalchemy.orm import sessionmaker
from setup_db import engine
from models import Gene, Protein, Phenotype, Disorder
from query_nedrex import needed_snomed_ids
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

#%%
Session = sessionmaker(bind=engine)
session = Session()

#%%
# check how many rows are in each table
genes = session.query(Gene).count()
proteins = session.query(Protein).count()
phenotypes = session.query(Phenotype).count()
disorders = session.query(Disorder).count()

f"Genes: {genes}, Proteins: {proteins}, Disorders: {disorders}, Phenotypes: {phenotypes}"

#%%
# get the overlap between phenotypes and disorders
phenotype_overlap = session.query(Phenotype, Disorder).filter(Phenotype.snomed_id == Disorder.snomed_id).all()
for phenotype, disorder in phenotype_overlap:
    print(f"Phenotype: {phenotype.hpo_id}, Disorder: {disorder.mondo_id}")
    print(f"Phenotype: {phenotype.snomed_id}, Disorder: {disorder.snomed_id}")

#%%
print(f"Phenotypes that are also considered disorders: {len(phenotype_overlap)}")

#%%
data_path = '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv'
needed_snomed = needed_snomed_ids(data_path)

#%%
# create venn diagram of phenotypes that are also disorders
venn2(subsets=(disorders - len(phenotype_overlap), phenotypes - len(phenotype_overlap), len(phenotype_overlap)),
      set_labels=('Disorders', 'Phenotypes'))
plt.savefig('phenotype_disorder_overlap_venn.png')
plt.show()

#%%
