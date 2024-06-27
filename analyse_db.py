#%%
from sqlalchemy.orm import sessionmaker
from setup_db import engine
from models import *

#%%
Session = sessionmaker(bind=engine)
session = Session()

all_models = [Gene, Protein, Phenotype, Disorder, Metabolite, MetaboliteAssocDisorder, ProteinAssocMetabolite,
              GeneAssocDisorder, DisorderAssocPhenotype, ProteinAssocProtein, EffectsDisorderDisorder,
              EffectsMetaboliteDisorder, EffectsMetaboliteMetabolite, EffectsMetabolitePhenotype,
              EffectsPhenotypeDisorder, EffectsPhenotypePhenotype, EffectsProteinDisorder, EffectsProteinMetabolite,
              EffectsProteinPhenotype, EffectsProteinProtein, Genomic_variant, Variant_affects_gene]
print(len(all_models))

#%%
# check how many rows are in each table
sum_all = 0
for model in all_models:
    count = session.query(model).count()
    print(f"{model.__name__}: {count}")
    sum_all += count

print(f"Total rows: {sum_all}")
