#%%
from sqlalchemy.orm import sessionmaker
from setup_db import engine
from models import *
import inspect
import sys

#%%
Session = sessionmaker(bind=engine)
session = Session()

all_models = [cls for name, cls in inspect.getmembers(sys.modules['models']) if inspect.isclass(cls)]
nodes = [Gene, Protein, Phenotype, Disorder, Metabolite, Genomic_variant, CohortProtein, CohortMetabolite,
         CohortPhenotype]
print(len(all_models))

#%%
# check how many rows are in each table
sum_all = 0
for model in all_models:
    if not hasattr(model, '__tablename__'):
        continue
    count = session.query(model).count()
    print(f"{model.__name__}: {count}")
    sum_all += count

print(f"Total rows: {sum_all}")

#%%
# check the two layers of the database, first layer is the cohort & calculated stuff, second layer is
# the external knowledge graph
layer_one = 0
layer_two = 0
names_layer_one = []
names_layer_two = []

for model in all_models:
    if not hasattr(model, '__tablename__'):
        continue
    if model.__name__.startswith('Cohort') or model.__name__.startswith('Effects'):
        layer_one += session.query(model).count()
        names_layer_one.append(model.__name__)
    else:
        layer_two += session.query(model).count()
        names_layer_two.append(model.__name__)

print(f"Cumulative rows in layer one: {layer_one:,}")
print(f"Cumulative rows in layer two: {layer_two:,}")
