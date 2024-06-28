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
nodes = [Gene, Protein, Phenotype, Disorder, Metabolite, Genomic_variant]
print(len(all_models))

#%%
# check how many rows are in each table
sum_all = 0
for model in nodes:
    if not hasattr(model, '__tablename__'):
        continue
    count = session.query(model).count()
    print(f"{model.__name__}: {count}")
    sum_all += count

print(f"Total rows: {sum_all}")
