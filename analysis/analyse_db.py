#%%
from sqlalchemy.orm import sessionmaker
from sqlalchemy import distinct
from setup_db import engine
from utils.models import *
import inspect
import sys

#%%
Session = sessionmaker(bind=engine)
session = Session()

all_models = [cls for name, cls in inspect.getmembers(sys.modules['utils.models']) if inspect.isclass(cls) and
              hasattr(cls, '__tablename__')]
nodes = [Gene, Protein, Phenotype, Disorder, Metabolite, GenomicVariant, CohortProtein, CohortMetabolite,
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


#%%
# calculate the coverage
def coverage(base_reference, model, model_2 = None):
    # retrieve the number of unique rows of cohort_id column
    total_cohort = session.query(base_reference.cohort_id).count()
    if not model_2:
        num_rows = session.query(distinct(model.cohort_id)).count()
        return num_rows / total_cohort

    rows_1 = session.query(distinct(model.cohort_id)).all()
    rows_2 = session.query(distinct(model_2.cohort_id)).all()
    combined_rows = set(rows_1) | (set(rows_2))

    return len(combined_rows) / total_cohort


variant_coverage = coverage(CohortVariant, CohortReferencesVariant)
print(f"Variant coverage: {variant_coverage:.2f}")
protein_coverage = coverage(CohortProtein, CohortReferencesProtein)
print(f"Protein coverage: {protein_coverage:.2f}")
metabolite_coverage = coverage(CohortMetabolite, CohortReferencesMetabolite)
print(f"Metabolite coverage: {metabolite_coverage:.2f}")
phenotype_coverage = coverage(CohortPhenotype, CohortReferencesPhenotype, CohortReferencesDisease)
print(f"Phenotype coverage: {phenotype_coverage:.2f}")

#%%
import matplotlib.pyplot as plt

coverage = [('Variant', variant_coverage), ('Protein', protein_coverage),
            ('Metabolite', metabolite_coverage), ('Phenotype', phenotype_coverage)]
coverage.sort(key=lambda x: x[1], reverse=True)

# plot the coverage
fig, ax = plt.subplots()

ax.bar([x[0] for x in coverage], [x[1] * 100 for x in coverage])
for i, v in enumerate([x[1] * 100 for x in coverage]):
    ax.text(i, v + 1, f"{v:.2f}%", ha='center', va='bottom')
ax.set_ylabel('Coverage in %')
ax.set_xlabel('Node type')
ax.set_title('Coverage of the different node types')
fig.tight_layout()
fig.savefig('coverage.png')
plt.show()
