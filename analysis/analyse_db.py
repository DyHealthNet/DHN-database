from sqlalchemy import distinct
from sqlalchemy.orm import Session
from sqlalchemy.ext.declarative import declarative_base
from utils.models import *
import inspect
import sys
import matplotlib.pyplot as plt
from utils.logger import get_logger

logger = get_logger(__name__)

all_models = [cls for name, cls in inspect.getmembers(sys.modules['utils.models']) if inspect.isclass(cls) and
              hasattr(cls, '__tablename__')]
nodes = [Gene, Protein, Phenotype, Disorder, Metabolite, GenomicVariant, CohortProtein, CohortMetabolite,
         CohortPhenotype]


def model_rows(session: Session):
    # check how many rows are in each table
    sum_all = 0
    model_counts = {}
    for model in all_models:
        if not hasattr(model, '__tablename__'):
            continue
        count = session.query(model).count()
        model_counts[model.__tablename__] = count
        sum_all += count

    model_counts['Total'] = sum_all
    return model_counts


def cumulative_rows(session: Session):
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

    return {'One': layer_one, 'Two': layer_two}


def write_csv(row_counts: dict, layer_counts: dict, filename: str = 'db_size.csv'):
    logger.debug(f"Writing the database size to {filename}")
    file = open(filename, 'w')
    file.write('Table, Count\n')
    for key, value in row_counts.items():
        if "cohort" in key or "effect" in key:
            file.write(f"{key}, {value}\n")
    file.write(f"Layer 1 total rows, {layer_counts['One']}\n")
    for key, value in row_counts.items():
        if "cohort" not in key and "effect" not in key:
            file.write(f"{key}, {value}\n")
    file.write(f"Layer 2 total rows, {layer_counts['Two']}\n")
    file.write(f"Total rows, {row_counts['Total']}")
    file.close()


def coverage(session: Session, base_reference: declarative_base, model: declarative_base,
             model_2: declarative_base = None):
    # retrieve the number of unique rows of cohort_id column
    total_cohort = session.query(base_reference.cohort_id).count()
    if not model_2:
        num_rows = session.query(distinct(model.cohort_id)).count()
        return num_rows / total_cohort

    rows_1 = session.query(distinct(model.cohort_id)).all()
    rows_2 = session.query(distinct(model_2.cohort_id)).all()
    combined_rows = set(rows_1) | (set(rows_2))

    coverage = len(combined_rows) / total_cohort
    logger.debug(f"Coverage of {model.__name__} and {model_2.__name__}: {coverage}")
    return coverage


def vis_coverage(node_coverages: dict, filename: str = 'coverage.png'):
    coverage = [('Variant', node_coverages['Variant']), ('Protein', node_coverages['Protein']),
                ('Metabolite', node_coverages['Metabolite']), ('Phenotype', node_coverages['Phenotype'])]
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
    fig.savefig(filename)
    plt.show()


def main(session: Session):
    logger.debug(f"Found {len(all_models)} models in the database")

    row_counts = model_rows(session)
    layer_counts = cumulative_rows(session)
    filename = 'db_size.csv'
    write_csv(row_counts, layer_counts, filename)
    logger.info(f"Wrote database size to {filename}")

    coverage_plot = 'coverage.png'
    node_coverages = {
        'Protein': coverage(session, CohortProtein, CohortReferencesProtein),
        'Variant': coverage(session, CohortVariant, CohortReferencesVariant),
        'Metabolite': coverage(session, CohortMetabolite, CohortReferencesMetabolite),
        'Phenotype': coverage(session, CohortPhenotype, CohortReferencesPhenotype, CohortReferencesDisease)
    }
    vis_coverage(node_coverages, coverage_plot)
    logger.info(f"Generated coverage plot: {coverage_plot}")
