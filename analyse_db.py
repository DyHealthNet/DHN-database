#%%
from sqlalchemy.orm import sessionmaker
from setup_db import engine
from models import *

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
# get the edges between each table
gene_disorder_edges = session.query(GeneAssocDisorder).count()
disorder_phenotype_edges = session.query(DisorderAssocPhenotype).count()
protein_assoc_protein_edges = session.query(ProteinAssocProtein).count()
protein_assoc_gene_edges = session.query(Protein).filter(Protein.gene_entrez_id is not None).count()

(f"Gene-Disorder edges: {gene_disorder_edges}, Disorder-Phenotype edges: {disorder_phenotype_edges}, "
 f"Protein-Protein edges: {protein_assoc_protein_edges}, Protein-Gene edges: {protein_assoc_gene_edges}")
