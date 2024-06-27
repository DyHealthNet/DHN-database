from nedrex.core import iter_edges

from models import *



from sqlalchemy import Column, String, Integer, ForeignKey, ARRAY, Float
from sqlalchemy.orm import declarative_base

# Create a declarative base
Base = declarative_base()



import unittest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from models import *

def test_gene():
    testgene = Gene(
        entrez_id="entrez.1",
        display_name="Test Gene",
        description="This is a test gene",
        synonyms=["Test1", "Test2"],
        chromosome="1",
        observation_source="Test Source"
    )
    return(testgene)

def test_disorder():
    disorder = Disorder(
        mondo_id="MONDO:0001",
        description="Test Disorder",
        xrefs=["xref1", "xref2"],
        observation_source="Test Source"
    )
    return(disorder)
def test_phenotype():
    phenotype = Phenotype(
        hpo_id="HP:0001",
        display_name="Test Phenotype",
        description="This is a test phenotype",
        xrefs=["xref1", "xref2"],
        synonyms="Test Synonym",
        observation_source="Test Source"
    )
    return(phenotype)

def test_protein():
    protein = Protein(
        uniprot_id="P12345",
        sequence="MSEQENCE",
        gene_entrez_id="entrez.1",
        description="Test Protein",
        observation_source="Test Source"
    )
    return(protein)

def test_metabolite():
    metabolite = Metabolite(
        hmdb_id="HMDB00001",
        display_name="Test Metabolite",
        description="This is a test metabolite",
        synonyms="Synonym1, Synonym2",
        xrefs=["xref1", "xref2"],
        observation_source="Test Source"
    )
    return(metabolite)

# Additional test cases for association and effects tables can follow the same pattern
def test_gene_assoc_disorder():
    gene = Gene(entrez_id="entrez.1", display_name="Test Gene")
    disorder = Disorder(mondo_id="MONDO:0001", description="Test Disorder")

def test_disorder_assoc_phenotype():
    disorder = Disorder(mondo_id="MONDO:0001", description="Test Disorder")
    phenotype = Phenotype(hpo_id="HP:0001", display_name="Test Phenotype")


def test_proteinAssocProtein(counter=5):
    counter = counter


def test_proteinAssocProtein(NoOfAssoc=5):
    counter = 0
    protein_interactions = []
    #edges = [e for e in iter_edges(edge_type) if e[first_node] in node_ids or e[second_node] in node_ids]
    for edge in iter_edges('protein_interacts_with_protein'):
        uniprot_id_memberOne = String(edge['memberOne'])
        uniprot_id_memberTwo = String(edge['memberTwo'])
        interaction = ProteinAssocProtein( uniprot_id_memberOne=uniprot_id_memberOne,
                                           uniprot_id_memberTwo=uniprot_id_memberTwo)
        protein_interactions.append(interaction)
        counter += 1
        if (counter == NoOfAssoc):
            break
    return (protein_interactions)
def test_variant():
    # Create a genomic variant
    variant = Genomic_variant(
        variant_primaryDomainId="clinvar.17735",
        alternativeSequence="T",
        chromosome="NW_009646201.1",
        created="2024-06-17T12:36:21.275000",
        dataSources="clinvar",
        domainIds="clinvar.17735,dbsnp.1556058284",
        position="83614",
        referenceSequence="TC",
        type="GenomicVariant",
        variantType="Deletion"
    )
    return(variant)


def test_variant_affects_gene():
    # Create a genomic variant
    variant = Genomic_variant(
        variant_primaryDomainId="clinvar.17735",
        alternativeSequence="T",
        chromosome="NW_009646201.1",
        created="2024-06-17T12:36:21.275000",
        dataSources="clinvar",
        domainIds="clinvar.17735,dbsnp.1556058284",
        position="83614",
        referenceSequence="TC",
        type="GenomicVariant",
        variantType="Deletion"
    )
    self.session.add(variant)

    # Create a gene
    gene = Gene(
        entrez_id="entrez.1",
        display_name="Test Gene",
        description="This is a test gene",
        synonyms=["Test1", "Test2"],
        chromosome="1",
        observation_source="Test Source"
    )



def test_protein_assoc_metabolite():
    # Create a protein
    protein = Protein(
        uniprot_id="P12345",
        sequence="MKT...",
        gene_entrez_id="entrez.1",
        description="Test protein",
        observation_source="Test source"
    )


    # Create a metabolite
    metabolite = Metabolite(
        hmdb_id="HMDB00001",
        display_name="Test Metabolite",
        description="This is a test metabolite",
        synonyms="Metab1, Metab2",
        xrefs=["xref1", "xref2"],
        observation_source="Test Source"
    )
def test_metabolite_assoc_disorder():
    # Create a metabolite
    metabolite = Metabolite(
        hmdb_id="HMDB00002",
        display_name="Test Metabolite 2",
        description="This is another test metabolite",
        synonyms="Metab3, Metab4",
        xrefs=["xref3", "xref4"],
        observation_source="Test Source 2"
    )

    # Create a disorder
    disorder = Disorder(
        mondo_id="MONDO:0002",
        description="Test Disorder 2",
        xrefs=["xref5", "xref6"],
        observation_source="Test Source 3"
    )
