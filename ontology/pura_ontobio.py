from ontobio.ontol_factory import OntologyFactory
from ontobio.assoc_factory import AssociationSetFactory

isoform = 'UniProtKB:Q00577'
protein_id = 'GO_REF:0000024'

ont = OntologyFactory().create("go")

[nucleus] = ont.search('nucleus')
ancestors = ont.ancestors(nucleus)

afactory = AssociationSetFactory()
aset = afactory.create(ontology=ont,
                       subject_category='gene',
                       object_category='function',
                       taxon='NCBITaxon:7955')