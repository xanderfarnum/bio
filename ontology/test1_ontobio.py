from ontobio.ontol_factory import OntologyFactory
from ontobio.assoc_factory import AssociationSetFactory


ont = OntologyFactory().create("go")
[nucleus] = ont.search('nucleus')
ancestors = ont.ancestors(nucleus)

afactory = AssociationSetFactory()
aset = afactory.create(ontology=ont,
                       subject_category='gene',
                       object_category='function',
                       taxon='NCBITaxon:7955') ## Zebrafish