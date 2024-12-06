import rdflib
from rdflib.extras.external_graph_libs import rdflib_to_networkx_digraph
import pydot

# Load the TTL file into an RDF graph
g = rdflib.Graph()
g.parse("/Users/alexanderfarnum/Documents/Bio/Data/A_Box/pura.ttl", format="turtle")

# Convert the RDF graph to a DOT graph
dot_string = rdflib_to_networkx_digraph(g)

# Create a pydot graph from the DOT string
graph = pydot.graph_from_dot_data(dot_string)

# Write the DOT graph to a file
graph.write_dot("output.dot")