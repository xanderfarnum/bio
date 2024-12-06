import rdflib
import networkx as nx
import matplotlib.pyplot as plt

# Create an RDF graph
g = rdflib.Graph()
g.parse("/Users/alexanderfarnum/Documents/Bio/Data/A_Box/pura.ttl", format="turtle")

# Convert RDF graph to NetworkX graph
G = nx.Graph()
for s, p, o in g:
    G.add_edge(s, o, label=p)

# Visualize the graph
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_size=1000)
nx.draw_networkx_edge_labels(G, pos)
plt.show()