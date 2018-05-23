#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(0)
G.add_node(7)
G.add_node(9)
G.add_node(10)
G.add_node(19)

G.add_edge(0,9)
G.add_edge(0,7)
G.add_edge(7,9)
G.add_edge(7,10)
G.add_edge(7,19)
G.add_edge(10,19)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/ZKZFPRUSWCYSGT-UHFFFAOYSA-N_unoriented.png')
