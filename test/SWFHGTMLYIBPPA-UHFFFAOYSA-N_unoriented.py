#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(0)
G.add_node(6)
G.add_node(13)
G.add_node(14)
G.add_node(27)

G.add_edge(0,13)
G.add_edge(0,6)
G.add_edge(6,27)
G.add_edge(6,13)
G.add_edge(6,14)
G.add_edge(14,27)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/SWFHGTMLYIBPPA-UHFFFAOYSA-N_unoriented.png')
