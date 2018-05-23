#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(0)
G.add_node(1)
G.add_node(3)
G.add_node(5)
G.add_node(7)

G.add_edge(0,1)
G.add_edge(0,3)
G.add_edge(0,5)
G.add_edge(0,7)
G.add_edge(1,3)
G.add_edge(5,7)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/LBJMHPIVBXUBLY-UHFFFAOYSA-N_unoriented.png')
