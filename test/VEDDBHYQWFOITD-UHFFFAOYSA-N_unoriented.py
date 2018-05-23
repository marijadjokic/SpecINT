#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(32)
G.add_node(34)
G.add_node(35)
G.add_node(39)
G.add_node(10)
G.add_node(12)
G.add_node(14)
G.add_node(15)
G.add_node(19)

G.add_edge(32,34)
G.add_edge(32,35)
G.add_edge(32,10)
G.add_edge(32,39)
G.add_edge(34,10)
G.add_edge(34,35)
G.add_edge(34,39)
G.add_edge(35,10)
G.add_edge(35,39)
G.add_edge(39,10)
G.add_edge(10,12)
G.add_edge(10,14)
G.add_edge(10,15)
G.add_edge(10,19)
G.add_edge(12,19)
G.add_edge(12,14)
G.add_edge(12,15)
G.add_edge(14,19)
G.add_edge(14,15)
G.add_edge(15,19)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/VEDDBHYQWFOITD-UHFFFAOYSA-N_unoriented.png')
