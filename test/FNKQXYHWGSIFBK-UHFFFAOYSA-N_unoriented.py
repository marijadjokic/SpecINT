#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(4)
G.add_node(10)
G.add_node(11)
G.add_node(12)
G.add_node(14)
G.add_node(15)
G.add_node(16)
G.add_node(21)
G.add_node(27)
G.add_node(28)

G.add_edge(4,10)
G.add_edge(4,11)
G.add_edge(4,12)
G.add_edge(4,14)
G.add_edge(4,15)
G.add_edge(4,16)
G.add_edge(10,11)
G.add_edge(10,12)
G.add_edge(10,14)
G.add_edge(10,15)
G.add_edge(10,16)
G.add_edge(11,12)
G.add_edge(11,14)
G.add_edge(11,15)
G.add_edge(11,16)
G.add_edge(12,14)
G.add_edge(12,15)
G.add_edge(12,16)
G.add_edge(12,21)
G.add_edge(12,27)
G.add_edge(12,28)
G.add_edge(14,15)
G.add_edge(14,16)
G.add_edge(15,16)
G.add_edge(21,27)
G.add_edge(21,28)
G.add_edge(27,28)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/FNKQXYHWGSIFBK-UHFFFAOYSA-N_unoriented.png')
