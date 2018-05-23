#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(8)
G.add_node(10)
G.add_node(11)
G.add_node(12)
G.add_node(13)
G.add_node(14)
G.add_node(15)
G.add_node(26)
G.add_node(27)

G.add_edge(8,10)
G.add_edge(8,11)
G.add_edge(8,12)
G.add_edge(8,13)
G.add_edge(8,14)
G.add_edge(8,15)
G.add_edge(8,26)
G.add_edge(8,27)
G.add_edge(10,11)
G.add_edge(10,12)
G.add_edge(10,13)
G.add_edge(10,14)
G.add_edge(10,15)
G.add_edge(11,12)
G.add_edge(11,13)
G.add_edge(11,14)
G.add_edge(11,15)
G.add_edge(12,13)
G.add_edge(12,14)
G.add_edge(12,15)
G.add_edge(13,14)
G.add_edge(13,15)
G.add_edge(14,15)
G.add_edge(26,27)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/VXIMPSPISRVBPZ-NWUMPJBXSA-N_unoriented.png')
