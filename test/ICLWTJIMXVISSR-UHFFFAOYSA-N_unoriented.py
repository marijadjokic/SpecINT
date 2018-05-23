#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(0)
G.add_node(3)
G.add_node(9)
G.add_node(11)
G.add_node(15)
G.add_node(21)

G.add_edge(0,11)
G.add_edge(0,9)
G.add_edge(0,3)
G.add_edge(0,21)
G.add_edge(0,15)
G.add_edge(3,9)
G.add_edge(3,11)
G.add_edge(9,11)
G.add_edge(15,21)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/ICLWTJIMXVISSR-UHFFFAOYSA-N_unoriented.png')
