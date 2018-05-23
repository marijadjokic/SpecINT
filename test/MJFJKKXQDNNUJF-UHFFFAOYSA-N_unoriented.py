#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(0)
G.add_node(1)
G.add_node(9)
G.add_node(10)
G.add_node(16)
G.add_node(17)
G.add_node(20)
G.add_node(21)
G.add_node(22)
G.add_node(23)
G.add_node(32)
G.add_node(38)
G.add_node(41)

G.add_edge(0,1)
G.add_edge(0,9)
G.add_edge(0,10)
G.add_edge(0,16)
G.add_edge(0,17)
G.add_edge(0,20)
G.add_edge(0,21)
G.add_edge(1,9)
G.add_edge(1,10)
G.add_edge(1,16)
G.add_edge(1,17)
G.add_edge(1,20)
G.add_edge(1,21)
G.add_edge(38,32)
G.add_edge(38,41)
G.add_edge(38,23)
G.add_edge(38,22)
G.add_edge(38,9)
G.add_edge(32,41)
G.add_edge(32,23)
G.add_edge(32,22)
G.add_edge(32,9)
G.add_edge(9,41)
G.add_edge(9,10)
G.add_edge(9,16)
G.add_edge(9,17)
G.add_edge(9,20)
G.add_edge(9,21)
G.add_edge(9,22)
G.add_edge(9,23)
G.add_edge(10,16)
G.add_edge(10,17)
G.add_edge(10,20)
G.add_edge(10,21)
G.add_edge(23,41)
G.add_edge(23,22)
G.add_edge(16,17)
G.add_edge(16,20)
G.add_edge(16,21)
G.add_edge(17,20)
G.add_edge(17,21)
G.add_edge(20,21)
G.add_edge(22,41)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/MJFJKKXQDNNUJF-UHFFFAOYSA-N_unoriented.png')
