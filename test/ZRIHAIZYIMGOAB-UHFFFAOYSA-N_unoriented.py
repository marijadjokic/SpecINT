#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(1)
G.add_node(5)
G.add_node(7)
G.add_node(8)
G.add_node(9)
G.add_node(14)
G.add_node(16)
G.add_node(17)
G.add_node(23)
G.add_node(25)
G.add_node(26)
G.add_node(27)
G.add_node(32)

G.add_edge(32,1)
G.add_edge(32,23)
G.add_edge(32,25)
G.add_edge(32,26)
G.add_edge(32,27)
G.add_edge(1,5)
G.add_edge(1,7)
G.add_edge(1,8)
G.add_edge(1,9)
G.add_edge(1,14)
G.add_edge(1,16)
G.add_edge(1,17)
G.add_edge(1,23)
G.add_edge(1,25)
G.add_edge(1,26)
G.add_edge(1,27)
G.add_edge(5,7)
G.add_edge(5,8)
G.add_edge(5,9)
G.add_edge(5,14)
G.add_edge(5,16)
G.add_edge(5,17)
G.add_edge(7,8)
G.add_edge(7,9)
G.add_edge(7,14)
G.add_edge(7,16)
G.add_edge(7,17)
G.add_edge(8,9)
G.add_edge(8,14)
G.add_edge(8,16)
G.add_edge(8,17)
G.add_edge(9,14)
G.add_edge(9,16)
G.add_edge(9,17)
G.add_edge(14,16)
G.add_edge(14,17)
G.add_edge(16,17)
G.add_edge(23,25)
G.add_edge(23,26)
G.add_edge(23,27)
G.add_edge(25,26)
G.add_edge(25,27)
G.add_edge(26,27)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/ZRIHAIZYIMGOAB-UHFFFAOYSA-N_unoriented.png')
