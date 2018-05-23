#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(1)
G.add_node(3)
G.add_node(5)
G.add_node(6)
G.add_node(7)
G.add_node(9)
G.add_node(13)
G.add_node(14)

G.add_edge(1,3)
G.add_edge(1,5)
G.add_edge(1,6)
G.add_edge(1,7)
G.add_edge(3,5)
G.add_edge(3,6)
G.add_edge(3,7)
G.add_edge(3,9)
G.add_edge(3,13)
G.add_edge(3,14)
G.add_edge(5,6)
G.add_edge(5,7)
G.add_edge(6,7)
G.add_edge(9,13)
G.add_edge(9,14)
G.add_edge(13,14)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/HKSZLNNOFSGOKW-ZGQXJOJZSA-N_unoriented.png')
