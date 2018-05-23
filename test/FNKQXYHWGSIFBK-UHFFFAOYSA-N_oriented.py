#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

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
#fixed node 100
G.add_node(100)
G.add_edge(4,100)
G.add_edge(10,100)
G.add_edge(11,100)
G.add_edge(12,100)
G.add_edge(14,100)
G.add_edge(15,100)
G.add_edge(16,100)
G.add_edge(21,100)
G.add_edge(27,100)
G.add_edge(28,100)

G.add_edge(10,16)
G.add_edge(10,12)
G.add_edge(10,4)
G.add_edge(10,14)
G.add_edge(10,15)
G.add_edge(14,4)
G.add_edge(14,10)
G.add_edge(14,12)
G.add_edge(14,15)
G.add_edge(14,16)
G.add_edge(15,4)
G.add_edge(21,12)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/FNKQXYHWGSIFBK-UHFFFAOYSA-N_oriented.png')
