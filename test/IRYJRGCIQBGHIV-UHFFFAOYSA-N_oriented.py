#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(0)
G.add_node(1)
G.add_node(15)
G.add_node(16)
G.add_node(21)
G.add_node(22)
G.add_node(24)
G.add_node(26)
G.add_node(27)
G.add_node(28)
G.add_node(42)
G.add_node(48)
G.add_node(50)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(1,100)
G.add_edge(15,100)
G.add_edge(16,100)
G.add_edge(21,100)
G.add_edge(22,100)
G.add_edge(24,100)
G.add_edge(26,100)
G.add_edge(27,100)
G.add_edge(28,100)
G.add_edge(42,100)
G.add_edge(48,100)
G.add_edge(50,100)

G.add_edge(0,1)
G.add_edge(1,16)
G.add_edge(1,26)
G.add_edge(1,15)
G.add_edge(15,16)
G.add_edge(15,1)
G.add_edge(15,26)
G.add_edge(26,1)
G.add_edge(27,28)
G.add_edge(28,16)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/IRYJRGCIQBGHIV-UHFFFAOYSA-N_oriented.png')
