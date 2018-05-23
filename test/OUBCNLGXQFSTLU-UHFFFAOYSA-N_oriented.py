#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(0)
G.add_node(1)
G.add_node(2)
G.add_node(7)
G.add_node(9)
G.add_node(17)
G.add_node(18)
G.add_node(22)
G.add_node(23)
G.add_node(24)
G.add_node(29)
G.add_node(31)
G.add_node(39)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(1,100)
G.add_edge(2,100)
G.add_edge(7,100)
G.add_edge(9,100)
G.add_edge(17,100)
G.add_edge(18,100)
G.add_edge(22,100)
G.add_edge(23,100)
G.add_edge(24,100)
G.add_edge(29,100)
G.add_edge(31,100)
G.add_edge(39,100)

G.add_edge(0,1)
G.add_edge(1,18)
G.add_edge(1,7)
G.add_edge(1,22)
G.add_edge(1,17)
G.add_edge(17,1)
G.add_edge(17,18)
G.add_edge(17,22)
G.add_edge(17,7)
G.add_edge(23,24)
G.add_edge(24,18)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/OUBCNLGXQFSTLU-UHFFFAOYSA-N_oriented.png')
