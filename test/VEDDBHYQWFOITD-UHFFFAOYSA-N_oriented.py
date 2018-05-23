#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(32)
G.add_node(34)
G.add_node(35)
G.add_node(39)
G.add_node(10)
G.add_node(12)
G.add_node(14)
G.add_node(15)
G.add_node(19)
#fixed node 100
G.add_node(100)
G.add_edge(32,100)
G.add_edge(34,100)
G.add_edge(35,100)
G.add_edge(39,100)
G.add_edge(10,100)
G.add_edge(12,100)
G.add_edge(14,100)
G.add_edge(15,100)
G.add_edge(19,100)

G.add_edge(32,10)
G.add_edge(12,10)
G.add_edge(12,15)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/VEDDBHYQWFOITD-UHFFFAOYSA-N_oriented.png')
