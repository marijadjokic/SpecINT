#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

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
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(1,100)
G.add_edge(9,100)
G.add_edge(10,100)
G.add_edge(16,100)
G.add_edge(17,100)
G.add_edge(20,100)
G.add_edge(21,100)
G.add_edge(22,100)
G.add_edge(23,100)
G.add_edge(32,100)
G.add_edge(38,100)
G.add_edge(41,100)

G.add_edge(0,1)
G.add_edge(1,17)
G.add_edge(1,21)
G.add_edge(1,20)
G.add_edge(1,9)
G.add_edge(23,9)
G.add_edge(20,1)
G.add_edge(20,9)
G.add_edge(20,21)
G.add_edge(20,17)
G.add_edge(21,1)
G.add_edge(22,23)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/MJFJKKXQDNNUJF-UHFFFAOYSA-N_oriented.png')
