#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(0)
G.add_node(1)
G.add_node(4)
G.add_node(5)
G.add_node(14)
G.add_node(17)
G.add_node(19)
G.add_node(20)
G.add_node(21)
G.add_node(22)
G.add_node(30)
G.add_node(32)
G.add_node(38)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(1,100)
G.add_edge(4,100)
G.add_edge(5,100)
G.add_edge(14,100)
G.add_edge(17,100)
G.add_edge(19,100)
G.add_edge(20,100)
G.add_edge(21,100)
G.add_edge(22,100)
G.add_edge(30,100)
G.add_edge(32,100)
G.add_edge(38,100)

G.add_edge(0,1)
G.add_edge(1,4)
G.add_edge(1,5)
G.add_edge(1,14)
G.add_edge(1,19)
G.add_edge(1,20)
G.add_edge(4,19)
G.add_edge(4,5)
G.add_edge(38,14)
G.add_edge(19,1)
G.add_edge(19,4)
G.add_edge(21,22)
G.add_edge(22,14)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/OGQICQVSFDPSEI-UHFFFAOYSA-N_oriented.png')
