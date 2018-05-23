#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(0)
G.add_node(4)
G.add_node(6)
G.add_node(9)
G.add_node(11)
G.add_node(16)
G.add_node(19)
G.add_node(24)
G.add_node(25)
G.add_node(26)
G.add_node(32)
G.add_node(35)
G.add_node(37)
G.add_node(42)
G.add_node(45)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(4,100)
G.add_edge(6,100)
G.add_edge(9,100)
G.add_edge(11,100)
G.add_edge(16,100)
G.add_edge(19,100)
G.add_edge(24,100)
G.add_edge(25,100)
G.add_edge(26,100)
G.add_edge(32,100)
G.add_edge(35,100)
G.add_edge(37,100)
G.add_edge(42,100)
G.add_edge(45,100)

G.add_edge(37,4)
G.add_edge(11,0)
G.add_edge(11,4)
G.add_edge(11,6)
G.add_edge(11,16)
G.add_edge(11,19)
G.add_edge(11,24)
G.add_edge(11,25)
G.add_edge(45,4)
G.add_edge(16,0)
G.add_edge(16,4)
G.add_edge(16,6)
G.add_edge(16,11)
G.add_edge(16,19)
G.add_edge(16,24)
G.add_edge(16,25)
G.add_edge(19,24)
G.add_edge(19,6)
G.add_edge(24,11)
G.add_edge(24,19)
G.add_edge(26,4)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/ROBVIMPUHSLWNV-UHFFFAOYSA-N_oriented.png')
