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
G.add_node(3)
G.add_node(7)
G.add_node(9)
G.add_node(12)
G.add_node(13)
G.add_node(21)
G.add_node(22)
G.add_node(23)
G.add_node(25)
G.add_node(26)
G.add_node(32)
G.add_node(40)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(1,100)
G.add_edge(2,100)
G.add_edge(3,100)
G.add_edge(7,100)
G.add_edge(9,100)
G.add_edge(12,100)
G.add_edge(13,100)
G.add_edge(21,100)
G.add_edge(22,100)
G.add_edge(23,100)
G.add_edge(25,100)
G.add_edge(26,100)
G.add_edge(32,100)
G.add_edge(40,100)

G.add_edge(0,1)
G.add_edge(1,3)
G.add_edge(1,9)
G.add_edge(1,7)
G.add_edge(3,1)
G.add_edge(3,9)
G.add_edge(3,7)
G.add_edge(7,21)
G.add_edge(7,13)
G.add_edge(21,1)
G.add_edge(21,7)
G.add_edge(23,9)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/UFLGIAIHIAPJJC-UHFFFAOYSA-N_oriented.png')
