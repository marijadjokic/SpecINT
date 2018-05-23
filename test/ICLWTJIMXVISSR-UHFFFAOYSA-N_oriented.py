#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(0)
G.add_node(3)
G.add_node(9)
G.add_node(11)
G.add_node(15)
G.add_node(21)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(3,100)
G.add_edge(9,100)
G.add_edge(11,100)
G.add_edge(15,100)
G.add_edge(21,100)

G.add_edge(9,0)
G.add_edge(9,11)
G.add_edge(11,9)
G.add_edge(15,0)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/ICLWTJIMXVISSR-UHFFFAOYSA-N_oriented.png')
