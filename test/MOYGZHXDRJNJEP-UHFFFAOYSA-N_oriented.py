#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(3)
G.add_node(4)
G.add_node(9)
G.add_node(10)
G.add_node(11)
G.add_node(13)
G.add_node(15)
G.add_node(16)
G.add_node(17)
G.add_node(21)
G.add_node(27)
G.add_node(28)
#fixed node 100
G.add_node(100)
G.add_edge(3,100)
G.add_edge(4,100)
G.add_edge(9,100)
G.add_edge(10,100)
G.add_edge(11,100)
G.add_edge(13,100)
G.add_edge(15,100)
G.add_edge(16,100)
G.add_edge(17,100)
G.add_edge(21,100)
G.add_edge(27,100)
G.add_edge(28,100)

G.add_edge(9,11)
G.add_edge(9,10)
G.add_edge(9,3)
G.add_edge(9,4)
G.add_edge(10,3)
G.add_edge(10,15)
G.add_edge(11,9)
G.add_edge(11,10)
G.add_edge(11,3)
G.add_edge(11,4)
G.add_edge(15,10)
G.add_edge(16,17)
G.add_edge(17,4)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/MOYGZHXDRJNJEP-UHFFFAOYSA-N_oriented.png')
