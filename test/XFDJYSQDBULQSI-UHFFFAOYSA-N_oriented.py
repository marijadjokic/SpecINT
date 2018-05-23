#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(0)
G.add_node(1)
G.add_node(3)
G.add_node(7)
G.add_node(10)
G.add_node(16)
G.add_node(17)
G.add_node(18)
G.add_node(19)
G.add_node(21)
G.add_node(25)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(1,100)
G.add_edge(3,100)
G.add_edge(7,100)
G.add_edge(10,100)
G.add_edge(16,100)
G.add_edge(17,100)
G.add_edge(18,100)
G.add_edge(19,100)
G.add_edge(21,100)
G.add_edge(25,100)

G.add_edge(0,1)
G.add_edge(1,3)
G.add_edge(1,7)
G.add_edge(1,10)
G.add_edge(1,16)
G.add_edge(1,17)
G.add_edge(7,16)
G.add_edge(7,1)
G.add_edge(7,10)
G.add_edge(7,3)
G.add_edge(7,17)
G.add_edge(16,1)
G.add_edge(18,19)
G.add_edge(19,10)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('D:/Users/Branko/Desktop/XFDJYSQDBULQSI-UHFFFAOYSA-N_oriented.png')
