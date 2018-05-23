#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(8)
G.add_node(10)
G.add_node(11)
G.add_node(12)
G.add_node(13)
G.add_node(14)
G.add_node(15)
G.add_node(26)
G.add_node(27)
#fixed node 100
G.add_node(100)
G.add_edge(8,100)
G.add_edge(10,100)
G.add_edge(11,100)
G.add_edge(12,100)
G.add_edge(13,100)
G.add_edge(14,100)
G.add_edge(15,100)
G.add_edge(26,100)
G.add_edge(27,100)

G.add_edge(10,8)
G.add_edge(10,11)
G.add_edge(10,12)
G.add_edge(10,13)
G.add_edge(10,14)
G.add_edge(11,12)
G.add_edge(11,15)
G.add_edge(13,10)
G.add_edge(13,14)
G.add_edge(15,11)
G.add_edge(15,13)
G.add_edge(27,8)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/VXIMPSPISRVBPZ-NWUMPJBXSA-N_oriented.png')
