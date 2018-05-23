#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(0)
G.add_node(3)
G.add_node(6)
G.add_node(12)
G.add_node(13)
G.add_node(15)
G.add_node(19)
G.add_node(20)
G.add_node(21)
G.add_node(22)
G.add_node(28)
G.add_node(34)
G.add_node(35)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(3,100)
G.add_edge(6,100)
G.add_edge(12,100)
G.add_edge(13,100)
G.add_edge(15,100)
G.add_edge(19,100)
G.add_edge(20,100)
G.add_edge(21,100)
G.add_edge(22,100)
G.add_edge(28,100)
G.add_edge(34,100)
G.add_edge(35,100)

G.add_edge(13,19)
G.add_edge(13,3)
G.add_edge(13,20)
G.add_edge(13,15)
G.add_edge(15,3)
G.add_edge(15,19)
G.add_edge(15,20)
G.add_edge(15,13)
G.add_edge(19,15)
G.add_edge(21,22)
G.add_edge(22,3)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/UNBRKDKAWYKMIV-QWQRMKEZSA-N_oriented.png')
