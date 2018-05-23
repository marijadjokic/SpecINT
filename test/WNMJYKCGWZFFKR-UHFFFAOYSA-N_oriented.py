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
G.add_node(8)
G.add_node(12)
G.add_node(19)
G.add_node(22)
G.add_node(23)
G.add_node(24)
G.add_node(25)
G.add_node(26)
G.add_node(27)
G.add_node(36)
G.add_node(46)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(1,100)
G.add_edge(2,100)
G.add_edge(8,100)
G.add_edge(12,100)
G.add_edge(19,100)
G.add_edge(22,100)
G.add_edge(23,100)
G.add_edge(24,100)
G.add_edge(25,100)
G.add_edge(26,100)
G.add_edge(27,100)
G.add_edge(36,100)
G.add_edge(46,100)

G.add_edge(0,1)
G.add_edge(1,2)
G.add_edge(1,19)
G.add_edge(1,22)
G.add_edge(1,23)
G.add_edge(1,24)
G.add_edge(22,24)
G.add_edge(22,1)
G.add_edge(22,2)
G.add_edge(22,19)
G.add_edge(22,23)
G.add_edge(23,1)
G.add_edge(25,26)
G.add_edge(26,19)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/WNMJYKCGWZFFKR-UHFFFAOYSA-N_oriented.png')
