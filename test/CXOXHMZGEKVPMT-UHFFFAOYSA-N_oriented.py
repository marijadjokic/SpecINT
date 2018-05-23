#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(0)
G.add_node(1)
G.add_node(36)
G.add_node(39)
G.add_node(11)
G.add_node(15)
G.add_node(17)
G.add_node(18)
G.add_node(20)
G.add_node(21)
G.add_node(22)
G.add_node(23)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(1,100)
G.add_edge(36,100)
G.add_edge(39,100)
G.add_edge(11,100)
G.add_edge(15,100)
G.add_edge(17,100)
G.add_edge(18,100)
G.add_edge(20,100)
G.add_edge(21,100)
G.add_edge(22,100)
G.add_edge(23,100)

G.add_edge(0,1)
G.add_edge(1,15)
G.add_edge(1,17)
G.add_edge(1,18)
G.add_edge(1,20)
G.add_edge(1,21)
G.add_edge(15,1)
G.add_edge(15,18)
G.add_edge(15,20)
G.add_edge(15,17)
G.add_edge(20,1)
G.add_edge(22,23)
G.add_edge(23,17)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/CXOXHMZGEKVPMT-UHFFFAOYSA-N_oriented.png')
