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
G.add_node(4)
#fixed node 100
G.add_node(100)
G.add_edge(0,100)
G.add_edge(1,100)
G.add_edge(2,100)
G.add_edge(3,100)
G.add_edge(4,100)

G.add_edge(0,1)
G.add_edge(0,2)
G.add_edge(3,2)
G.add_edge(4,2)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/YKFRUALXTUDSBW-UHFFFAOYSA-N_oriented.png')
