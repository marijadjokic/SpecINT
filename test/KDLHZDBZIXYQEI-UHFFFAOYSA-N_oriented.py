#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(6)
G.add_node(23)
G.add_node(11)
#fixed node 100
G.add_node(100)
G.add_edge(6,100)
G.add_edge(23,100)
G.add_edge(11,100)

G.add_edge(23,6)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/KDLHZDBZIXYQEI-UHFFFAOYSA-N_oriented.png')
