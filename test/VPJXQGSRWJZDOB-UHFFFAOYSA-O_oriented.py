#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node(1)
G.add_node(4)
G.add_node(6)
G.add_node(13)
G.add_node(15)
G.add_node(17)
G.add_node(18)
G.add_node(19)
G.add_node(20)
G.add_node(22)
G.add_node(25)
G.add_node(27)
G.add_node(34)
#fixed node 100
G.add_node(100)
G.add_edge(1,100)
G.add_edge(4,100)
G.add_edge(6,100)
G.add_edge(13,100)
G.add_edge(15,100)
G.add_edge(17,100)
G.add_edge(18,100)
G.add_edge(19,100)
G.add_edge(20,100)
G.add_edge(22,100)
G.add_edge(25,100)
G.add_edge(27,100)
G.add_edge(34,100)

G.add_edge(13,1)
G.add_edge(13,15)
G.add_edge(13,17)
G.add_edge(13,18)
G.add_edge(13,19)
G.add_edge(13,20)
G.add_edge(18,1)
G.add_edge(18,13)
G.add_edge(18,17)
G.add_edge(18,19)
G.add_edge(18,20)
G.add_edge(19,18)
G.add_edge(22,15)

print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/VPJXQGSRWJZDOB-UHFFFAOYSA-O_oriented.png')
