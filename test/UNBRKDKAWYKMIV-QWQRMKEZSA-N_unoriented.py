#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

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

G.add_edge(0,3)
G.add_edge(0,6)
G.add_edge(0,12)
G.add_edge(0,13)
G.add_edge(0,15)
G.add_edge(0,19)
G.add_edge(0,20)
G.add_edge(35,34)
G.add_edge(35,3)
G.add_edge(35,21)
G.add_edge(35,22)
G.add_edge(35,28)
G.add_edge(34,3)
G.add_edge(34,21)
G.add_edge(34,22)
G.add_edge(34,28)
G.add_edge(3,6)
G.add_edge(3,12)
G.add_edge(3,13)
G.add_edge(3,15)
G.add_edge(3,19)
G.add_edge(3,20)
G.add_edge(3,21)
G.add_edge(3,22)
G.add_edge(3,28)
G.add_edge(6,12)
G.add_edge(6,13)
G.add_edge(6,15)
G.add_edge(6,19)
G.add_edge(6,20)
G.add_edge(12,13)
G.add_edge(12,15)
G.add_edge(12,19)
G.add_edge(12,20)
G.add_edge(13,15)
G.add_edge(13,19)
G.add_edge(13,20)
G.add_edge(15,19)
G.add_edge(15,20)
G.add_edge(19,20)
G.add_edge(21,22)
G.add_edge(21,28)
G.add_edge(22,28)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/UNBRKDKAWYKMIV-QWQRMKEZSA-N_unoriented.png')