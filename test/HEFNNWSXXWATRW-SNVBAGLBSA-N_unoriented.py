#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(0)
G.add_node(1)
G.add_node(2)
G.add_node(8)
G.add_node(12)
G.add_node(13)
G.add_node(14)
G.add_node(15)
G.add_node(25)

G.add_edge(0,8)
G.add_edge(0,1)
G.add_edge(0,2)
G.add_edge(0,12)
G.add_edge(1,8)
G.add_edge(1,2)
G.add_edge(1,12)
G.add_edge(2,12)
G.add_edge(2,8)
G.add_edge(8,12)
G.add_edge(8,13)
G.add_edge(8,14)
G.add_edge(8,15)
G.add_edge(8,25)
G.add_edge(13,25)
G.add_edge(13,14)
G.add_edge(13,15)
G.add_edge(14,25)
G.add_edge(14,15)
G.add_edge(15,25)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/HEFNNWSXXWATRW-SNVBAGLBSA-N_unoriented.png')
