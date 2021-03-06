#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(0)
G.add_node(3)
G.add_node(8)
G.add_node(10)
G.add_node(12)
G.add_node(13)
G.add_node(14)
G.add_node(15)
G.add_node(16)
G.add_node(17)
G.add_node(20)
G.add_node(27)

G.add_edge(0,3)
G.add_edge(0,8)
G.add_edge(0,10)
G.add_edge(0,12)
G.add_edge(0,13)
G.add_edge(0,14)
G.add_edge(0,15)
G.add_edge(0,16)
G.add_edge(3,8)
G.add_edge(3,10)
G.add_edge(3,12)
G.add_edge(3,13)
G.add_edge(3,14)
G.add_edge(3,15)
G.add_edge(3,16)
G.add_edge(3,17)
G.add_edge(3,20)
G.add_edge(3,27)
G.add_edge(8,10)
G.add_edge(8,12)
G.add_edge(8,13)
G.add_edge(8,14)
G.add_edge(8,15)
G.add_edge(8,16)
G.add_edge(10,12)
G.add_edge(10,13)
G.add_edge(10,14)
G.add_edge(10,15)
G.add_edge(10,16)
G.add_edge(12,13)
G.add_edge(12,14)
G.add_edge(12,15)
G.add_edge(12,16)
G.add_edge(13,14)
G.add_edge(13,15)
G.add_edge(13,16)
G.add_edge(14,15)
G.add_edge(14,16)
G.add_edge(15,16)
G.add_edge(17,27)
G.add_edge(17,20)
G.add_edge(20,27)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/WZPBZJONDBGPKJ-IYZXUIDESA-N_unoriented.png')
