#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(1)
G.add_node(10)
G.add_node(11)
G.add_node(13)
G.add_node(14)
G.add_node(15)
G.add_node(17)
G.add_node(26)

G.add_edge(1,10)
G.add_edge(1,11)
G.add_edge(1,13)
G.add_edge(1,14)
G.add_edge(1,15)
G.add_edge(10,11)
G.add_edge(10,13)
G.add_edge(10,14)
G.add_edge(10,15)
G.add_edge(11,13)
G.add_edge(11,14)
G.add_edge(11,15)
G.add_edge(11,17)
G.add_edge(11,26)
G.add_edge(13,14)
G.add_edge(13,15)
G.add_edge(14,15)
G.add_edge(17,26)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/URXWVWVPMJSAJD-KOORYGTMSA-N_unoriented.png')
