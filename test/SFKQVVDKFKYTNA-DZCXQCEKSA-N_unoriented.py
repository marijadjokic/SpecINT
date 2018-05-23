#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(0)
G.add_node(1)
G.add_node(4)
G.add_node(7)
G.add_node(8)
G.add_node(9)
G.add_node(10)
G.add_node(15)

G.add_edge(0,8)
G.add_edge(0,1)
G.add_edge(0,4)
G.add_edge(0,7)
G.add_edge(1,8)
G.add_edge(1,4)
G.add_edge(1,7)
G.add_edge(4,8)
G.add_edge(4,7)
G.add_edge(7,8)
G.add_edge(8,9)
G.add_edge(8,10)
G.add_edge(8,15)
G.add_edge(9,10)
G.add_edge(9,15)
G.add_edge(10,15)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/SFKQVVDKFKYTNA-DZCXQCEKSA-N_unoriented.png')
