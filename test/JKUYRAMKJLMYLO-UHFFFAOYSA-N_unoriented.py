#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(3)
G.add_node(8)
G.add_node(12)

G.add_edge(8,3)
G.add_edge(8,12)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/JKUYRAMKJLMYLO-UHFFFAOYSA-N_unoriented.png')
