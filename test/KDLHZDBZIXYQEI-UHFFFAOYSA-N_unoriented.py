#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(6)
G.add_node(23)
G.add_node(11)

G.add_edge(11,6)
G.add_edge(6,23)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/KDLHZDBZIXYQEI-UHFFFAOYSA-N_unoriented.png')
