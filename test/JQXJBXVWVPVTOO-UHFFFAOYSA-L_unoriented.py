#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(1)
G.add_node(3)
G.add_node(5)

G.add_edge(1,3)
G.add_edge(1,5)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/JQXJBXVWVPVTOO-UHFFFAOYSA-L_unoriented.png')
