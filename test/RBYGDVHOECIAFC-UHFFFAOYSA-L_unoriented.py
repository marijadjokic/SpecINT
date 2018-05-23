#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(5)
G.add_node(8)
G.add_node(9)
G.add_node(15)
G.add_node(18)

G.add_edge(5,8)
G.add_edge(5,9)
G.add_edge(5,18)
G.add_edge(5,15)
G.add_edge(8,9)
G.add_edge(15,18)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/RBYGDVHOECIAFC-UHFFFAOYSA-L_unoriented.png')
