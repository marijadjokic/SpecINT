#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node(0)
G.add_node(1)
G.add_node(6)
G.add_node(7)
G.add_node(13)

G.add_edge(0,1)
G.add_edge(0,6)
G.add_edge(1,13)
G.add_edge(1,6)
G.add_edge(1,7)
G.add_edge(7,13)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/JBVNBBXAMBZTMQ-CEGNMAFCSA-N_unoriented.png')
