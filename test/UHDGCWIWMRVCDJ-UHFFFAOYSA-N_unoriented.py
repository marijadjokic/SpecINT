#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.Graph()

G.add_node('0')
G.add_node('1')
G.add_node('2')
G.add_node('9')
G.add_node('10')
G.add_node('11')

G.add_edge('0','1')
G.add_edge('0','2')
G.add_edge('0','9')
G.add_edge('1','9')
G.add_edge('1','2')
G.add_edge('2','9')
G.add_edge('9','11')
G.add_edge('9','10')
G.add_edge('10','11')

 print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/UHDGCWIWMRVCDJ-UHFFFAOYSA-N_unoriented.png')
