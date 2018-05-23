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
G.add_node(7)
G.add_node(11)
G.add_node(12)
G.add_node(17)
G.add_node(20)
G.add_node(25)
G.add_node(29)
G.add_node(30)
G.add_node(31)
G.add_node(32)
G.add_node(33)
G.add_node(38)
G.add_node(42)
G.add_node(47)
G.add_node(50)

G.add_edge(0,1)
G.add_edge(0,2)
G.add_edge(0,7)
G.add_edge(0,11)
G.add_edge(0,12)
G.add_edge(0,17)
G.add_edge(0,20)
G.add_edge(0,25)
G.add_edge(0,29)
G.add_edge(0,30)
G.add_edge(1,2)
G.add_edge(1,7)
G.add_edge(1,11)
G.add_edge(1,12)
G.add_edge(1,17)
G.add_edge(1,20)
G.add_edge(1,25)
G.add_edge(1,29)
G.add_edge(1,30)
G.add_edge(2,7)
G.add_edge(2,11)
G.add_edge(2,12)
G.add_edge(2,17)
G.add_edge(2,20)
G.add_edge(2,25)
G.add_edge(2,29)
G.add_edge(2,30)
G.add_edge(32,33)
G.add_edge(32,38)
G.add_edge(32,42)
G.add_edge(32,12)
G.add_edge(32,47)
G.add_edge(32,50)
G.add_edge(32,31)
G.add_edge(33,38)
G.add_edge(33,42)
G.add_edge(33,12)
G.add_edge(33,47)
G.add_edge(33,50)
G.add_edge(33,31)
G.add_edge(38,42)
G.add_edge(38,12)
G.add_edge(38,47)
G.add_edge(38,50)
G.add_edge(38,31)
G.add_edge(7,11)
G.add_edge(7,12)
G.add_edge(7,17)
G.add_edge(7,20)
G.add_edge(7,25)
G.add_edge(7,29)
G.add_edge(7,30)
G.add_edge(42,12)
G.add_edge(42,47)
G.add_edge(42,50)
G.add_edge(42,31)
G.add_edge(11,12)
G.add_edge(11,17)
G.add_edge(11,20)
G.add_edge(11,25)
G.add_edge(11,29)
G.add_edge(11,30)
G.add_edge(12,47)
G.add_edge(12,17)
G.add_edge(12,50)
G.add_edge(12,20)
G.add_edge(12,25)
G.add_edge(12,29)
G.add_edge(12,30)
G.add_edge(12,31)
G.add_edge(47,50)
G.add_edge(47,31)
G.add_edge(17,20)
G.add_edge(17,25)
G.add_edge(17,29)
G.add_edge(17,30)
G.add_edge(50,31)
G.add_edge(20,25)
G.add_edge(20,29)
G.add_edge(20,30)
G.add_edge(25,29)
G.add_edge(25,30)
G.add_edge(29,30)

print G.nodes()
print nx.fiedler_vector(G,method='lobpcg')

pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/GSDSWSVVBLHKDQ-UHFFFAOYSA-N_unoriented.png')
