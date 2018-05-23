#!/usr/bin/env python
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G = nx.DiGraph()

G.add_node('0')
G.add_node('1')
G.add_node('3')
G.add_node('4')
G.add_node('7')
G.add_node('12')
G.add_node('14')
G.add_node('18')
G.add_node('19')
G.add_node('20')
G.add_node('25')
G.add_node('27')
G.add_node('31')
G.add_node('36')
#fixed node 100
G.add_node('100')
G.add_edge('0','100')
G.add_edge('1','100')
G.add_edge('3','100')
G.add_edge('4','100')
G.add_edge('7','100')
G.add_edge('12','100')
G.add_edge('14','100')
G.add_edge('18','100')
G.add_edge('19','100')
G.add_edge('20','100')
G.add_edge('25','100')
G.add_edge('27','100')
G.add_edge('31','100')
G.add_edge('36','100')

G.add_edge('0','1')
G.add_edge('1','18')
G.add_edge('4','1')
G.add_edge('4','3')
G.add_edge('4','7')
G.add_edge('4','12')
G.add_edge('4','18')
G.add_edge('19','20')
G.add_edge('20','3')
G.add_edge('25','3')
G.add_edge('27','3')
G.add_edge('36','3')
G.add_edge('31','3')

 print 'Degree: '
print G.degree(G.nodes())
print '
PG:
 'print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/IKMDFBPHZNJCSN-UHFFFAOYSA-N_oriented.png')
