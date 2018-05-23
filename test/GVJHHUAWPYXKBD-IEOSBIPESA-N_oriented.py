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
G.add_node('10')
G.add_node('18')
G.add_node('23')
G.add_node('24')
G.add_node('26')
G.add_node('27')
G.add_node('28')
G.add_node('30')
G.add_node('37')
G.add_node('45')
G.add_node('50')
#fixed node 100
G.add_node('100')
G.add_edge('0','100')
G.add_edge('1','100')
G.add_edge('3','100')
G.add_edge('10','100')
G.add_edge('18','100')
G.add_edge('23','100')
G.add_edge('24','100')
G.add_edge('26','100')
G.add_edge('27','100')
G.add_edge('28','100')
G.add_edge('30','100')
G.add_edge('37','100')
G.add_edge('45','100')
G.add_edge('50','100')

G.add_edge('0','1')
G.add_edge('1','24')
G.add_edge('1','18')
G.add_edge('1','3')
G.add_edge('1','26')
G.add_edge('3','18')
G.add_edge('45','26')
G.add_edge('24','1')
G.add_edge('24','18')
G.add_edge('24','3')
G.add_edge('24','26')
G.add_edge('27','28')
G.add_edge('28','26')
G.add_edge('30','26')

 print 'Degree: '
print G.degree(G.nodes())
print 'PG:'
print nx.pagerank(G)
pos=nx.circular_layout(G,dim=50, scale=100)
plt.clf()
nx.draw(G,with_labels=True)
plt.savefig('C:/Users/Branko/Desktop/GVJHHUAWPYXKBD-IEOSBIPESA-N_oriented.png')
