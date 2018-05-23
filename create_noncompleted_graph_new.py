#!/usr/bin/env python
import rdflib
import urlparse
import networkx as nx
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from StringIO import StringIO    
import pycurl
import json
import sys
import os
from itertools import combinations
import re
from array import *
import chemspipy
from SPARQLWrapper import SPARQLWrapper, JSON
from datetime import datetime
import time
import random
from random import randint
import numpy.linalg
import operator
from operator import itemgetter


#function create_python_file
def create_python_file(nodes,edges,inchikey):
    file = open("/var/www/specint.org/public_html/specint/test/"+inchikey+"_oriented.py", "w");
    #add preambula
    file.write("#!/usr/bin/env python\n"+
    "import networkx as nx\n" +
    "import numpy as np\n" + 
    "import matplotlib\n" +
    "matplotlib.use('Agg')\n" +
    "import matplotlib.pyplot as plt\n\n")

    file.write("G = nx.DiGraph()\n\n")
   
    #add nodes
    for n in nodes:
     #print n
     file.write("G.add_node("+ str(n) + ")\n")
    
	#fixed node 100
    file.write("#fixed node 100\n")
    file.write("G.add_node(100)\n")
    for n in nodes:
     #print n
     file.write("G.add_edge(" + str(n) + ",100)\n")
	
    #add edges
    file.write("\n")	
    for u,v in edges:
     file.write("G.add_edge(" + str(u) + "," + str(v) + ")\n")
    file.write("\n")	
    file.write("print 'Degree: '\n")
    file.write("print G.degree(G.nodes())\n") 
    file.write("print 'PG:'\n")
    file.write("print nx.pagerank(G)")

    #draw graph
    file.write("\n")
    file.write("pos=nx.circular_layout(G,dim=50, scale=100)\n")
    file.write("plt.clf()\n")
    file.write("nx.draw(G,with_labels=True)\n")
    file.write("plt.savefig('C:/Users/Branko/Desktop/"+inchikey+"_oriented.png')\n")

    file.close()


##########Bio2RDF and Chem2Bio2RDF################################################



#***************Using Bio2RDF********************************************************

def non_completed_graph(inchikey,connection_node):

    #graphs generated using PIBAS dataset, UniChem and Bio2RDF
    G6 = nx.DiGraph()
    G8 = nx.DiGraph()
    #print connection_node

    #we start from PIBAS local ontology and  for a given InChiKey we found compound acronym 

    sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

 

    url = 'https://www.ebi.ac.uk/unichem/rest/inchikey/'+inchikey
    #print url
    j=0
    storage = StringIO()
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.WRITEFUNCTION, storage.write)
    c.perform()
    c.close()
    content = storage.getvalue()
    unichem_content = json.loads(content)
    for rs in unichem_content: 
     src_compound_id= rs['src_compound_id']
     url = 'https://www.ebi.ac.uk/unichem/rest/sources/'+rs['src_id']
     storage = StringIO()
     c = pycurl.Curl()
     c.setopt(c.URL, url)
     c.setopt(c.WRITEFUNCTION, storage.write)
     c.perform()
     c.close()
     content = storage.getvalue()
     unichem_src_name = json.loads(content)
     for rs in unichem_src_name:
      name=rs['name']
      G6.add_node(j,name_label=name+'/'+src_compound_id)
      j=j+1



    bio2rdf_dataset=['bindingdb','pubchem','pharmgkb','chebi','kegg_ligand','pdb','drugbank','chembl','pibas','ndc']

    num1=G6.number_of_nodes()
    nodes1=G6.nodes()
    for x in range(0, num1):
     if(((G6.node[x]['name_label']).split('/')[0] in bio2rdf_dataset) and ("pibas" not in (G6.node[x]['name_label']).split('/')[0])):
      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
      if((G6.node[x]['name_label']).split('/')[0]=='kegg_ligand'):
       sparql.setQuery(
       """  
        PREFIX bindingdb: <http://bio2rdf.org/bindingdb:>
        PREFIX pubchem: <http://bio2rdf.org/pubchem:>
        PREFIX pharmgkb: <http://bio2rdf.org/pharmgkb:>
        PREFIX chebi: <http://bio2rdf.org/chebi:>
        PREFIX kegg_ligand: <http://bio2rdf.org/kegg:>
        PREFIX pdb: <http://bio2rdf.org/pdb:>
        PREFIX drugbank: <http://bio2rdf.org/drugbank:>
        PREFIX chembl: <http://bio2rdf.org/chembl:> 
        PREFIX ndc: <http://bio2rdf.org/ndc:> 
        select ?p ?o
        WHERE
        {

         SERVICE SILENT<http://kegg.bio2rdf.org/sparql> { 
         OPTIONAL{
          %s:%s ?p ?o.
         FILTER(CONTAINS(str(?p),"http://bio2rdf.org/kegg_vocabulary:x-drugbank") || CONTAINS(str(?p),"http://bio2rdf.org/kegg_vocabulary:x-pubchem.compound") || CONTAINS(str(?p),"http://bio2rdf.org/kegg_vocabulary:x-kegg")  || CONTAINS(str(?p),"http://bio2rdf.org/kegg_vocabulary:x-chembl") || CONTAINS(str(?p),"http://bio2rdf.org/kegg_vocabulary:x-chebi") || CONTAINS(str(?p),"http://bio2rdf.org/kegg_vocabulary:same-as")).
        }
        }

        }

       """ % (((G6.node[x]['name_label']).split('/')[0]),((G6.node[x]['name_label']).split('/')[1])))

      else:
       sparql.setQuery(
       """  
        PREFIX bindingdb: <http://bio2rdf.org/bindingdb:>
        PREFIX pubchem: <http://bio2rdf.org/pubchem:>
        PREFIX pharmgkb: <http://bio2rdf.org/pharmgkb:>
        PREFIX chebi: <http://bio2rdf.org/chebi:>
        PREFIX kegg_ligand: <http://bio2rdf.org/kegg:>
        PREFIX pdb: <http://bio2rdf.org/pdb:>
        PREFIX drugbank: <http://bio2rdf.org/drugbank:>
        PREFIX chembl: <http://bio2rdf.org/chembl:>
        PREFIX ndc: <http://bio2rdf.org/ndc:> 
        select ?p ?o
        WHERE
        {

         SERVICE SILENT<http://%s.bio2rdf.org/sparql> { 
         OPTIONAL{
          %s:%s ?p ?o.
         FILTER(CONTAINS(str(?p),"http://bio2rdf.org/%s_vocabulary:x-drugbank") || CONTAINS(str(?p),"http://bio2rdf.org/%s_vocabulary:x-pubchem.compound") || CONTAINS(str(?p),"http://bio2rdf.org/%s_vocabulary:x-kegg")  || CONTAINS(str(?p),"http://bio2rdf.org/%s_vocabulary:x-chembl") || CONTAINS(str(?p),"http://bio2rdf.org/%s_vocabulary:x-chebi") || CONTAINS(str(?p),"http://bio2rdf.org/%s_vocabulary:x-ndc") || CONTAINS(str(?p),"http://bio2rdf.org/%s_vocabulary:same-as")).
        }
        }

        }

       """ % (((G6.node[x]['name_label']).split('/')[0]),((G6.node[x]['name_label']).split('/')[0]),((G6.node[x]['name_label']).split('/')[1]),((G6.node[x]['name_label']).split('/')[0]),((G6.node[x]['name_label']).split('/')[0]),((G6.node[x]['name_label']).split('/')[0]),((G6.node[x]['name_label']).split('/')[0]),((G6.node[x]['name_label']).split('/')[0]),((G6.node[x]['name_label']).split('/')[0]),((G6.node[x]['name_label']).split('/')[0])))

      sparql.setReturnFormat(JSON)
      final_results1 = sparql.query().convert()
      try:
       v=G6.number_of_nodes()
       for result1 in final_results1["results"]["bindings"]:
        if(result1["o"]["value"]!=""):  
         node_for_add=result1["o"]["value"] 
         node_of_second_level=node_for_add.split("/")[-1]
         if(node_of_second_level.split(':')[0]=='kegg'):
          node_for_add='kegg_ligand/'+node_of_second_level.split(':')[1]
         else:
          node_for_add=node_of_second_level.split(':')[0]+'/'+node_of_second_level.split(':')[1]
          #print node_for_add
         compare_node=[]
         for x in range(0,G6.number_of_nodes()):
           compare_node.append(str((G6.node[x]['name_label']).split('/')[0])+'/'+str((G6.node[x]['name_label']).split('/')[1]))

         if(node_for_add not in compare_node):
          G6.add_node(v,name_label=node_for_add,predicate=result1["p"]["value"])
          v=v+1
          #print v
      except:     
        continue



    nodes=G6.nodes()

    #print nodes

    edges = combinations(nodes, 2)
    G6.add_nodes_from(nodes)
    G6.add_edges_from(edges)

    num1=G6.number_of_nodes()
    nodes1=G6.nodes()
    left_graph_for_remove=[]
    for x in range(0, num1):
     if(((G6.node[x]['name_label']).split('/')[0] not in bio2rdf_dataset)):
      G6.remove_node(x)


    num1=G6.number_of_nodes()
    nodes1=G6.nodes()
    left_graph={}
    for x in nodes1:
     left_graph[x]=G6.nodes()

    #print left_graph

    ebunch=G6.edges()
    G6.remove_edges_from(ebunch)
    #G6.add_edge(0,1)

    custom_labels1={}
    for x in nodes1:
     custom_labels1[x] = str(x)+':'+G6.node[x]['name_label']


    for key,value in left_graph.iteritems():
     if(((G6.node[key]['name_label']).split('/')[0]!='pibas') and ((G6.node[key]['name_label']).split('/')[0]!='kegg_ligand')):
      for y in range(0, len(value)):
       if(("pibas" not in (G6.node[value[y]]['name_label']).split('/')[0]) and ((G6.node[value[y]]['name_label']).split('/')[0]!='kegg_ligand')):
        sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
        sparql.setQuery(
                """  
                PREFIX bindingdb: <http://bio2rdf.org/bindingdb:>
                PREFIX pubchem: <http://bio2rdf.org/pubchem:>
                PREFIX pharmgkb: <http://bio2rdf.org/pharmgkb:>
                PREFIX chebi: <http://bio2rdf.org/chebi:>
                PREFIX kegg_ligand: <http://bio2rdf.org/kegg:>
                PREFIX pdb: <http://bio2rdf.org/pdb:>
                PREFIX drugbank: <http://bio2rdf.org/drugbank:>
                PREFIX chembl: <http://bio2rdf.org/chembl:>
                PREFIX ndc: <http://bio2rdf.org/ndc:> 
                select ?p
                WHERE
                {

                 SERVICE SILENT<http://%s.bio2rdf.org/sparql> { 
                 OPTIONAL{
                   %s:%s ?p ?o.
                  FILTER(CONTAINS(str(?p),"%s") && CONTAINS(str(?o),"%s")).
                 }
                }

                }

                """ % (((G6.node[key]['name_label']).split('/')[0]),((G6.node[key]['name_label']).split('/')[0]),((G6.node[key]['name_label']).split('/')[1]),((G6.node[value[y]]['name_label']).split('/')[0]),((G6.node[value[y]]['name_label']).split('/')[1])))
        sparql.setReturnFormat(JSON)
        final_results3 = sparql.query().convert()
        try:
          for result3 in final_results3["results"]["bindings"]:
            if(result3["p"]["value"]!=""): 
             #print 'result:'+key+'-'+result3["p"]["value"]+'-'+value[y]
             G6.add_edge(key,value[y])   

        except:     
          continue

    for key,value in left_graph.iteritems():
        if(((G6.node[key]['name_label']).split('/')[0]!='pibas') and ((G6.node[key]['name_label']).split('/')[0]!='kegg_ligand')):
         for y in range(0, len(value)):
          if(("pibas" not in (G6.node[value[y]]['name_label'])) and ((G6.node[value[y]]['name_label']).split('/')[0]=='kegg_ligand')):
            sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
            sparql.setQuery(
                """  
                PREFIX bindingdb: <http://bio2rdf.org/bindingdb:>
                PREFIX pubchem: <http://bio2rdf.org/pubchem:>
                PREFIX pharmgkb: <http://bio2rdf.org/pharmgkb:>
                PREFIX chebi: <http://bio2rdf.org/chebi:>
                PREFIX kegg_ligand: <http://bio2rdf.org/kegg:>
                PREFIX pdb: <http://bio2rdf.org/pdb:>
                PREFIX drugbank: <http://bio2rdf.org/drugbank:>
                PREFIX chembl: <http://bio2rdf.org/chembl:>
                PREFIX ndc: <http://bio2rdf.org/ndc:> 
                select ?p
                WHERE
                {

                 SERVICE SILENT<http://%s.bio2rdf.org/sparql> { 
                 OPTIONAL{
                   %s:%s ?p ?o.
                  FILTER(CONTAINS(str(?p),"kegg") && CONTAINS(str(?o),"%s")).
                 }
                }

                }

                """ % (((G6.node[key]['name_label']).split('/')[0]),((G6.node[key]['name_label']).split('/')[0]),((G6.node[key]['name_label']).split('/')[1]),((G6.node[value[y]]['name_label']).split('/')[1])))
            sparql.setReturnFormat(JSON)
            final_results4 = sparql.query().convert()
            try:
             for result4 in final_results4["results"]["bindings"]:
              if(result4["p"]["value"]!=""): 
                #print 'result:'+key+'-'+result4["p"]["value"]+'-'+value[y]
                G6.add_edge(key,value[y])
            except:      
              continue




    for key,value in left_graph.iteritems():
        if(((G6.node[key]['name_label']).split('/')[0]!='pibas') and ((G6.node[key]['name_label']).split('/')[0]=='kegg_ligand')):
         for y in range(0, len(value)):
          if(("pibas" not in (G6.node[value[y]]['name_label'])) and ((G6.node[value[y]]['name_label']).split('/')[0]!='kegg_ligand')):
            #print (G6.node[key]['name_label']).split('/')[0] +' '+(G6.node[value[y]]['name_label']).split('/')[0]
            sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
            sparql.setQuery(
                """  
                PREFIX bindingdb: <http://bio2rdf.org/bindingdb:>
                PREFIX pubchem: <http://bio2rdf.org/pubchem:>
                PREFIX pharmgkb: <http://bio2rdf.org/pharmgkb:>
                PREFIX chebi: <http://bio2rdf.org/chebi:>
                PREFIX kegg_ligand: <http://bio2rdf.org/kegg:>
                PREFIX pdb: <http://bio2rdf.org/pdb:>
                PREFIX drugbank: <http://bio2rdf.org/drugbank:>
                PREFIX chembl: <http://bio2rdf.org/chembl:> 
                PREFIX ndc: <http://bio2rdf.org/ndc:>
                select ?p
                WHERE
                {

                 SERVICE SILENT<http://kegg.bio2rdf.org/sparql> { 
                 OPTIONAL{
                   %s:%s ?p ?o.
                  FILTER(CONTAINS(str(?p),"%s") && CONTAINS(str(?o),"%s")).
                 }
                }

                }

                """ % (((G6.node[key]['name_label']).split('/')[0]),((G6.node[key]['name_label']).split('/')[1]),((G6.node[value[y]]['name_label']).split('/')[0]),((G6.node[value[y]]['name_label']).split('/')[1])))
            sparql.setReturnFormat(JSON)
            final_results5 = sparql.query().convert()
            try:
             for result5 in final_results5["results"]["bindings"]:
              if(result5["p"]["value"]!=""): 
                #print 'result:'+key+'-'+result5["p"]["value"]+'-'+value[y]
                G6.add_edge(key,value[y])
            except:          
              continue



    for key,value in left_graph.iteritems():
        if(((G6.node[key]['name_label']).split('/')[0]!='pibas') and ((G6.node[key]['name_label']).split('/')[0]=='kegg_ligand')):
         for y in range(0, len(value)):
          if(("pibas" not in (G6.node[value[y]]['name_label'])) and ((G6.node[value[y]]['name_label']).split('/')[0]=='kegg_ligand')):
            sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
            sparql.setQuery(
                """  
                PREFIX bindingdb: <http://bio2rdf.org/bindingdb:>
                PREFIX pubchem: <http://bio2rdf.org/pubchem:>
                PREFIX pharmgkb: <http://bio2rdf.org/pharmgkb:>
                PREFIX chebi: <http://bio2rdf.org/chebi:>
                PREFIX kegg_ligand: <http://bio2rdf.org/kegg:>
                PREFIX pdb: <http://bio2rdf.org/pdb:>
                PREFIX drugbank: <http://bio2rdf.org/drugbank:>
                PREFIX chembl: <http://bio2rdf.org/chembl:>
                PREFIX ndc: <http://bio2rdf.org/ndc:> 
                select ?p
                WHERE
                {

                 SERVICE SILENT<http://kegg.bio2rdf.org/sparql> { 
                 OPTIONAL{
                   %s:%s ?p ?o.
                  FILTER(CONTAINS(str(?p),"kegg") && CONTAINS(str(?o),"%s")).
                 }
                }

                }

                """ % (((G6.node[key]['name_label']).split('/')[0]),((G6.node[key]['name_label']).split('/')[1]),((G6.node[value[y]]['name_label']).split('/')[1])))
            sparql.setReturnFormat(JSON)
            final_results6 = sparql.query().convert()
            try:
             for result6 in final_results6["results"]["bindings"]:
              if(result6["p"]["value"]!=""): 
                #print 'result:'+key+'-'+result6["p"]["value"]+'-'+value[y]
                G6.add_edge(key,value[y])
            except:         
              continue

    nodes1=G6.nodes()

    custom_labels1={}
    for x in nodes1:
     custom_labels1[x] = str(x)+':'+G6.node[x]['name_label']

    #create and draw graph G6
    pos=nx.circular_layout(G6,dim=2, scale=100)
    plt.clf()
    nx.draw(G6, labels=custom_labels1, with_labels=True)
    plt.savefig('/var/www/specint.org/public_html/specint/img/noncompleted_grpah_bio2rdf_'+inchikey+'.png')


    #******************Using Chem2Bio2RDF***********************************************************************************

    sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

    #graphs generated using Chem2Bio2RDF and Chembl
    G7 = nx.DiGraph()

    index_for_new_graph=max(G6.nodes())

   

    url = 'https://www.ebi.ac.uk/unichem/rest/inchikey/'+inchikey
    #print url
    j=index_for_new_graph+1
    storage = StringIO()
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.WRITEFUNCTION, storage.write)
    c.perform()
    c.close()
    content = storage.getvalue()
    unichem_content = json.loads(content)
    for rs in unichem_content: 
     src_compound_id= rs['src_compound_id']
     url = 'https://www.ebi.ac.uk/unichem/rest/sources/'+rs['src_id']
     storage = StringIO()
     c = pycurl.Curl()
     c.setopt(c.URL, url)
     c.setopt(c.WRITEFUNCTION, storage.write)
     c.perform()
     c.close()
     content = storage.getvalue()
     unichem_src_name = json.loads(content)
     for rs in unichem_src_name:
      name=rs['name']
      G7.add_node(j,name_label=name+'/'+src_compound_id)
      j=j+1


    chem2bio2rdf_dataset=['bindingdb','pubchem','chebi','kegg_ligand','kegg','pdb','drugbank','chembl','uniprot','matador','ctd','dcdb','hgnc','pharmgkb','hprd','chemogenomics']


    num1=G7.number_of_nodes()
    nodes1=G7.nodes()
    for x in range(index_for_new_graph+1,index_for_new_graph+num1):
     if(((G7.node[x]['name_label']).split('/')[0] in bio2rdf_dataset) and ("chembl" not in (G7.node[x]['name_label']).split('/')[0])):
      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
      if((G7.node[x]['name_label']).split('/')[0]=='kegg_ligand'):
       sparql.setQuery(
       """  
        PREFIX bindingdb: <http://chem2bio2rdf.org/bindingdb/resource/>
        PREFIX pubchem: <http://chem2bio2rdf.org/pubchem/resource/>
        PREFIX uniprot: <http://chem2bio2rdf.org/uniprot/resource/>
        PREFIX chebi: <http://chem2bio2rdf.org/chebi/resource/chebi/CHEBI~>
        PREFIX kegg_ligand: <http://chem2bio2rdf.org/kegg/resource/kegg_ligand/>
        PREFIX pdb: <http://chem2bio2rdf.org/pdb/resource/pdb_ligand/>
        PREFIX drugbank: <http://chem2bio2rdf.org/drugbank/resource/>
        PREFIX matador: <http://chem2bio2rdf.org/matador/resource/>
        PREFIX chembl: <http://chem2bio2rdf.org/chembl/resource/> 
        PREFIX uniprot: <http://chem2bio2rdf.org/uniprot/resource/>
        PREFIX db: <http://chem2bio2rdf.org/kegg/resource/>
        select ?o
        WHERE
        {

         SERVICE SILENT<http://cheminfov.informatics.indiana.edu:8080/kegg/sparql> { 
         OPTIONAL{
          %s:%s db:CID ?o.

            }

        }

        }

       """ % (((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[1])))

      elif((G7.node[x]['name_label']).split('/')[0]=='drugbank'):
       
       sparql.setQuery(
       """
        PREFIX db:<http://chem2bio2rdf.org/drugbank/resource/>
        PREFIX drugbank:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
        select ?o
        WHERE
        {

         
         SERVICE SILENT<http://147.91.203.161:8890/sparql>{
          drugbank:%s db:CID ?o.

            }

        }
       """ % (G7.node[x]['name_label']).split('/')[1])  
      else:
       sparql.setQuery(
       """  
        PREFIX bindingdb: <http://bio2rdf.org/bindingdb:>
        PREFIX pubchem: <http://bio2rdf.org/pubchem:>
        PREFIX pharmgkb: <http://bio2rdf.org/pharmgkb:>
        PREFIX chebi: <http://bio2rdf.org/chebi:>
        PREFIX kegg_ligand: <http://bio2rdf.org/kegg:>
        PREFIX pdb: <http://bio2rdf.org/pdb:>
        PREFIX drugbank: <http://bio2rdf.org/drugbank:>
        PREFIX chembl: <http://bio2rdf.org/chembl:>
        PREFIX ndc: <http://bio2rdf.org/ndc:> 
        PREFIX db: <http://chem2bio2rdf.org/%s/resource/>
        select ?o
        WHERE
        {

         SERVICE SILENT<http://cheminfov.informatics.indiana.edu:8080/%s/sparql> { 
         OPTIONAL{
          %s:%s db:CID ?o.
            }
        }

        }

       """ % (((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[1])))

      sparql.setReturnFormat(JSON)
      final_results1 = sparql.query().convert()
      try:
       v=G7.number_of_nodes()
       for result1 in final_results1["results"]["bindings"]:
        if(result1["o"]["value"]!=""):  
         node_for_add=result1["o"]["value"] 
         node_of_second_level=node_for_add.split("/")[-1]
         if(node_of_second_level.split(':')[0]=='kegg'):
          node_for_add='kegg_ligand/'+node_of_second_level.split(':')[1]
         else:
          node_for_add=node_of_second_level.split(':')[0]+'/'+node_of_second_level.split(':')[1]
          #print node_for_add
         compare_node=[]
         for x in range(index_for_new_graph+1,index_for_new_graph+num1):
           compare_node.append(str((G7.node[x]['name_label']).split('/')[0])+'/'+str((G7.node[x]['name_label']).split('/')[1]))
         if(node_for_add not in compare_node):
          G7.add_node(v,name_label=node_for_add,predicate=result1["p"]["value"])
          v=v+1
          #print v
      except:     
        continue



    nodes=G7.nodes()

    #print nodes

    edges = combinations(nodes, 2)
    G7.add_nodes_from(nodes)
    G7.add_edges_from(edges)

    num1=G7.number_of_nodes()
    nodes1=G7.nodes()
    left_graph_for_remove=[]
    for x in nodes1:
     if(((G7.node[x]['name_label']).split('/')[0] not in chem2bio2rdf_dataset)):
      G7.remove_node(x)


    #print left_graph
    ebunch=G7.edges()
    G7.remove_edges_from(ebunch)
    #G7.add_edge(index_for_new_graph+1,index_for_new_graph+2)


    num1=G7.number_of_nodes()
    nodes1=G7.nodes()
    left_graph={}
    for x in nodes1:
     left_graph[x]=G7.nodes()

    #print left_graph
    
    for key,value in left_graph.iteritems():
     if(((G7.node[key]['name_label']).split('/')[0]=='drugbank')):
      for y in range(0, len(value)):
       sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
       sparql.setQuery(
        """
         PREFIX drugbank:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         SELECT ?p
         WHERE {
          SERVICE SILENT<http://147.91.203.161:8890/sparql>{
           %s:%s ?p ?o.
           FILTER(regex(str(?o), "%s") && regex(str(?o), "%s")).
          }
         }
        """ % (((G7.node[key]['name_label']).split('/')[0]),((G7.node[key]['name_label']).split('/')[1]),((G7.node[value[y]]['name_label']).split('/')[0]),((G7.node[value[y]]['name_label']).split('/')[1])))
       sparql.setReturnFormat(JSON)
       final_results7 = sparql.query().convert()
       try:
        for result7 in final_results7["results"]["bindings"]:
         if(result7["p"]["value"]!=""):
          #print 'result:'+key+'-'+result7["p"]["value"]+'-'+value[y]
          G7.add_edge(key,value[y])
       except:
        continue      
      

    for key,value in left_graph.iteritems():
     if(((G7.node[key]['name_label']).split('/')[0]!='kegg_ligand')):
      for y in range(0, len(value)):
       if(((G7.node[value[y]]['name_label']).split('/')[0]!='kegg_ligand')):
        sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
        sparql.setQuery(
                """  
                PREFIX bindingdb: <http://chem2bio2rdf.org/bindingdb/resource/bindingdb_ligand/>
                PREFIX pubchem: <http://chem2bio2rdf.org/pubchem/resource/pubchem_compound/>
                PREFIX uniprot: <http://chem2bio2rdf.org/uniprot/resource/>
                PREFIX chebi: <http://chem2bio2rdf.org/chebi/resource/chebi/CHEBI~>
                PREFIX kegg_ligand: <http://chem2bio2rdf.org/kegg/resource/kegg_ligand/>
                PREFIX pdb: <http://chem2bio2rdf.org/pdb/resource/pdb_ligand/>
                PREFIX drugbank: <http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
                PREFIX pharmgkb:<http://chem2bio2rdf.org/pharmgkb/resource/pharmgkb_drug/>
                PREFIX matador: <http://chem2bio2rdf.org/matador/resource/matador/>
                PREFIX chembl: <http://chem2bio2rdf.org/chembl/resource/chembl_compounds/> 
                PREFIX uniprot: <http://chem2bio2rdf.org/uniprot/resource/uniprot/>
                select ?p
                WHERE
                {

                 SERVICE SILENT<http://cheminfov.informatics.indiana.edu:8080/%s/sparql> { 
                 OPTIONAL{
                   %s:%s ?p ?o.
                  FILTER(regex(str(?o), "%s") && regex(str(?o), "%s")).
                 }
                }

                }

                """ % (((G7.node[key]['name_label']).split('/')[0]),((G7.node[key]['name_label']).split('/')[0]),((G7.node[key]['name_label']).split('/')[1]),((G7.node[value[y]]['name_label']).split('/')[0]),((G7.node[value[y]]['name_label']).split('/')[1])))
        sparql.setReturnFormat(JSON)
        final_results3 = sparql.query().convert()
        try:
          for result3 in final_results3["results"]["bindings"]:
            if(result3["p"]["value"]!=""): 
             #print 'result:'+key+'-'+result3["p"]["value"]+'-'+value[y]
             G7.add_edge(key,value[y])   

        except:     
          continue

    for key,value in left_graph.iteritems():
        if(((G7.node[key]['name_label']).split('/')[0]!='kegg_ligand')):
         for y in range(0, len(value)):
          if(((G7.node[value[y]]['name_label']).split('/')[0]=='kegg_ligand')):
            sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
            sparql.setQuery(
                """  
                PREFIX bindingdb: <http://chem2bio2rdf.org/bindingdb/resource/bindingdb_ligand/>
                PREFIX pubchem: <http://chem2bio2rdf.org/pubchem/resource/pubchem_compound/>
                PREFIX uniprot: <http://chem2bio2rdf.org/uniprot/resource/>
                PREFIX chebi: <http://chem2bio2rdf.org/chebi/resource/chebi/CHEBI~>
                PREFIX kegg_ligand: <http://chem2bio2rdf.org/kegg/resource/kegg_ligand/>
                PREFIX pdb: <http://chem2bio2rdf.org/pdb/resource/pdb_ligand/>
                PREFIX drugbank: <http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
                PREFIX matador: <http://chem2bio2rdf.org/matador/resource/matador/>
                PREFIX chembl: <http://chem2bio2rdf.org/chembl/resource/chembl_compounds/> 
                PREFIX uniprot: <http://chem2bio2rdf.org/uniprot/resource/uniprot/> 
                PREFIX pharmgkb:<http://chem2bio2rdf.org/pharmgkb/resource/pharmgkb_drug/>
                select ?p
                WHERE
                {

                 SERVICE SILENT<http://cheminfov.informatics.indiana.edu:8080/%s/sparql> { 
                 OPTIONAL{
                   %s:%s ?p ?o.
                  FILTER(regex(str(?o), "%s") && regex(str(?o), "%s")).
                 }
                }

                }

                """ % (((G7.node[key]['name_label']).split('/')[0]),((G7.node[key]['name_label']).split('/')[0]),((G7.node[key]['name_label']).split('/')[1]),((G7.node[value[y]]['name_label']).split('/')[0]),((G7.node[value[y]]['name_label']).split('/')[1])))
            sparql.setReturnFormat(JSON)
            final_results4 = sparql.query().convert()
            try:
             for result4 in final_results4["results"]["bindings"]:
              if(result4["p"]["value"]!=""): 
                #print 'result:'+key+'-'+result4["p"]["value"]+'-'+value[y]
                G7.add_edge(key,value[y])
            except:      
              continue




    for key,value in left_graph.iteritems():
        if(((G7.node[key]['name_label']).split('/')[0]=='kegg_ligand')):
         for y in range(0, len(value)):
          if(((G7.node[value[y]]['name_label']).split('/')[0]!='kegg_ligand')):
            #print (G7.node[key]['name_label']).split('/')[0] +' '+(G7.node[value[y]]['name_label']).split('/')[0]
            sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
            sparql.setQuery(
                """  
                PREFIX bindingdb: <http://chem2bio2rdf.org/bindingdb/resource/bindingdb_ligand/>
                PREFIX pubchem: <http://chem2bio2rdf.org/pubchem/resource/pubchem_compound/>
                PREFIX uniprot: <http://chem2bio2rdf.org/uniprot/resource/>
                PREFIX chebi: <http://chem2bio2rdf.org/chebi/resource/chebi/CHEBI~>
                PREFIX kegg_ligand: <http://chem2bio2rdf.org/kegg/resource/kegg_ligand/>
                PREFIX pdb: <http://chem2bio2rdf.org/pdb/resource/pdb_ligand/>
                PREFIX drugbank: <http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
                PREFIX matador: <http://chem2bio2rdf.org/matador/resource/matador/>
                PREFIX chembl: <http://chem2bio2rdf.org/chembl/resource/chembl_compounds/> 
                PREFIX uniprot: <http://chem2bio2rdf.org/uniprot/resource/uniprot/>
                PREFIX pharmgkb:<http://chem2bio2rdf.org/pharmgkb/resource/pharmgkb_drug/>
                select ?p
                WHERE
                {

                 SERVICE SILENT<http://cheminfov.informatics.indiana.edu:8080/kegg/sparql> { 
                 OPTIONAL{
                   %s:%s ?p ?o.
                  FILTER(regex(str(?o), "%s") && regex(str(?o), "%s")).
                 }
                }

                }

                """ % (((G7.node[key]['name_label']).split('/')[0]),((G7.node[key]['name_label']).split('/')[1]),((G7.node[value[y]]['name_label']).split('/')[0]),((G7.node[value[y]]['name_label']).split('/')[1])))
            sparql.setReturnFormat(JSON)
            final_results5 = sparql.query().convert()
            try:
             for result5 in final_results5["results"]["bindings"]:
              if(result5["p"]["value"]!=""): 
                #print 'result:'+key+'-'+result5["p"]["value"]+'-'+value[y]
                G7.add_edge(key,value[y])
            except:          
              continue



    for key,value in left_graph.iteritems():
        if(((G7.node[key]['name_label']).split('/')[0]=='kegg_ligand')):
         for y in range(0, len(value)):
          if(((G7.node[value[y]]['name_label']).split('/')[0]=='kegg_ligand')):
            sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
            sparql.setQuery(
                """  
                PREFIX bindingdb: <http://chem2bio2rdf.org/bindingdb/resource/bindingdb_ligand/>
                PREFIX pubchem: <http://chem2bio2rdf.org/pubchem/resource/pubchem_compound/>
                PREFIX uniprot: <http://chem2bio2rdf.org/uniprot/resource/>
                PREFIX chebi: <http://chem2bio2rdf.org/chebi/resource/chebi/CHEBI~>
                PREFIX kegg_ligand: <http://chem2bio2rdf.org/kegg/resource/kegg_ligand/>
                PREFIX pdb: <http://chem2bio2rdf.org/pdb/resource/pdb_ligand/>
                PREFIX drugbank: <http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
                PREFIX matador: <http://chem2bio2rdf.org/matador/resource/matador/>
                PREFIX chembl: <http://chem2bio2rdf.org/chembl/resource/chembl_compounds/> 
                PREFIX uniprot: <http://chem2bio2rdf.org/uniprot/resource/uniprot/>
                PREFIX pharmgkb:<http://chem2bio2rdf.org/pharmgkb/resource/pharmgkb_drug/> 
                select ?p
                WHERE
                {

                 SERVICE SILENT<http://cheminfov.informatics.indiana.edu:8080/kegg/sparql> { 
                 OPTIONAL{
                   %s:%s ?p ?o.
                  FILTER(regex(str(?o), "%s") && regex(str(?o), "%s")).
                 }
                }

                }

                """ % (((G7.node[key]['name_label']).split('/')[0]),((G7.node[key]['name_label']).split('/')[1]),((G7.node[value[y]]['name_label']).split('/')[0]),((G7.node[value[y]]['name_label']).split('/')[1])))
            sparql.setReturnFormat(JSON)
            final_results6 = sparql.query().convert()
            try:
             for result6 in final_results6["results"]["bindings"]:
              if(result6["p"]["value"]!=""): 
                #print 'result:'+key+'-'+result6["p"]["value"]+'-'+value[y]
                G7.add_edge(key,value[y])
            except:         
              continue
    

    nodes1=G7.nodes()
   
    custom_labels1={}
    for x in nodes1:
     custom_labels1[x] = str(x)+':'+G7.node[x]['name_label']
    
   
    

    #create and draw graph G7
    pos=nx.circular_layout(G7,dim=2, scale=100)
    plt.clf()
    nx.draw(G7, labels=custom_labels1, with_labels=True)
    plt.savefig('/var/www/specint.org/public_html/specint/img/noncompleted_graph_chem2bio2rdf_'+inchikey+'.png')


    #*************Join grpahs*************************
    G10=nx.union(G6,G7)


    nodes1=G10.nodes()
    custom_labels1={}
    for x in nodes1:
     custom_labels1[x] = str(x)+':'+G10.node[x]['name_label']

    #create and draw graph G10
    pos=nx.circular_layout(G10,dim=2, scale=100)
    plt.clf()
    nx.draw(G10, labels=custom_labels1, with_labels=True)
    plt.savefig('/var/www/specint.org/public_html/specint/img/noncompleted_graph_'+inchikey+'.png')
    plt.savefig('/var/www/specint.org/public_html/specint/img/noncompleted_graph_'+inchikey+'.eps',format='eps', dpi=300)


    #************************Removing node************************
    nodes_for_connection=[]
    nodes_value_for_connection=[]

    #connection_node='kegg_ligand'
    nodes5=G10.nodes()
    for x in nodes5:
     if(((G10.node[x]['name_label']).split('/')[0]==connection_node) and (G10.neighbors(x)>=0)):
      nodes_for_connection.append(x)
      nodes_value_for_connection.append((G10.node[x]['name_label']).split('/')[1])

    
    #print nodes_for_connection
    my_array_for_node_connected=[]
    for x in range(0,len(nodes_value_for_connection)-1):
      for y in range(x+1,len(nodes_value_for_connection)):
        if(nodes_value_for_connection[x]==nodes_value_for_connection[y]):
         #print nodes_value_for_connection[x]
         #print nodes_value_for_connection[y]
         my_array_for_node_connected.append(nodes_for_connection[x])
         my_array_for_node_connected.append(nodes_for_connection[y])
    
    #print my_array_for_node_connected
    #print max(my_array_for_node_connected)
    #print min(my_array_for_node_connected)
    if((len(my_array_for_node_connected)>1) and (G10.degree(max(my_array_for_node_connected))>0)):
     out_edge=G10.out_edges(max(my_array_for_node_connected))
     in_edge=G10.in_edges(max(my_array_for_node_connected))
     
     out_edge=[item for item in out_edge if item[0] != item[1]]
     in_edge=[item for item in in_edge if item[0] != item[1]]
     if(len(out_edge)==0): neighbors=in_edge
     if(len(in_edge)==0): neighbors=out_edge
     if((len(in_edge)>0) and (len(out_edge)>0)): neighbors=zip(out_edge,in_edge)


     #print neighbors
    
     neighbors_new = []
     for n,i in enumerate(neighbors):
         if i[0] == max(my_array_for_node_connected):
             t = (min(my_array_for_node_connected), i[1])
             neighbors_new.append(t)
         if i[1] == max(my_array_for_node_connected):
             t = (i[0], min(my_array_for_node_connected))
             neighbors_new.append(t)

     #print neighbors_new

     G10.add_edges_from(neighbors_new)

     G10.remove_node(max(my_array_for_node_connected))
     
     nodes1=G10.nodes()
     custom_labels1={}
     for x in nodes1:
      custom_labels1[x] = str(x)+':'+G10.node[x]['name_label']


     #create and draw graph G10
     #pos=nx.circular_layout(G10,dim=2, scale=100)
     plt.clf()
     nx.draw_circular(G10, labels=custom_labels1, with_labels=True)
     plt.savefig('/var/www/specint.org/public_html/specint/img/noncompleted_graph_'+inchikey+'.png')
     plt.savefig('/var/www/specint.org/public_html/specint/img/noncompleted_graph_'+inchikey+'.eps',format='eps', dpi=300)


     #*************Pagerank***************
     #print number_connected_components(G10)
     #print number_of_nodes(G10)
     
     G11=G10.to_undirected()
     remove = [node for node,degree in G11.degree().items() if degree == 0]
     G11.remove_nodes_from(remove)
     #print nx.is_connected(G11)
     if(nx.is_connected(G11)==True):
      G10.remove_edges_from(G10.selfloop_edges())
      pr = nx.pagerank(G10)      
      create_python_file(G10.nodes(),G10.edges(),inchikey)
      H = nx.DiGraph()
      H.add_nodes_from(G10.nodes(data=True))
      H.add_edges_from(G10.edges())	
      return {'pr':pr, 'di_graph':H,'node':index_for_new_graph+1}
     else:
      return {'pr':'Directed graph can not be connected for a given connected node! Please, try again!'}
    else:
     #return non_completed_graph(inchikey,connection_node)
     return {'pr':'Directed graph can not be connected for a given connected node! Please, try again!'}


#non_completed_graph("IHUNBGSDBOWDMA-AQFIFDHZSA-N","pubchem")
