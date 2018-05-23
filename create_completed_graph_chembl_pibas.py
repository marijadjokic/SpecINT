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
import json



#function create_python_file
def create_python_file(nodes,edges,inchikey):
    file = open("/var/www/specint.org/public_html/specint/test/"+inchikey+"_unoriented.py", "w");
    #add preambula
    file.write("#!/usr/bin/env python\n"+
    "import networkx as nx\n" +
    "import numpy as np\n" + 
    "import matplotlib\n" +
    "matplotlib.use('Agg')\n" +
    "import matplotlib.pyplot as plt\n\n")

    file.write("G = nx.Graph()\n\n")
    
   
    
    #add nodes
    for n in nodes:
     #print n
     file.write("G.add_node("+ str(n) + ")\n")

    #add edges
    file.write("\n")	
    for u,v in edges:
     file.write("G.add_edge(" + str(u) + "," + str(v) + ")\n")
    file.write("\n")
    file.write("print G.nodes()\n")
    file.write("print nx.fiedler_vector(G,method='lobpcg')\n")
    #draw graph
    file.write("\n")
    file.write("pos=nx.circular_layout(G,dim=50, scale=100)\n")
    file.write("plt.clf()\n")
    file.write("nx.draw(G,with_labels=True)\n")
    file.write("plt.savefig('C:/Users/Branko/Desktop/"+inchikey+"_unoriented.png')\n")

    file.close()


##########Bio2RDF and Chem2Bio2RDF################################################



#***************Using Bio2RDF********************************************************

#image='myimage_for_completed_graph.png'

def completed_graph(inchikey):


    #graphs generated using PIBAS dataset, UniChem and Bio2RDF
    G6 = nx.Graph()
    G8 = nx.Graph()


    #we start from PIBAS local ontology and  for a given InChiKey we found compound acronym 


    

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
      #print name
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

    nodes1=G6.nodes()
    #print nodes1
    custom_labels1={}
    for x in nodes1:
     custom_labels1[x] = str(x)+':'+G6.node[x]['name_label']
   
    
    #print custom_labels1
    #create and draw graph G6
    pos=nx.circular_layout(G6,dim=2, scale=100)
    plt.clf()
    nx.draw(G6, labels=custom_labels1, with_labels=True)
    plt.savefig('/var/www/specint.org/public_html/specint/img/completed_graph_bio2rdf_'+inchikey+'.png')


    #******************Using Chem2Bio2RDF***********************************************************************************

    sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

    #graphs generated using Chem2Bio2RDF and Chembl
    G7 = nx.Graph()

    sparql.setQuery(
     """  
        PREFIX chembl: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
        SELECT ?chembl_compound_acronym
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblontology.owl>
        WHERE
        {
          ?chembl_compound chembl:InChIkey ?inchikey;
                           chembl:acronym ?chembl_compound_acronym.				
          FILTER(?inchikey="%s").
        }
     """% inchikey)  


    sparql.setReturnFormat(JSON)
    final_results = sparql.query().convert()
    for result in final_results["results"]["bindings"]: 
      chembl_acronym=result["chembl_compound_acronym"]["value"]

    #print chembl_acronym

    index_for_new_graph=max(G6.nodes())
    #print index_for_new_graph

    #add chembl acronym in corresponding form to G7 graph
    chembl_node=str(chembl_acronym.split(':')[0]+'/'+chembl_acronym.split(':')[1])
    #print chembl_node

    G7.add_node(index_for_new_graph+1,name_label=chembl_node)

    #for a given acronym and Chembl mapp we found mapping substance
    sparql.setQuery(
     """  
        PREFIX chembl:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
        SELECT ?mapping_node ?sourceNumber
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl>
        WHERE
        {
          %s chembl:sameAS ?mapping_node;
             chembl:sourceNumber ?sourceNumber.
        }
     """ % chembl_acronym) 

    sparql.setReturnFormat(JSON)
    final_results = sparql.query().convert()
    for result in final_results["results"]["bindings"]: 
      chembl_mapping=str(result["mapping_node"]["value"].split(':')[0]+'/'+result["mapping_node"]["value"].split(':')[1])
      chembl_mapping_source_number=result["sourceNumber"]["value"]

    #print chembl_mapping
    #print '----'

    #add chembl mapping node to graph G7
    G7.add_node(index_for_new_graph+2,name_label=chembl_mapping)


    url = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/'+str(chembl_mapping.split('/')[1])+'/'+str(chembl_mapping_source_number)
    #print url
    j=index_for_new_graph+3
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
      if (((name+'/'+src_compound_id)!=chembl_node)):
        if (((name+'/'+src_compound_id)!=chembl_mapping)):
          G7.add_node(j,name_label=name+'/'+src_compound_id)
          j=j+1


    chem2bio2rdf_dataset=['bindingdb','pubchem','chebi','kegg_ligand','kegg','pdb','drugbank','chembl','uniprot','matador','ctd','dcdb','hgnc','pharmgkb','hprd']


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


    nodes1=G7.nodes()
    custom_labels1={}
    for x in nodes1:
     custom_labels1[x] = str(x)+':'+G7.node[x]['name_label']


    #create and draw graph G7
    pos=nx.circular_layout(G7,dim=2, scale=100)
    plt.clf()
    nx.draw(G7, labels=custom_labels1, with_labels=True)
    plt.savefig('/var/www/specint.org/public_html/specint/img/completed_graph_chem2bio2rdf_'+inchikey+'.png')


    #*************Join grpahs*************************
    G10=nx.union(G6,G7)

    nodes1=G10.nodes()
    custom_labels1={}
    for x in nodes1:
     custom_labels1[x] = str(x)+':'+G10.node[x]['name_label']

    #create and draw graph G10
    pos=nx.circular_layout(G10,dim=2, scale=300)
    plt.clf()
    nx.draw(G10,  labels=custom_labels1, with_labels=True)
    plt.savefig('/var/www/specint.org/public_html/specint/img/completed_graph_'+inchikey+'.eps',format='eps', dpi=300)
    plt.savefig('/var/www/specint.org/public_html/specint/img/completed_graph_'+inchikey+'.png')
    

    #************************Removing node************************
    nodes_for_connection=[]
    nodes_value_for_connection=[]

    values_for_connected_nodes=['pubchem','chebi','drugbank','kegg_ligand']
    for x in values_for_connected_nodes:
     if(x==chembl_mapping.split('/')[0]):
       values_for_connected_nodes.remove(x)
   
    #print values_for_connected_nodes
    for x in values_for_connected_nodes:
        nodes_for_connection=[]
        #connection_node=random.choice(values_for_connected_nodes)
        connection_node=x
        #print connection_node

        nodes5=G10.nodes()
		
        for x in nodes5:
         if((G10.node[x]['name_label']).split('/')[0]==connection_node):
          nodes_for_connection.append(x)
          nodes_value_for_connection.append((G10.node[x]['name_label']).split('/')[1])


        #print nodes_for_connection

        if(len(nodes_for_connection)>=2):
            my_array_for_node_connected=[]
            for x in range(0,len(nodes_value_for_connection)-1):
              for y in range(x+1,len(nodes_value_for_connection)):
                if(nodes_value_for_connection[x]==nodes_value_for_connection[y]):
                 my_array_for_node_connected.append(nodes_for_connection[x])
                 my_array_for_node_connected.append(nodes_for_connection[y])

            #print my_array_for_node_connected
            #print max(my_array_for_node_connected)
            #print min(my_array_for_node_connected)
            if(len(my_array_for_node_connected)>0):
             G10.remove_node(max(my_array_for_node_connected))

             for_delete=G10.nodes()
            #print for_delete

             position=(G10.nodes()).index(min(my_array_for_node_connected))
            #print position

             for x in for_delete:
              if((x>min(my_array_for_node_connected)) and (x not in G6.nodes())):
               G10.add_edge(min(my_array_for_node_connected),x)

             nodes1=G10.nodes()
             custom_labels1={}
             for x in nodes1:
              custom_labels1[x] = str(x)+':'+G10.node[x]['name_label']
             #create and draw graph G10
             pos=nx.circular_layout(G10,dim=2, scale=300)
             plt.clf()
             nx.draw(G10,  labels=custom_labels1, with_labels=True)
             plt.savefig('/var/www/specint.org/public_html/specint/img/completed_graph_'+inchikey+'.eps',format='eps', dpi=300)
             plt.savefig('/var/www/specint.org/public_html/specint/img/completed_graph_'+inchikey+'.png')
             #H = nx.Graph()
             #H.add_nodes_from(G10.nodes())
             #H.add_edges_from(G10.edges())			
             image='myimage_for_completed_graph.png'
             if(nx.is_connected(G10)==True):
              fv=nx.fiedler_vector(G10,method='lobpcg')
              create_python_file(G10.nodes(),G10.edges(), inchikey)
              return {'fv':fv, 'connection_node':connection_node}
              break
             else:
              #print 'Undirected graph is not connected! Please, try again!';
              return {'fv':'Undirected graph is not connected! Please, try again!'}
            else:
             #print 'Undirected graph is not connected! Please, try again!'
             return {'fv':'Undirected graph is not connected! Please, try again!'}  
        else:
         #return completed_graph(inchikey)
         #continue
         #print 'Undirected graph is not connected! Please, try again!'
         return {'fv':'Undirected graph is not connected! Please, try again!'}
      

#completed_graph("OGQICQVSFDPSEI-UHFFFAOYSA-N")