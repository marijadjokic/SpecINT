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

##########Bio2RDF and Chem2Bio2RDF################################################

#***************Using Bio2RDF********************************************************

def create_graphs_pibas_chembl(inchikey):

    #graphs generated using PIBAS dataset, UniChem and Bio2RDF
    G6 = nx.Graph()
    

    #we start from PIBAS local ontology and  for a given InChiKey we found compound acronym 

    sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

    sparql.setQuery(
     """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?pibas_compound_acronym ?inchikey
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasOntology.owl>
        WHERE
        {
          ?pibas_compound pibas:InChIkey ?inchikey;
                          pibas:acronym ?pibas_compound_acronym.
          FILTER(?inchikey="%s").	

        }
     """% inchikey)  


    sparql.setReturnFormat(JSON)
    final_results = sparql.query().convert()
    for result in final_results["results"]["bindings"]: 
      pibas_acronym=result["pibas_compound_acronym"]["value"]
    #print pibas_acronym

    #add pibas acronym in corresponding form to G6 graph
    pibas_node=pibas_acronym.split(':')[0]+'/'+pibas_acronym.split(':')[1]
    G6.add_node(0,name_label=pibas_node)

    #for a given acronym and PIBAS mapp we found mapping substance
    sparql.setQuery(
     """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?mapping_node ?sourceNumber
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>
        WHERE
        {
          %s pibas:sameAs ?mapping_node;
             pibas:sourceNumber ?sourceNumber.
        }
     """% pibas_acronym) 

    sparql.setReturnFormat(JSON)
    final_results = sparql.query().convert()
    for result in final_results["results"]["bindings"]: 
      pibas_mapping=result["mapping_node"]["value"].split(':')[0]+'/'+result["mapping_node"]["value"].split(':')[1]
      pibas_mapping_source_number=result["sourceNumber"]["value"]

    #print pibas_mapping

    #add pibas mapping node to graph G6
    G6.add_node(1,name_label=pibas_mapping)



    url = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/'+str(pibas_mapping.split('/')[1])+'/'+str(pibas_mapping_source_number)
    #print url
    j=2
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
      if ((name+'/'+src_compound_id)!=pibas_mapping):
        if((name+'/'+src_compound_id)!=pibas_node):
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
    custom_labels1={}
    for x in nodes1:
     custom_labels1[x] = str(x)+':'+G6.node[x]['name_label']


    #create and draw graph G6
    pos=nx.circular_layout(G6,dim=2, scale=100)
    plt.clf()
    nx.draw(G6, labels=custom_labels1, with_labels=True)
    plt.savefig('/var/www/specint.org/public_html/specint/img/completed_graph_'+inchikey+'.png')



    fv=nx.fiedler_vector(G6,method='lobpcg')
    

    #*****************************************************************************************

    #graphs generated using PIBAS dataset, UniChem and Bio2RDF
    G7 = nx.DiGraph()
    


    #we start from PIBAS local ontology and  for a given InChiKey we found compound acronym 

    sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

    sparql.setQuery(
     """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?pibas_compound_acronym ?inchikey
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasOntology.owl>
        WHERE
        {
          ?pibas_compound pibas:InChIkey ?inchikey;
                          pibas:acronym ?pibas_compound_acronym.
          FILTER(?inchikey="%s").	

        }
     """% inchikey)  


    sparql.setReturnFormat(JSON)
    final_results = sparql.query().convert()
    for result in final_results["results"]["bindings"]: 
      pibas_acronym=result["pibas_compound_acronym"]["value"]
    #print pibas_acronym

    #add pibas acronym in corresponding form to G7 graph
    pibas_node=pibas_acronym.split(':')[0]+'/'+pibas_acronym.split(':')[1]
    G7.add_node(0,name_label=pibas_node)
    
    #for a given acronym and PIBAS mapp we found mapping substance
    sparql.setQuery(
     """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?mapping_node ?sourceNumber
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>
        WHERE
        {
          %s pibas:sameAs ?mapping_node;
             pibas:sourceNumber ?sourceNumber.
        }
     """% pibas_acronym) 

    sparql.setReturnFormat(JSON)
    final_results = sparql.query().convert()
    for result in final_results["results"]["bindings"]: 
      pibas_mapping=result["mapping_node"]["value"].split(':')[0]+'/'+result["mapping_node"]["value"].split(':')[1]
      pibas_mapping_source_number=result["sourceNumber"]["value"]

    #print pibas_mapping
 
    #add pibas mapping node to graph G7
    G7.add_node(1,name_label=pibas_mapping)



    url = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/'+str(pibas_mapping.split('/')[1])+'/'+str(pibas_mapping_source_number)
    #print url
    j=2
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
      if ((name+'/'+src_compound_id)!=pibas_mapping):
        if((name+'/'+src_compound_id)!=pibas_node):
         G7.add_node(j,name_label=name+'/'+src_compound_id)
         j=j+1



    bio2rdf_dataset=['bindingdb','pubchem','pharmgkb','chebi','kegg_ligand','pdb','drugbank','pibas','ndc']

    num1=G7.number_of_nodes()
    nodes1=G7.nodes()
    for x in range(0, num1):
     if(((G7.node[x]['name_label']).split('/')[0] in bio2rdf_dataset) and ("pibas" not in (G7.node[x]['name_label']).split('/')[0])):
      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
      if((G7.node[x]['name_label']).split('/')[0]=='kegg_ligand'):
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

       """ % (((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[1])))

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

       """ % (((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[1]),((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[0]),((G7.node[x]['name_label']).split('/')[0])))

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
         for x in range(0,G7.number_of_nodes()):
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
    for x in range(0, num1):
     if(((G7.node[x]['name_label']).split('/')[0] not in bio2rdf_dataset)):
      G7.remove_node(x)


    num1=G7.number_of_nodes()
    nodes1=G7.nodes()
    left_graph={}
    for x in nodes1:
     left_graph[x]=G7.nodes()

    #print left_graph

    ebunch=G7.edges()
    G7.remove_edges_from(ebunch)
    G7.add_edge(0,1)

    custom_labels1={}
    for x in nodes1:
     custom_labels1[x] = str(x)+':'+G7.node[x]['name_label']


    for key,value in left_graph.iteritems():
     if(((G7.node[key]['name_label']).split('/')[0]!='pibas') and ((G7.node[key]['name_label']).split('/')[0]!='kegg_ligand')):
      for y in range(0, len(value)):
       if(("pibas" not in (G7.node[value[y]]['name_label']).split('/')[0]) and ((G7.node[value[y]]['name_label']).split('/')[0]!='kegg_ligand')):
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
        if(((G7.node[key]['name_label']).split('/')[0]!='pibas') and ((G7.node[key]['name_label']).split('/')[0]!='kegg_ligand')):
         for y in range(0, len(value)):
          if(("pibas" not in (G7.node[value[y]]['name_label'])) and ((G7.node[value[y]]['name_label']).split('/')[0]=='kegg_ligand')):
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

                """ % (((G7.node[key]['name_label']).split('/')[0]),((G7.node[key]['name_label']).split('/')[0]),((G7.node[key]['name_label']).split('/')[1]),((G7.node[value[y]]['name_label']).split('/')[1])))
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
        if(((G7.node[key]['name_label']).split('/')[0]!='pibas') and ((G7.node[key]['name_label']).split('/')[0]=='kegg_ligand')):
         for y in range(0, len(value)):
          if(("pibas" not in (G7.node[value[y]]['name_label'])) and ((G7.node[value[y]]['name_label']).split('/')[0]!='kegg_ligand')):
            #print (G7.node[key]['name_label']).split('/')[0] +' '+(G7.node[value[y]]['name_label']).split('/')[0]
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
        if(((G7.node[key]['name_label']).split('/')[0]!='pibas') and ((G7.node[key]['name_label']).split('/')[0]=='kegg_ligand')):
         for y in range(0, len(value)):
          if(("pibas" not in (G7.node[value[y]]['name_label'])) and ((G7.node[value[y]]['name_label']).split('/')[0]=='kegg_ligand')):
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

                """ % (((G7.node[key]['name_label']).split('/')[0]),((G7.node[key]['name_label']).split('/')[1]),((G7.node[value[y]]['name_label']).split('/')[1])))
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
    plt.savefig('/var/www/specint.org/public_html/specint/img/noncompleted_graph_'+inchikey+'.png')


    pr = nx.pagerank(G7)

    return {'fv':fv, 'pr':pr,'di_graph':G7}
