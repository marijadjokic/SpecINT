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


#******************Using Chem2Bio2RDF***********************************************************************************



def create_graphs_chembl_pibas(inchikey):
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


    #add chembl acronym in corresponding form to G7 graph
    chembl_node=str(chembl_acronym.split(':')[0]+'/'+chembl_acronym.split(':')[1])
    #print chembl_node

    G7.add_node(0,name_label=chembl_node)

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
    G7.add_node(1,name_label=chembl_mapping)


    url = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/'+str(chembl_mapping.split('/')[1])+'/'+str(chembl_mapping_source_number)
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
      if (((name+'/'+src_compound_id)!=chembl_node)):
        if (((name+'/'+src_compound_id)!=chembl_mapping)):
          G7.add_node(j,name_label=name+'/'+src_compound_id)
          j=j+1


    chem2bio2rdf_dataset=['bindingdb','pubchem','chebi','kegg_ligand','kegg','pdb','drugbank','chembl','uniprot','matador','ctd','dcdb','hgnc','pharmgkb','hprd']


    num1=G7.number_of_nodes()
    nodes1=G7.nodes()
    for x in range(0,G7.number_of_nodes()):
     if(((G7.node[x]['name_label']).split('/')[0] in chem2bio2rdf_dataset) and ("chembl" not in (G7.node[x]['name_label']).split('/')[0])):
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
    plt.savefig('/var/www/specint.org/public_html/specint/img/completed_graph_'+inchikey+'.png')

    fv=nx.fiedler_vector(G7,method='lobpcg')
    #print fv



    #******************Using Chem2Bio2RDF***********************************************************************************

    G7=nx.DiGraph()
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


    #add chembl acronym in corresponding form to G7 graph
    chembl_node=str(chembl_acronym.split(':')[0]+'/'+chembl_acronym.split(':')[1])
    #print chembl_node

    G7.add_node(0,name_label=chembl_node)

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
    G7.add_node(1,name_label=chembl_mapping)


    url = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/'+str(chembl_mapping.split('/')[1])+'/'+str(chembl_mapping_source_number)
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
      if (((name+'/'+src_compound_id)!=chembl_node)):
        if (((name+'/'+src_compound_id)!=chembl_mapping)):
          G7.add_node(j,name_label=name+'/'+src_compound_id)
          j=j+1


    chem2bio2rdf_dataset=['bindingdb','pubchem','chebi','kegg_ligand','kegg','pdb','drugbank','chembl','uniprot','matador','ctd','dcdb','hgnc','pharmgkb','hprd']


    num1=G7.number_of_nodes()
    nodes1=G7.nodes()
    for x in range(0,G7.number_of_nodes()):
     if(((G7.node[x]['name_label']).split('/')[0] in chem2bio2rdf_dataset) and ("chembl" not in (G7.node[x]['name_label']).split('/')[0])):
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
          # print node_for_add
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
    G7.add_edge(0,1)


    num1=G7.number_of_nodes()
    nodes1=G7.nodes()
    left_graph={}
    for x in nodes1:
     left_graph[x]=G7.nodes()

    #print left_graph

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
    plt.savefig('/var/www/specint.org/public_html/specint/img/noncompleted_graph_'+inchikey+'.png')


    #print G7.nodes()

    pr = nx.pagerank(G7)
    return {'fv':fv, 'pr':pr,'di_graph':G7}


