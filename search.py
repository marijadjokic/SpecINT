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
from array import *
import urllib
from urlparse import urlparse

#InChiKey get from search.php page
inchikey=sys.argv[1]
information=sys.argv[2]
#print inchikey
#information="target"


#graphs generated using PIBAS dataset, UniChem and Bio2RDF
G6 = nx.Graph()
G8 = nx.Graph()


#we start from PIBAS local ontology and  for a given InChiKey we found compound acronym 

sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

sparql.setQuery(
 """  
    PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
    PREFIX chembl: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
    SELECT ?pibas_compound_acronym ?chembl_compound_acronym
    FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasOntology.owl>
    FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblontology.owl>
    WHERE
    {
      OPTIONAL{?pibas_compound pibas:InChIkey ?inchikey_pibas;
                      pibas:acronym ?pibas_compound_acronym.
     FILTER(?inchikey_pibas="%s").} 
     OPTIONAL{?chembl_compound chembl:InChIkey ?inchikey_chembl;
                      chembl:acronym ?chembl_compound_acronym.
      FILTER(?inchikey_chembl="%s").	}
      
    }
 """% (inchikey,inchikey))  
 

sparql.setReturnFormat(JSON)
final_results = sparql.query().convert()
for result in final_results["results"]["bindings"]:
  if(("pibas_compound_acronym" not in result) and ("chembl_compound_acronym" not in result)):
   #print 'Entered InChiKey dosen not exist in PIBAS and CHEMBL mapping ontologies!'
   import create_completed_graph_new
   #create_completed_graph.completed_graph(inchikey)   
   result=create_completed_graph_new.completed_graph(inchikey)
   fv=result['fv']
   if(type(fv) == str):
    print 'fv-'+str(fv)
   else:
    fv_around=np.around(fv,decimals=3)
    print 'fv-'+str(fv_around)
    connection_node=result['connection_node']
    import create_noncompleted_graph_new
    #create_noncompleted_graph.non_completed_graph(inchikey,connection_node)
   
    result1=create_noncompleted_graph_new.non_completed_graph(inchikey,connection_node)
    pr=result1['pr']
    if(type(pr) == str):
     print 'pr-'+str(pr)
    else:
     pr_around = {k:round(v,2) for k, v in pr.items()}
     print 'pr-'+str(pr_around)
     fv = fv_around   
     pg=pr

     G = result1['di_graph']
     nodes = G.nodes()	 
     print result1['node']
     positive_path = [max((t for t in pg if fv[nodes.index(t)] > 0 and abs(fv[nodes.index(t)]) > 0.0000001), key=pg.get)]
     negative_path = [max((t for t in pg if fv[nodes.index(t)] < 0 and abs(fv[nodes.index(t)]) > 0.0000001), key=pg.get)] 
     neutral_node = []
     visited = []

     # cvorovi odredjuju raspored u vektorima
     nodes = G.nodes()
     # positive vertices
     if (len(positive_path) > 0) :
      while positive_path[-1] not in visited:
        visited.append(positive_path[-1])
        neigh = G.neighbors(positive_path[-1])
        #print neigh		
        if not neigh:
            break
        all_ranks = []
        for nvalue in neigh:
         all_ranks.append(pg[nvalue])
        najveci = max(all_ranks)	
	flag = 0
	v_pos = -1
	
	for nvalue in neigh:
	    v_pos = nodes.index(nvalue)
	    if pg[nvalue] == najveci and abs(fv[v_pos]) < 0.0000001:
		positive_path.append(nvalue)
		flag = 1
		break
	    elif pg[nvalue] == najveci:
		v_pos = nvalue
			
	if flag == 1:
	    break
	
	if v_pos != -1:
		positive_path.append(v_pos)

     #print positive_path       
     

     # negative vertices
     if (len(negative_path) > 0) :
      while negative_path[-1] not in visited:
        visited.append(negative_path[-1])
        neigh = G.neighbors(negative_path[-1])
        #print neigh
        if not neigh:
            break
        all_ranks = []
        for nvalue in neigh:
           all_ranks.append(pg[nvalue])
        najveci = max(all_ranks)	
	flag = 0
	v_pos = -1
        
        for nvalue in neigh:
          v_pos = nodes.index(nvalue)
          if pg[nvalue] == najveci and abs(fv[v_pos]) < 0.0000001:
            negative_path.append(nvalue)
	    flag = 1
	    break
          elif pg[nvalue] == najveci:
	    v_pos = nvalue

        if flag == 1:
		break
	
	if v_pos != -1:
		negative_path.append(v_pos)
     
     
     #print negative_path
     
    
     value_for_positive_path=[]
     for x in positive_path:
       name=nx.get_node_attributes(G,'name_label')
       value_for_positive_path.append(name[x])

     #print value_for_positive_path

     value_for_negative_path=[]
     for x in negative_path:
       name=nx.get_node_attributes(G,'name_label')
       value_for_negative_path.append(name[x])
        
     #print value_for_negative_path
     
     if(information=='target'):
      array_for_target_bio2rdf=[]
      for x in range(0,len(value_for_positive_path)):
       #print value_for_positive_path[x].split('/')[0];
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?target_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasTarget ?targettemplate.
        ?targettemplate pibas:targetSchema ?target_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        target_schema=result["target_schema"]["value"]
        array_for_target_bio2rdf.append(target_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_target_bio2rdf
      array_for_target_chem2bio2rdf=[]
      for x in range(0,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?target_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasTarget ?targettemplate.
        ?targettemplate pibas:targetSchema ?target_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        target_schema=result["target_schema"]["value"]
        array_for_target_chem2bio2rdf.append(target_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_target_chem2bio2rdf

      bio2rdf=''
      if len(array_for_target_bio2rdf)>0:
        for x in range(0,len(array_for_target_bio2rdf)):
         bio2rdf=bio2rdf+'{'+array_for_target_bio2rdf[x]+'} UNION'
        #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_target_chem2bio2rdf)>0:
        for x in range(0,len(array_for_target_chem2bio2rdf)):
         chem2bio2rdf=chem2bio2rdf+'{'+array_for_target_chem2bio2rdf[x]+'} UNION'
        #print chem2bio2rdf
      
      if (bio2rdf!="" and chem2bio2rdf!=""):
         final_query= bio2rdf.rsplit(' ', 1)[0] + ' UNION ' + chem2bio2rdf.rsplit(' ', 1)[0]
      elif (bio2rdf=="" and chem2bio2rdf!=""):
          final_query= chem2bio2rdf.rsplit(' ', 1)[0]
      elif (bio2rdf!="" and chem2bio2rdf==""):
          final_query= bio2rdf.rsplit(' ', 1)[0]
      else:
          final_query= ""
          
      #print final_query
      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")
      
      sparql.setQuery(
      """
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?target
         { {   

         %s
         
         
         }
         }
      """ % (final_query)) 
      
      my_result='target-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["target"]["value"]
      string_for_query="PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?target \n { \n %s \n }" %(final_query)
      my_result=my_result+' Query^'+string_for_query
      #print string_for_query
      print my_result




     if(information=='cellline'):
      array_for_cellline_bio2rdf=[]
      for x in range(0,len(value_for_positive_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?cellline_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasCellLine ?celllinetemplate.
        ?celllinetemplate pibas:cellLineSchema ?cellline_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        cellline_schema=result["cellline_schema"]["value"]
        array_for_target_bio2rdf.append(cellline_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_cellline_bio2rdf
      array_for_cellline_chem2bio2rdf=[]
      for x in range(0,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?cellline_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasCellLine ?celllinetemplate.
        ?celllinetemplate pibas:cellLineSchema ?cellline_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        cellline_schema=result["cellline_schema"]["value"]
        array_for_cellline_chem2bio2rdf.append(cellline_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_cellline_chem2bio2rdf
      
      bio2rdf=''
      if len(array_for_cellline_bio2rdf)>0:
       for x in range(0,len(array_for_cellline_bio2rdf)):
        bio2rdf=bio2rdf+'{'+array_for_cellline_bio2rdf[x]+'} UNION'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_cellline_chem2bio2rdf)>0:
       for x in range(0,len(array_for_cellline_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'{'+array_for_cellline_chem2bio2rdf[x]+'} UNION'
       #print chem2bio2rdf
      
      if (bio2rdf!="" and chem2bio2rdf!=""):
         final_query= bio2rdf.rsplit(' ', 1)[0] + ' UNION ' + chem2bio2rdf.rsplit(' ', 1)[0]
      elif (bio2rdf=="" and chem2bio2rdf!=""):
          final_query= chem2bio2rdf.rsplit(' ', 1)[0]
      elif (bio2rdf!="" and chem2bio2rdf==""):
          final_query= bio2rdf.rsplit(' ', 1)[0]
      else:
          final_query= ""
          
      #print final_query
      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?cellline
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl>
         { {  

         %s
        
         }
         }
      """ % (final_query))  
      
      my_result='cellline-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["cellline"]["value"]
      string_for_query="PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?target \n { \n %s \n }" %(final_query)
      my_result=my_result+' Query^'+string_for_query
      print my_result

     if(information=='ic50'):
      array_for_ic50_bio2rdf=[]
      for x in range(0,len(value_for_positive_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?ic50_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasIC50value ?ic50template.
        ?ic50template pibas:IC50valueSchema ?ic50_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        ic50_schema=result["ic50_schema"]["value"]
        array_for_ic50_bio2rdf.append(ic50_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_ic50_bio2rdf
      array_for_ic50_chem2bio2rdf=[]
      for x in range(0,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?ic50_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasIC50value ?ic50template.
        ?ic50template pibas:IC50valueSchema ?ic50_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        ic50_schema=result["ic50_schema"]["value"]
        array_for_ic50_chem2bio2rdf.append(ic5_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_ic50_chem2bio2rdf

      bio2rdf=''
      if len(array_for_ic50_bio2rdf)>0:
       for x in range(0,len(array_for_ic50_bio2rdf)):
        bio2rdf=bio2rdf+'{'+array_for_ic50_bio2rdf[x]+'} UNION'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_ic50_chem2bio2rdf)>0:
       for x in range(0,len(array_for_ic50_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'{'+array_for_ic50_chem2bio2rdf[x]+'} UNION'
       #print chem2bio2rdf
      
       if (bio2rdf!="" and chem2bio2rdf!=""):
         final_query= bio2rdf.rsplit(' ', 1)[0] + ' UNION ' + chem2bio2rdf.rsplit(' ', 1)[0]
      elif (bio2rdf=="" and chem2bio2rdf!=""):
          final_query= chem2bio2rdf.rsplit(' ', 1)[0]
      elif (bio2rdf!="" and chem2bio2rdf==""):
          final_query= bio2rdf.rsplit(' ', 1)[0]
      else:
          final_query= ""
          
      #print final_query
      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?cellline ?IC50value
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl>
         { {  
         %s
         
         }
         }
      """ % (final_query))  
      
      my_result='ic50-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["cellline"]["value"]
        my_result=my_result+' '+result["IC50value"]["value"]
      string_for_query="PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?target \n { \n %s \n }" %(final_query)
      my_result=my_result+' Query^'+string_for_query
      print my_result

  if (("pibas_compound_acronym" in result) and ("chembl_compound_acronym" in result)):
   pibas_acronym=result["pibas_compound_acronym"]["value"]
   chembl_acronym=result["chembl_compound_acronym"]["value"]
     
   import create_completed_graph
   #create_completed_graph.completed_graph(inchikey)   
   result=create_completed_graph.completed_graph(inchikey)
   fv=result['fv']
   if(type(fv) == str):
    print 'fv-'+str(fv)
   else:
    fv_around=np.around(fv,decimals=3)
    print 'fv-'+str(fv_around)
    connection_node=result['connection_node']
    import create_noncompleted_graph
    #create_noncompleted_graph.non_completed_graph(inchikey,connection_node)
   
    result1=create_noncompleted_graph.non_completed_graph(inchikey,connection_node)
    pr=result1['pr']
    if(type(pr) == str):
     print 'pr-'+str(pr)
    else:
     pr_around = {k:round(v,2) for k, v in pr.items()}
     print 'pr-'+str(pr_around)
     fv = fv_around   
     pg=pr

     G = result1['di_graph']
     #print G.nodes()	 
     print result1['node']
     positive_path = [0]
     negative_path = [result1['node']] 
     neutral_node = []
     visited = []

     # cvorovi odredjuju raspored u vektorima
     nodes = G.nodes()
     # positive vertices
     if (len(positive_path) > 0) :
      while positive_path[-1] not in visited:
        visited.append(positive_path[-1])
        neigh = G.neighbors(positive_path[-1])
        #print neigh		
        if not neigh:
            break
        all_ranks = []
        for nvalue in neigh:
         all_ranks.append(pg[nvalue])
        najveci = max(all_ranks)	
	flag = 0
	v_pos = -1
	
	for nvalue in neigh:
	    v_pos = nodes.index(nvalue)
	    if pg[nvalue] == najveci and abs(fv[v_pos]) < 0.0000001:
		positive_path.append(nvalue)
		flag = 1
		break
	    elif pg[nvalue] == najveci:
		v_pos = nvalue
			
	if flag == 1:
	    break
	
	if v_pos != -1:
		positive_path.append(v_pos)
      #print positive_path       
     

     # negative vertices
     if (len(negative_path) > 0) :
      while negative_path[-1] not in visited:
        visited.append(negative_path[-1])
        neigh = G.neighbors(negative_path[-1])
        #print neigh
        if not neigh:
            break
        all_ranks = []
        for nvalue in neigh:
           all_ranks.append(pg[nvalue])
        najveci = max(all_ranks)	
	flag = 0
	v_pos = -1
        
        for nvalue in neigh:
          v_pos = nodes.index(nvalue)
          if pg[nvalue] == najveci and abs(fv[v_pos]) < 0.0000001:
            negative_path.append(nvalue)
	    flag = 1
	    break
          elif pg[nvalue] == najveci:
	    v_pos = nvalue

        if flag == 1:
		break
	
	if v_pos != -1:
		negative_path.append(v_pos)
     
     #print negative_path
     
     


     value_for_positive_path=[]
     for x in positive_path:
       name=nx.get_node_attributes(G,'name_label')
       value_for_positive_path.append(name[x])

     #print value_for_positive_path

     value_for_negative_path=[]
     for x in negative_path:
       name=nx.get_node_attributes(G,'name_label')
       value_for_negative_path.append(name[x])
     #print value_for_negative_path
     
    
     
     if(information=='target'):
      array_for_target_bio2rdf=[]
      for x in range(1,len(value_for_positive_path)):
       #print value_for_positive_path[x].split('/')[0];
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?target_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasTarget ?targettemplate.
        ?targettemplate pibas:targetSchema ?target_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        target_schema=result["target_schema"]["value"]
        array_for_target_bio2rdf.append(target_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_target_bio2rdf
      array_for_target_chem2bio2rdf=[]
      for x in range(1,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?target_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasTarget ?targettemplate.
        ?targettemplate pibas:targetSchema ?target_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        target_schema=result["target_schema"]["value"]
        array_for_target_chem2bio2rdf.append(target_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_target_chem2bio2rdf

      bio2rdf=''
      if len(array_for_target_bio2rdf)>0:
       for x in range(0,len(array_for_target_bio2rdf)):
        bio2rdf=bio2rdf+'UNION{'+array_for_target_bio2rdf[x]+'}'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_target_chem2bio2rdf)>0:
       for x in range(0,len(array_for_target_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'UNION{'+array_for_target_chem2bio2rdf[x]+'}'
       #print chem2bio2rdf

      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?target
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl>
         { {   {  pibas:%s pibas:sameAs ?mapping_node;
                           pibas:hasTarget ?target.
         }

         %s %s
         UNION
         {
         ?result chembl_mapp:sameAS "%s:%s".
         SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/>
           { 
             ?activity a cco:Activity ;
             cco:hasMolecule chembl_molecule:%s;
                             cco:hasAssay ?assay.
             ?assay cco:hasTarget ?target.
           }
         }
         }
         }
      """ % (value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf,value_for_negative_path[1].split('/')[0],value_for_negative_path[1].split('/')[1],value_for_negative_path[0].split('/')[1]))  
      
      my_result='target-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["target"]["value"]
      string_for_query="PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#> \n PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?target \n FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl> \n FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl> \n { \n { \n  { \n  pibas:%s pibas:sameAs ?mapping_node; \n pibas:hasTarget ?target.\n } \n %s  \n %s \n UNION \n { \n ?result chembl_mapp:sameAS '%s:%s'. \n SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/> { \n ?activity a cco:Activity ; \n cco:hasMolecule chembl_molecule:%s; \n cco:hasAssay ?assay. \n ?assay cco:hasTarget ?target. \n } \n } \n } \n }" %(value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf,value_for_negative_path[1].split('/')[0],value_for_negative_path[1].split('/')[1],value_for_negative_path[0].split('/')[1])
      my_result=my_result+' Query^'+string_for_query
      print my_result




     if(information=='cellline'):
      array_for_cellline_bio2rdf=[]
      for x in range(1,len(value_for_positive_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?cellline_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasCellLine ?celllinetemplate.
        ?celllinetemplate pibas:cellLineSchema ?cellline_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        cellline_schema=result["cellline_schema"]["value"]
        array_for_target_bio2rdf.append(cellline_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_cellline_bio2rdf
      array_for_cellline_chem2bio2rdf=[]
      for x in range(1,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?cellline_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasCellLine ?celllinetemplate.
        ?celllinetemplate pibas:cellLineSchema ?cellline_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        cellline_schema=result["cellline_schema"]["value"]
        array_for_cellline_chem2bio2rdf.append(cellline_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_cellline_chem2bio2rdf
      
      bio2rdf=''
      if len(array_for_cellline_bio2rdf)>0:
       for x in range(0,len(array_for_cellline_bio2rdf)):
        bio2rdf=bio2rdf+'UNION{'+array_for_cellline_bio2rdf[x]+'}'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_cellline_chem2bio2rdf)>0:
       for x in range(0,len(array_for_cellline_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'UNION{'+array_for_cellline_chem2bio2rdf[x]+'}'
       #print chem2bio2rdf
      
      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?cellline
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl>
         { {   {  pibas:%s pibas:sameAs ?mapping_node;
                           pibas:hasCellLine ?cellline.
         }

         %s %s
         UNION
         {
         ?result chembl_mapp:sameAS "%s:%s".
         SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/>
           { 
              ?act ?o chembl_molecule:%s .
              ?act cco:standardUnits "nM" .
              ?act cco:standardType "IC50" .
              ?act cco:standardValue ?IC50value .
              ?assay cco:hasActivity ?act .
              ?assay cco:hasCellLine ?cellline.
           }
         }
         }
         }
      """ % (value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf,value_for_negative_path[1].split('/')[0],value_for_negative_path[1].split('/')[1],value_for_negative_path[0].split('/')[1]))  
      
      my_result='cellline-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["cellline"]["value"]
      string_for_query="PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#> \n PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?cellline \n FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl> \n FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl> \n { \n { \n {  pibas:%s pibas:sameAs ?mapping_node; \n pibas:hasCellLine ?cellline. \n } \n %s \n  %s \n UNION \n { \n ?result chembl_mapp:sameAS '%s:%s'. \n SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/> { \n ?act ?o chembl_molecule:%s . \n ?act cco:standardUnits 'nM' . \n ?act cco:standardType 'IC50' . \n ?act cco:standardValue ?IC50value . \n ?assay cco:hasActivity ?act . \n ?assay cco:hasCellLine ?cellline. \n } \n } \n } \n}" %(value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf,value_for_negative_path[1].split('/')[0],value_for_negative_path[1].split('/')[1],value_for_negative_path[0].split('/')[1])
      my_result=my_result+' Query^'+string_for_query
      print my_result

     if(information=='ic50'):
      array_for_ic50_bio2rdf=[]
      for x in range(1,len(value_for_positive_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?ic50_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasIC50value ?ic50template.
        ?ic50template pibas:IC50valueSchema ?ic50_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        ic50_schema=result["ic50_schema"]["value"]
        array_for_ic50_bio2rdf.append(ic50_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_ic50_bio2rdf
      array_for_ic50_chem2bio2rdf=[]
      for x in range(1,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?ic50_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasIC50value ?ic50template.
        ?ic50template pibas:IC50valueSchema ?ic50_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        ic50_schema=result["ic50_schema"]["value"]
        array_for_ic50_chem2bio2rdf.append(ic5_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_ic50_chem2bio2rdf

      bio2rdf=''
      if len(array_for_ic50_bio2rdf)>0:
       for x in range(0,len(array_for_ic50_bio2rdf)):
        bio2rdf=bio2rdf+'UNION{'+array_for_ic50_bio2rdf[x]+'}'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_ic50_chem2bio2rdf)>0:
       for x in range(0,len(array_for_ic50_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'UNION{'+array_for_ic50_chem2bio2rdf[x]+'}'
       #print chem2bio2rdf

      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?cellline ?IC50value
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl>
         { {   {  pibas:%s pibas:sameAs ?mapping_node;
                           pibas:IC50value ?IC50value;
                           pibas:hasCellLine ?cellline.
         }

         %s %s
         UNION
         {
         ?result chembl_mapp:sameAS "%s:%s".
         SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/>
           { 
              ?act ?o chembl_molecule:%s.
              ?act cco:standardUnits "nM" .
              ?act cco:standardType "IC50" .
              ?act cco:standardValue ?IC50value .
              ?assay cco:hasActivity ?act.
              ?assay cco:hasCellLine ?cellline. 
           }
         }
         }
         }
      """ % (value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf,value_for_negative_path[1].split('/')[0],value_for_negative_path[1].split('/')[1],value_for_negative_path[0].split('/')[1]))  
      
      my_result='ic50-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["cellline"]["value"]
        my_result=my_result+' '+result["IC50value"]["value"]
      string_for_query="PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#> \n PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?cellline ?IC50value \n FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl> \n FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl> \n { \n { \n {  pibas:%s pibas:sameAs ?mapping_node; \n pibas:IC50value ?IC50value; pibas:hasCellLine ?cellline. \n } %s \n %s \n UNION \n { \n ?result chembl_mapp:sameAS '%s:%s'. \n SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/> { \n ?act ?o chembl_molecule:%s.?act cco:standardUnits 'nM' . \n ?act cco:standardType 'IC50' . \n ?act cco:standardValue ?IC50value . \n ?assay cco:hasActivity ?act. \n ?assay cco:hasCellLine ?cellline. \n } \n } \n } \n }" %(value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf,value_for_negative_path[1].split('/')[0],value_for_negative_path[1].split('/')[1],value_for_negative_path[0].split('/')[1])
      my_result=my_result+' Query^'+string_for_query
      print my_result

  if (("pibas_compound_acronym" in result) and ("chembl_compound_acronym" not in result)):
   pibas_acronym=result["pibas_compound_acronym"]["value"]
   

   import create_completed_graph_pibas_chembl
   #create_completed_graph.completed_graph(inchikey)   
   result=create_completed_graph_pibas_chembl.completed_graph(inchikey)
   fv=result['fv']
   if(type(fv) == str):
    print 'fv-'+str(fv)
   else:
    fv_around=np.around(fv,decimals=3)
    print 'fv-'+str(fv_around)
    connection_node=result['connection_node']
    import create_noncompleted_graph_pibas_chembl
    #create_noncompleted_graph.non_completed_graph(inchikey,connection_node)
   
    result1=create_noncompleted_graph_pibas_chembl.non_completed_graph(inchikey,connection_node)
    pr=result1['pr']
    if(type(pr) == str):
     print 'pr-'+str(pr)
    else:
     pr_around = {k:round(v,2) for k, v in pr.items()}
     print 'pr-'+str(pr_around)
     fv = fv_around   
     pg=pr

     G = result1['di_graph']
     nodes = G.nodes()	 
     print result1['node']
     positive_path = [0]
     negative_path = [max((t for t in pg if fv[nodes.index(t)] < 0 and abs(fv[nodes.index(t)]) > 0.0000001), key=pg.get)] 
     neutral_node = []
     visited = []

     # cvorovi odredjuju raspored u vektorima
     nodes = G.nodes()
     # positive vertices
     if (len(positive_path) > 0) :
      while positive_path[-1] not in visited:
        visited.append(positive_path[-1])
        neigh = G.neighbors(positive_path[-1])
        #print neigh		
        if not neigh:
            break
        all_ranks = []
        for nvalue in neigh:
         all_ranks.append(pg[nvalue])
        najveci = max(all_ranks)	
	flag = 0
	v_pos = -1
	
	for nvalue in neigh:
	    v_pos = nodes.index(nvalue)
	    if pg[nvalue] == najveci and abs(fv[v_pos]) < 0.0000001:
		positive_path.append(nvalue)
		flag = 1
		break
	    elif pg[nvalue] == najveci:
		v_pos = nvalue
			
	if flag == 1:
	    break
	
	if v_pos != -1:
		positive_path.append(v_pos)

     #print positive_path       
     

     # negative vertices
     if (len(negative_path) > 0) :
      while negative_path[-1] not in visited:
        visited.append(negative_path[-1])
        neigh = G.neighbors(negative_path[-1])
        #print neigh
        if not neigh:
            break
        all_ranks = []
        for nvalue in neigh:
           all_ranks.append(pg[nvalue])
        najveci = max(all_ranks)	
	flag = 0
	v_pos = -1
        
        for nvalue in neigh:
          v_pos = nodes.index(nvalue)
          if pg[nvalue] == najveci and abs(fv[v_pos]) < 0.0000001:
            negative_path.append(nvalue)
	    flag = 1
	    break
          elif pg[nvalue] == najveci:
	    v_pos = nvalue

        if flag == 1:
		break
	
	if v_pos != -1:
		negative_path.append(v_pos)
     
     
     #print negative_path
     
     


     value_for_positive_path=[]
     for x in positive_path:
       name=nx.get_node_attributes(G,'name_label')
       value_for_positive_path.append(name[x])

     #print value_for_positive_path

     value_for_negative_path=[]
     for x in negative_path:
       name=nx.get_node_attributes(G,'name_label')
       value_for_negative_path.append(name[x])
     #print value_for_negative_path
     
    
     
     if(information=='target'):
      array_for_target_bio2rdf=[]
      for x in range(1,len(value_for_positive_path)):
       #print value_for_positive_path[x].split('/')[0];
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?target_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasTarget ?targettemplate.
        ?targettemplate pibas:targetSchema ?target_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        target_schema=result["target_schema"]["value"]
        array_for_target_bio2rdf.append(target_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_target_bio2rdf
      array_for_target_chem2bio2rdf=[]
      for x in range(1,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?target_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasTarget ?targettemplate.
        ?targettemplate pibas:targetSchema ?target_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        target_schema=result["target_schema"]["value"]
        array_for_target_chem2bio2rdf.append(target_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_target_chem2bio2rdf

      bio2rdf=''
      if len(array_for_target_bio2rdf)>0:
       for x in range(0,len(array_for_target_bio2rdf)):
        bio2rdf=bio2rdf+'UNION{'+array_for_target_bio2rdf[x]+'}'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_target_chem2bio2rdf)>0:
       for x in range(0,len(array_for_target_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'UNION{'+array_for_target_chem2bio2rdf[x]+'}'
       #print chem2bio2rdf

      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?target
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>         
         { {   {  pibas:%s pibas:sameAs ?mapping_node;
                           pibas:hasTarget ?target.
         }

         %s %s
         
        
         }
         }
      """ % (value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf))  
      
      my_result='target-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["target"]["value"]
      string_for_query="PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#> \n PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?target \n FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl> \n { \n { \n  { \n  pibas:%s pibas:sameAs ?mapping_node; \n pibas:hasTarget ?target.\n } \n %s  \n %s \n  } \n }" %(value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf)
      my_result=my_result+' Query^'+string_for_query
      print my_result




     if(information=='cellline'):
      array_for_cellline_bio2rdf=[]
      for x in range(1,len(value_for_positive_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?cellline_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasCellLine ?celllinetemplate.
        ?celllinetemplate pibas:cellLineSchema ?cellline_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        cellline_schema=result["cellline_schema"]["value"]
        array_for_target_bio2rdf.append(cellline_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_cellline_bio2rdf
      array_for_cellline_chem2bio2rdf=[]
      for x in range(1,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?cellline_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasCellLine ?celllinetemplate.
        ?celllinetemplate pibas:cellLineSchema ?cellline_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        cellline_schema=result["cellline_schema"]["value"]
        array_for_cellline_chem2bio2rdf.append(cellline_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_cellline_chem2bio2rdf
      
      bio2rdf=''
      if len(array_for_cellline_bio2rdf)>0:
       for x in range(0,len(array_for_cellline_bio2rdf)):
        bio2rdf=bio2rdf+'UNION{'+array_for_cellline_bio2rdf[x]+'}'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_cellline_chem2bio2rdf)>0:
       for x in range(0,len(array_for_cellline_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'UNION{'+array_for_cellline_chem2bio2rdf[x]+'}'
       #print chem2bio2rdf
      
      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?cellline
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>         
         { {   {  pibas:%s pibas:sameAs ?mapping_node;
                           pibas:hasCellLine ?cellline.
         }

         %s %s
         
         }
         }
      """ % (value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf))  
      
      my_result='cellline-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["cellline"]["value"]
      string_for_query="PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#> \n PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?cellline \n FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl> \n { \n { \n {  pibas:%s pibas:sameAs ?mapping_node; \n pibas:hasCellLine ?cellline. \n } \n %s \n  %s \n  } \n}" %(value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf)
      my_result=my_result+' Query^'+string_for_query
      print my_result

     if(information=='ic50'):
      array_for_ic50_bio2rdf=[]
      for x in range(1,len(value_for_positive_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?ic50_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasIC50value ?ic50template.
        ?ic50template pibas:IC50valueSchema ?ic50_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        ic50_schema=result["ic50_schema"]["value"]
        array_for_ic50_bio2rdf.append(ic50_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_ic50_bio2rdf
      array_for_ic50_chem2bio2rdf=[]
      for x in range(1,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?ic50_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasIC50value ?ic50template.
        ?ic50template pibas:IC50valueSchema ?ic50_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        ic50_schema=result["ic50_schema"]["value"]
        array_for_ic50_chem2bio2rdf.append(ic5_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_ic50_chem2bio2rdf

      bio2rdf=''
      if len(array_for_ic50_bio2rdf)>0:
       for x in range(0,len(array_for_ic50_bio2rdf)):
        bio2rdf=bio2rdf+'UNION{'+array_for_ic50_bio2rdf[x]+'}'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_ic50_chem2bio2rdf)>0:
       for x in range(0,len(array_for_ic50_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'UNION{'+array_for_ic50_chem2bio2rdf[x]+'}'
       #print chem2bio2rdf

      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?cellline ?IC50value
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>         
         { {   {  pibas:%s pibas:sameAs ?mapping_node;
                           pibas:IC50value ?IC50value;
                           pibas:hasCellLine ?cellline.
         }

         %s %s
        
         }
         }
      """ % (value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf))  
      
      my_result='ic50-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["cellline"]["value"]
        my_result=my_result+' '+result["IC50value"]["value"]
      string_for_query="PREFIX pibas:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#> \n PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?cellline ?IC50value \n FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl> \n { \n { \n {  pibas:%s pibas:sameAs ?mapping_node; \n pibas:IC50value ?IC50value; pibas:hasCellLine ?cellline. \n } %s \n %s \n  } \n }" %(value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf)
      my_result=my_result+' Query^'+string_for_query
      print my_result

  if (("pibas_compound_acronym" not in result) and ("chembl_compound_acronym" in result)):
   chembl_acronym=result["chembl_compound_acronym"]["value"]
  
   import create_completed_graph_chembl_pibas
   #create_completed_graph.completed_graph(inchikey)   
   result=create_completed_graph_chembl_pibas.completed_graph(inchikey)
   fv=result['fv']
   if(type(fv) == str):
    print 'fv-'+str(fv)
   else:
    fv_around=np.around(fv,decimals=3)
    print 'fv-'+str(fv_around)
    connection_node=result['connection_node']
    import create_noncompleted_graph_chembl_pibas
    #create_noncompleted_graph.non_completed_graph(inchikey,connection_node)
   
    result1=create_noncompleted_graph_chembl_pibas.non_completed_graph(inchikey,connection_node)
    pr=result1['pr']
    if(type(pr) == str):
     print 'pr-'+str(pr)
    else:
     pr_around = {k:round(v,2) for k, v in pr.items()}
     print 'pr-'+str(pr_around)
     fv = fv_around   
     pg=pr

     G = result1['di_graph']
     nodes = G.nodes()	 
     print result1['node']
     positive_path = [max((t for t in pg if fv[nodes.index(t)] > 0 and abs(fv[nodes.index(t)]) > 0.0000001), key=pg.get)]
     negative_path = [result1['node']] 
     neutral_node = []
     visited = []

     # cvorovi odredjuju raspored u vektorima
     nodes = G.nodes()
     # positive vertices
     if (len(positive_path) > 0) :
      while positive_path[-1] not in visited:
        visited.append(positive_path[-1])
        neigh = G.neighbors(positive_path[-1])
        #print neigh		
        if not neigh:
            break
        all_ranks = []
        for nvalue in neigh:
         all_ranks.append(pg[nvalue])
        najveci = max(all_ranks)	
	flag = 0
	v_pos = -1
	
	for nvalue in neigh:
	    v_pos = nodes.index(nvalue)
	    if pg[nvalue] == najveci and abs(fv[v_pos]) < 0.0000001:
		positive_path.append(nvalue)
		flag = 1
		break
	    elif pg[nvalue] == najveci:
		v_pos = nvalue
			
	if flag == 1:
	    break
	
	if v_pos != -1:
		positive_path.append(v_pos)

     #print positive_path       
     

     # negative vertices
     if (len(negative_path) > 0) :
      while negative_path[-1] not in visited:
        visited.append(negative_path[-1])
        neigh = G.neighbors(negative_path[-1])
        #print neigh
        if not neigh:
            break
        all_ranks = []
        for nvalue in neigh:
           all_ranks.append(pg[nvalue])
        najveci = max(all_ranks)	
	flag = 0
	v_pos = -1
        
        for nvalue in neigh:
          v_pos = nodes.index(nvalue)
          if pg[nvalue] == najveci and abs(fv[v_pos]) < 0.0000001:
            negative_path.append(nvalue)
	    flag = 1
	    break
          elif pg[nvalue] == najveci:
	    v_pos = nvalue

        if flag == 1:
		break
	
	if v_pos != -1:
		negative_path.append(v_pos)
     
     
     #print negative_path
     
     


     value_for_positive_path=[]
     for x in positive_path:
       name=nx.get_node_attributes(G,'name_label')
       value_for_positive_path.append(name[x])

     #print value_for_positive_path

     value_for_negative_path=[]
     for x in negative_path:
       name=nx.get_node_attributes(G,'name_label')
       value_for_negative_path.append(name[x])
     #print value_for_negative_path
     
    
     
     if(information=='target'):
      array_for_target_bio2rdf=[]
      for x in range(1,len(value_for_positive_path)):
       #print value_for_positive_path[x].split('/')[0];
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?target_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasTarget ?targettemplate.
        ?targettemplate pibas:targetSchema ?target_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        target_schema=result["target_schema"]["value"]
        array_for_target_bio2rdf.append(target_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_target_bio2rdf
      array_for_target_chem2bio2rdf=[]
      for x in range(1,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?target_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasTarget ?targettemplate.
        ?targettemplate pibas:targetSchema ?target_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        target_schema=result["target_schema"]["value"]
        array_for_target_chem2bio2rdf.append(target_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_target_chem2bio2rdf

      bio2rdf=''
      if len(array_for_target_bio2rdf)>0:
       for x in range(0,len(array_for_target_bio2rdf)):
        bio2rdf=bio2rdf+'UNION{'+array_for_target_bio2rdf[x]+'}'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_target_chem2bio2rdf)>0:
       for x in range(0,len(array_for_target_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'UNION{'+array_for_target_chem2bio2rdf[x]+'}'
       #print chem2bio2rdf

      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?target
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl>
         { {   
         SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/>
           { 
             ?activity a cco:Activity ;
             cco:hasMolecule chembl_molecule:%s;
                             cco:hasAssay ?assay.
             ?assay cco:hasTarget ?target.
           }
          }
         %s %s
         
         
         }
      """ % (value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf))  
      
      my_result='target-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["target"]["value"]
      string_for_query="PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?target \n  \n { \n { \n  { \n SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/> { ?activity a cco:Activity ; cco:hasMolecule chembl_molecule:%s; cco:hasAssay ?assay. ?assay cco:hasTarget ?target.  \n } \n } \n %s  \n %s \n  } \n }" %(value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf)
      my_result=my_result+' Query^'+string_for_query
      print my_result




     if(information=='cellline'):
      array_for_cellline_bio2rdf=[]
      for x in range(0,len(value_for_positive_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?cellline_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasCellLine ?celllinetemplate.
        ?celllinetemplate pibas:cellLineSchema ?cellline_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        cellline_schema=result["cellline_schema"]["value"]
        array_for_target_bio2rdf.append(cellline_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_cellline_bio2rdf
      array_for_cellline_chem2bio2rdf=[]
      for x in range(0,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?cellline_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasCellLine ?celllinetemplate.
        ?celllinetemplate pibas:cellLineSchema ?cellline_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        cellline_schema=result["cellline_schema"]["value"]
        array_for_cellline_chem2bio2rdf.append(cellline_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_cellline_chem2bio2rdf
      
      bio2rdf=''
      if len(array_for_cellline_bio2rdf)>0:
       for x in range(0,len(array_for_cellline_bio2rdf)):
        bio2rdf=bio2rdf+'UNION{'+array_for_cellline_bio2rdf[x]+'}'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_cellline_chem2bio2rdf)>0:
       for x in range(0,len(array_for_cellline_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'UNION{'+array_for_cellline_chem2bio2rdf[x]+'}'
       #print chem2bio2rdf
      
      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?target
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl>
         { {   
         SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/>
           { 
              ?act ?o chembl_molecule:%s .
              ?act cco:standardUnits "nM" .
              ?act cco:standardType "IC50" .
              ?act cco:standardValue ?IC50value .
              ?assay cco:hasActivity ?act .
              ?assay cco:hasCellLine ?cellline.
           }
          }
         %s %s
         
         
         }
      """ % (value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf))  
      
      my_result='cellline-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["cellline"]["value"]
      string_for_query="PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?celline \n  \n { \n { \n  { \n SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/> { ?act ?o chembl_molecule:%s . ?act cco:standardUnits 'nM'. ?act cco:standardType 'IC50'. ?act cco:standardValue ?IC50value. ?assay cco:hasActivity ?act. ?assay cco:hasCellLine ?cellline.  \n } \n } \n %s  \n %s \n  } \n }" %(value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf)
      my_result=my_result+' Query^'+string_for_query
      print my_result

     if(information=='ic50'):
      array_for_ic50_bio2rdf=[]
      for x in range(0,len(value_for_positive_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?ic50_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-bio2rdf'.
        ?dataset pibas:hasIC50value ?ic50template.
        ?ic50template pibas:IC50valueSchema ?ic50_schema.

        } 
 
        """% value_for_positive_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        ic50_schema=result["ic50_schema"]["value"]
        array_for_ic50_bio2rdf.append(ic50_schema % value_for_positive_path[x].split('/')[1])
      #print array_for_ic50_bio2rdf
      array_for_ic50_chem2bio2rdf=[]
      for x in range(0,len(value_for_negative_path)):
       sparql.setQuery(
        """  
        PREFIX pibas: <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS#>
        SELECT ?ic50_schema
        FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/OntologyIntegration.owl>
        WHERE{
        ?dataset pibas:datasetName '%s-chem2bio2rdf'.
        ?dataset pibas:hasIC50value ?ic50template.
        ?ic50template pibas:IC50valueSchema ?ic50_schema.

        } 
 
        """% value_for_negative_path[x].split('/')[0])
       sparql.setReturnFormat(JSON)
       final_results = sparql.query().convert()
       for result in final_results["results"]["bindings"]:
        ic50_schema=result["ic50_schema"]["value"]
        array_for_ic50_chem2bio2rdf.append(ic5_schema % value_for_negative_path[x].split('/')[1])
       #print array_for_ic50_chem2bio2rdf

      bio2rdf=''
      if len(array_for_ic50_bio2rdf)>0:
       for x in range(0,len(array_for_ic50_bio2rdf)):
        bio2rdf=bio2rdf+'UNION{'+array_for_ic50_bio2rdf[x]+'}'
       #print bio2rdf
      
      chem2bio2rdf=''
      if len(array_for_ic50_chem2bio2rdf)>0:
       for x in range(0,len(array_for_ic50_chem2bio2rdf)):
        chem2bio2rdf=chem2bio2rdf+'UNION{'+array_for_ic50_chem2bio2rdf[x]+'}'
       #print chem2bio2rdf

      sparql = SPARQLWrapper("http://cpctas-lcmb.pmf.kg.ac.rs:2020/sparql")

      sparql.setQuery(
      """
         PREFIX drugbank:<http://bio2rdf.org/drugbank:>
         PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/>
         PREFIX kegg_ligand:<http://bio2rdf.org/kegg:>
         PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#>
         PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/>
         PREFIX cco:  <http://rdf.ebi.ac.uk/terms/chembl#>
         SELECT  DISTINCT ?target
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/pibasmapping.owl>
         FROM <http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/PIBAS/chemblmapping.owl>
         { {   
         SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/>
           { 
              ?act ?o chembl_molecule:%s .
              ?act cco:standardUnits "nM" .
              ?act cco:standardType "IC50" .
              ?act cco:standardValue ?IC50value .
              ?assay cco:hasActivity ?act .
              ?assay cco:hasCellLine ?cellline.
           }
          }
         %s %s
         
         
         }
      """ % (value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf))
      my_result='ic50-'
      sparql.setReturnFormat(JSON)
      final_results = sparql.query().convert()
      for result in final_results["results"]["bindings"]: 
        my_result=my_result+' '+result["cellline"]["value"]
        my_result=my_result+' '+result["IC50value"]["value"]
      string_for_query="PREFIX drugbank:<http://bio2rdf.org/drugbank:> \n PREFIX drugbank1:<http://chem2bio2rdf.org/drugbank/resource/drugbank_drug/> \n PREFIX kegg_ligand:<http://bio2rdf.org/kegg:> \n PREFIX chembl_mapp:<http://cpctas-lcmb.pmf.kg.ac.rs/2012/3/chembl#> \n PREFIX chembl_molecule:<http://rdf.ebi.ac.uk/resource/chembl/molecule/> \n PREFIX cco:<http://rdf.ebi.ac.uk/terms/chembl#> \n SELECT  DISTINCT ?cellline ?IC50value \n  \n { \n { \n  { \n SERVICE SILENT<https://www.ebi.ac.uk/rdf/services/sparql/> { ?act ?o chembl_molecule:%s . ?act cco:standardUnits 'nM'. ?act cco:standardType 'IC50'. ?act cco:standardValue ?IC50value. ?assay cco:hasActivity ?act. ?assay cco:hasCellLine ?cellline.  \n } \n } \n %s  \n %s \n  } \n }" %(value_for_positive_path[0].split('/')[1],bio2rdf,chem2bio2rdf)
      my_result=my_result+' Query^'+string_for_query
      print my_result

    
    
  
 
   
  
   
 