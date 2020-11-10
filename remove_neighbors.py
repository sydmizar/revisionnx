# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 17:04:55 2020

@author: BALAMLAPTOP2
"""
import networkx as nx
from networkx.algorithms import bipartite
import collections
import pandas as pd
import multiprocessing
import time
import random

def suffle_edges_lc(G, iteration, degY):
    print("Getting largest component ...")
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    largest_component = components[0]
    C = G.subgraph(largest_component)
        
    degX,degY=bipartite.degrees(C,nodes_0)
    degATC = dict(degX).values()
    degCIE = dict(degY).values()
    counterATC = collections.Counter(degATC)
    counterCIE = collections.Counter(degCIE)
    
    nodes_0_c = []
    nodes_1_c = []
    for n in C.nodes(data=True):
        if n[1]['bipartite'] == 0:
            nodes_0_c.append(n[0])
        if n[1]['bipartite'] == 1:
            nodes_1_c.append(n[0])
            
    print("Shuffling edges ... ")
    unfrozen_graph = nx.Graph(C)
    
    k=0
    iter=2*unfrozen_graph.size()
    while k<iter:
        r1=random.choice(sorted(dict(degY).keys()))
        d1=random.choice(list(unfrozen_graph.neighbors(r1)))
        
        r2=random.choice(sorted(dict(degY).keys()))
        d2=random.choice(list(unfrozen_graph.neighbors(r2)))
        
        if (unfrozen_graph.has_edge(r1,d2)==False) & (unfrozen_graph.has_edge(r2,d1)==False):
            unfrozen_graph.add_edge(r1,d2)
            unfrozen_graph.remove_edge(r1,d1)       
            unfrozen_graph.add_edge(r2,d1)
            unfrozen_graph.remove_edge(r2,d2)
#            print(k)
            k=k+1
            
    nx.write_graphml(unfrozen_graph,'shuffled_network_'+str(iteration)+'.graphml')
            
if __name__ == '__main__':
    # p = multiprocessing.Pool()
        #timing it...
    # start = time.time()
    
    vdmdata = pd.read_csv('vdmdata_reduce.csv', encoding = 'utf-8-sig')
    vdmdata.columns = ['icd_code', 'atc_code','atc_name','nrows']
    
    nodes_0 = []
    nodes_1 = []
    for m in vdmdata.iterrows():
        nodes_0.append(m[1][0]) #ICD
        nodes_1.append(m[1][1]) #ATC
        
    nodes_0 = list(dict.fromkeys(nodes_0))
    nodes_1 = list(dict.fromkeys(nodes_1))
    
    # Build a bipartite graph:
    G = nx.Graph()
    G.add_nodes_from(nodes_0, bipartite=0) # Add the node attribute “bipartite” disease
    G.add_nodes_from(nodes_1, bipartite=1) # active substance
    
    
    for m in vdmdata.iterrows():
        enfermedad = m[1][0];
        #peso = m[1][3];
        sustancia = m[1][1];
        G.add_edge(enfermedad, sustancia)
        
    print("Getting largest component ...")
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    largest_component = components[0]
    C = G.subgraph(largest_component)
        
    degX,degY=bipartite.degrees(C,nodes_0)
    degATC = dict(degX).values()
    degCIE = dict(degY).values()
    counterATC = collections.Counter(degATC)
    counterCIE = collections.Counter(degCIE)
    
    nodes_0_c = []
    nodes_1_c = []
    for n in C.nodes(data=True):
        if n[1]['bipartite'] == 0:
            nodes_0_c.append(n[0])
        if n[1]['bipartite'] == 1:
            nodes_1_c.append(n[0])
        
    df_icd = pd.DataFrame(dict(degY).items(), columns=['node', 'degree'])
    df_atc = pd.DataFrame(dict(degX).items(), columns=['node', 'degree'])
    
    df_icd = df_icd.sort_values(by=['degree', 'node'], ascending=False)
    df_atc = df_atc.sort_values(by=['degree', 'node'], ascending=False)
    
    unfrozen_graph = nx.Graph(C)
    
    for i in range(0,100,1):
        n1 = df_atc.iloc[i][0]
        #n2 = list(unfrozen_graph.neighbors(n1))
        n2 = df_icd[df_icd['node'].isin(list(unfrozen_graph.neighbors(n1)))].iloc[0][0]
        print(n1+' '+n2)
        unfrozen_graph.remove_edge(n1,n2)
    
    k=0
    iter=2*unfrozen_graph.size()
    while k<iter:
        r1=random.choice(sorted(dict(degY).keys()))
        d1=random.choice(list(unfrozen_graph.neighbors(r1)))
        
        r2=random.choice(sorted(dict(degY).keys()))
        d2=random.choice(list(unfrozen_graph.neighbors(r2)))
        
        if (unfrozen_graph.has_edge(r1,d2)==False) & (unfrozen_graph.has_edge(r2,d1)==False):
            unfrozen_graph.add_edge(r1,d2)
            unfrozen_graph.remove_edge(r1,d1)       
            unfrozen_graph.add_edge(r2,d1)
            unfrozen_graph.remove_edge(r2,d2)
#            print(k)
            k=k+1
            
    nx.write_graphml(unfrozen_graph,'shuffled_network_remove.graphml')
    
    
        