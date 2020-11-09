# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 16:48:54 2020

@author: BALAMLAPTOP2
"""

import networkx as nx
from networkx.algorithms import bipartite
import collections
import pandas as pd
import sys
import multiprocessing
import time

def threshold_analysis_ws(C, th, type_proj, nn, iteration, attr):
    print("Creating new graph by given threshold ... "+str(th))
    H = nx.Graph()
    index_source = 1
    index_target = 1
    for n in C.nodes(data=True):
        if n[1][attr] == type_proj:
            #print(index_source)
            sourceNode = n[0]
            s_neighbors = set(C.neighbors(n[0]))
            for m in C.nodes(data = True):
                if m[1][attr] == type_proj: #### Change to 1 to change the projection to active ingredient
                    targetNode = m[0]
                    t_neighbors = set(C.neighbors(m[0]))
                    if sourceNode != targetNode and index_target > index_source:
                        if len(s_neighbors & t_neighbors) >= th:
                            H.add_node(sourceNode)
                            H.add_node(targetNode)
                            H.add_edge(sourceNode,targetNode)
                    index_target += 1
            index_target = 1
            index_source += 1        					
    components = sorted(nx.connected_components(H), key=len, reverse=True)
    nodes_connected = sum(list(map(lambda c: len(c), components)))
    mean_size_components = nodes_connected / len(components)
    nodes_unconnected = nn - nodes_connected
    lcs = len(components[0])
    degrees = H.degree()
    sum_of_edges = sum(list(dict(degrees).values()))
    avg_degree = sum_of_edges / H.number_of_nodes()
    
    print("Saving values for the given threshold ..."+str(th))
    if type_proj == 0:
        with open("threshold_shuffle_icd_ws_"+str(iteration)+".txt", "a+") as f:
            f.write(str(th)+","+str(len(components))+","+str(nodes_connected)+","+str(nodes_unconnected)+","+str(lcs)+","+str(avg_degree)+","+str(mean_size_components)+"\n")
        # nx.write_graphml(H,'ICD/projICD_th_'+str(th)+'.graphml')
    elif type_proj == 1:
        with open("threshold_shuffle_atc_ws_"+str(iteration)+".txt", "a+") as f:
            f.write(str(th)+","+str(len(components))+","+str(nodes_connected)+","+str(nodes_unconnected)+","+str(lcs)+","+str(avg_degree)+","+str(mean_size_components)+"\n")
        # nx.write_graphml(H,'ATC/projATC_th_'+str(th)+'.graphml')
    else:
        print("The option doesn't exist. Try again.")

if __name__ == '__main__':
    
    # Remover al azar 0 o dirigido 1
    # iteration number of the shuffled file 
    # type_proj ICD 0 o ATC 1 
    # attr d0 read from shuffle data and bipartite read from original data
    # type_method 0 original and 1 shuffle
    type_method = int(sys.argv[1])
    # type_proj = int(sys.argv[2])
    
    p = multiprocessing.Pool()
        #timing it...
    start = time.time()
    
    if type_method == 0: # ORIGINAL DATA
        vdmdata = pd.read_csv('vdmdata_reduce.csv', encoding = 'utf-8-sig')
        #vdmdata.columns = ['icd_code', 'atc_code','atc_name','nrows']
        
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
            
        
        components = sorted(nx.connected_components(G), key=len, reverse=True)
        largest_component = components[0]
        C = G.subgraph(largest_component)
            
        degX,degY=bipartite.degrees(C,nodes_0)
        degATC = dict(degX).values()
        degCIE = dict(degY).values()
        counterATC = collections.Counter(degATC)
        counterCIE = collections.Counter(degCIE)
        
        print("Apply threshold analysis to original bipartite ... ")
        for th_icd in sorted(list(counterCIE.keys())):
            p.apply_async(threshold_analysis_ws, [C, th_icd, 0, len(degY), 'original', 'bipartite'])
            
        for th_atc in sorted(list(counterATC.keys())):
            p.apply_async(threshold_analysis_ws, [C, th_atc, 1, len(degX), 'original', 'bipartite'])
        
        
    else: # SHUFFLE DATA
        
        for i in range(1,4,1):
            print("Read data ... iteration "+str(i))
            edges = pd.read_csv('edges_'+str(i)+'.csv', sep=',')
            G = nx.from_pandas_edgelist(edges, 'Source', 'Target', edge_attr=True)
            
            nodes = pd.read_csv('nodes_'+str(i)+'.csv', sep=',')
            data = nodes.set_index('Id').to_dict('index').items()
            G.add_nodes_from(data)
            # print(G.nodes(data=True))
            # print(G.edges(data=True))
            
            nodes_0 = [x for x,y in G.nodes(data=True) if y['d0']==0]
            nodes_1 = [x for x,y in G.nodes(data=True) if y['d0']==1]
            
            degX,degY=bipartite.degrees(G,nodes_0)
            degATC = dict(degX).values()
            degCIE = dict(degY).values()
            counterATC = collections.Counter(degATC)
            counterCIE = collections.Counter(degCIE)
            
            print("Apply threshold analysis to original bipartite ... ")
            for th_icd in sorted(list(counterCIE.keys())):
                p.apply_async(threshold_analysis_ws, [G, th_icd, 0, len(degY), i, 'd0'])
                
            for th_atc in sorted(list(counterATC.keys())):
                p.apply_async(threshold_analysis_ws, [G, th_atc, 1, len(degX), i, 'd0'])
            
    p.close()
    p.join()
    print("Complete")
    end = time.time()
    print('total time (s)= ' + str(end-start))
    
########################################################################### funcional
    # type_proj = 1
    # iteration = 0
    
    
    # for th in sorted(list(counterATC.keys())):
    #     # threshold_analysis_ws(G, th_atc, 1, len(degX), 0)
    #     print("Creating new graph by given threshold ... "+str(th))
    #     H = nx.Graph()
    #     # index_source = 1
    #     # index_target = 1
    #     for n in C.nodes(data=True):
    #         if n[1]['bipartite'] == 1: #d0
    #             #print(index_source)
    #             sourceNode = n[0]
    #             s_neighbors = set(C.neighbors(n[0]))
    #             for m in C.nodes(data = True):
    #                 if m[1]['bipartite'] == 1: #### Change to 1 to change the projection to active ingredient
    #                     targetNode = m[0]
    #                     t_neighbors = set(C.neighbors(m[0]))
    #                     if sourceNode != targetNode:
    #                     # if sourceNode != targetNode and index_target > index_source:
    #                         if len(s_neighbors & t_neighbors) >= th:
    #                             H.add_node(sourceNode)
    #                             H.add_node(targetNode)
    #                             H.add_edge(sourceNode,targetNode)
    #             #         index_target += 1
    #             # index_target = 1
    #             # index_source += 1        					
    #     components = sorted(nx.connected_components(H), key=len, reverse=True)
    #     nodes_connected = sum(list(map(lambda c: len(c), components)))
    #     mean_size_components = nodes_connected / len(components)
    #     nodes_unconnected = len(degX) - nodes_connected
    #     lcs = len(components[0])
    #     degrees = H.degree()
    #     sum_of_edges = sum(list(dict(degrees).values()))
    #     avg_degree = sum_of_edges / H.number_of_nodes()
        
    #     print("Saving values for the given threshold ..."+str(th))
    #     if type_proj == 0:
    #         with open("threshold_shuffle_icd_ws_"+str(iteration)+".txt", "a+") as f:
    #             f.write(str(th)+","+str(len(components))+","+str(nodes_connected)+","+str(nodes_unconnected)+","+str(lcs)+","+str(avg_degree)+","+str(mean_size_components)+"\n")
    #         # nx.write_graphml(H,'ICD/projICD_th_'+str(th)+'.graphml')
    #     elif type_proj == 1:
    #         with open("threshold_shuffle_atc_ws_"+str(iteration)+".txt", "a+") as f:
    #             f.write(str(th)+","+str(len(components))+","+str(nodes_connected)+","+str(nodes_unconnected)+","+str(lcs)+","+str(avg_degree)+","+str(mean_size_components)+"\n")
    #         # nx.write_graphml(H,'ATC/projATC_th_'+str(th)+'.graphml')
    #     else:
    #         print("The option doesn't exist. Try again.")
###########################################################################
    # degX,degY=bipartite.degrees(C,nodes_0)
    # degATC = dict(degX).values()
    # degCIE = dict(degY).values()
    # counterATC = collections.Counter(degATC)
    # counterCIE = collections.Counter(degCIE)
    # c_list_cie = []
    # nc_list_cie = []
    # nuc_list_cie = []
    # for th in sorted(list(counterCIE.keys())):
    #     #th = 1
    #     H = nx.Graph()
    #     #for v in G.nodes(data = True):
    #     #    if v[1]['bipartite'] == 0:
    #     #        H.add_node(v[0])
    
    #     for n in C.nodes(data=True):
    #         if n[1]['bipartite'] == 0:
    #             sourceNode = n[0]
    #             s_neighbors = set(C.neighbors(n[0]))
    #             for m in C.nodes(data = True):
    #                 if m[1]['bipartite'] == 0: #### Change to 1 to change the projection to active ingredient
    #                     targetNode = m[0]
    #                     t_neighbors = set(C.neighbors(m[0]))
    #                     if sourceNode != targetNode:
    #                         if len(s_neighbors & t_neighbors) >= th:
    #                             H.add_node(sourceNode)
    #                             H.add_node(targetNode)
    #                             H.add_edge(sourceNode,targetNode)                  
    #     components = sorted(nx.connected_components(H), key=len, reverse=True)
    #     #sum(list(map(lambda c: len(c), components)))i
    #     c_list_cie.append(len(components))
    #     nodes_connected = sum(list(map(lambda c: len(c), components)))
    #     nc_list_cie.append(nodes_connected)
    #     nuc_list_cie.append(len(nodes_0) - nodes_connected)
    #     nx.write_graphml(H,'ICD/projICD_th_'+str(th)+'.graphml')
    
              #H = add_and_remove_edges(C, type_proj, dict(degX), dict(degY))
    #         C, counterATC, counterCIE, degX, degY, nodes_0_c, nodes_1_c = suffle_edges_lc(G)
    # #        suffle_edges(C, sorted(dict(degX).keys()), sorted(dict(degY).keys()))
    #         degX_sh,degY_sh=bipartite.degrees(C,nodes_0_c)
    #         degATC_sh = dict(degX_sh).values()
    #         degCIE_sh = dict(degY_sh).values()
    #         counterATC_sh = collections.Counter(degATC_sh)
    #         counterCIE_sh = collections.Counter(degCIE_sh)
    #         nx.write_graphml(C,'networks/bipartite_sh_'+str(i)+'.graphml')
        
    
    
    # G = nx.read_graphml('graph_shufflenx/bipartite_sh_1.graphml', force_multigraph=True)
    
        
        # nx.write_graphml(C,'networks/bipartite_sh_'+str(i)+'.graphml')
        
        # print("Apply threshold analysis to shuffled graph ... "+str(i))


