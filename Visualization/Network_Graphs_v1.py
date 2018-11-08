# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 20:23:55 2018

@author: Eric
"""


import agglomcluster
from hac import GreedyAgglomerativeClusterer
import datetime
import numpy
from scipy.cluster import hierarchy
from scipy.spatial import distance
from collections import defaultdict
import sys
from subprocess import call
import csv
from contextlib import contextmanager
import os
import random
import scipy.stats as stats
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
#Simple network viewer based off a known network database 
#

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

#convert tab-sep interaction database to simple edges
def BIOGRID_to_edges(fh):
    edges = []
    for line in fh:
        array1 = line.rstrip().split("\t")
        edge1 = (array1[7],array1[8])
        edges.append(edge1)
    return edges   

#Parse network DB file line by line and save to edges


#get only relevant interactions for our Data
def extract_edges(bank, DataSet):
    bank2 = []
    for edge in bank:
        if edge[0] in DataSet and edge[1] in DataSet:
            bank2.append(edge)
    return bank2
#get direct interaction and 1st neighbors        
def extract_edges_first_neighbor(bank, DataSet):
    bank2 = []
    for edge in bank:
        if edge[0] in DataSet or edge[1] in DataSet:
            bank2.append(edge)
    return bank2


#Need to add some sort of clustering algorithm, and integrate path analysis.
#re-draw network after clustering (perhaps node attribute = cluster)    
    
#Create a network graph with only the genes in our dataset and 1st neighbors. 
def Network_graph(list_of_edges, outname):
    G = nx.Graph()
    G.add_edges_from(list_of_edges)
    pos = nx.spring_layout(G,k=0.5)
    dis = nx.degree(G)
    plt1 = plt.figure(figsize=(12,12),dpi=600)
    #here set labels to True or False to not overplot big networks
    nx.draw_networkx(G, with_labels = False, node_size = [v * 2 for v in dis.values()],node_color="yellow", font_size=7,edge_color = 'grey',alpha=0.5)
    plt1.savefig(outname,format="PNG")





def Main():
    assert (len(sys.argv)==4),"Usage: Python This_script Path_DB.tsv path_to_geneSet OutFileName.png"
    script = sys.argv[0]
    db = str(sys.argv[1])
    Exp = str(sys.argv[2])
    outname = str(sys.argv[3])
    
#    #testing
#    db = "C:\\Users\\Eric\\Documents\\Programming\\Test\\PWAT\\BIOGRID-ALL-3.4.160.tab2.txt"
#    Exp = "C:\\Users\\Eric\\Documents\\Programming\\Test\\PWAT\\Test_Set.txt"
#    outname = "OutFile_gaph.png"
#    outfile2 = "Blocks.png"
    
    bank = BIOGRID_to_edges(open(db))
    
    with open(Exp) as fh:
        result = fh.readlines()
    result = [x.strip() for x in result]
    #Edges = extract_edges(bank,result)
    Edges = extract_edges_first_neighbor(bank, result)

    Network_graph(Edges, outname)
    #testing
    #clustering(Edges,outfile2)
    
Main()