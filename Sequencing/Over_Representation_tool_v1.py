# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 16:49:14 2018

@author: Eric
"""


import datetime
import sys
from subprocess import call
import csv
from contextlib import contextmanager
import os
import random
import scipy.stats as stats

#This tool has limits, for instance it does not take hierarchy or dependence of genes sets into account
#advantages: it is data driven, based on permutation tests, so tissue types don't get certain gene sets automatically overrepresented
#it would be interesting to correct for genes in a same pathway with possibly antagonistic functions ie: a biosynthesis enzyme in same pathway as the eliminating enzyme
#this should be amended to include more current methods such as tools in the Self-Contained algorithms, N-statistic, or others...
#

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


#Parse GMT file by line, and create DICT with pathway as index with list of genes as value
def read_gmt(fh):
    bank = {}
    for line in fh:
        array1 = line.rstrip().split("\t")
        bank[array1[0]] = True
        bank[array1[0]] = list(array1[2:len(array1)])
    return bank

#function determining the expected value of random sampling of our set
#Monte-carlo type re-sampling to determine median based on our dataset
#Re-sample a set from experimental data of the tested set size, and compare to initial sample of same size
#10000 iterations, return expected value from random sampling. 
def expected(experimental, gene_sets_len):
    samp = random.sample(experimental, gene_sets_len)
    sum_of_hits = 0
    Theor_max = 10000*gene_sets_len
    x=0
    while x<10000:
        #can I get mean + error from montecarlo sim?
        #If so, I could base it off the error
        simul = random.sample(experimental, gene_sets_len)
        for item in simul:
            if item in samp:
                sum_of_hits += 1                
        x+=1
    return float(float(sum_of_hits)/float(Theor_max))
    
def enrichments(experimental, bank, outname):
    #determine discrepancy from expected
    wf = open(outname, 'w')    
    for sets in bank:
        N = len(sets)
        M = 0
        expect = expected(experimental, N)
        for gene in bank[sets]:
            if gene in experimental:
                M +=1
        #calculated score is simply the observed number of genes - computed expected value        
        score = float(float(M)-expect)
        Od,Pv = stats.fisher_exact([[M,int(round(expect))],[N,N]])
        #return the pathway name, expected score, real score, Pvalue of fishers exact test, and newline
        #Needs a BH correction?
        combi = (sets,expect, score, Pv, "\n")
        combi = '\t'.join(str(i) for i in combi)
        wf.write(combi)

        
        
def Main():
    assert (len(sys.argv)==4),"Usage: Python This_script Path_DB path_to_geneSet OutFileName    DB should be .GMT format    geneset one per line"
    script = sys.argv[0]
    db = str(sys.argv[1])
    Exp = str(sys.argv[2])
    outname = str(sys.argv[3])
    
    #testing
#    db = "C:\\Users\\Eric\\Documents\\Programming\\Test\\PWAT\\ReactomePathways.gmt"
#    Exp = "C:\\Users\\Eric\\Documents\\Programming\\Test\\PWAT\\Test_Set.txt"
#    outname = "OutFile_test1.txt"
    
    bank = read_gmt(open(db))
    
    with open(Exp) as fh:
        result = fh.readlines()
    result = [x.strip() for x in result]
    
#    call(['mkdir','HGT_out'])
    
#    with cd('HGT_output'):
#        #The body of the code goes here
#        #some sort of writer
#        enrich_result = enrichments(result, bank)
#        wf = open(outname, 'w')
#        for item in enrich_result:
#            wf.write("%s\n" % item)


    #testing
    enrichments(result, bank, outname)
    
Main()