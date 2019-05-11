# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 09:29:23 2017

@author: Eric
"""

## read csv file containing SRR numbers, one per line
import sys
from subprocess import call
from csv import reader
from contextlib import contextmanager
import os

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
## loop through the csv file and call trimmomatic on each one

def makeConfigFile(direct):
    with cd(direct):
        call('ls *.fastq.gz > filenames.csv', shell=True)

def Trim():
    assert (len(sys.argv)==3),"Usage: PreProcess.py fastq_file_path adapters_to_clip.fa \
    This script takes the three arguments listed above. Clipping is done with trimmomatic, \
    the .jar of trimmomatic and adapter fasta must be in the same directory as this script. \
    Output will be sent to Trim_output directory, which will be a sub-directory of the current working directory. \
    If any Trimmomatic params must be modified, it should be modifies directly from the source code. \
    Defaults are currently: paired-end, phred33, sliding window of 4 with min quality of 20. \
    Please do not enter full path to adapter file"
    script = sys.argv[0]
    file_path = str(sys.argv[1])
    
    makeConfigFile(file_path)
    files = file_path+"\\filenames.csv"

    adapt_file = str(sys.argv[2])

    csv_reader = reader(open(files,"r"), quotechar="\"")

    list_files =[]
    for row in csv_reader:
        list_files.append(row)

    list_files.sort

    
    call(['mkdir', 'Trim_output'])
    with cd("Trim_output"): 
        x=0
        while x<len(list_files):
            if list_files[x][0][-5:] == "fastq" or list_files[x][0][-8:] == "fastq.gz":
                FW=list_files[x][0]
                RV=list_files[x+1][0]
                pFW=file_path+"\\"+FW #Path to forward read
                print(pFW) #test
                pRV=file_path+"\\"+RV #Path to Reverse read
                FUP='%s_FW_UP.fq.gz' % FW #output file names (unpaired and paired)
                FP='%s_FW_P.fq.gz' % FW
                RUP='%s_RV_UP.fq.gz' % RV
                RP='%s_RV_P.fq.gz' % RV
                Adapt = 'ILLUMINACLIP:'+"..\\"+adapt_file+':1:30:10'
                print(Adapt) #test
                call(['java','-jar','..\\trimmomatic-0.36.jar','PE','-phred33','-threads','8',pFW,pRV,FUP,FP,RUP,RP,Adapt,'SLIDINGWINDOW:4:20','MINLEN:20'])
            x+=2

Trim()