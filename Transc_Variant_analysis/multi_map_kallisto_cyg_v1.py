# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 10:34:20 2017

@author: Eric
"""


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

def makeConfigFile(direct):
    with cd(direct):
        call('ls *_UP.fq.gz > Trimmed_filenames.csv', shell=True)

def Kallisto():
    assert (len(sys.argv)==3),"Usage: Python This_script path_to_trimmed_files path_to_gff_index \
    This script does mapping of many fastq files to the provided annotation and generates estimated counts with kallisto\
    You should therefore have kallisto installed on your machine for this script to work.\
    please see https://pachterlab.github.io/kallisto/"
    script = sys.argv[0]
    file_path = str(sys.argv[1])
    gff = str(sys.argv[2])

    makeConfigFile(str(file_path))
    ConfigFile=file_path+"\\Trimmed_filenames.csv"
    csv_reader = reader(open(ConfigFile,"r"), quotechar="\"")

    list_files =[]
    for row in csv_reader:
        list_files.append(row)

    list_files.sort
    
    call(['mkdir','Kallisto_out'])
    with cd('Kallisto_out'):
        x=0
        while x<len(list_files):
            if list_files[x][0][-11:] == "FW_UP.fq.gz":
                FW=list_files[x][0]
                RV=list_files[x+1][0]
                pFW=file_path+"\\"+FW
                pRV=file_path+"\\"+RV
                dirname=list_files[x][0][0:10]
                call(['kallisto','quant','-i',gff,'--threads','4','-o','.\\%s' % dirname,pFW,pRV])
                x+=2
            else:
                break
    
    with open("headers_kallisto.csv", "wb") as csvfile:
        writer = csv.writer(csvfile, delimiter = " ", quotechar="\"")
        for item in list_files:
            writer.writerow([item])      
    call(['mv','headers_kallisto.csv','Kallisto_out'])

Kallisto()
