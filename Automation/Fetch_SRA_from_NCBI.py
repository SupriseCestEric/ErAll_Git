# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 10:34:20 2017

@author: Eric
"""
#This script requires a download_config.txt file with one line, comma 
seperated SRR identifiers. It will download all SRA files and dump fastq 
files. 

from subprocess import call
from csv import reader
csv_reader = reader(open("download_config.txt","r"), quotechar="\"")

list1 =[]
for row in csv_reader:
    list1.append(row)

for item in list1[0]:
    y="ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra" %(item[0:3],item[0:6],item,item) 
    print y
    call(['wget',y])

#Change the path below to your path to fastq-dump file.

for item in list1[0]:
    x="%s.sra" %item
    call(['/Path_to/sratoolkit.2.8.1-3-mac64/bin/fastq-dump','--split-3',x])
