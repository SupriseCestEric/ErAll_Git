#!/usr/bin/env python3
## Created on 19-10-2020

#THIS PYTHON SCRIPT REINITIALIZES THE HD200 STANDARDS DATABASE
#RUNNING IT WILL ERASE ALL PREVIOUS DATA IN THE HD200 DB

#read-in vcf header
import os
import re
import numpy
import pandas
import vcf
import json
import sqlite3
from io import StringIO
import ast

conn = sqlite3.connect('HD200_analysis.sqlite')
cur = conn.cursor()

# Do some setup
# 7 Metadata tables with reference genomes, workflow info, disease info etc. will be created
#The Other variant info metadata tables with gene and chromosome info will also be created. 
# The main tables are CallData (genotypes), VarData (variants), and RunInfo (Samples)
cur.executescript('''

DROP TABLE IF EXISTS ref_genome;
DROP TABLE IF EXISTS AnnotationVersion;
DROP TABLE IF EXISTS WorkFlowName;
DROP TABLE IF EXISTS WorkFlowVer;
DROP TABLE IF EXISTS Chromosomes;
DROP TABLE IF EXISTS BaseCallVer;
DROP TABLE IF EXISTS VarType;
DROP TABLE IF EXISTS DisType;
DROP TABLE IF EXISTS RunInfo;
DROP TABLE IF EXISTS CallData;
DROP TABLE IF EXISTS VarData;
DROP TABLE IF EXISTS Genes;

CREATE TABLE ref_genome (
    id   INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name    VARCHAR UNIQUE
);

CREATE TABLE AnnotationVersion (
    id   INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name    VARCHAR UNIQUE
);

CREATE TABLE WorkFlowName (
    id   INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name    VARCHAR UNIQUE
);

CREATE TABLE WorkFlowVer (
    id   INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name    VARCHAR UNIQUE
);        

CREATE TABLE BaseCallVer (
    id   INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name    VARCHAR UNIQUE
);        

CREATE TABLE Chromosomes (
    id   INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name    VARCHAR UNIQUE
);        

CREATE TABLE VarType (
    id   INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name    VARCHAR UNIQUE
);           

 
CREATE TABLE DisType (
    id   INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name    VARCHAR UNIQUE
);        

CREATE TABLE Genes (
    id   INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name    VARCHAR UNIQUE
);        
        
        
CREATE TABLE RunInfo (
    id     INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name   VARCHAR UNIQUE,
    sex    VARCHAR,
    pheno    VARCHAR,
    fileformat    VARCHAR,
    filedate    INTEGER,
    reference    VARCHAR,
    OM_anno_version    VARCHAR,
    IonWF    VARCHAR,
    IonWF_version    VARCHAR,
    Cellularity    FLOAT,
    fusionOC    VARCHAR,
    fusionQC    VARCHAR,
    fusionReads    INTEGER,
    basecall_ver    VARCHAR,
    perc_mapped    FLOAT,
    FOREIGN KEY (reference)
        REFERENCES ref_genome (id),
    FOREIGN KEY (OM_anno_version)
        REFERENCES AnnotationVersion (id),
    FOREIGN KEY (IonWF)
        REFERENCES WorkFlowName (id),
    FOREIGN KEY (IonWF_version)
        REFERENCES WorkFlowVer (id),
    FOREIGN KEY (basecall_ver)
        REFERENCES BaseCallVer (id)
    FOREIGN KEY (pheno)
        REFERENCES DisType (id)
);

CREATE TABLE VarData (
    id    INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    name    TEXT UNIQUE,
    chr    INTEGER,
    start    INTEGER,
    len    VARCHAR,
    ref    VARCHAR,
    alt    VARCHAR,
    type    VARCHAR,
    annotation    VARCHAR,
    gene    VARCHAR,
    FOREIGN KEY (chr)
        REFERENCES Chromosomes (id)
    FOREIGN KEY (type)
        REFERENCES VarType (id)
    FOREIGN KEY (gene)
        REFERENCES Genes (id)    
);

CREATE TABLE CallData (
    genotype    VARCHAR,
    geno_qual    VARCHAR,
    pass_filter    VARCHAR,
    afreq    VARCHAR,
    coverage    INTEGER,
    ref_count    VARCHAR,
    alt_count    VARCHAR,
    norm_count    FLOAT,
    cn    FLOAT,
    variant    INTEGER,
    sample   INTEGER,
    FOREIGN KEY (variant)
        REFERENCES VarData (id),
    FOREIGN KEY (sample)
        REFERENCES Run (name),
    PRIMARY KEY (variant, sample)
)
''')


conn.close()

