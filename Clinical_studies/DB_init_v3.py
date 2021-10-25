#!/usr/bin/env python3
## Created on 19-10-2020

#THIS PYTHON SCRIPT REINITIALIZES THE HD200 STANDARDS DATABASE
#RUNNING IT WILL ERASE ALL PREVIOUS DATA IN THE HD200 DB

#read-in vcf header
import mysql.connector
import os
import re
import numpy
import pandas
import vcf
import json
from io import StringIO
import ast

conn = mysql.connector.connect(host='localhost',user='##',password='##',database='##')
cur = conn.cursor()

# Do some setup
# 7 Metadata tables with reference genomes, workflow info, disease info etc. will be created
#The Other variant info metadata tables with gene and chromosome info will also be created. 
# The main tables are CallData (genotypes), VarData (variants), and RunInfo (Samples)
cur.execute('''DROP TABLE IF EXISTS CallData''')
cur.execute('''DROP TABLE IF EXISTS VarData''')
cur.execute('''DROP TABLE IF EXISTS RunInfo''')
cur.execute('''DROP TABLE IF EXISTS WorkFlowVer''')
cur.execute('''DROP TABLE IF EXISTS Chromosomes''')
cur.execute('''DROP TABLE IF EXISTS BaseCallVer''')
cur.execute('''DROP TABLE IF EXISTS VarType''')
cur.execute('''DROP TABLE IF EXISTS DisType''')
cur.execute('''DROP TABLE IF EXISTS ref_genome''')
cur.execute('''DROP TABLE IF EXISTS AnnotationVersion''')
cur.execute('''DROP TABLE IF EXISTS WorkFlowName''')
cur.execute('''DROP TABLE IF EXISTS Genes''')
cur.execute('''CREATE TABLE ref_genome (id INTEGER NOT NULL AUTO_INCREMENT, name VARCHAR(100), UNIQUE (id), UNIQUE (name), PRIMARY KEY (id))''')
cur.execute('''CREATE TABLE AnnotationVersion (id   INTEGER NOT NULL AUTO_INCREMENT, name VARCHAR(100), UNIQUE (id), UNIQUE (name), PRIMARY KEY (id))''')
cur.execute('''CREATE TABLE WorkFlowName (id   INTEGER NOT NULL AUTO_INCREMENT, name VARCHAR(100), UNIQUE (id), UNIQUE (name), PRIMARY KEY (id))''')
cur.execute('''CREATE TABLE WorkFlowVer (id   INTEGER NOT NULL AUTO_INCREMENT, name VARCHAR(100), UNIQUE (id), UNIQUE (name), PRIMARY KEY (id))''')
cur.execute('''CREATE TABLE BaseCallVer (id   INTEGER NOT NULL AUTO_INCREMENT, name VARCHAR(255), UNIQUE (id), UNIQUE (name), PRIMARY KEY (id))''')
cur.execute('''CREATE TABLE Chromosomes (id   INTEGER NOT NULL AUTO_INCREMENT, name VARCHAR(255), UNIQUE (id), UNIQUE (name), PRIMARY KEY (id))''')
cur.execute('''CREATE TABLE VarType (id   INTEGER NOT NULL AUTO_INCREMENT, name VARCHAR(255), UNIQUE (id), UNIQUE (name), PRIMARY KEY (id))''')
cur.execute('''CREATE TABLE DisType (id   INTEGER NOT NULL AUTO_INCREMENT,name    VARCHAR(255), UNIQUE (id), UNIQUE (name), PRIMARY KEY (id))''')
cur.execute('''CREATE TABLE Genes (id   INTEGER NOT NULL AUTO_INCREMENT, name    VARCHAR(255), UNIQUE (id), UNIQUE (name), PRIMARY KEY (id))''')

cur.execute('''CREATE TABLE RunInfo (
    id     INTEGER NOT NULL AUTO_INCREMENT,
    name   VARCHAR(255),
    sex    VARCHAR(50),
    pheno    INTEGER,
    fileformat    VARCHAR(255),
    filedate    INTEGER,
    reference    INTEGER,
    OM_anno_version    INTEGER,
    IonWF    INTEGER,
    IonWF_version    INTEGER,
    Cellularity    FLOAT,
    fusionOC    VARCHAR(255),
    fusionQC    VARCHAR(255),
    fusionReads    INTEGER,
    basecall_ver    INTEGER,
    perc_mapped    FLOAT,
    CONSTRAINT fk_reference
    FOREIGN KEY (reference)
        REFERENCES ref_genome (id)
        ON DELETE SET NULL
        ON UPDATE CASCADE,
    CONSTRAINT fk_annoversion
    FOREIGN KEY (OM_anno_version)
        REFERENCES AnnotationVersion (id)
	ON DELETE SET NULL
	ON UPDATE CASCADE,
    CONSTRAINT fk_IonWF
    FOREIGN KEY (IonWF)
        REFERENCES WorkFlowName (id)
	ON DELETE SET NULL
        ON UPDATE CASCADE,
    CONSTRAINT fk_IonWFversion
    FOREIGN KEY (IonWF_version)
        REFERENCES WorkFlowVer (id)
        ON DELETE SET NULL
        ON UPDATE CASCADE,
    CONSTRAINT fk_BaseCallVer
    FOREIGN KEY (basecall_ver)
        REFERENCES BaseCallVer (id)
        ON DELETE SET NULL
        ON UPDATE CASCADE,
    CONSTRAINT fk_pheno
    FOREIGN KEY (pheno)
        REFERENCES DisType (id)
        ON DELETE SET NULL
        ON UPDATE CASCADE,
    UNIQUE (id),
    UNIQUE (name),
    PRIMARY KEY (id)
) ENGINE INNODB''')
cur.execute('''CREATE TABLE VarData (
    id    INTEGER NOT NULL AUTO_INCREMENT,
    name    VARCHAR(255),
    chr    INTEGER,
    start    INTEGER,
    len    VARCHAR(50),
    ref    VARCHAR(255),
    alt    VARCHAR(255),
    type    INTEGER,
    annotation    VARCHAR(255),
    gene    INTEGER,
    FOREIGN KEY (chr)
        REFERENCES Chromosomes (id)
        ON DELETE SET NULL
        ON UPDATE CASCADE,
    FOREIGN KEY (type)
        REFERENCES VarType (id)
        ON DELETE SET NULL
        ON UPDATE CASCADE,
    FOREIGN KEY (gene)
        REFERENCES Genes (id)
        ON DELETE SET NULL
        ON UPDATE CASCADE,
    UNIQUE (id),
    UNIQUE (name),
    PRIMARY KEY (id)
) ENGINE INNODB''')
cur.execute('''CREATE TABLE CallData (
    id INTEGER NOT NULL AUTO_INCREMENT,
    genotype    VARCHAR(50),
    geno_qual    VARCHAR(255),
    pass_filter    VARCHAR(50),
    afreq    VARCHAR(50),
    coverage    INTEGER,
    ref_count    VARCHAR(255),
    alt_count    VARCHAR(255),
    norm_count    FLOAT,
    cn    FLOAT,
    variant    INTEGER,
    sample   INTEGER,
    PRIMARY KEY (id),
    FOREIGN KEY (variant)
        REFERENCES VarData (id)
        ON DELETE SET NULL
        ON UPDATE CASCADE,
    FOREIGN KEY (sample)
        REFERENCES RunInfo (id)
        ON DELETE SET NULL
        ON UPDATE CASCADE,
    UNIQUE INDEX variant_sample (variant, sample)
 ) ENGINE INNODB''')

conn.commit()
conn.close()
