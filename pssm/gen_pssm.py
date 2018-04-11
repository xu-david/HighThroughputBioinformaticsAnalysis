#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
# Filename   : gen_pssm.py
# Updated    : 2014/03/27
# Author     : David Xu
# Class      : I590 Computational Methods for High-Throughput Data
# Assignment : Assignment 4
# Term       : IUPUI Spring 2014
# Usage      : ./gen_pssm.py
# Description: This script requires the default library installed. It was
#              written and tested in a Linux/UNIX-like environment using Python
#              2.7.
#------------------------------------------------------------------------------

from bz2 import BZ2File
from math import log
from itertools import islice

def read_count_matrix(filen):
    '''
    Read the initial count matrix
    '''
    dmatrix = {}
    with open(filen) as matrixf:
        for matrixline in matrixf:
            matrixline = matrixline.strip().split()
            dmatrix[matrixline[0]] = [int(x) for x in matrixline[2:]]
    return dmatrix
    
def generate_freq_weight_matrices(dcountmatrix):
    '''
    Generates the frequency and weight matrices from the count matrix
    '''
    dfreqmatrix = {}
    dweightmatrix = {}
    for k,v in dcountmatrix.iteritems():
        dfreqmatrix[k] = [(float(x)+1.0)/31.0 for x in v]
        dweightmatrix[k] = [log(x/0.25) for x in dfreqmatrix.get(k)]
    return dfreqmatrix, dweightmatrix

def sliding_window(seq, n=18):
    '''
    Sliding window generator object from the itertools module docs
    '''
    it = iter(seq)
    result = tuple(islice(it,n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
    
def get_max_score(geneseq,dweightmatrix):
    '''
    Returns the max score of the gene sequence based on the PSSM
    '''
    maxscore = None
    for swindow in sliding_window(geneseq):
        score = 0
        for i,j in enumerate(swindow):
            score += dweightmatrix[j][i]
        if score > maxscore:
            maxscore = score
    return maxscore
    
def scan_binding_site(dweightmatrix,inputfilen):
    '''
    Reads the compressed bz2 input file and returns the max score per gene
    '''
    dbindingsitematrix = {}
    with BZ2File(inputfilen,'r') as f:
        for line in f:
            line = line.strip().split('\\')
            geneid = line[0].strip()
            geneseq = line[1].strip()
            genescore = get_max_score(geneseq,dweightmatrix)
            dbindingsitematrix[geneid] = genescore
    return dbindingsitematrix

def write_outputs(dcountmatrix,dfreqmatrix,dweightmatrix,dbindingsitematrix):
    '''
    Write frequency matrix, weight matrix, and top 30 genes to output files
    '''
    with open('output.files/top30_genes.txt','w') as fout, open('output.files/all_genes.txt','w') as foutall:
        print >> fout, '\t'.join(('Rank','GeneID','Score'))
        print >> foutall, '\t'.join(('Rank','GeneID','Score'))
        print '\t'.join(('Rank','GeneID','Score'))
        for i,(k,v) in enumerate(sorted(dbindingsitematrix.iteritems(),key=lambda x:x[1],reverse=True),start=1):
            if i <= 30:
                print >> fout, '\t'.join(str(x) for x in (i,k,v))
                print '\t'.join(str(x) for x in (i,k,v))
            print >> foutall, '\t'.join(str(x) for x in (i,k,v))                        
    for dmatrix,filen in ((dfreqmatrix,'output.files/freq_matrix.txt'),(dweightmatrix,'output.files/weight_matrix.txt')):
        with open(filen,'w') as fout:
            for k,v in dmatrix.iteritems():
                print >> fout, '\t'.join([k,'|'] + [str(x) for x in (v)])
    return
    
def main():
    print 'Reading the counts matrix...'
    dcountmatrix = read_count_matrix('input.files/argR-counts-matrix.txt')
    print 'Generating the frequency and weight matrices...'
    dfreqmatrix, dweightmatrix = generate_freq_weight_matrices(dcountmatrix)
    print 'Scanning the binding sites...'
    dbindingsitematrix = scan_binding_site(dweightmatrix, 'input.files/E_coli_K12_MG1655.400_50.bz2')
    print 'Writing outputs...'
    print '-' * 80
    write_outputs(dcountmatrix,dfreqmatrix,dweightmatrix,dbindingsitematrix)
    return
    
if __name__ == '__main__':
    main()
