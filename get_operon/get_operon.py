#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
# Filename   : get_operon.py
# Updated    : 2014/03/17
# Author     : David Xu
# Class      : I590 Computational Methods for High-Throughput Data
# Assignment : Assignment 3
# Term       : IUPUI Spring 2014
# Usage      : ./get_operon.py
# Description: This script requires no addition packages installed other than
#              the default library. It was written and tested in a
#              Linux/UNIX-like environment using Python 2.7.
#------------------------------------------------------------------------------

import gzip

def read_ptt_gz(inputgzfile,startn):
    '''
    Read a PTT.gz file
    '''
    loperons = []
    with gzip.open(inputgzfile,'rb') as gzipfile:
        gzipflines = gzipfile.readlines()
        old_gene = (None,None,None,None,None)
        loperon = []
        for gzipline in gzipflines[startn:]:
            gzipline = gzipline.strip().split('\t')
            sdirection = gzipline[1]
            pid = gzipline[3]
            sstart, send = gzipline[0].split('..')
            sstart = int(sstart)
            send = int(send)
            genename = gzipline[4]
            if old_gene != (None,None,None,None,None):
                dist_ij = sstart - old_gene[2]
                if old_gene[3] == sdirection and dist_ij < 50:
                    old_gene = (pid,sstart,send,sdirection,genename)
                    loperon.append(old_gene)
                else:
                    loperons.append(loperon)
                    loperon = []
                    old_gene = (pid,sstart,send,sdirection,genename)
                    loperon.append(old_gene)
            else:
                old_gene = (pid,sstart,send,sdirection,genename)
                loperon.append(old_gene)
        else:
            loperons.append(loperon)
    return loperons

def read_gff(inputfile,startn):
    '''
    Read a GFF file
    '''
    loperons = []
    with open(inputfile) as infile:
        infilelines = infile.readlines()
        old_gene = (None,None,None,None,None)
        loperon = []
        for infileline in infilelines[startn:]:
            infileline = infileline.strip().split('\t')
            sdirection = infileline[6]
            pid = infileline[0]
            sstart = int(infileline[3])
            send = int(infileline[4])
            genename = infileline[8].split(';')[0][3:]
            if old_gene != (None,None,None,None,None):
                dist_ij = sstart - old_gene[2]
                if old_gene[3] == sdirection and dist_ij < 50:
                    old_gene = (pid,sstart,send,sdirection,genename)
                    loperon.append(old_gene)
                else:
                    loperons.append(loperon)
                    loperon = []
                    old_gene = (pid,sstart,send,sdirection,genename)
                    loperon.append(old_gene)
            else:
                old_gene = (pid,sstart,send,sdirection,genename)
                loperon.append(old_gene)
        else:
            loperons.append(loperon)
    return loperons
    
def write_output(loperons,outputn):
    '''
    Write the list of operons to the output file
    '''
    with open(outputn,'w') as fout:
        print >> fout, '#Operon','\t','Genes in Operon'
        count = 0
        for operon in loperons:
            #Skip single genes
            if len(operon) > 1:
                count += 1
                #print count,'\t',operon
                print >> fout, count,'\t',operon
    return
    
def main():
    '''
    Wrapper function to read and write each transcript
    '''
    #E. coli
    print 'Reading E. coli...'
    ecoli_operons = read_ptt_gz('input.files/E_coli_K12_MG1655.ptt.gz',3)
    write_output(ecoli_operons,'output.files/E_coli_operons.txt')
    
    #B. subtilis
    print 'Reading B. subtilis...'
    bsubtilis_operons = read_ptt_gz('input.files/B_subtilis_168.ptt.gz',3)
    write_output(bsubtilis_operons,'output.files/B_subtilis_operons.txt')
    
    #Halobacterium
    print 'Reading Halobacterium...'
    halobacterium_operons = read_ptt_gz('input.files/Halobacterium_NRC1.ptt.gz',3)
    write_output(halobacterium_operons,'output.files/Halobacterium_operons.txt')
    
    #Synechocystis
    print 'Reading Synechocystis...'
    synechocystis_operons = read_ptt_gz('input.files/Synechocystis_PCC6803_uid159873.ptt.gz',3)
    write_output(synechocystis_operons,'output.files/Synechocystis_operons.txt')
    
    #Hoatzin
    print 'Reading Hoatzin...'
    hoatzin_operons = read_gff('input.files/2088090036.gff',1)
    write_output(hoatzin_operons,'output.files/Hoatzin_operons.txt')
    return

if __name__ == '__main__':
    main()
