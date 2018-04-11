#!/usr/bin/env python2.7
#------------------------------------------------------------------------------
# Filename   : gencorrelation.py
# Updated    : 2014/02/12
# Author     : David Xu
# Class      : I590 Computational Methods for High-Throughput Data
# Assignment : Assignment 1
# Term       : IUPUI Spring 2014
# Usage      : ./gencorrelation.py
# Description: This script requires NumPy, SciPy, and MatPlotLib packages
#              installed. It was written and tested in a Linux/UNIX-like
#              environment using Python 2.7.
#------------------------------------------------------------------------------

import os
import numpy as np
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt

def get_pearson(x,y):
    '''
    Returns the Pearson correlation coefficient for 2 ordered lists
    '''
    assert len(x) == len(y)
    return pearsonr(x,y)[0]

def visualize_matrix(dmatrix):
    '''
    Visualize a matrix
    '''
    lmatrixvizualize = [ '\t'.join([''] + [ i for i in sorted(dmatrix.iterkeys())]) ]
    for k,v in sorted(dmatrix.iteritems()):
        lmatrixvizualize.append('\t'.join(str(x) for x in [k] + [ round(j,3) for i,j in sorted(v.iteritems())]))
    return '\n'.join(lmatrixvizualize)

def read_matrix(pathtomatrix):
    '''
    Reads the tab delimited matrix with column names and row names
    '''
    with open(pathtomatrix) as matrixfile:
        dmatrix = {}
        dmatrixa = {}
        for i, cancertype in enumerate(matrixfile.readline().strip().split('\t')[1:],start=1):
            dmatrix[i] = cancertype
            dmatrixa[cancertype] = {}
        for line in matrixfile:
            line = line.strip().split('\t')
            mirnaname = line[0]
            for i,j in enumerate(line[1:],start=1):
                cancertype = dmatrix.get(i)
                dmatrixa[cancertype][mirnaname] = float(j)
    return dmatrixa

def correlation_matrix(dmatrixa,dmatrixb):
    '''
    Generate the correlation matrix
    '''
    dmatrixc = {}
    for cancertype1, j1 in sorted(dmatrixa.iteritems()):
        dmatrixc[cancertype1] = {}
        for cancertype2, j2 in sorted(dmatrixb.iteritems()):
            i1 = [ y for x,y in sorted(j1.iteritems())]
            i2 = [ y for x,y in sorted(j2.iteritems())]
            dmatrixc[cancertype1][cancertype2] = get_pearson(i1,i2)
    return dmatrixc

def generate_heatmap(dmatrix,outname):
    '''
    Formats and displays a matrix as a heatmap
    '''
    labels = sorted(dmatrix.iterkeys())
    nparray = np.array([[j for i,j in sorted(v.iteritems())] for k,v in sorted(dmatrix.iteritems())])
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(nparray, cmap=plt.cm.Blues,alpha=0.8)
    ax.set_yticks(np.arange(nparray.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(nparray.shape[1]) + 0.5, minor=False)
    ax.invert_xaxis()
    ax.set_xticklabels(labels, minor=False)
    ax.set_yticklabels(labels, minor=False)
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.25,left=0.2)
    #Display the plot
    plt.show()
    #Saves the figure at 300 DPI if working directory is writable
    if os.access(os.getcwd(), os.W_OK):
        fig.savefig(outname,dpi=300)
    return

if __name__ == '__main__':
    '''
    Matrix 1
    '''
    #Read the first matrix and generate the correlation matrix
    part1_matrix = read_matrix('matrix1.txt')
    part1a_matrix = part1_matrix.copy()
    part1_corr_matrix = correlation_matrix(part1_matrix,part1a_matrix)

    #Visualize the correlation matrix from part 1
    print '#####\n#Correlation matrix for part 1\n#####'
    print visualize_matrix(part1_corr_matrix)

    #Heatmap generation of matrix1
    generate_heatmap(part1_corr_matrix,'figure1.png')
    
    print '-' * 50

    '''
    Matrix 2
    '''

    #Read the second matrix and generate the correlation matrix
    part2_matrix = read_matrix('matrix2.txt')
    part2a_matrix = part2_matrix.copy()
    part2_corr_matrix = correlation_matrix(part2_matrix,part2a_matrix)

    #Visualize the correlation matrix from part 2
    print '#####\n#Correlation matrix for part 2\n#####'
    print visualize_matrix(part2_corr_matrix)
    
    #Heatmap generation of matrix2
    generate_heatmap(part2_corr_matrix,'figure2.png')

    print '-' * 50
    
    '''
    Correlation matrix between first and second matrix
    '''
    #correlation matrix and visualization
    combined_corr_matrix = correlation_matrix(part1_matrix, part2_matrix)
    print '#####\n#Correlation matrix for part 2, where rows are matrix1, columns are matrix2\n#####'
    print visualize_matrix(combined_corr_matrix)

    #Heatmap generation of combined correlation matrix
    generate_heatmap(combined_corr_matrix,'figure3.png')
