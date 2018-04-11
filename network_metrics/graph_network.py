#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
# Filename   : graph_network.py
# Updated    : 2014/04/15
# Author     : David Xu
# Class      : I590 Computational Methods for High-Throughput Data
# Assignment : Assignment 5
# Term       : IUPUI Spring 2014
# Usage      : ./gen_pssm.py
# Description: This script requires the NumPy, SciPy, MatPlotLib, and NetworkX
#              packages installed. It was written and tested in both a
#              Linux/UNIX-like and Windows environment using Python 2.7.
#------------------------------------------------------------------------------

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import ranksums
import powerlaw

#Ignore Numpy errors in powerlaw function
np.seterr(invalid='ignore')

def generate_network(infile):
    '''
    Generate the graph object, removing duplicate edges
    '''
    return nx.read_edgelist(infile, nodetype=str)

def calc_node_stats(G, outfile):
    '''
    Calculate each node's degree and clustering coefficient into output file.
    Also calculates average clustering coefficient of the network.
    '''
    print 'Writing node statistics to \'{}\''.format(outfile)
    with open(outfile, 'w') as fout:
        print >> fout, '\t'.join(str(x) for x in ('Symbol', 'Degree', 'Clustering Coefficient'))
        for node in sorted(G.nodes()):
            node_cc = nx.clustering(G, node)
            node_degree = G.degree(node)
            print >> fout, '\t'.join(str(x) for x in (node, node_degree, node_cc))
    print '-' * 80
    print 'Network Average Clustering Coefficient:', nx.average_clustering(G)
    print '-' * 80
    return

def plot_degree_distribution(G):
    '''
    Plot the degree distribution and calculate power law fit
    '''
    degree_sequence = nx.degree(G).values()
    lcounts = [(x, degree_sequence.count(x)) for x in set(degree_sequence)]
    plt.loglog([x[0] for x in lcounts], [x[1] for x in lcounts], 'ko')
    plt.title('Log-Log Degree Distribution')
    plt.xlabel('Node Degree')
    plt.ylabel('Number of Nodes')
    plt.show()
    results = powerlaw.Fit(degree_sequence)
    print 'Power law alpha:', results.power_law.alpha
    print 'Power law sigma:', results.power_law.sigma
    print 'Power law x-min:', results.power_law.xmin
    print '-' * 80
    return

def read_prot_list(infile):
    '''
    Read initial lists of proteins
    '''
    with open(infile) as f:
        l = []
        for line in f:
            if '\r' in line:
                l = line.strip().split('\r')
            else:
                l.append(line.strip())
    return l

def calc_shortest_path(G, plist, outfile):
    '''
    Calculate the shortest paths for a list of genes in a network
    '''
    dlength = {}
    with open(outfile, 'w') as fout:
        print >> fout, '\t'.join(('Source', 'Target', 'Path', 'Path Length'))
        for i in plist:
            for j in plist:
                #Ignore proteins not in network, duplicate pairs, and same pairs
                if i not in G or j not in G or i >= j:
                    continue
                spl_ij = nx.shortest_path_length(G, i, j)
                sp_ij = nx.shortest_path(G, source=i, target=j)
                print >> fout, '{}\t{}\t{}\t{}'.format(i, j, sp_ij, spl_ij)
                dlength[(i, j)] = spl_ij
    return dlength

def compare_lists_wilcox(dshortpath1, dshortpath2):
    '''
    Calculates the Wilcoxon rank-sums z-statistic and p-value for two samples
    '''
    zstat, pvalue = ranksums(dshortpath1.values(), dshortpath2.values())
    print 'Wilcoxon rank-sum test z-statistic:', zstat
    print 'Wilcoxon rank-sum test two-tailed p-value:', pvalue
    return

def main():
    '''
    Main wrapper
    '''
    dnetwork = generate_network('input.files/Human-PPI.txt')
    calc_node_stats(dnetwork, 'output.files/node_properties.txt')
    plot_degree_distribution(dnetwork)
    plist1 = read_prot_list('input.files/protein-list1.txt')
    plist2 = read_prot_list('input.files/protein-list2.txt')
    dshortpath1 = calc_shortest_path(dnetwork, plist1, 'output.files/list1_shortestpaths.txt')
    dshortpath2 = calc_shortest_path(dnetwork, plist2, 'output.files/list2_shortestpaths.txt')
    compare_lists_wilcox(dshortpath1, dshortpath2)
    return

if __name__ == '__main__':
    main()
