#!/usr/bin/env python2.7
#------------------------------------------------------------------------------
# Filename   : gethalflife.py
# Updated    : 2014/02/28
# Author     : David Xu
# Class      : I590 Computational Methods for High-Throughput Data
# Assignment : Assignment 2
# Term       : IUPUI Spring 2014
# Usage      : ./gethalflife.py
# Description: This script requires the NumPy and SciPy packages
#              installed. It was written and tested in a Linux/UNIX-like
#              environment using Python 2.7.
#------------------------------------------------------------------------------

import numpy
import csv
from scipy import stats
from math import log

#Supress errors from linear regression calculations for remove later
numpy.seterr(all='ignore')

def linear_regression(ltimecourse):
    '''
    Returns slope from the linear regression
    '''
    return stats.linregress(ltimecourse)[0]

def list_mean(inlist):
    '''
    Returns the mean of a list, otherwise returns None
    '''
    if len(inlist):
        return float(sum(inlist))/len(inlist)
    else:
        return

def read_timecourse(timedata):
    '''
    Reads the initial timecourse file
    '''
    dtimedata = {}
    timecourset = [0.,5.,10.,15.,20.,30.,40.,50.,60.]
    with open(timedata) as timedatafile:
        for i in range(2):
            timedatafile.next()
        csvr = csv.reader(timedatafile,delimiter='\t',quotechar='"')
        for line in csvr:
            gene = line[0]
            timecourse1 = line[1:10]
            timecourse2 = line[10:19]
            timecourse3 = line[19:28]
            halflifes = ['','','']
            for timecoursen,timecourse in enumerate((timecourse1,timecourse2,timecourse3)):
                #Ignore timecourse if no timepoints
                if all(x=='' for x in timecourse):
                    continue
                ltimecourse = []
                for i,j in zip(timecourset,timecourse):
                    try:
                        #Natural log transform value if possible
                        j = log(float(j))
                        ltimecourse.append((i,j))
                    except:
                        continue
                #Ignore timecourse with less than 2 points
                if len(ltimecourse) > 1:
                    #Transform into two-dimensional array
                    ltimecourse2 = [[ltimecourse[i][0] for i in range(len(ltimecourse))], [ltimecourse[i][1] for i in range(len(ltimecourse))]]
                    #Absolute value of slope from linear regression
                    slope = abs(linear_regression(ltimecourse2))
                    halflife = log(2.0)/slope
                    halflifes[timecoursen] = halflife
            #Remove empty and inf values from average half-life calculation
            listin = [x for x in halflifes if x != '' and not any((numpy.isnan(x),numpy.isinf(x)))]
            dtimedata[gene] = halflifes + [list_mean(listin)]
    return dtimedata

def percentage_of_dict(dtimedata,fractioncut,orderby,outfilen):
    '''
    Returns top or bottom x percent of the time course data
    '''
    cutoff = int(len(dtimedata.keys()) * fractioncut)
    count = 0
    with open(outfilen,'w') as outfile:
        if outfilen == 'results/sorted_timecourse.txt':
            print >> outfile, '\t'.join(['Rank','Gene','HalfLife1','HalfLife2','HalfLife3','AverageHalfLife'])
        for gene,halflife in sorted(dtimedata.iteritems(),key=lambda x: x[1][3],reverse=orderby):
            #Prints the full list of genes and halflifes
            if outfilen == 'results/sorted_timecourse.txt':
                count += 1
                print >> outfile, '\t'.join(str(x) for x in [count, gene, halflife[0], halflife[1], halflife[2], halflife[3]])
            #Otherwise print the top/bottom list
            elif halflife[3]:
                count += 1
                print >> outfile, gene
                if count >= cutoff:
                    break
    return

def main():
    '''
    Main wrapper
    '''
    timedata = 'DecayTimecourse.txt'
    print 'Reading the initial time course and calculating the average half-lives...'
    dtimedata = read_timecourse(timedata)

    print 'Writing the sorted average half-lives into results/sorted_timecourse.txt...'
    percentage_of_dict(dtimedata,1.0,True,'results/sorted_timecourse.txt')
    
    print 'Writing the top 10% into results/top10percent.txt...'
    percentage_of_dict(dtimedata,0.1,True,'results/top10percent.txt')

    print 'Writing the bottom 10% into results/top10percent.txt...'
    percentage_of_dict(dtimedata,0.1,False,'results/bot10percent.txt')
    return

if __name__ == '__main__':
    main()
