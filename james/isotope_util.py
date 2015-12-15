# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:32:42 2015

@author: noore
"""
import numpy as np
from scipy.misc import comb # comb(N,k) = The number of combinations of N things taken k at a time

THETA = 0.011 # the natural abundance of 13C among the two isotopes (13C and 12C).

def compute_fractions(counts):
    """
        Calculates the isotope fractions of a compound, given the list of 
        counts (assuming it starts from M+0).
        
        Usage:
            counts    - a list of positive values representing the counts of each
                        isotope starting from M+0
        
        Returns:
            fractions - a list of values between 0..1 that represent the fraction
                        of each isotope from the total pool, after correcting for the
                        natural abundance of 13C
    """
    N = len(counts)-1
    F = np.matrix(np.zeros((N+1, N+1)))
    for i in xrange(N+1):
        for j in xrange(i+1):
            F[i,j] = comb(N-j, i-j) * THETA**(i-j) * (1-THETA)**(N-j)
            
    X = np.matrix(counts, dtype=float).T
    corrected_counds = list((F.I * X).flat)
    
    return corrected_counds

if __name__ == '__main__':
    counts = [900, 100, 5, 900, 5000]
    
    Y = compute_fractions(counts)
    
    print "The corrected isotope relative abundances are:"
    print '-'*50
    print '  |  '.join(map(lambda d: ' M + %d ' % d, range(len(counts))))
    print '  |  '.join(map(lambda s: '%4.1e' % (s*100), Y))
    print '-'*50
    
    print compute_fractions([1]*7)