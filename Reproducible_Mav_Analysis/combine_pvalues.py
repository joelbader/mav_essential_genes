'''
Truncated Product Method: (Zaykin, Genet Epidemiol. 2002)
'''

from tpm import tpm
import numpy as np
import time

def combine_pvalues(pvalues, tau=0.05, seed=time.time(),
                    nperms=1000):
    '''
        returns: pvalue

        pvalues: list of floats
        tau: (used for tpm) between 0 and 1
        seed: (used for tpm) random number generator seed
        nperms: (used for tpm) number of permutations to do
    '''
    
    #check if all values are NA
    if all(p == "NA" for p in pvalues):
        print("all values are NA")
        return np.nan
    #remove NAs from list
    pvalues = [p for p in pvalues if p != "NA"]
    
    combined_pval = tpm(tau, nperms, seed, pvalues)
        
    return combined_pval
