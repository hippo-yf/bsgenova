import time
import random
import os
import io
import math
import numpy as np
from numba import jit
from multiprocessing import Pool

from scipy.stats import multinomial
dmultinom = multinomial.pmf

 

def shrink_depth(depth, threshold = 60):

    k = threshold - np.sqrt(threshold)

    depth[depth > threshold] = np.round(np.sqrt(depth[depth > threshold]) + k)
    return depth


# mutation rate
pm = 1/1000/3

# error rate
# in oocyte samples, error rate is set triple
# pe = 1/100/3
pe = 3/100/3

# total mis rate
p = pm + pe

# methylation rate/proportion 
pr_cg = 0.6        # CG content
pr_ncg = 1/100     # non-CG content


# transition prob of haploidy

def transA(pr):
    return (1-3*pm-3*pe, 2*pm-pm*pr+pe, pm*pr+pe, pm+pe)
def transT(pr):
    return (pm+pe, 1-2*pm-pm*pr-3*pe, pm*pr+pe, pm+pe)
def transC(pr):
    return (pm+pe, pm+pe+(1-3*pm-3*pe)*(1-pr), (1-3*pm-3*pe)*pr, pm+pe)
def transG(pr):
    return (pm+pe, 2*pm-pm*pr+pe, pm*pr+pe, 1-3*pm-3*pe)

PAs = {'CG': np.array(transA(pr_cg)), 'CH': np.array(transA(pr_ncg))}
PTs = {'CG': np.array(transT(pr_cg)), 'CH': np.array(transT(pr_ncg))}
PCs = {'CG': np.array(transC(pr_cg)), 'CH': np.array(transC(pr_ncg))}
PGs = {'CG': np.array(transG(pr_cg)), 'CH': np.array(transG(pr_ncg))}


# STATUS
HOMO = ('A', 'T', 'C', 'G')
HETER = ('AC', 'AG', 'AT', 'CG', 'CT', 'GT')
STATUS = HOMO + HETER


# prior

ps = np.array(((1-3*p)**2, p**2, 2*p*(1-3*p)))

# 0-based
priA = ps[np.array([1,2,2,2,3,3,3,2,2,2]) -1]
priT = ps[np.array([2,1,2,2,2,2,3,2,3,3]) -1]
priC = ps[np.array([2,2,1,2,3,2,2,3,3,2]) -1]
priG = ps[np.array([2,2,2,1,2,3,2,3,2,3]) -1]

pris = {'A': priA, 'T': priT, 'C': priC, 'G': priG}

## allele frequencies

allele_weights = np.array(
    (1, 0, 0, 0, 0.5, 0.5, 0.5, 0  , 0  , 0  ,
    0, 1, 0, 0, 0,   0  , 0.5, 0  , 0.5, 0.5,
    0, 0, 1, 0, 0.5, 0  , 0  , 0.5, 0.5, 0  ,
    0, 0, 0, 1, 0  , 0.5, 0  , 0.5, 0  , 0.5
    ), 
    dtype='float32').reshape(4, 10)

# @jit(nopython=True)
def postp(line: str, args: dict):

    chr, ref, pos, pattern, dinuc, AW, TW, CW, GW, _, TC, AC, GC, CC, _, _ = line.split('\t')
    
    if(ref == 'N'):
        return(None)
  
    coverage = shrink_depth(np.array([AW, TW, CW, GW, AC, TC, CC, GC], dtype='int'))

  
    # prior
  
    theta = pris[ref]
  
    # conditional prob
    
    if(pattern != 'CG'):
        pattern = 'CH'
  
    PA = PAs[pattern]
    PT = PTs[pattern]
    PC = PCs[pattern]
    PG = PGs[pattern]

    watson = coverage[:4]
    crick = coverage[4:]
    DP = coverage.sum()

    prob_cond = np.array(
        [dmultinom(coverage, DP, np.append(PA, PT)/2), # A
             dmultinom(coverage, DP, np.append(PT, PA)/2), # T
             dmultinom(coverage, DP, np.append(PC, PG)/2), # C
             dmultinom(coverage, DP, np.append(PG, PC)/2), # G
             dmultinom(coverage, DP, np.append(PA+PC, PT+PG)/4), # AC
             dmultinom(coverage, DP, np.append(PA+PG, PT+PC)/4), # AG
             dmultinom(coverage, DP, np.append(PA+PT, PT+PA)/4), # AT
             dmultinom(coverage, DP, np.append(PC+PG, PG+PC)/4), # CG
             dmultinom(coverage, DP, np.append(PC+PT, PG+PA)/4), # CT
             dmultinom(coverage, DP, np.append(PG+PT, PC+PA)/4)  # GT
    ])
  
    # posterior prob
    post_unnorm = prob_cond*theta
    post = post_unnorm/post_unnorm.sum()
  
    # prob not mutation (same with ref)
    # regarded as p.value
  
    p_value = float(post[:4][np.array(['A', 'T', 'C', 'G']) == ref])

    if p_value > args['p_value']:
        return None
    
    # allele frequencies
    allele_freq = allele_weights @ post

    DP_watson = watson.sum()
    DP_crick = crick.sum()
    p_homozyte = np.sum(post[1:4])

    line_res = '%s\t%s\t%s\t%e\t%e\t%e\t%e\t%e\t%d\t%d\t%e\n' % (
        chr, pos, ref,
        p_value,
        allele_freq[0], allele_freq[1], allele_freq[2], allele_freq[3],
        DP_watson, DP_crick,
        p_homozyte
    )
    return line_res

def calculate(func, args):
    return func(*args)

# def writeLine(args, OUT):
#     # global TASKS_IN_QUEUE 

#     p_value, allele_freq, DP_watson, DP_crick, p_homozyte = args
#     OUT.write('%e\t%e\t%e\t%e\t%e\t%d\t%d\t%e\n' % 
#               (p_value,
#                allele_freq[0], allele_freq[1], allele_freq[2], allele_freq[3],
#                DP_watson, DP_crick,
#                p_homozyte
#                ))
#     # TASKS_IN_QUEUE -= 1


def writeLine(line):
    OUT.write(line)



if __name__ == '__main__':

    ##
    TASKS_IN_QUEUE = 1 
    BATCH_SIZE = 100000

    # start 4 worker processes
    
    infile = 'D:/Documents/GitHub/BS-SNV-Caller/data/atcg.example'
    IN = io.open(infile, 'r')
    outfile = 'D:/Documents/GitHub/BS-SNV-Caller/data/atcg.out'
    OUT = io.open(outfile, 'w+')
    # OUT = 'data/out'
    
    args = {'p_value': 0.01}

    with Pool(processes=4) as pool:
        

        N = 1
        LAST_BATCH = False
        while True:
            
            k_batch = 0
            TASKS = []
            while (k_batch < BATCH_SIZE):

                if line:= IN.readline().strip():
                    # print( line.strip().split("\t"))

                    TASKS.append((postp, (line, args)))
                    k_batch += 1
                else:
                    LAST_BATCH = True
                    break

            # pool.apply_async(postp, line, callback=writeLine)
            # pool.apply_async(f, (int(id), name, float(x), float(y)), callback=writeLine)
            results = [pool.apply_async(calculate, task) for task in TASKS]
            
            # results = pool.imap(calculate, lines_batch)
            for r in results:
                R = r.get()
                if R:= r.get():
                    OUT.write(R)
                    
            
            if LAST_BATCH:
                break

        pool.close()
        pool.join()
