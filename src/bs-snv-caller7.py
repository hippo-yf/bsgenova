import time
import random
import os
import io
import math
import gzip
import numpy as np
# from numba import jit
from multiprocessing import Pool

# from scipy.stats import multinomial
# dmultinom = multinomial.pmf

 

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


# transition prob of haploidy, likelihood

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

zeros4 = np.array([0]*4)

# Wastson strand

P_AWs = {'CG': np.append(transA(pr_cg), zeros4), 
        'CH': np.append(transA(pr_ncg), zeros4)}
P_TWs = {'CG': np.append(transT(pr_cg), zeros4), 
        'CH': np.append(transT(pr_ncg), zeros4)}
P_CWs = {'CG': np.append(transC(pr_cg), zeros4), 
        'CH': np.append(transC(pr_ncg), zeros4)}
P_GWs = {'CG': np.append(transG(pr_cg), zeros4), 
        'CH': np.append(transG(pr_ncg), zeros4)}

# Crick strand

P_ACs = {'CG': np.append(zeros4, transA(pr_cg)), 
        'CH': np.append(zeros4, transA(pr_ncg))}
P_TCs = {'CG': np.append(zeros4, transT(pr_cg)), 
        'CH': np.append(zeros4, transT(pr_ncg))}
P_CCs = {'CG': np.append(zeros4, transC(pr_cg)), 
        'CH': np.append(zeros4, transC(pr_ncg))}
P_GCs = {'CG': np.append(zeros4, transG(pr_cg)), 
        'CH': np.append(zeros4, transG(pr_ncg))}


# log p_i in multinomial distribution
# likelihood, i.e. conditional

def log_likelihood(pattern):
    P_AW = P_AWs[pattern]
    P_TW = P_TWs[pattern]
    P_CW = P_CWs[pattern]
    P_GW = P_GWs[pattern]
    P_AC = P_ACs[pattern]
    P_TC = P_TCs[pattern]
    P_CC = P_CCs[pattern]
    P_GC = P_GCs[pattern]
    
    mx = np.log(np.vstack([
        (P_AW+P_TC)/2, (P_TW+P_AC)/2, (P_CW+P_GC)/2, (P_GW+P_CC)/2, # A T C G
        (P_AW+P_CW+P_TC+P_GC)/4, # AC
        (P_AW+P_GW+P_TC+P_CC)/4, # AG
        (P_AW+P_TW+P_AC+P_TC)/4, # AT
        (P_CW+P_GW+P_CC+P_GC)/4, # CG
        (P_CW+P_TW+P_AC+P_GC)/4, # CT
        (P_GW+P_TW+P_AC+P_CC)/4  # GT
        ], dtype=float).T)
    return mx

log_p_i = {'CG': log_likelihood('CG'), 'CH': log_likelihood('CH')}



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

## allele frequencies of each genotype

allele_weights = np.array(
    (1, 0, 0, 0, 0.5, 0.5, 0.5, 0  , 0  , 0  ,
    0, 1, 0, 0, 0,   0  , 0.5, 0  , 0.5, 0.5,
    0, 0, 1, 0, 0.5, 0  , 0  , 0.5, 0.5, 0  ,
    0, 0, 0, 1, 0  , 0.5, 0  , 0.5, 0  , 0.5
    ), 
    dtype='float32').reshape(4, 10)


# @jit(nopython=True)
def BS_SNV_Caller(lines: list, args: dict):

    if len(lines) == 0:
        return None
    
    line_res = [l.split('\t') for l in lines]

    array = np.array(line_res)

    ## exclude Ns
    array = array[array[:,1] != 'N',:]
    
    # diff bases in the Crick strand
    reads = array[:,(5,6,7,8, 11,10,13,12)].astype(int) 

    # sqrt-transform read counts
    
    i = reads > args['shrink.reads']
    if i.sum() > 0:
        reads[i] = shrink_depth(reads[i], args['shrink.reads'])

    # exclude sites of low coverage

    i = np.sum(reads, axis=1) < args['min.DP']
    if i.sum() > 0:
        reads = reads[np.logical_not(i), :]
        array = array[np.logical_not(i), :]

    ## basic vars
    
    BASES = ['A', 'T', 'C', 'G']
    N_rows, _ = array.shape

    refs = array[:,1]
    patterns = array[:,3]

    is_CG = patterns == 'CG'

    ## posterior \prop likelihood (multinomial) \times prior

    post_mx = np.zeros((N_rows, 10))

    post_mx[is_CG, :] = reads[is_CG, :] @ log_p_i['CG']
    post_mx[np.logical_not(is_CG), :] = reads[np.logical_not(is_CG), :] @ log_p_i['CH']

    post_mx = np.exp(post_mx - np.max(post_mx, axis=1, keepdims=True))

    # priors

    prior_mx = np.zeros((N_rows, 10))

    for ref in BASES:
        prior_mx[refs == ref, :] = pris[ref]

    # post prob and normalization

    post_mx *= prior_mx
    post_mx = post_mx / np.sum(post_mx, axis=1, keepdims=True)

        
    # prob of unmutation (same with ref) was
    # regarded as p.value

    p_value = np.zeros(N_rows, dtype=float)
    for i in range(len(BASES)):
        I = refs == BASES[i]
        p_value[I] = post_mx[I, i]

    # significant sites

    i_sig = p_value < args['p.value']

    if i_sig.sum() == 0:
        return(None)

    #  only return singnificant sites

    p_values = p_value[i_sig]
    chrs = array[i_sig, 0]
    poss = array[i_sig, 2]
    reffs = refs[i_sig]


    # allele frequencies
    allele_freq = post_mx[i_sig,:] @ allele_weights.T

    DP_watson = np.sum(reads[i_sig, :4], axis=1)
    DP_crick = np.sum(reads[i_sig, 4:8], axis=1)
    p_homozyte = np.sum(post_mx[i_sig, :4], axis=1)

    lines_res = []

    for i in range(i_sig.sum()):
        lines_res.append('%s\t%s\t%s\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d\t%d\n' % (
        chrs[i], poss[i], reffs[i],
        p_values[i], p_homozyte[i],
        allele_freq[i,0], allele_freq[i,1], allele_freq[i,2], allele_freq[i,3],
        DP_watson[i], DP_crick[i]
        ))

    return(lines_res)




def calculate(func, args):
    return func(*args)


def writeLine(lines):
    global TASKS_IN_QUEUE
    TASKS_IN_QUEUE -= 1

    for l in lines:
        OUT.write(l)

def readBatch(IN, BATACH_SIZE = 1000):
    i = 0
    line_batch = []
    while i < BATACH_SIZE:
        line = IN.readline().strip()
        if line:
            line_batch.append(line)
        else:
            break
        i += 1
    return line_batch


class LineFile:
    def __init__(self, filename: str, batchSize: int):
        if filename.endswith(".gz") :
            self.input = gzip.open(filename, 'rt')
        else :
            self.input = io.open(filename, 'r')
        self.batchSize = batchSize
        self.exhausted = False
    # def __iter__(self):
    #     return self
    def __next__(self):
        if self.exhausted:
            return None
        
        # number if readed lines 
        i = 0
        lines = []
        while (l := self.input.readline().strip()) and (i < self.batchSize):
            lines.append(l)
            i += 1
        if i < self.batchSize:
            self.exhausted = True

        if i > 0:
            return lines
        else:
            return None

    def close(self):
        if not self.input.closed:
            self.input.close()

class WaitTimeSchimitter:
    def __init__(self, thres_u:int, thres_l:int, wtime: float, FLAG_WAIT: str):
        self.thres_u = thres_u
        self.thres_l = thres_l
        self.wtime = wtime
        self.FLAG_WAIT = FLAG_WAIT

    def setWaitTimeFlag(self, k:int):
        if self.FLAG_WAIT == 'upper':
            if k < self.thres_l:
                self.FLAG_WAIT = 'lower'
        elif self.FLAG_WAIT == 'lower':
            if k > self.thres_u:
                self.FLAG_WAIT = 'upper'

    def waitTime(self, k: int):
        self.setWaitTimeFlag(k)
        if self.FLAG_WAIT == 'upper':
            time.sleep(self.wtime)
            
            print(f'Waiting: {self.FLAG_WAIT}, {k}')
    

if __name__ == '__main__':

    ##
    TASKS_IN_QUEUE = 0
    BATCH_SIZE = 100000

    # FLAG_WAIT = 'upper'

    # start 4 worker processes
    
    infile = 'D:/Documents/GitHub/BS-SNV-Caller/data/atcg.example'
    # IN = io.open(infile, 'r')
    outfile = 'D:/Documents/GitHub/BS-SNV-Caller/data/atcg.example.out7.gz'
    OUT = gzip.open(outfile, 'wt')
    # OUT = io.open(outfile, 'w+')
    # OUT = 'data/out'
    


    ATCGfile = LineFile(infile, BATCH_SIZE)

    args = {'p.value': 0.01, 'min.DP': 10, 'shrink.reads': 60}

    num_workers = 1

    wait = WaitTimeSchimitter(num_workers*20, num_workers*50, 0.01, 'lower')


    TASKS = []
    with Pool(processes=num_workers) as pool:
        

        while True: 
            if wait.FLAG_WAIT == 'lower':
                line_batch = next(ATCGfile)
                if not line_batch: break

                TASKS_IN_QUEUE += 1
                # pool.apply_async(calculate, (postp, (line_batch, args)), callback=writeLine)
                pool.apply_async(BS_SNV_Caller, (line_batch, args), callback=writeLine)

                # pool.apply_async(f, (int(id), name, float(x), float(y)), callback=writeLine)
                # results = [pool.apply_async(calculate, task) for task in TASKS]

            wait.waitTime(TASKS_IN_QUEUE)



        ATCGfile.close()
        pool.close()
        pool.join()

