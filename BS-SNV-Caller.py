import time
import os
import io
import gzip
import re
import numpy as np
from collections import namedtuple
# from numba import jit
from multiprocessing import Pool

# from scipy.stats import multinomial
# dmultinom = multinomial.pmf

from argparse import ArgumentParser
from argparse import Namespace


class SNVparams:
    def __init__(self, args: Namespace):
        self.infile = args.infile
        self.outprefix = args.outprefix
        self.mutation_rate = args.mutation_rate/3
        self.error_rate = args.error_rate/3
        # total mismatch rate
        self.mis_rate = self.error_rate + self.mutation_rate
        self.methy_cg = args.methy_cg
        self.methy_ncg = args.methy_ncg
        self.min_depth = args.min_depth

        
        self.pvalue = args.pvalue
        self.shrink_depth = args.shrink_depth

        self.batch_size = args.batch_size
        self.num_process = args.num_process

        self.pool_lower_num = args.pool_lower_num
        self.pool_upper_num = args.pool_upper_num


    # def reset_param(self):
    #     self.mutation_rate = self.mutation_rate/3
    #     self.error_rate = self.error_rate/3
    #     self.mis_rate = self.error_rate + self.mutation_rate

    # transition prob of haploidy, likelihood

    def transA(self, pr):
        return (1-3*self.mutation_rate-3*self.error_rate, 2*self.mutation_rate-self.mutation_rate*pr+self.error_rate, self.mutation_rate*pr+self.error_rate, self.mutation_rate+self.error_rate)
    def transT(self, pr):
        return (self.mutation_rate+self.error_rate, 1-2*self.mutation_rate-self.mutation_rate*pr-3*self.error_rate, self.mutation_rate*pr+self.error_rate, self.mutation_rate+self.error_rate)
    def transC(self, pr):
        return (self.mutation_rate+self.error_rate, self.mutation_rate+self.error_rate+(1-3*self.mutation_rate-3*self.error_rate)*(1-pr), (1-3*self.mutation_rate-3*self.error_rate)*pr, self.mutation_rate+self.error_rate)
    def transG(self, pr):
        return (self.mutation_rate+self.error_rate, 2*self.mutation_rate-self.mutation_rate*pr+self.error_rate, self.mutation_rate*pr+self.error_rate, 1-3*self.mutation_rate-3*self.error_rate)

    def set_trans_prob(self):
            

        zeros4 = np.array([0]*4)

        pr_cg = self.methy_cg
        pr_ncg = self.methy_ncg

        # Wastson strand

        self.P_AWs = {'CG': np.append(self.transA(pr_cg), zeros4), 
                'CH': np.append(self.transA(pr_ncg), zeros4)}
        self.P_TWs = {'CG': np.append(self.transT(pr_cg), zeros4), 
                'CH': np.append(self.transT(pr_ncg), zeros4)}
        self.P_CWs = {'CG': np.append(self.transC(pr_cg), zeros4), 
                'CH': np.append(self.transC(pr_ncg), zeros4)}
        self.P_GWs = {'CG': np.append(self.transG(pr_cg), zeros4), 
                'CH': np.append(self.transG(pr_ncg), zeros4)}

        # Crick strand

        self.P_ACs = {'CG': np.append(zeros4, self.transA(pr_cg)), 
                'CH': np.append(zeros4, self.transA(pr_ncg))}
        self.P_TCs = {'CG': np.append(zeros4, self.transT(pr_cg)), 
                'CH': np.append(zeros4, self.transT(pr_ncg))}
        self.P_CCs = {'CG': np.append(zeros4, self.transC(pr_cg)), 
                'CH': np.append(zeros4, self.transC(pr_ncg))}
        self.P_GCs = {'CG': np.append(zeros4, self.transG(pr_cg)), 
                'CH': np.append(zeros4, self.transG(pr_ncg))}

    # log p_i in multinomial distribution
    # likelihood, i.e. conditional probability

    def log_likelihood(self, pattern):
        P_AW = self.P_AWs[pattern]
        P_TW = self.P_TWs[pattern]
        P_CW = self.P_CWs[pattern]
        P_GW = self.P_GWs[pattern]
        P_AC = self.P_ACs[pattern]
        P_TC = self.P_TCs[pattern]
        P_CC = self.P_CCs[pattern]
        P_GC = self.P_GCs[pattern]
        
        mx = np.log(np.vstack([
            (P_AW+P_TC)/2, (P_TW+P_AC)/2, (P_CW+P_GC)/2, (P_GW+P_CC)/2, # A T C G
            (P_AW+P_CW+P_TC+P_GC)/4, # AC
            (P_AW+P_GW+P_TC+P_CC)/4, # AG
            (P_AW+P_TW+P_AC+P_TC)/4, # AT
            (P_CW+P_GW+P_CC+P_GC)/4, # CG
            (P_CW+P_TW+P_AC+P_GC)/4, # CT
            (P_GW+P_TW+P_AC+P_CC)/4  # GT
            # ], dtype=float).T) # dtype is introduced in Ver 1.24
            ]).T)
        return mx

    def set_likelihood(self):
        self.loglik= {'CG': self.log_likelihood('CG'), 'CH': self.log_likelihood('CH')}

    def set_prior(self):

        # # STATUS
        # HOMO = ('A', 'T', 'C', 'G')
        # HETER = ('AC', 'AG', 'AT', 'CG', 'CT', 'GT')
        # STATUS = HOMO + HETER

        # prior
        p = self.mis_rate
        ps = np.array(((1-3*p)**2, p**2, 2*p*(1-3*p)))

        # 0-based
        priA = ps[np.array([1,2,2,2,3,3,3,2,2,2]) -1]
        priT = ps[np.array([2,1,2,2,2,2,3,2,3,3]) -1]
        priC = ps[np.array([2,2,1,2,3,2,2,3,3,2]) -1]
        priG = ps[np.array([2,2,2,1,2,3,2,3,2,3]) -1]

        self.priors = {'A': priA, 'T': priT, 'C': priC, 'G': priG}

    def set_allele_weights(self):
            
        ## allele frequencies of each genotype

        self.allele_weights = np.array(
            (1, 0, 0, 0, 0.5, 0.5, 0.5, 0  , 0  , 0  ,
            0, 1, 0, 0, 0,   0  , 0.5, 0  , 0.5, 0.5,
            0, 0, 1, 0, 0.5, 0  , 0  , 0.5, 0.5, 0  ,
            0, 0, 0, 1, 0  , 0.5, 0  , 0.5, 0  , 0.5
            ), 
            dtype='float32').reshape(4, 10)

    def set_out_files(self):

        # set default sampel name
        if self.outprefix is None or self.outprefix == "":
            dir = os.path.dirname(self.infile)
            filename = os.path.basename(self.infile)
            fn = re.split(r'\.', filename)[0]
            self.sampleID = fn # sample name
            self.outprefix = os.path.join(dir, fn)
        else:
            fn = re.split(r'/', self.outprefix)[-1]
            self.sampleID = fn # sample name

        self.out_snv = self.outprefix + ".snv.gz"
        self.out_vcf = self.outprefix + ".vcf.gz"

                
        ## vcf header
        self.VCF_HEADER = [
            '##fileformat=VCFv4.4\n',
            f'##fileDate=%s\n' % time.strftime("%Y%m%d", time.localtime()),
            '##source=BS-SNV-Caller\n',
            '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n',
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n',
            '##INFO=<ID=DPW,Number=1,Type=Integer,Description="Total Depth of Wastson Strand">\n',
            '##INFO=<ID=DPC,Number=1,Type=Integer,Description="Total Depth of Crick Strand">\n',
            '##FILTER=<ID=q30,Description="Quality < 30">\n',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality.In some cases of single-stranded coverge, we are sure there is a SNV, but we can not determine the alternative variant. So, we express the GQ as the Phred score (-10log10 (p-value)) of posterior probability of homozygote/heterozygote, namely, Prob(heterozygote) for homozygous sites and Prob(homozygote) for heterozygous sites. This is somewhat different with SNV calling from WGS data.">\n',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
            '##FORMAT=<ID=DPW,Number=1,Type=Integer,Description="Read Depth of Wastson Strand">\n',
            '##FORMAT=<ID=DPC,Number=1,Type=Integer,Description="Read Depth of Crick Strand">\n',
            f'#CHROM\tPOS\tID\tREF\tALT\t\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sampleID}\n'
            ]

    def set_model_params(self):
        self.set_out_files()
        self.set_trans_prob()
        self.set_likelihood()
        self.set_prior()
        self.set_allele_weights()
        


def shrink_depth(depth, threshold = 60):

    k = threshold - np.sqrt(threshold)

    depth[depth > threshold] = np.round(np.sqrt(depth[depth > threshold]) + k)
    return depth

class ControlWriteFile:

    def __init__(self, params: SNVparams):
        self.first = True
        # self.current_snv = None
        # self.current_vcf = None
        self.TASKS_IN_QUEUE = 0

        self.out_snv = gzip.open(params.out_snv, 'wt')
        self.out_vcf = gzip.open(params.out_vcf, 'wt')

        # Schimitter waiting
        self.thres_u = params.pool_upper_num*params.num_process
        self.thres_l = params.pool_lower_num*params.num_process
        self.wtime = 0.01
        self.FLAG_WAIT = 'lower'

    def setWaitTimeFlag(self):
        if self.FLAG_WAIT == 'upper':
            if self.TASKS_IN_QUEUE < self.thres_l:
                self.FLAG_WAIT = 'lower'
        elif self.FLAG_WAIT == 'lower':
            if self.TASKS_IN_QUEUE > self.thres_u:
                self.FLAG_WAIT = 'upper'

    def add_task_num(self):
        self.TASKS_IN_QUEUE += 1
        self.setWaitTimeFlag()

    def minus_task_num(self):
        self.TASKS_IN_QUEUE -= 1
        self.setWaitTimeFlag()

    def waitTime(self):
        if self.FLAG_WAIT == 'upper':
            time.sleep(self.wtime)

    def close(self):
        if not self.out_snv.closed:
            self.out_snv.close()
        if not self.out_vcf.closed:
            self.out_vcf.close()

SNVResultsBatch = namedtuple('SNVResultsBatch', ['snv', 'vcf'])

# class SNVResultsBatch:

#     def __init__(self, chrs, poss, reffs, p_values, p_homozyte, allele_freq, DP_watson, DP_crick):
#         self.chrs = chrs
#         self.poss = poss
#         self.reffs = reffs
#         self.p_values = p_values
#         self.p_homozyte = p_homozyte
#         self.allele_freq = allele_freq
#         self.DP_watson = DP_watson
#         self.DP_crick = DP_crick

#         self.len = len(chrs)
#         self.res_snv = None
#         self.res_vcf = None
    
#     def format_snv(self):
#         # .snv format output
#         res_snv = []
#         for i in range(self.len):
#             res_snv.append('%s\t%s\t%s\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d\t%d\n' % (
#             self.chrs[i], self.poss[i], self.reffs[i],
#             self.p_values[i], self.p_homozyte[i],
#             self.allele_freq[i,0], self.allele_freq[i,1], self.allele_freq[i,2], self.allele_freq[i,3],
#             self.DP_watson[i], self.DP_crick[i]
#             ))
#         self.res_snv = res_snv

#     def format_vcf(self):
        
#         return None
    
#     def format_ouptput(self):
#         self.format_snv()
#         self.format_vcf()
            
#     def write_files(self, wrt_ctl: ControlWriteFile):
#         wrt_ctl.minus_task_num()
#         if self.res_snv is not None:
#             for l in self.res_snv:
#                 wrt_ctl.out_snv.write(l)
#         if self.res_vcf is not None:
#             for l in self.res_vcf:
#                 wrt_ctl.out_vcf.write(l)


# .snv format output
def format_snv(chrs, poss, reffs, p_values, p_homozyte, allele_freq, DP_watson, DP_crick):

    res_snv = []
    for i in range(len(chrs)):
        res_snv.append('%s\t%s\t%s\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d\t%d\n' % (
        chrs[i], poss[i], reffs[i],
        p_values[i], p_homozyte[i],
        allele_freq[i,0], allele_freq[i,1], allele_freq[i,2], allele_freq[i,3],
        DP_watson[i], DP_crick[i]
        ))
    return res_snv


def format_vcf(chrs, poss, reffs, p_values, p_homozyte, allele_freq, DP_watson, DP_crick):
    return None


def write_lines(res: SNVResultsBatch):
    global wrt_ctl

    wrt_ctl.minus_task_num()
    if res is None: return

    if res.snv is not None:
        wrt_ctl.out_snv.writelines(res.snv)
    if res.vcf is not None:
        wrt_ctl.out_vcf.writelines(res.vcf)


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
        
        # number of readed lines 
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


# @jit(nopython=True)
def BS_SNV_Caller_batch(lines: list, params: SNVparams):

    if len(lines) == 0:
        return None
    
    line_res = [l.split('\t') for l in lines]

    array = np.array(line_res)

    ## exclude Ns
    i = array[:,1] != 'N'
    if sum(i) < len(lines):
        array = array[i,:]
    
    # base representation of Crick strand in this model
    # is different from .ATCG file
    reads = array[:,(5,6,7,8, 11,10,13,12)].astype(int) 

    # sqrt-transform read counts
    
    i = reads > params.shrink_depth
    if i.sum() > 0:
        reads[i] = shrink_depth(reads[i], params.shrink_depth)

    # exclude sites of low coverage

    depth = np.sum(reads, axis=1)
    i = depth >= params.min_depth
    if i.sum() < len(lines):
        reads = reads[i, :]
        array = array[i, :]
        depth = depth[i]

    ## basic vars
    
    BASES = np.array(['A', 'T', 'C', 'G'])
    N_rows, _ = array.shape

    refs = array[:,1]
    patterns = array[:,3]

    is_CG = patterns == 'CG'

    ## posterior \prop likelihood (multinomial) \times prior

    post_mx = np.zeros((N_rows, 10))

    post_mx[is_CG, :] = reads[is_CG, :] @ params.loglik['CG']
    post_mx[np.logical_not(is_CG), :] = reads[np.logical_not(is_CG), :] @ params.loglik['CH']

    post_mx = np.exp(post_mx - np.max(post_mx, axis=1, keepdims=True))

    # priors

    prior_mx = np.zeros((N_rows, 10))

    for ref in BASES:
        prior_mx[refs == ref, :] = params.priors[ref]

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

    i_sig = p_value < params.pvalue

    # number of mutative sites
    n_sig = i_sig.sum()
    if n_sig == 0: return(None)

    # return singnificant sites only

    p_values = p_value[i_sig]
    chrs = array[i_sig, 0]
    poss = array[i_sig, 2]
    reffs = refs[i_sig]

    # allele frequencies
    allele_freq = post_mx[i_sig,:] @ params.allele_weights.T

    # strand coverage
    DP_watson = np.sum(reads[i_sig, :4], axis=1)
    DP_crick = np.sum(reads[i_sig, 4:8], axis=1)
    depths = depth[i_sig]

    # probability of homozygote
    p_homozyte = np.sum(post_mx[i_sig, :4], axis=1)


    # # .snv format output
    # res_snv = format_snv(chrs, poss, reffs, p_values, p_homozyte, allele_freq, DP_watson, DP_crick)

    # # .vcf format output
    # res_vcf = format_vcf(chrs, poss, reffs, p_values, p_homozyte, allele_freq, DP_watson, DP_crick)

    # res = SNVResultsBatch(chrs, poss, reffs, p_values, p_homozyte, allele_freq, DP_watson, DP_crick)
    # res.format_ouptput()

    res_snv = []
    for i in range(n_sig):
        res_snv.append('%s\t%s\t%s\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t%d\n' % (
        chrs[i], poss[i], reffs[i],
        p_values[i], p_homozyte[i],
        allele_freq[i,0], allele_freq[i,1], allele_freq[i,2], allele_freq[i,3],
        DP_watson[i], DP_crick[i]
        ))

    ## vcf results

    # BASES_dict = {'A':0, 'T':1, 'C':2, 'G':3}

    vcf_QUAL = np.fmin(512, np.int32(np.around(-10*np.log10(p_values))))

    res_vcf = []
    for i in range(n_sig):
        # except the ref
        a = BASES == reffs[i]
        # alt = np.logical_and([True]*4, ~a)
        alt = ~a
        order = np.argsort(-allele_freq[i, alt])
        BASES_ord = BASES[alt][order]
        PROP_ALELLE = allele_freq[i, alt][order]

        vcf_ID = '.'
        # if p_homozyte[i] < 0.8:
        b = PROP_ALELLE > 0.05
        if b.sum() == 0:
            vcf_ALT = '.'
        else:
            vcf_ALT = ','.join(BASES_ord[b])

        # vcf_FILTER = []
        # if p_values[i] < 0.01: vcf_FILTER.append('q30')
        # if depths[i] <= 6: vcf_FILTER.append('dp6')
        vcf_FILTER = 'PASS'

        # allele frequency
        vcf_af = ','.join(['%.2f' % x for x in PROP_ALELLE[b]])
        vcf_INFO = f'NS=1,DP={depths[i]},DPW={DP_watson[i]},DPC={DP_crick[i]}'
        if b.sum() != 0:
            vcf_INFO += f',AF={vcf_af}'

        vcf_FORMAT = 'GT:GQ:DP:DPW:DPC'
        
        # score of homozygote
        vcf_QUAL_HM = np.fmin(512, np.int32(np.around(-10*np.log10(1-p_homozyte))))
        # score of heterzygote
        vcf_QUAL_HT = np.fmin(512, np.int32(np.around(-10*np.log10(p_homozyte))))

        # genotype 
        if p_homozyte[i] > 0.5:
            vcf_GQ = vcf_QUAL_HM[i]
            if allele_freq[i, a] > PROP_ALELLE[0]:
                vcf_GT = '0'
                vcf_ALT = '.'
            else:
                vcf_GT = '1'
                vcf_ALT = BASES_ord[0]
        else:
            vcf_GQ = vcf_QUAL_HT[i]
            if allele_freq[i, a] > PROP_ALELLE[1]:
                vcf_GT = '0/1'
            else:
                vcf_GT = '1/2'

        vcf_SAMPLE = f'{vcf_GT}:{vcf_GQ}:{depths[i]}:{DP_watson[i]}:{DP_crick[i]}'

        res_vcf.append(f'{chrs[i]}\t{poss[i]}\t{vcf_ID}\t{reffs[i]}\t{vcf_ALT}\t{vcf_QUAL[i]}\t{vcf_FILTER}\t{vcf_INFO}\t{vcf_FORMAT}\t{vcf_SAMPLE}\n')

    return SNVResultsBatch(snv=res_snv, vcf=res_vcf)
    # return (res_snv, res_vcf)

def BS_SNV_CAller(options: Namespace):
    
    ############################
    ##

    # model params
    params = SNVparams(options)
    params.set_model_params()

    # control writing results
    global wrt_ctl 
    wrt_ctl = ControlWriteFile(params)
    ATCGfile = LineFile(params.infile, params.batch_size)

    # write header lines of vcf file
    wrt_ctl.out_vcf.writelines(params.VCF_HEADER)

    # multi-process parallelization
    with Pool(processes=params.num_process) as pool:
        while True: 
            if wrt_ctl.FLAG_WAIT == 'lower':
                line_batch = next(ATCGfile)
                if not line_batch: break

                wrt_ctl.add_task_num()
                pool.apply_async(BS_SNV_Caller_batch, (line_batch, params), callback=write_lines)

                # pool.apply_async(calculate, (postp, (line_batch, args)), callback=writeLine)
                # pool.apply_async(f, (int(id), name, float(x), float(y)), callback=writeLine)
                # results = [pool.apply_async(calculate, task) for task in TASKS]

            wrt_ctl.waitTime()

        pool.close()
        pool.join()
        ATCGfile.close()
        wrt_ctl.close()
    

def BS_SNV_CAller_keep_order(options: Namespace):
    
    ############################
    ##

    # model params
    params = SNVparams(options)
    params.set_model_params()

    # control writing results
    global wrt_ctl 
    wrt_ctl = ControlWriteFile(params)
    ATCGfile = LineFile(params.infile, params.batch_size)

    # write header lines of vcf file
    wrt_ctl.out_vcf.writelines(params.VCF_HEADER)

    FLAG_LAST_BATCH = False
    
    # multi-process parallelization
    with Pool(processes=params.num_process) as pool:
        while True: 
            if FLAG_LAST_BATCH: break
            line_batches = []
            i = 0
            for _ in range(params.pool_upper_num*params.num_process):
                line_batch = next(ATCGfile)
                if not line_batch: 
                    FLAG_LAST_BATCH = True
                    break
                else:
                    line_batches.append(line_batch)
                    i += 1

            if i == 0: break
            # pool.apply_async(BS_SNV_Caller_batch, (line_batch, params), callback=write_lines)

            # pool.apply_async(calculate, (postp, (line_batch, args)), callback=writeLine)
            # results = pool.apply_async(BS_SNV_Caller_batch, (params,))
            results = [pool.apply_async(BS_SNV_Caller_batch, (lines, params)) for lines in line_batches]

            for res in results:
                write_lines(res.get())


        pool.close()
        pool.join()
        ATCGfile.close()
        wrt_ctl.close()

def as_bool(x: str):
    x = x.upper()
    if x == 'TRUE': return True
    if x == 'False': return False
    return None

if __name__ == '__main__':

    # parse command line
    
    usage = 'Usage: BS-SNA-Caller.py -i sample.atcg.gz [options]'
    desc = 'SNV caller with bisulite-converted sequencing data'

    parser = ArgumentParser(description=desc)
    parser.add_argument('-i', '--atcg-file', dest='infile', help='an input .atcg[.gz] file', type=str, required=True)
    parser.add_argument('-o', '--output-prefix', dest='outprefix', help='prefix of output files, a prefix.snv.gz and a prefix.vcf.gz will be returned, by default, same with input filename except suffix, say for input of path/sample.atcg.gz, the output is path/sample.snv.gz and path/sample.vcf.gz which is equilant to setting -o path/sample', type=str)
    parser.add_argument('-m', '--mutation-rate', dest='mutation_rate', help='mutation rate a hyploid base is different with reference base', type=float, default=0.001)
    parser.add_argument('-e', '--error-rate', dest='error_rate', help='error rate a base is misdetected due to sequencing or mapping', type=float, default=0.03)
    parser.add_argument('-c', '--methy-cg', dest='methy_cg', help='Cytosine methylation rate of CpG-context', type=float, default=0.6)
    parser.add_argument('-n', '--methy-ch', dest='methy_ncg', help='Cytosine methylation rate of non-CpG-context', type=float, default=0.01)
    parser.add_argument('-d', '--min-depth', dest='min_depth', help='sites with coverage depth less than min DP will be skipped', type=int, default=10)
    parser.add_argument('-p', '--pvalue', dest='pvalue', help='p-value threshold', type=float, default=0.01)
    parser.add_argument('--shrink-depth', dest='shrink_depth', help='sites with coverage larger than this value will be shrinked by a square-root transform', type=int, default=60)
    parser.add_argument('--batch-size', dest='batch_size', help='a batch of sites will be processed at the same time', type=int, default=100000)
    parser.add_argument('-P', '--num-process', dest='num_process', help='number of processes in parallel', type=int, default=4)
    parser.add_argument('--pool-lower-num', dest='pool_lower_num', help='lower number of bacthes in memory pool per process', type=int, default=10)
    parser.add_argument('--pool-upper-num', dest='pool_upper_num', help='upper number of bacthes in memory pool per process', type=int, default=30)
    parser.add_argument('--keep-order', dest='sort', help='keep the results same order with input', type=as_bool, default=True)

    options = parser.parse_args()

    if options.sort:
        BS_SNV_CAller_keep_order(options)
    else:
        BS_SNV_CAller(options)
