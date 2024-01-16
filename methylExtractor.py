#!/usr/bin/env python

import pysam # ver 0.22.0
import math
import numpy as np
from typing import NamedTuple
# import array
import gzip
import io
import sys
from argparse import ArgumentParser

##

# 3-nucleotide context, CG/CHG/CHH etc.

CG_CONTEXT_FORWARD_HASH = {'CGA':'CG', 'CGT':'CG', 'CGC':'CG', 'CGG':'CG', # 4 CG
                           'CAG':'CHG', 'CTG':'CHG', 'CCG':'CHG',          # 3 CHG
                           'CAA':'CHH', 'CAT':'CHH', 'CAC':'CHH',          # 9 CHH
                           'CTA':'CHH', 'CTT':'CHH', 'CTC':'CHH',
                           'CCA':'CHH', 'CCT':'CHH', 'CCC':'CHH'
                           }

CG_CONTEXT_REVERSE_HASH = {'ACG':'CG', 'TCG':'CG', 'CCG':'CG', 'GCG':'CG', # 4 CG
                           'CAG':'CHG', 'CTG':'CHG', 'CGG':'CHG',          # 3 CHG
                           'AAG':'CHH', 'ATG':'CHH', 'AGG':'CHH',          # 9 CHH
                           'TAG':'CHH', 'TTG':'CHH', 'TGG':'CHH',
                           'GAG':'CHH', 'GTG':'CHH', 'GGG':'CHH'
                           }

# 2-nucleotide context of reverse strand

DI_CONTEXT_REVERSE_HASH = {'AG':'CT', 'TG':'CA', 'CG':'CG', 'GG':'CC'}


def as_bool(x: str):
    x = x.upper()
    if x in ["TRUE", "T", "YES", "Y"] : return True
    if x in ["FALSE", "F", "NO", "N"] : return False
    return None


# interval must be involved in single chr
class GenomicInterval(NamedTuple):
    chr: str
    chr_length: int
    start: int
    end: int

# `ref` includes ref bases from `start`-2 to `end`-2
class FaGenomicInterval(NamedTuple):
    chr: str
    start: int
    end: int
    bases: str


class GenomicIntervalGenerator:
    
    def __init__(self, 
                 fa: pysam.FastaFile, 
                 chrs,
                 start: int,
                 end: int,
                 step: int,
                 ) -> None:
        
        self.chrs = fa.references
        self.lens = fa.lengths
        if chrs == "all":
            self.chrs_selected = self.chrs
        else:
            self.chrs_selected = chrs
        
        self.start = start
        self.end = end
        self.step = step

        assert(step > 0 and start < end)
        assert(len(self.chrs) > 0 and len(self.chrs) == len(self.lens))


    def __repr__(self):
        return f'GenomicIntervalGenerator({len(self.chrs)} contig(s) with step {self.step})'
    
    def __iter__(self) -> GenomicInterval:
        for (chr, len) in zip(self.chrs, self.lens):
            if chr not in self.chrs_selected:
                continue

            end2 = min(self.end, len)
            start = self.start
            end = start + self.step
            while start < end2:
                # if start < end2:
                end = min(end, end2)
                yield GenomicInterval(chr=chr, chr_length=len, start=start, end=end)
                start = end
                end += self.step
                # else:
                #     break

class MyFastaFile(pysam.FastaFile):

    def rich_fetch(self, intrv: GenomicInterval, padding:int) -> FaGenomicInterval:

        bases = self.fetch(reference=intrv.chr, 
                           start=max(0, intrv.start-padding), 
                           end=max(intrv.chr_length, intrv.end+padding)
                           ).upper()
        
        # padding N
        if intrv.start < padding:
            bases = "N"*(padding-intrv.start) + bases
        if intrv.end+padding > intrv.chr_length:
            bases = bases + "N"*(intrv.end+padding-intrv.chr_length)

        return bases


class Coverage(NamedTuple):
    watson: np.array
    crick: np.array


class Parameters(NamedTuple):
    fafile: str
    bamfile: str
    out_atcg: str
    chr: str
    start: int
    end: int
    step: int
    quality_threshold: int
    context_size: int
    coordinate_base: int  # 0/1-based
    swap_strand: bool # swap the counts of forward and reverse strands
    read_quality: int

# def reverse_read(read) -> bool:
#     return read.is_reverse

# def forward_read(read) -> bool:
#     return not read.is_reverse

def check_read(forward_read: bool, read_quality: int):
    def valid_read(read: pysam.AlignedSegment):
        return (forward_read ^ read.is_reverse) and (not read.is_unmapped) and (not read.is_duplicate) and (not read.is_secondary) and (not read.is_qcfail) and (read.mapping_quality>=read_quality)
    return valid_read


class MyAlignmentFile(pysam.AlignmentFile):

    # def __init__(self, )

    def Watson_Crick_coverage(self, intrv: GenomicInterval, 
                              params: Parameters
                              ) -> Coverage:
        cov_watson = self.count_coverage(contig=intrv.chr, 
                                         start=intrv.start, 
                                         stop=intrv.end, 
                                         quality_threshold=params.quality_threshold,
                                         read_callback=check_read(forward_read=True, read_quality=params.read_quality)
                                         )
        cov_crick = self.count_coverage(contig=intrv.chr, 
                                         start=intrv.start, 
                                         stop=intrv.end, 
                                         quality_threshold=params.quality_threshold,
                                         read_callback=check_read(forward_read=False, read_quality=params.read_quality)
                                         )
        
        # in some bams, the read strandness seem be reversly flaged in `FLAG`
        # parsed in read.is_reverse
        # for example gemBS

        if params.swap_strand:
            return Coverage(crick=np.array(cov_watson), watson=np.array(cov_crick))
        else:
            return Coverage(watson=np.array(cov_watson), crick=np.array(cov_crick))



def methylExtractor(params: Parameters) -> None:
    
    fa = MyFastaFile(params.fafile)
    bam = MyAlignmentFile(params.bamfile, 'rb')


    # outputs
    if params.out_atcg == '-':
        outfile_atcg = sys.stdout
    elif params.out_atcg.endswith('.gz'):
        outfile_atcg = gzip.open(params.out_atcg, 'wt')
    else:
        outfile_atcg = io.open(params.out_atcg, 'wt')

    intervals = iter(GenomicIntervalGenerator(fa, 
                                              chrs= params.chr, 
                                              start = params.start, 
                                              end = params.end,
                                              step=params.step
                                              ))

    for intrv in intervals:

        # context size
        consize = params.context_size
        # ref sequences
        bases = fa.rich_fetch(intrv, padding=consize-1).upper()


        # bam coverages
        try:
            covs = bam.Watson_Crick_coverage(intrv, 
                                            quality_threshold=params.quality_threshold, 
                                            swap_strand=params.swap_strand
                                            )
            cov_sum_W = np.sum(covs.watson, axis=0)
            cov_sum_C = np.sum(covs.crick, axis=0)
        except KeyError:
            continue
        else:
            pass

        # chr = intrv.chr
        # pos
        for i in range(intrv.end-intrv.start):
            if cov_sum_W[i]+cov_sum_C[i] == 0:
                continue
            

            j = i+2
            base = bases[j]
            if base == 'C':

                # CG/CHG/CHH
                bases_con = bases[j:(j+consize)]
                CG_context = '--' if 'N' in bases_con else CG_CONTEXT_FORWARD_HASH[bases_con]

                # dinucleatide context CA/CT/...
                dicontext = bases[j:(j+2)]
                nCT = covs.watson[1,i]+covs.watson[3,i]
                meth_ratio = covs.watson[1,i]/nCT if nCT>0 else math.nan

            elif base == 'G':
                bases_con = bases[(j-consize+1):(j+1)]
                CG_context = '--' if 'N' in bases_con else CG_CONTEXT_REVERSE_HASH[bases_con]

                bases2 = bases[(j-1):(j+1)]
                dicontext = '--' if 'N' in bases2 else DI_CONTEXT_REVERSE_HASH[bases2]

                nGA = covs.crick[2,i]+covs.crick[0,i]
                meth_ratio = covs.crick[2,i]/nGA if nGA>0 else math.nan
            else:
                CG_context = dicontext = "--"
                meth_ratio = math.nan

            # write files
            # compatiable with ATCGmap file foramt

            outfile_atcg.write(f'{intrv.chr}\t{base}\t{intrv.start+i+params.coordinate_base}\t{CG_context}\t{dicontext}\t{covs.watson[0,i]}\t{covs.watson[3,i]}\t{covs.watson[1,i]}\t{covs.watson[2,i]}\t0\t{covs.crick[0,i]}\t{covs.crick[3,i]}\t{covs.crick[1,i]}\t{covs.crick[2,i]}\t0\t{meth_ratio:.2f}\n')

    
    # close file handles

    fa.close()
    bam.close()
    outfile_atcg.close()


if __name__ == '__main__':


    # parse command line
    
    usage = 'Usage: methylExtrator -i sample.bam -g genome.fa -o sample.ATCGmap.gz [options]'
    desc = 'Extract ATCG (ATCGmap) and CG (CGmap) profiles from bam file'

    parser = ArgumentParser(description=desc)
    parser.add_argument('-b', '--bam-file', dest='in_bam', help='an input .bam file', type=str, required=True)
    parser.add_argument('-g', '--reference-genome', dest='in_fa', help='genome reference file .fa with index (.fai) in the same path', type=str, required=True)
    parser.add_argument('-a', '--output-atcgmap', dest='out_atcg', help='ATCGmap file', type=str, required=False, default='-')
    parser.add_argument('-c', '--chr', dest='chr', help='chromosomes/contigs', type=str, default='all')
    parser.add_argument('-s', '--start', dest='start', help='start coordinate of chromosomes/contigs', type=int, default=0)
    parser.add_argument('-e', '--end', dest='end', help='end coordinate of chromosomes/contigs', type=int, default=math.inf)
    parser.add_argument('--batch-size', dest='step', help='batch size of genomic intervals', type=int, default=2_000_000)

    parser.add_argument('--swap-strand', dest='swap_strand', help='swap read counts on two strands, true/false, or yes/no', type=as_bool, required=False, default='no')
    parser.add_argument('--base-quality', dest='base_quality', help='base sequencing quality threshold', type=int, default=15)
    parser.add_argument('--read-quality', dest='read_quality', help='read mapping quality threshold', type=int, default=20)
    parser.add_argument('--coordinate-base', dest='coordinate_base', help='0/1-based coordinate of output', type=int, default=1)

    options = parser.parse_args()

    params = Parameters(fafile=options.in_fa, 
                        bamfile=options.in_bam,
                        out_atcg=options.out_atcg,
                        chr=options.chr,  # 'all' for all chrs
                        start=options.start,
                        end=options.end,
                        # end=50_000_000,
                        step=options.step,
                        quality_threshold=options.base_quality, # base seq quality
                        read_quality=options.read_quality, # read mapping quality threshold
                        context_size=3, # size of CHG/...
                        coordinate_base=options.coordinate_base, # 0/1-based
                        swap_strand=options.swap_strand  # swap counts of two strands
                        )
    methylExtractor(params)



