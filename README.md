# BS-SNV-Caller
 
An **ultra-fast** and **accurate** Single-Nucleotide Variant (SNV) Caller with **bisulfite-converted** sequencing data for both single-cell and bulk samples.

- Enable synergistic analysis of genetic and epigenetic questions using BS-seq, RRBS, WGBS data and so on.

- Finish whole-genome SNV calling within an hour in a 4-core CPU and 1GB memory laptop

## Example

`python bs-snv-caller.py -i data/example.atcg.gz -o output/example`

It will generate the results: `output/example.snv.gz` and `output/example.vcf.gz`

## Installation and dependencies

- python >= 3.8
- numpy >= 1.13

Just git clone / download this repository, run in an environment with python3 and numpy.

## Usage and parameters

`python BS-SNV-Caller.py -i /path/to/atcgfile`

|**parameter** | **type** | **description**| **defaults** |
|  ----  | ----  | ----  | ----  |
| `-i`, `--atcg-file` | `string` | an input `.atcg[.gz]` file||
| `-o`, `--output-prefix` | `string` | prefix of output files, a prefix.snv.gz and a prefix.vcf.gz will be returned, by default, same with input filename except suffix, say for input of path/sample.atcg.gz, the output is `path/sample`.snv.gz and path/sample.vcf.gz which is equilant to setting -o `path/sample` | prefix of input file, this program will not create directoy automatically|
|`-m`, `--mutation-rate`| `float` | mutation rate a hyploid base is different with reference base |0.001 |
|`-e`, `--error-rate` | `float` |error rate a base is misdetected due to sequencing or mapping | 0.03|
|`-c`, `--methy-cg` |`float` | Cytosine methylation rate of CpG-context | 0.6|
|`-n`, `--methy-ch` |`float` |Cytosine methylation rate of non-CpG-context | 0.01|
|`-d`, `--min-depth` |`float` | sites with coverage depth less than min DP will be skipped | 10|
|`-p`, `--pvalue` |`float` | *p*-value threshold |0.01|
|`--shrink-depth` |`integer` | sites with coverage larger than this value will be shrinked by a square-root transform | 60|
|`--batch-size` |`integer` | a batch of genomic sites will be processed at the same time |100000|
|`-P`, `--num-process` |`integer` | number of processes in parallel |4|
|`--pool-lower-num` |`integer` | lower number of bacthes in memory pool per process |10|
|`--pool-upper-num` |`integer` | upper number of bacthes in memory pool per process |30|
|`--keep-order` |`logical` | keep the results same order with input, if the order of sites makes no difference, set `False` to enable faster non-blocking asynchronous IO |True|
|`-h`, `--help` | | show this help message and exit ||

## Input

a .ATCG file, returned by `bsseeker2/3` or `cgmaptools`

## Output

### .snv file

|**column** | description|
|  ----  | ----  |
|1| chromosome|
|2| position|
|3| reference base|
|4| posterior probability of not a SNV (different from reference)|
|5| posterior probability of a homozygote |
|6-9| posterior allele frequencies of A,T,C, and G respectively |
|10-11| coverage depths of Watson and Crick strands respectively|

an example

```
1	1023917	G	4.56e-27	9.98e-14	1.88e-20	6.71e-09	5.00e-01	5.00e-01	20	13
1	1023921	C	3.26e-53	9.97e-01	1.32e-03	9.31e-09	7.28e-07	9.99e-01	20	13
1	1024083	A	1.43e-16	8.29e-01	8.56e-02	3.86e-06	3.86e-06	9.14e-01	10	8
1	1024085	C	2.30e-09	1.00e-04	4.59e-10	5.00e-01	5.00e-01	4.59e-10	12	8
1	1024093	C	1.79e-18	1.26e-10	5.00e-01	1.99e-06	5.00e-01	4.72e-14	13	8
1	1024131	C	8.62e-06	7.54e-02	8.84e-08	5.38e-01	4.62e-01	8.84e-08	15	4
1	1024399	G	2.80e-32	9.99e-01	1.00e+00	7.61e-08	7.61e-08	4.33e-04	18	5
1	1024462	A	3.17e-36	1.00e+00	1.10e-04	4.97e-09	4.97e-09	1.00e+00	20	7
1	1024652	G	8.20e-03	8.20e-03	4.96e-01	5.53e-11	5.53e-11	5.04e-01	21	5
1	1024664	C	5.64e-10	4.81e-01	7.25e-08	7.40e-01	2.60e-01	7.25e-08	18	5
```

### .vcf file

Standard .vcf file.

an example (header lines are ommited)

```
##fileformat=VCFv4.4
##fileDate=20230717
##source=BS-SNV-Caller
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=DPW,Number=1,Type=Integer,Description="Total Depth of Wastson Strand">
##INFO=<ID=DPC,Number=1,Type=Integer,Description="Total Depth of Crick Strand">
##FILTER=<ID=q30,Description="Quality < 30">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality.In some cases of single-stranded coverge, we are sure there is a SNV, but we can not determine the alternative variant. So, we express the GQ as the Phred score (-10log10 (p-value)) of posterior probability of homozygote/heterozygote, namely, Prob(heterozygote) for homozygous sites and Prob(homozygote) for heterozygous sites. This is somewhat different with SNV calling from WGS data.">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=DPW,Number=1,Type=Integer,Description="Read Depth of Wastson Strand">
##FORMAT=<ID=DPC,Number=1,Type=Integer,Description="Read Depth of Crick Strand">
#CHROM	POS	ID	REF	ALT		QUAL	FILTER	INFO	FORMAT	example
1	1023917	.	G	C	263	PASS	NS=1,DP=33,DPW=20,DPC=13,AF=0.50	GT:GQ:DP:DPW:DPC	0/1:130:33:20:13
1	1023921	.	C	G	512	PASS	NS=1,DP=33,DPW=20,DPC=13,AF=1.00	GT:GQ:DP:DPW:DPC	1:26:33:20:13
1	1024083	.	A	G	158	PASS	NS=1,DP=18,DPW=10,DPC=8,AF=0.91	GT:GQ:DP:DPW:DPC	1:8:18:10:8
1	1024085	.	C	T	86	PASS	NS=1,DP=20,DPW=12,DPC=8,AF=0.50	GT:GQ:DP:DPW:DPC	0/1:40:20:12:8
1	1024093	.	C	A	177	PASS	NS=1,DP=21,DPW=13,DPC=8,AF=0.50	GT:GQ:DP:DPW:DPC	0/1:99:21:13:8
1	1024131	.	C	T	51	PASS	NS=1,DP=19,DPW=15,DPC=4,AF=0.54	GT:GQ:DP:DPW:DPC	0/1:11:19:15:4
1	1024399	.	G	A	316	PASS	NS=1,DP=23,DPW=18,DPC=5,AF=1.00	GT:GQ:DP:DPW:DPC	1:31:23:18:5
1	1024462	.	A	G	355	PASS	NS=1,DP=27,DPW=20,DPC=7,AF=1.00	GT:GQ:DP:DPW:DPC	1:37:27:20:7
1	1024652	.	G	A	21	PASS	NS=1,DP=26,DPW=21,DPC=5,AF=0.50	GT:GQ:DP:DPW:DPC	0/1:21:26:21:5
1	1024664	.	C	T	92	PASS	NS=1,DP=23,DPW=18,DPC=5,AF=0.74	GT:GQ:DP:DPW:DPC	0/1:3:23:18:5
```


## Performance

*BS-SNV-Caller* implementes parallelism with `multiprocessing`, the default number of processes is 4, at which the acceratation is about 3.8 times faster than single-process run.

By maintaining a pool in memory, the maximum memory usage can be limited by parameter `--pool-upper-num`. The maxumun lines/sites is  `num-process` $\times$ `pool-upper-num` $\times$ `batch-size`. With defaults: `--pool-upper-num 30`, `--num-process 4`, and `--batch-size 100000`  the memory usage is ~1GB.

