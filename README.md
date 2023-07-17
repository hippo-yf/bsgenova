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
|`-p`, `--pvalue` |`float` | p-value threshodl |0.01|
|`--shrink-depth` |`integer` | sites with coverage larger than this value will be shrinked by a square-root transform | 60|
|`--batch-size` |`integer` | a batch of genomic sites will be processed at the same time |100000|
|`-P`, `--num-process` |`integer` | number of processes in parallel |4|
|`--pool-lower-num` |`integer` | lower number of bacthes in memory pool per process |10|
|`--pool-upper-num` |`integer` | upper number of bacthes in memory pool per process |30|
|`--keep-order` |`logical` | keep the results same order with input, if the order of sites makes no difference, set `False` to enable faster non-blocking asynchronous IO |True|
|`-h`, `--help` | | show this help message and exit ||

## Performance

*BS-SNV-Caller* implementes parallelism with `multiprocessing`, the default number of processes is 4, at which the acceratation is about 3.8 times than single-process run.

By maintaining a pool in memory, the maximum memory usage can be limited by parameter `--pool-upper-num`. The maxumun lines/sites is  `num-process` $\times$ `pool-upper-num` $\times$ `batch-size`. With defaults: `--pool-upper-num 30`, `--num-process 4`, and `--batch-size 100000`  the memory usage is ~1GB.

