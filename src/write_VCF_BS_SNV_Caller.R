

library(argparser, quietly=TRUE)
library(readr)
library(dplyr)
library(magrittr)

## parse arguments

arg <- arg_parser("transform the results of BS-SNV-Caller to .vcf file")

arg <- add_argument(arg, "--input", help="input", type="character", default = 'STDIN')
# arg <- add_argument(arg, "--header", help="files include headers ?", type="logic", default = TRUE)
# arg <- add_argument(arg, "--remove_suffix", help="remove suffix in the additional column of sample", type="logic", default = TRUE)
# arg <- add_argument(arg, "--remove_suffix_greedy", help="remove suffix as much as possible", type="logic", default = FALSE)
# 
# arg <- add_argument(arg, "--delim", help="delimiter", type="character", default = ',')

arg <- add_argument(arg, "--output", help="output", type="character", default = 'STDOUT')

argv <- parse_args(arg)

if (argv$input == 'STDIN'){
  d = read_delim(stdin(), col_names = F, 
                 col_types = cols(
                   X1 = col_character())
  )
}else{
  d = read_delim(argv$input, col_names = F, 
                 col_types = cols(
                   X1 = col_character())
  )
}

# d = read_delim('AF5_12.bssnv.gz', col_names = F, 
#                col_types = cols(
#                  X1 = col_character())
#                )


# t = c('A', 'T', 'C', 'G')[freq]

# a = d[,7:10] %>% as.matrix
# t = rowRanks(-as.matrix(d[,7:10]), ties.method = 'first')

ALT = apply(d[,7:10] >= 0.1, 1, function(x){
  paste(c('A', 'T', 'C', 'G')[x], collapse = ',')
})

d2 = d %>% transmute(
  `#CHROM` = X1,
  POS = X3,
  ID = '.',
  REF = X2,
  ALT = ALT,
  QUAL = pmin(255, ceiling(-log10(X6))),
  FILTER = 'PASS',
  INFO = paste0('NS=1:DP=', X11 + X12),
  FORMAT = 'GT:GQ:DP',
  sample = paste(if_else(X13 > 0.5, '1/1', '1/2'),
                 0,
                 X11 + X12,
                 sep = ':')
)

regexp = '\\..+$'
sample = sub(regexp, '', basename(argv$input))
colnames(d2)[10] = sample

header = sprintf(
  '##fileformat=VCFv4.2
##fileDate=%s
##source=BS-SNV-Caller
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FILTER=<ID=q10,Description="Quality < 3 and DP < 10">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
Sys.Date()
)


if(argv$output == 'STDOUT'){
  write(header, stdout())
  write_delim(d2, stdout(), append = T, col_names = T, delim = '\t')
}else{
  write(header, argv$output)
  write_delim(d2, argv$output, append = T, col_names = T, delim = '\t')
}

# argv$output = 'AF5_12.bssnv.vcf'
