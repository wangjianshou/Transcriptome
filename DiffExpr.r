#usage: Rscript --input readCountMatrix.xls --group group1:sample1,sample2,sample3 group2:sample4,sample5,sample6 --compare group1_vs_group2
library(argparser, lib.loc="/USCIMD/usr/wangjianshou/software/R3.2.2_lib")
library(DESeq2)
p <- arg_parser("This script is used for differential expression analysis")
p <- add_argument(p, "--input", help="input file, reads count matrix", nargs=1, short='-i')
p <- add_argument(p, "--group", help="groupName:sample1Name,sample2Name ...", nargs=Inf, short="-g")
p <- add_argument(p, "--compare", help="group1Name_vs_group2Name ...", nargs=Inf, short='-c')
argv <- parse_args(p)

readCountMatrix <- read.table(argv$input, header=T, sep='\t', check.names=FALSE)
dds <- DESeqDataSetFromMatrix(readCountMatrix, )
