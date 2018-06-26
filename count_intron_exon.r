#usage:
#Rscript count_intron_exon.r transcript.count exon.count outfile
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
outfile <- args[3]
transcript <- fread(args[1], header=F)
exon <- fread(args[2], header=F)
transcript_reads <- sum(transcript$V2[1:(dim(transcript)[1]-5)])
exon_reads <- sum(exon$V2[1:(dim(exon)[1]-5)])
exon_ambiguous <- exon[dim(exon)[1]-3]$V2
transcript_ambiguous <- transcript[dim(transcript)[1]-3]$V2
intron_reads <- transcript_reads + transcript_ambiguous - exon_ambiguous - exon_reads
no_feature <- transcript[dim(transcript)[1]-4]$V2
cat("exon_reads", "intron_reads", "__no_feature", sep="\t", end="\n", file=outfile)
cat(exon_reads, intron_reads, no_feature, file=outfile, sep='\t', end="\n", append=T)
all_reads <- exon_reads + intron_reads + no_feature
cat(exon_reads/all_reads, intron_reads/all_reads, no_feature/all_reads, file=outfile, sep='\t', end="\n", append=T)

