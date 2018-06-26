#usage: Rscript getReadCountMatrix.r outdir sample1:readCountFile sample2:readCountFile ...
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
transcript_length_mean <- fread("/USCIMD/usr/wangjianshou/tr/database/transcript_length_mean", header=T, sep='\t')
for(sample in args[-1])
{
  sampleName <- unlist(strsplit(sample, split=":"))[1]
  readcount <- fread(unlist(strsplit(sample, split=":"))[2], header=F, sep='\t')
  setnames(readcount, c("Gene", sampleName))
  ReadCountMatrix <- if(exists("ReadCountMatrix")) ReadCountMatrix[readcount, on="Gene"] else readcount
}
rr <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
ReadCountMatrix <- ReadCountMatrix[!(Gene %in% rr)]
write.table(ReadCountMatrix, file=file.path(args[1], "ReadCountMatrix.xls"), sep='\t', row.names=F, quote=F, col.names=T)

ReadCountMatrix <- ReadCountMatrix[rowSums(ReadCountMatrix[,2:length(ReadCountMatrix)]) != 0,]
setkey(transcript_length_mean, gene_id)
setkey(ReadCountMatrix, Gene)

#TPM
TPM <- transcript_length_mean[ReadCountMatrix]
TPM[, names(TPM)[3:dim(TPM)[2]] := lapply(.SD, function(x) x/transcript_length_mean), .SDcols=names(TPM)[3:dim(TPM)[2]]]
TPM[, names(TPM)[3:dim(TPM)[2]] := lapply(.SD, function(x) x/sum(x)*10^6), .SDcols=names(TPM)[3:dim(TPM)[2]]]
TPM[, transcript_length_mean:=NULL]
write.table(TPM, file=file.path(args[1], "TPM.xls"), sep='\t', row.names=F, col.names=T, quote=F)

#FPKM
FPKM <- transcript_length_mean[ReadCountMatrix]
FPKM[, names(FPKM)[3:dim(FPKM)[2]] := lapply(.SD, function(x) x/transcript_length_mean/sum(x)*10^9), .SDcols=names(FPKM)[3:dim(FPKM)[2]]]
FPKM[, transcript_length_mean:=NULL]
write.table(FPKM, file=file.path(args[1], "FPKM.xls"), sep='\t', row.names=F, col.names=T, quote=F)
