#uasge
#Rscript diff_A_vs_B.r ReadCountMatrix.xls group1:sample1,sample2,... group2:sample3,sample4,...
library(DESeq2)
args <- commandArgs(TRUE)
ReadCountMatrix <- read.table(args[1], header=T, sep='\t', quote='', check.names=F)
ReadCountMatrix <- ReadCountMatrix[rowSums(ReadCountMatrix[,2:length(ReadCountMatrix)]) != 0,]
groupA <- unlist(strsplit(args[2], split="[:,]"))
groupB <- unlist(strsplit(args[3], split="[:,]"))
datacol <- ifelse(colnames(ReadCountMatrix)[-1] %in% groupA[-1], groupA[1], groupB[1])
datacol <- data.frame(condition=datacol)
dds <- DESeqDataSetFromMatrix(ReadCountMatrix, colData=datacol, design=~condition, tidy=TRUE)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", groupA[1], groupB[1]), alpha=0.05)
res <- res[!is.na(res$padj), ]
gene_info <- read.table("/USCIMD/usr/wangjianshou/tr/database/gene_info.xls", sep='\t', header=T)
res_table  <- as.data.frame(res[res$padj < 0.05, ])
gene_info_diff <- gene_info[gene_info$gene_id %in% rownames(res_table),]
rownames(gene_info_diff) <- gene_info_diff$gene_id
result <- cbind(gene_info_diff,res_table[rownames(gene_info_diff),])
write.table(result, file=args[4], sep='\t', col.names=T, row.names=F, quote=F)
