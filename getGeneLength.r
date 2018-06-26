#usage:
#Rscript getGeneLength.r *gtf
args <- commandArgs(TRUE)
gtf <- fread(args[1], header=F, sep='\t')
gtf <- gtf[V3=="exon"]
gtf[, `:=`(gene_id=sub(".*gene_id \"(.*?)\";.*", "\\1", V9),
           transcript_id=sub(".*transcript_id \"(.*?)\";.*", "\\1", V9),
           gene_biotype=sub(".*gene_biotype \"(.*?)\";.*", "\\1", V9),
           gene_name=sub(".*gene_name \"(.*?)\";.*", "\\1", V9),
           V9=NULL)
   ]
gtf[, transcript_length:=sum(abs(V5-V4)+1), by=.(gene_id, transcript_id)]
gtf <- gtf[, .SD[1], by=.(gene_id, transcript_id)]
gtf[, transcript_length_mean:=mean(transcript_length), by=gene_id]
gtf <- gtf[, .SD[1], by=gene_id]
write.table(gtf[, .(gene_id, transcript_length_mean)], file="transcript_length_mean", col.names=T, row.names=F, sep='\t', quote=F)

