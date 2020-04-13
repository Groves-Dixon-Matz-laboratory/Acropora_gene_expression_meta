#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
x0=read.table(infile, sep='\t', header =T)
x=x0[,7:ncol(x0)]
rownames(x) = x0$Geneid
sums = data.frame(apply(x, 2, sum))
sums$stat = 'geneCounted'
write.table(sums, file='gene_count_sumsFC.tsv', sep='\t', quote=F, row.names=T, col.names=F)
