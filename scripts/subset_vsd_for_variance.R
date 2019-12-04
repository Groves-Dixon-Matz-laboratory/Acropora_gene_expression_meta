#!/usr/bin/env Rscript
#subset_vsd_for_variance.R


#PARSE ARGUMENTS
library(optparse)
option_list = list(
  
  make_option(c("--i"), type="character", default=NULL, 
              help="The R oject with datExpr and datTraits"),
  
  make_option(c("--v"), type="double", default=0.8, 
              help="The amount of the data you want based on variacne. 
              Eg: 0.8 means you want 80% of the data, 
              with all but the lowest 20% variance genes."),
  make_option(c("--e"), type="double", default=0.8, 
              help="The percentage of data you want based on expression. 
              Eg: 0.8 means the top 80% of the data, all but the lowest 20% expression genes.")
)



print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
infile = opt$i
varCut=opt$v
expCut=opt$e
prefix = sub('.Rdata', '', infile)
varString = paste('_varCut', varCut, sep='')
expString = paste('_expCut', expCut, sep='')
outfile = paste(c(prefix, varString, expString, '.Rdata'), collapse='')
ll=load(infile)
ll





library(DESeq2)

print('Subsetting for variacne...')
rv <- rowVars(t(datExpr))
genes=colnames(datExpr)
ntopV=round(varCut*length(genes), digits=0)
selectV <- order(rv, decreasing = TRUE)[seq_len(min(ntopV,length(rv)))]
svGenes=genes[selectV]

print('Subsetting for expression level...')
medLvl = apply(datExpr, 2, median)
ntopE = round(expCut*length(genes), digits=0)
selectE <- order(medLvl, decreasing = TRUE)[seq_len(min(ntopE,length(rv)))]
seGenes = genes[selectE]


#get the genes that are both
select = intersect(svGenes, seGenes)
nsel = length(select)
psel = round(nsel/length(genes)*100, digits=2)
print(paste(nsel, 'total genes passed both cutoffs'))
print(paste(psel, '% of total'))
datExpr=datExpr[,select]


print(paste('Saving results as', outfile))
save(datExpr, datTraits, file=outfile)