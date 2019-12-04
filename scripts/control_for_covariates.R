#!/usr/bin/env Rscript
library(optparse)



#upload output from initialize_raw_featureCounts

library(limma)
print('Reading in data...')

#parse input file name
option_list = list(
  
  make_option(c("--i"), type="character", default=NULL, 
              help="infile"),
  make_option(c("--o"), type="character", default=NULL, 
              help="out prefix")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
infileName = opt$i
outPrefix = opt$o

#load
ll=load(infileName)
ll
toFix = t(datExpr)



# control for project -----------------------------------------------------
if (length(unique(datTraits$my_title)) > 1){
	print('Controlling for project...')
	contRes <- removeBatchEffect(toFix, batch=datTraits$my_title)
	datExpr=t(contRes)
	print('saving project controlled data as project_controlled.Rdata...')
	save(datExpr,datTraits, file=paste(outPrefix, 'project_controlled.Rdata', sep='_'))
}

# control for read counts -------------------------------------------------
print('Controlling for counts...')
countFix <- removeBatchEffect(toFix, covariates=log(datTraits$geneCounted, 10))
datExpr=t(countFix)
print('saving project controlled data as project_controlled.Rdata...')
save(datExpr,datTraits, file=paste(outPrefix, 'count_controlled.Rdata', sep='_'))

# control for both -------------------------------------------------
if (length(unique(datTraits$my_title)) > 1){
	print('Controlling for read counts on top of controlling for project...')
	contRes2 <- removeBatchEffect(contRes, covariates=log(datTraits$geneCounted, 10))
	datExpr=t(contRes2)
	print('saving project controlled data as project_and_count_controlled.Rdata...')
	save(datExpr,datTraits, file=paste(outPrefix, 'project_and_count_controlled.Rdata', sep='_'))
}

