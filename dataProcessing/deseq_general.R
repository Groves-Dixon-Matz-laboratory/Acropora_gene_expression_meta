#!/usr/bin/env Rscript

library(optparse)
options(stringsAsFactors = FALSE) #important for WGCNA

option_list = list(
  
  make_option(c("--i"), type="character", default='deseqBaselineInput.Rdata', 
              help="Input file"),
  make_option(c("--v"), type="character",
              help="variable to compare between."),
  make_option(c("--t"), type="character",
              help="variable value indicating treated samples"),
  make_option(c("--c"), type="character",
              help="variable value indicating control samples"),
  make_option(c("--pt"), type="character", default='using_all_types',
              help="Indicate which project type in coldata to use. Leave blank to use all types in coldata table"),
  make_option(c("--runInd"), type="logical", default=FALSE,
              help="Logical for whether to run individual projects. Only need to both doing with for the full stress dataset once")

)

print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
infile = opt$i
variable = opt$v
treatString = opt$t
controlString = opt$c
projType = opt$pt
runIndividual = opt$runInd


#load the data
library(DESeq2)
library(tidyverse)
ll = load(infile)
ll
coldata
head(counts)
sum(colnames(counts) == rownames(coldata))==ncol(counts)


#----- format coldata
print('Unique treatments:')
print(unique(coldata$treat))


if (projType != 'using_all_types'){
  print(paste('Only using samples with project type', projType))
  before = nrow(coldata)
  coldata = coldata[coldata$projType==projType,]
  counts = counts[,colnames(counts) %in% coldata$Run]
  after = nrow(coldata)
  diff = before - after
  print('total samples removed:')
  print(diff)
  print('double-checking. Tables still line up?')
  print(sum(colnames(counts)==coldata$Run)==ncol(counts))
  coldata[,variable]=factor(coldata[,variable]) #reset factor
  print('Unique treatments after subsetting:')
  print(unique(coldata$treat))
} else{
  print("--pt left blank. Will not subset the dataset before running deseq.")
}



# run for full dataset ----------------------------------------------------
formula = paste(c('~', variable, '+my_title'), collapse='')
dds<-DESeqDataSetFromMatrix(counts,
                            colData = coldata,
                            design = formula(formula))
dds <- DESeq(dds)
res = results(dds, contrast = c(variable, treatString, controlString))
save(res, file='all_deseqResults.Rdata')


# loop through projects individually --------------------------------------
if (runIndividual==TRUE){
  print('Running individual projects...')
  formula = paste('~', variable,sep='')
  projects = unique(coldata$my_title)
  for (proj in projects){
    print('----------')
    print(paste(c('Running project', proj,'...'), collapse=''))
    subColdata = coldata %>% 
      filter(my_title==proj)
    subCounts = counts %>% 
      select(subColdata$Run)
    dds<-DESeqDataSetFromMatrix(subCounts,
      colData = subColdata,
      design = formula(formula))
    dds <- DESeq(dds)
    res = results(dds, contrast = c(variable, treatString, controlString))
    outName = paste(proj,'deseqResults.Rdata', sep='_')
    save(res, file=outName)
  }
}