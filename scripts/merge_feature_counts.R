#!/usr/bin/env Rscript
#merge_feature_counts.R


#PARSE ARUGMENTS
library(optparse)
option_list = list(
  
  make_option(c("--i"), type="character", default="^feature_counts_out_.*.txt$", 
              help="glob to infiles"),
  
  make_option(c("--pth"), type="character", default='.', 
              help="Path to data file directory"),
  
  make_option(c("--o"), type="character", default='all_featureCounts_geneCounts.tsv', 
              help="Name for output file")
  
)



print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
datGlob = opt$i
datPath = opt$pth
outName = opt$o


suppressMessages(library(tidyverse))

#get files names
fileList = list.files(path = datPath, pattern = datGlob)
print('Merging the following feature counts files:')
print(paste(length(fileList), 'total'))
for (f in fileList){
  print(f)
}


#read in function
read_dat = function(x){
  print(paste(x,'...',sep=''))
  read.table(x, sep="\t", header = TRUE) %>% 
    as_tibble()
}


#read in and join the datasets
print('Reading in and joining...')
dat <- fileList %>%
  map(function(x) read_dat(x)) %>% 
  reduce(full_join, by = c('Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'))
print('Done.')

#print summary
nGenes = nrow(dat)
nSample = ncol(dat)-6
print(paste(nGenes, 'total genes'))
print(paste(nSample, 'total samples'))

#get gene count totals
counts = dat %>% 
  select(-seq(1,6))
tots = data.frame(apply(counts, 2, function(x) sum(x, na.rm=T)))
tots$stat='geneCountedFCs'
print('Writing out total reads counted on genes as featureCounts_totalGeneCounted.tsv...')
write.table(tots, file='featureCounts_totalGeneCounted.tsv', sep="\t", row.names=T, col.names=F, quote=F)


#write out
print(paste(c('Writing out results to file ', outName, '...'), collapse=''))
write_tsv(dat, path=outName)