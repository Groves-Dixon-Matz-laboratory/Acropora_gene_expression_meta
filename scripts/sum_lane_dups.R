#!/usr/bin/env Rscript
#sum_lane_dups.R

#currently very specific for PRJNA319662


library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
dupsIn = args[1]
fcIn = args[2]
outName = args[3]

#read in the samples with lane duplicates (taken from SRA table, column 1 = runs column2 = sample names with lane dup info)
dups = read_csv(dupsIn, col_names=c('run', 'id')) %>% 
  separate(id, into=c('id1', 'id2', 'lane'), extra='drop') %>% 
  mutate(idBar = paste(id1, id2,sep='_'))

  
samples = unique(dups$idBar)
fcdat = read_tsv(fcIn)
fcOriginal = fcdat

ndups = c()
keptNames = c()
for (s in samples){
  sub = dups %>% 
    filter(idBar==s)
  runs = sub %>% 
    pull(run)
  ndups = append(ndups, length(runs))
  subfc = fcdat %>% 
    select(runs)
  fcdat = fcdat %>% 
    select(-runs)
  sums = apply(subfc, 1, sum)
  newcol = runs[1]
  fcdat[[newcol]]=sums
  keptNames = append(keptNames, newcol)
}


#check that everything worked
resSum = data.frame('sample'=samples, 'nlaneDups'=ndups, 'keptName'=keptNames)
oldTot = ncol(fcOriginal)
newTot = ncol(fcdat)
expectedLoss = sum(resSum$nlaneDups-1)
print(paste('Original number of samples included in table=', oldTot-6))
print(paste('After merging lane dups, new total =', newTot-6))
print('Checking that correct number of columns remain after summing lane duplicates...')
right = newTot+expectedLoss == oldTot
if (right){
  print('All looks good')
} else{
  print('No! Does not match!')
}

#write out
write_tsv(resSum, path='./laneDupSummary.tsv')
write_tsv(fcdat, path=outName)

