#!/usr/bin/env Rscript
#deseq_BEWW.R
#this script runs differential expression on each of 
#the experiments from this study individually
#The experiments were:
#cold
#hyposalinity
#heat1
#heat2
#three-way

#load the data
library(DESeq2)
library(tidyverse)
ll=load('beww_experiment/bewwDeseq/deseqBaselineInput.Rdata')
ll

#make general formats to coldata
coldata=coldata %>% 
  separate(treatDescription, into=c('experiment','genotype', 'treatment', 'timepoint', 'replicate'))
head(coldata)

# run for cold ----------------------------------------------------
#cold experiment has 12 treated and 12 control
#4 genotypes, single timepoint
sub.coldata = coldata %>% 
  filter(grepl('^Ch', experiment))
table(sub.coldata$treatment)
table(sub.coldata$replicate)
dim(sub.coldata)
head(sub.coldata)
sub.counts = counts[,sub.coldata$Run]
sum(colnames(sub.counts)==sub.coldata$Run)==ncol(sub.counts)

dds<-DESeqDataSetFromMatrix(sub.counts,
                            colData = sub.coldata,
                            design = ~treatment +  genotype)
dds <- DESeq(dds)
res = results(dds, contrast = c('treatment', 'T', 'C'))
save(res, file='deseqResults/individual_projects/bewwCold_deseqResults.Rdata')


# run for hyposalinity ----------------------------------------------------
#4 genotypes, two timepoints
#some are missing
sub.coldata = coldata %>% 
  filter(experiment=='Ls')
table(sub.coldata$treatment)
table(sub.coldata$replicate)
table(sub.coldata$timepoint)
dim(sub.coldata)
head(sub.coldata)
sub.counts = counts[,sub.coldata$Run]
sum(colnames(sub.counts)==sub.coldata$Run)==ncol(sub.counts)

dds<-DESeqDataSetFromMatrix(sub.counts,
                            colData = sub.coldata,
                            design = ~treatment + genotype + timepoint)
dds <- DESeq(dds)
res = results(dds, contrast = c('treatment', 'T', 'C'))
save(res, file='deseqResults/individual_projects/bewwHyposalinity_deseqResults.Rdata')

# run for heat ----------------------------------------------------
#2 experimental rounds, 4 genotypes, two timepoints each
#some are missing
sub.coldata = coldata %>% 
  filter(experiment %in% c('Ht1', 'Ht2') | (experiment=='Tw' & treatment=='C' & timepoint=='tp1')) %>% 
  mutate(experiment=if_else(experiment=='Tw',
                            'Ht2',
                            experiment))
#check stats
sub.coldata %>% 
  group_by(experiment, genotype, treatment, timepoint) %>% 
  summarize(N=n()) %>% 
  data.frame()


sub.counts = counts[,sub.coldata$Run]
sum(colnames(sub.counts)==sub.coldata$Run)==ncol(sub.counts)

dds<-DESeqDataSetFromMatrix(sub.counts,
                            colData = sub.coldata,
                            design = ~treatment + experiment + genotype + timepoint)
dds <- DESeq(dds)
res = results(dds, contrast = c('treatment', 'T', 'C'))
save(res, file='deseqResults/individual_projects/bewwHeat_deseqResults.Rdata')


# run for multi ----------------------------------------------------
#2 experimental rounds, 4 genotypes, two timepoints each
#some are missing

sub.coldata = coldata %>% 
  filter(experiment %in% c('Tw', 'HtC', 'CtH'))

#check stats
sub.coldata %>% 
  group_by(experiment, genotype, treatment, timepoint) %>% 
  summarize(N=n()) %>% 
  data.frame()


sub.counts = counts[,sub.coldata$Run]
sum(colnames(sub.counts)==sub.coldata$Run)==ncol(sub.counts)

dds<-DESeqDataSetFromMatrix(sub.counts,
                            colData = sub.coldata,
                            design = ~treatment + genotype + timepoint)
dds <- DESeq(dds)
res = results(dds, contrast = c('treatment', 'T', 'C'))
save(res, file='deseqResults/individual_projects/bewwMulti_deseqResults.Rdata')
