#!/usr/bin/env Rscript
#load the data
library(DESeq2)
library(tidyverse)
ll = load('deseqBaselineInput.Rdata')
ll
coldata
head(counts)
sum(colnames(counts) == rownames(coldata))==ncol(counts)


#----- format coldata
unique(coldata$treat)

#look at totals
coldata %>% 
  group_by(treat,my_title,my_stage) %>% 
  summarize(N=n())


#filter to just adults and larvae
keep = coldata %>%
  mutate(Run=rownames(coldata)) %>%
  filter(my_stage %in% c('adult', 'larval')) %>%
  pull(Run)

coldata = coldata[keep,]
counts = counts[,keep]


### get variance stabilized counts for timepoint 2
dds<-DESeqDataSetFromMatrix(counts,
	colData = coldata, 
	design = formula(~my_stage))
dds <- DESeq(dds)
res = results(dds, contrast = c('my_stage', 'adult', 'larval'))
save(res, file='all_deseqResults.Rdata')



