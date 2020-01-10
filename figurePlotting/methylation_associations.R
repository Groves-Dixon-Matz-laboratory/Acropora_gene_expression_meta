#methylation_associations.R
library(tidyverse)
source('figurePlotting/rna_functions.R')


# LOAD DNA METHYLATION DATA -----------------------------------------------

ll=load('figurePlotting/gbmLvl.Rdata')
ll
gbm.dat = gbm.dat %>% 
  mutate(logFracMeth = log(fracMeth, 2)) %>% 
  dplyr::select(name, logFracMeth, mbd.score, mrB)


# LOAD LEVEL --------------------------------------------------------------




# LOAD DESEQ RESULTS ------------------------------------------------------

ll=load('deseqResults/adultVlarva_deseqResults.Rdata')


dat = res %>% 
  data.frame() %>% 
  mutate(name=rownames(res)) %>% 
  dplyr::select(name, log2FoldChange, pvalue) %>% 
  left_join(gbm.dat, by = 'name') %>% 
  pivot_longer(logFracMeth:mrB, names_to = 'gbmMeasure', values_to = 'methLvl') %>% 
  as_tibble()

meth_scatter = function(df){

  df %>% 
    ggplot(aes(x=methLvl, y=log2FoldChange)) + 
    geom_point(alpha=0.2) +
    geom_smooth(method='lm')
}

dat %>% 
  pivot_longer(logFracMeth:mrB, names_to = 'gbmMeasure', values_to = 'methLvl') %>% 
  meth_scatter() + facet_grid(~gbmMeasure, scales='free')
