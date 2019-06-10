setwd("/Users/grovesdixon/lab_files/projects/gene_expression_meta/by_project/A_moya_acid_PRJNA149513")
setwd("/Users/grovesdixon/lab_files/projects/gene_expression_meta/pipeline_stats")


library(tidyverse)
library(cowplot)
sdat0 = read_tsv('allStats.txt', col_names=c('run', 'value', 'stat'))
unique(sdat0$stat)



cdat = read_tsv('heatTolLat_rnaCounts.tsv')
colnames(cdat) = sub('_dupsRemoved.counts.txt', '', colnames(cdat))
sdat <- cdat %>% 
  filter(!grepl('__', geneID)) %>% 
  select(2:ncol(cdat)) %>% 
  summarise_all(funs(sum)) %>% 
  gather(key='run') %>% 
  mutate(stat='counted') %>% 
  rbind(sdat0) %>% 
  filter(!stat %in% c('dedupEff', 'dedupPropPair', 'dupRemProp', 'predupEff', 'predupPropPaired')) %>% 
  mutate(stat = factor(stat, levels=c('rawCounts', 'trimmedCounts', 'predupMapped', 'dedupMapped', 'counted')))



#plot barplot
bp<-sdat %>% 
  mutate(value=as.numeric(value)) %>% 
  ggplot(aes(x=stat, y=value, color=run, fill=run)) +
    geom_bar(stat='identity', position='dodge') +
    labs(y='Read count', x='Pipeline step', title='Moya experiment') +
    theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.5))


#plot scatter abs
lp<-sdat %>% 
  mutate(value=as.numeric(value)) %>% 
  ggplot(aes(x=stat, y=value, color=run)) +
  geom_point() +
  geom_line(aes(group=run)) +
  theme(legend.position='none') +
  labs(y='Read count', x='Pipeline step', subtitle='read counts')
  

#plot_grid(bp, lp)

  
#plot scatter prop raw
pp<-sdat %>% 
  mutate(value=as.numeric(value)) %>% 
  group_by(run) %>% 
  mutate(prop = value/max(value) ) %>% 
  ggplot(aes(x=stat, y=prop, color=run)) +
    geom_point() +
    geom_line(aes(group=run)) +
    labs(y='Proportion raw reads', x='Pipeline step', subtitle='read proportions') +
    theme(legend.position='none')

plot_grid(lp, pp)

cdat %>% 
  filter(!grepl('__', geneID))
