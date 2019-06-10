#pipeline_read_counts_all.R
#organize and plot the counts from each file from throughout data processing pipeline


library(tidyverse)
library(cowplot)

input='./pipeline_counts/all_pipeline_counts.txt'
sdat0 = read.table(input, header=F, col.names=c('Run', 'value', 'stat'), stringsAsFactors=F) %>% 
  mutate(statRun = paste(stat, Run, sep='_')) %>% 
  filter(!duplicated(statRun))  %>% 
  mutate(Run=as.character(Run),
         value=as.double(value)) %>% 
  select(-statRun) %>% 
  as_tibble()


j = sdat %>% 
  filter(my_title=="j1_BLEACH_EWW") %>% 
  pull(Run)

sdat0 %>% 
  filter(Run %in% j) %>% 
  view()


#UPLOAD SAMPLE DATA AND REMOVE UNWANTED INFORMATION

#raw runInfo from NCBI
cdat0 = read_csv('./metadata/all_acropora_sra_runInfoTable.csv') %>% 
  select(Run, LibraryLayout)

#my modified table with arbitrary titles to organize projects and lots of the junk trimmed
cdat1 = read_csv('./metadata/ALL_Coldata.csv') %>% 
  select(Run, my_title)

#merge and handle the two additional projects that weren't in all_acropora_sra_runInfoTable.csv
cdat = cdat1 %>% 
  left_join(cdat0, by = 'Run') %>% 
  mutate(LibraryLayout = if_else(my_title %in% c('j1_BLEACH_EWW', "k1_Palumbi_lab_heat_resilience_PRJNA274410"),
                                 'SINGLE',
                                 LibraryLayout))


#merge and divide paired end read counts from samtools flagstat by 2, since it counts all reads
#remove uneeded columns
sdat = cdat %>% 
  left_join(sdat0, by = 'Run') %>% 
  mutate(fixValue=if_else( !stat %in% c('rawCounts', 'trimmedCounts', 'geneCounted') & LibraryLayout=='PAIRED',
                       value/2,
                       value)
  )



#get project means
mdat = sdat %>% 
  group_by(my_title, stat) %>% 
  summarize(m = median(fixValue)) %>% 
  filter(stat %in% c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted")) %>% 
  mutate(stat=factor(stat, levels=c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted")))




#order the stats in pipline order
order = mdat %>% 
  filter(stat=='rawCounts') %>% 
  arrange(by=m) %>% 
  pull(my_title) %>% 
  rev() 
mdat$my_title = factor(mdat$my_title, levels=order)





#plot full barplot
bp<-mdat %>% 
  ggplot(aes(x=stat, y=m, color=my_title, fill=my_title)) +
    geom_bar(stat='identity', position='dodge') +
    labs(y='Read count', x='Pipeline step', title='Median Counts') +
    theme(axis.text.x=element_text(angle=20, vjust=0.75))
bp


#plot gene counted only
gc = mdat %>% 
  filter(stat=='geneCounted') %>% 
  arrange(by=m)
gc$my_title = factor(gc$my_title, levels = rev(pull(gc, my_title)))

histMeans = gc %>% 
  ggplot(aes(x=m)) +
  geom_histogram(bins=12) +
  labs(x='Mean gene counted reads accross projects')

rankedMeans = gc %>% 
  rename('Project'=my_title) %>% 
  ggplot(aes(x=stat, y=m, color=Project, fill=Project)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y='Read count', x='Pipeline step', title='Mean Counts')

plot_grid(histMeans, rankedMeans, nrow=1, rel_widths=c(.5,1))

#plot scatter abs
lp<-mdat %>% 
  ggplot(aes(x=stat, y=m, color=my_title)) +
  geom_point() +
  geom_line(aes(group=my_title)) +
  theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75)) +
  labs(y='Read count', x='Pipeline step', subtitle='read counts')
  

#plot scatter prop raw
pp<-mdat %>% 
  mutate(prop = m/max(m) ) %>% 
  ggplot(aes(x=stat, y=prop, color=my_title)) +
    geom_point() +
    geom_line(aes(group=my_title)) +
    labs(y='Proportion raw reads', x='Pipeline step', subtitle='read proportions') +
    theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75))

plot_grid(lp, pp)


lp = sdat %>% 
  filter(stat %in% c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted")) %>% 
  mutate(stat=factor(stat, levels=c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted"))) %>% 
  ggplot(aes(x=stat, y=fixValue, color=my_title)) +
  geom_point() +
  geom_line(aes(group=Run)) +
  theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75)) +
  labs(y='Read count', x='Pipeline step', subtitle='read counts')


pp<-sdat %>% 
  filter(stat %in% c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted")) %>% 
  mutate(stat=factor(stat, levels=c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted"))) %>% 
  group_by(Run) %>% 
  mutate(prop = fixValue/max(fixValue) ) %>% 
  ggplot(aes(x=stat, y=prop, color=my_title)) +
    geom_point() +
    geom_line(aes(group=Run)) +
    labs(y='Proportion raw reads', x='Pipeline step', subtitle='read proportions') +
    theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75))
plot_grid(lp, pp)


# stats and writing out ---------------------------------------------------

#get stats
library(plotrix)
statTable = sdat %>% 
  filter(stat %in% c('rawCounts', 'trimmedCounts', 'predupMapped', 'dedupMapped', 'geneCounted')) %>% 
  group_by(stat) %>% 
  summarize(mean=mean(fixValue),
            sd = sd(fixValue),
            stdErr = std.error(fixValue)) 
statTable
statTable %>% 
  write_tsv('./pipeline_counts/summary_table.tsv')


#Write out formatted table
head(sdat)
long = sdat %>% 
  select(-value) %>% 
  spread(stat, fixValue) %>% 
  select('Run', 'my_title', 'LibraryLayout', 'rawCounts', 'trimmedCounts', 'predupMapped', 'dedupMapped', 'geneCounted') %>% 
  write_tsv(path='./pipeline_counts/pipeline_counts_table.tsv')

