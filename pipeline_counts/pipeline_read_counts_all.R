#pipeline_read_counts_all.R
#organize and plot the counts from each file from throughout data processing pipeline


library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

input='./pipeline_counts/all_pipeline_counts.txt'
sdat0 = read.table(input, header=F, col.names=c('Run', 'value', 'stat'), stringsAsFactors=F) %>% 
  mutate(statRun = paste(stat, Run, sep='_')) %>% 
  filter(!duplicated(statRun))  %>% 
  mutate(Run=as.character(Run),
         value=as.double(value)) %>% 
  select(-statRun) %>% 
  as_tibble()


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


#gather bioProjects
titleSplits = sapply(as.character(mdat$my_title), function(x) strsplit(x, '_'))
bpList = list()
for (ts in titleSplits){
  bp=ts[[length(ts)]][1]
  bpList = append(bpList, bp)
}
bps=unlist(bpList)
bps[bps=='EWW']<-'This study'
mdat$BioProject = bps


#gather refs
mr = mdat %>% 
  left_join(read_csv('metadata/detailed_tables/BioProject_Publication_Table.csv'), by = 'BioProject') %>% 
  mutate(Reference = if_else(Reference=='none',
                             '(none)',
                             Reference))
sum(mr$my_title==mdat$my_title)==nrow(mdat)
mdat$Bioproject = paste(mdat$BioProject, mr$Reference)

#plot full barplot
bp<-mdat %>% 
  ggplot(aes(x=stat, y=m, color=my_title, fill=my_title)) +
    geom_bar(stat='identity', position='dodge') +
    labs(y='read count', x='pipeline step', title='Median Counts') +
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
  labs(y='read count', x='pipeline step', title='mean Counts')

plot_grid(histMeans, rankedMeans, nrow=1, rel_widths=c(.5,1))

#plot scatter abs
lp<-mdat %>% 
  ggplot(aes(x=stat, y=m, color=Bioproject)) +
  geom_point() +
  geom_line(aes(group=my_title)) +
  theme(legend.position='none',
        axis.text.x=element_text(angle=20, vjust=0.75),
        axis.title.x = element_blank()) +
  labs(y='read count', x='pipeline step', subtitle='read counts')
  

#plot scatter prop raw
pp<-mdat %>% 
  mutate(prop = m/max(m) ) %>% 
  ggplot(aes(x=stat, y=prop, color=Bioproject)) +
    geom_point() +
    geom_line(aes(group=my_title)) +
    labs(y='proportion raw reads', x='pipeline step', subtitle='read proportions') +
    theme(legend.position='none',
          axis.text.x=element_text(angle=20, vjust=0.75),
          axis.title.x = element_blank())



#build the plot with multiple calls of plot_grid
xlab = ggdraw() + draw_label("pipeline step")
plts = plot_grid(lp, pp)
sharedLegend <- cowplot::get_legend(pp + theme(legend.position='right'))
final=plot_grid(plts, xlab, sharedLegend, nrow=3, rel_heights=c(1, 0.05, 1))
final




# COMPARE FEATURE COUNTING EFFICIENCY FOR TAG-SEQ -------------------------

view(mdat)
mdat %>% 
  filter(BioProject=='PRJNA559404')
tag_seq_projects = c('H_matz_heatTolLat_PRJNA279192',
                     'K_strader_redDiapause_PRJNA292574',
                     'S_rachel_immune_PRJNA319662',
                     'V_kenkel_co2seep_PRJNA362652',
                     'j1_thisStudy_PRJNA559404',
                     'X_strader_larvalCompetance_PRJNA379147')
mdat %>% 
  mutate(tagseq = my_title %in% tag_seq_projects) %>% 
  filter(stat %in% c('dedupMapped', 'geneCounted')) %>% 
  pivot_wider(stat,
              id_cols = -stat,
              values_from = m) %>% 
  mutate(counting_efficiency = geneCounted/dedupMapped) %>% 
  ggplot(aes(x=tagseq, y = counting_efficiency, color=my_title, shape=tagseq)) +
  geom_jitter(size=5,
              width=0.1) +
  labs(x='Tag-seq',
       y='feature counting efficiency',
       shape = 'Tag-seq',
       color='Bioproject')
  


# 
# 
# 
# lp = sdat %>% 
#   filter(stat %in% c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted")) %>% 
#   mutate(stat=factor(stat, levels=c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted"))) %>% 
#   ggplot(aes(x=stat, y=fixValue, color=my_title)) +
#   geom_point() +
#   geom_line(aes(group=Run)) +
#   theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75)) +
#   labs(y='Read count', x='Pipeline step', subtitle='read counts')
# 
# 
# pp<-sdat %>% 
#   filter(stat %in% c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted")) %>% 
#   mutate(stat=factor(stat, levels=c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted"))) %>% 
#   group_by(Run) %>% 
#   mutate(prop = fixValue/max(fixValue) ) %>% 
#   ggplot(aes(x=stat, y=prop, color=my_title)) +
#     geom_point() +
#     geom_line(aes(group=Run)) +
#     labs(y='Proportion raw reads', x='Pipeline step', subtitle='read proportions') +
#     theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75))
# plot_grid(lp, pp)


# stats and writing out ---------------------------------------------------

#gather bioProjects
titleSplits = sapply(as.character(sdat$my_title), function(x) strsplit(x, '_'))
bpList = list()
for (ts in titleSplits){
  bp=ts[[length(ts)]][1]
  bpList = append(bpList, bp)
}
bps=unlist(bpList)
bps[bps=='EWW']<-'This study'
sdat$Bioproject = bps


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


#WRITE OUT FORMATTED TABLE

#gather bioProjects
titleSplits = sapply(as.character(sdat$my_title), function(x) strsplit(x, '_'))
bpList = list()
for (ts in titleSplits){
  bp=ts[[length(ts)]][1]
  bpList = append(bpList, bp)
}
bps=unlist(bpList)
bps[bps=='EWW']<-'This study'
sdat$Bioproject = bps


long = sdat %>% 
  select(-value) %>% 
  spread(stat, fixValue) %>% 
  select('Run', 'Bioproject', 'LibraryLayout', 'rawCounts', 'trimmedCounts', 'predupMapped', 'dedupMapped', 'geneCounted') %>% 
  write_tsv(path='./pipeline_counts/pipeline_counts_table.tsv')

