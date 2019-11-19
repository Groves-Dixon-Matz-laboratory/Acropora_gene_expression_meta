#core_stress_go_correlationsV2.R
rm(list=ls())
library(tidyverse)
library(DESeq2)
source('figurePlotting/rna_functions.R')


# PLOT CORRELATION OF TWO-TAILED DELTA RANKS FROM GO MWU FOR COMBINATIONS ------------------
ALPHA=0.3


#COMBINATIONS OF PROJECTS
#combinations of bioprojects by stress type
ll=load('./go_mwu/groupedInputs.Rdata')
groupedNames = c('bleached', 'heatNOBLEACH', 'immune', 'salinityNOBLEACH', 'ph')
inputNames = groupedNames


# READ IN AND BIND INTO SINGLE DATAFRAME 
goDivision='BP'
cdatList = lapply(groupedNames, function(x) read_in_two_tailed_go_results(x))
cdat = cdatList %>% 
  purrr::reduce(rbind)


#READ IN TRAIT DATA
tdat = read_csv('metadata/subset_tables/allStress_Coldata.csv') %>% 
  filter(stress=='stressed') %>% 
  select(my_title, projType, stress, bleached) %>% 
  mutate(bleached = if_else(is.na(bleached),
                            'no',
                            bleached),
         stressType = projType,
         stressType = if_else(bleached == 'yes',
                              'bleached',
                              stressType),
         stressType = if_else(stressType=='temp' & bleached=='no',
                              'heat',
                              stressType))
unique(tdat$stressType)
  

# READ IN DATA TO BE X-AXIS 

xInputName = 'corStress'
xdat = read_in_two_tailed_go_results(xInputName) %>% 
  select(delta.rank, p.adj, name)

# MERGE 

pdat = xdat %>% 
  left_join(cdat, by = 'name')
  


# PLOT  THE GROUPED RESULTS

grouped = pdat %>% 
  filter(!is.na(inputName)) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point(alpha=ALPHA) +
  geom_smooth(method='lm') +
  labs(x='All high-stress delta rank',
       y='Specific stress type delta rank',
       color='stress type')
grouped



# DELTA RANK CORRELATION WITH INDIVIDUAL PROJECTS AGAINST CLUSTER A ------------------

#INDIVIDUAL PROJECTS
ll=load('./go_mwu/individualInputs.Rdata')
ll
ll=load('metadata/corStressProjs.Rdata')
ll

# READ IN AND BIND INTO SINGLE DATAFRAME 

icdatList = lapply(individualNames, function(x) read_in_two_tailed_go_results(x))
icdat = icdatList %>% 
  purrr::reduce(rbind)



# READ IN DATA TO BE X-AXIS 

xInputName = 'clusterAstress'
xdat = read_in_two_tailed_go_results(xInputName) %>% 
  select(delta.rank, p.adj, term, name)


# MERGE 

ipdat = xdat %>% 
  left_join(icdat, by = 'name')

# PLOT
legendPos = 'right'
highstress = ipdat %>% 
  filter(!is.na(inputName),
         inputName %in% corStressProjs) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point(alpha=ALPHA) +
  geom_smooth(method='lm', se=FALSE) +
  theme(legend.position=legendPos)

lowstress = ipdat %>% 
  filter(!is.na(inputName),
         inputName %in% lowStressProjs) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point(alpha=ALPHA) +
  geom_smooth(method='lm', se=FALSE) +
  theme(legend.position=legendPos)


# DELTA RANK CORRELATION WITH INDIVIDUAL PROJECTS AGAINST CLUSTER B -------
# 
# READ IN DATA TO BE X-AXIS

B.xInputName = 'clusterBstress'
B.xdat = read_in_two_tailed_go_results(B.xInputName) %>%
  select(delta.rank, p.adj, term, name)



# MERGE

B.ipdat = B.xdat %>%
  left_join(icdat, by = 'name')

# PLOT
legendPos = 'right'
highstress = B.ipdat %>%
  filter(!is.na(inputName),
         inputName %in% corStressProjs) %>%
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point(alpha=ALPHA) +
  geom_smooth(method='lm', se=FALSE) +
  theme(legend.position=legendPos)

lowstress = B.ipdat %>%
  filter(!is.na(inputName),
         inputName %in% lowStressProjs) %>%
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point(alpha=ALPHA) +
  geom_smooth(method='lm', se=FALSE) +
  theme(legend.position=legendPos)




# CREATE CUSTOM CORAL STRESS GROUPS ---------------------------------------

#long the GO annotations in long format
ll=load('metadata/Amil.v2.eggnogWebsite.gene2go.long.Rdata')
ll
lgo

#get annotations
adat = read_tsv('metadata/Amillepora_euk.emapper.annotations.tsv') %>% 
  select(c('query_name',
           'GO_terms',
           'COG cat',
           'eggNOG.annot')) %>% 
  dplyr::rename('gene'= query_name,
         'cogCat'=`COG cat`) 


# PLOT OVERLAY PLOT FOR CLUSTER A AGAINST ALL CLUSTER A ------------------------------

#LOAD MODIFIED NAMES FOR LABELING
mnames = read_tsv('metadata/detailed_tables/stressNamesModified.txt') %>%
  mutate('inputName'=my_title) %>%
  filter(!grepl("^All", treat2)) %>%
  mutate(treat2=factor(treat2, levels=c('bleached', 'heat', 'hyposalinity', 'immune', 'pH')))

#LOAD CUSTOM GO GROUPINGS FROM creat_custom_go_groups.R
ll=load('figurePlotting/selectGoGroups.Rdata')
ll


YLIMS = c(-2600, 2000)
XLIMS = c(-1500, 1400)
STICKADD=200


corStressProjsSepBEWWs = corStressProjs[corStressProjs != 'j1_thisStudy_PRJNA559404'] %>% 
  append(c('bewwCold', 'bewwHeat', 'bewwHyposalinity', 'bewwMulti'))

AonA.dat = ipdat %>% 
  filter(inputName %in% corStressProjs)


upcoords = build_annotation_coords(upSelect,
                                   AonA.dat,
                                   sigOnly=FALSE,
                                   xlabel_left_bound=10,
                                   xlabel_right_bound=1100,
                                   ylabel_bottom_bound=YLIMS[1],
                                   ylabel_top_bound = -600,
                                   stickAdd=STICKADD,
                                   pointTop=FALSE,
                                   center=FALSE)
downcoords = build_annotation_coords(downSelect,
                                     AonA.dat,
                                     sigOnly=FALSE,
                                     xlabel_left_bound=-1500,
                                     xlabel_right_bound=-700,
                                     ylabel_bottom_bound=700,
                                     ylabel_top_bound = 1800,
                                     stickAdd=STICKADD,
                                     pointTop=TRUE,
                                     center=TRUE)

xcoords = rbind(upcoords[[1]], downcoords[[1]])
sub = rbind(upcoords[[2]], downcoords[[2]])
more.sub = rbind(upcoords[[3]], downcoords[[3]])

#check inputs
sub %>% 
  group_by(treat2) %>% 
  summarize(length(unique(inputName)))



#re-build original
textSize=6
overlayPointSize = 3
ALPHA=0.1
legendPos='bottom'


AonA.0 = AonA.dat %>% 
  filter(!is.na(inputName)) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 0, lty=2, color='grey') +
  geom_point(alpha=ALPHA,
             color='grey') +
  # geom_smooth(method='lm', se=FALSE) +
  lims(y=YLIMS,
       x=XLIMS)

AonA.1 = add_custom_go_overlay(AonA.0, sub, more.sub, hjust=0) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  scale_fill_manual(values=c('white', 'firebrick', 'dodgerblue', 'olivedrab1', 'goldenrod')) +
  labs(x='',
       y='',
       subtitle='cluster A',
       fill='',
       shape='') +
  theme(legend.position=legendPos,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 
# highPlot

high.legend <- cowplot::get_legend(AonA.1) %>% 
  ggdraw()

AonA.1


# PLOT OVERLAY PLOT FOR TYPE B AGAINST B ----------------------------------

#filter for type B projects using overall type B as X axis
BonB.dat = B.ipdat %>% 
  filter(inputName %in% lowStressProjs)


#set coordinates for processes exepected to be upgregulated by stress
upcoords = build_annotation_coords(upSelect,
                                   BonB.dat,
                                   sigOnly=FALSE,
                                   xlabel_left_bound=-1300,
                                   xlabel_right_bound=800,
                                   ylabel_bottom_bound=YLIMS[1],
                                   ylabel_top_bound = -1500,
                                   stickAdd=STICKADD,
                                   pointTop=FALSE,
                                   center=FALSE)

#set coordinates for processes exepected to be downregulated by stress
downcoords = build_annotation_coords(downSelect,
                                     BonB.dat,
                                     sigOnly=FALSE,
                                     xlabel_left_bound=-300,
                                     xlabel_right_bound=1200,
                                     ylabel_bottom_bound=1600,
                                     ylabel_top_bound = 2000,
                                     stickAdd=STICKADD,
                                     pointTop=TRUE,
                                     center=TRUE)






xcoords = rbind(upcoords[[1]], downcoords[[1]])
sub = rbind(upcoords[[2]], downcoords[[2]])
more.sub = rbind(upcoords[[3]], downcoords[[3]])

#plot
BonB.0 = BonB.dat %>% 
  filter(!is.na(inputName),
         inputName %in% lowStressProjs) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 0, lty=2, color='grey') +
  geom_point(alpha=ALPHA,
             color='grey') +
  # geom_smooth(method='lm', se=FALSE) +
  theme(legend.position=legendPos,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(x='',
       y='',
       subtitle='cluster B',
       fill='',
       shape='') +
  lims(y=YLIMS)


# lowPlot = add_custom_go_overlay(lowstress, sub, more.sub)
BonB.1 = add_custom_go_overlay(BonB.0, sub, more.sub) +
  scale_shape_manual(values=c(21, 22, 24, 25)) +
  scale_fill_manual(values=c('white', 'firebrick', 'olivedrab1', 'goldenrod'))
BonB.1



# PLOT OVERLAY FOR TYPE B AGAINST A ---------------------------------------

#filter for type B projects using type A as X axis
BonA.dat = ipdat %>% 
  filter(inputName %in% lowStressProjs)


#set coordinates for up and down regulated process labeles
upcoords = build_annotation_coords(upSelect,
                                   BonA.dat,
                                   sigOnly=FALSE,
                                   xlabel_left_bound=100,
                                   xlabel_right_bound=1600,
                                   ylabel_bottom_bound=YLIMS[1],
                                   ylabel_top_bound = -1500,
                                   stickAdd=STICKADD,
                                   pointTop=FALSE,
                                   center=FALSE)
downcoords = build_annotation_coords(downSelect,
                                     BonA.dat,
                                     sigOnly=FALSE,
                                     xlabel_left_bound=-1500,
                                     xlabel_right_bound=-700,
                                     ylabel_bottom_bound=1600,
                                     ylabel_top_bound = 2000,
                                     stickAdd=STICKADD,
                                     pointTop=TRUE,
                                     center=TRUE)

xcoords = rbind(upcoords[[1]], downcoords[[1]])
sub = rbind(upcoords[[2]], downcoords[[2]])
more.sub = rbind(upcoords[[3]], downcoords[[3]])

#plot
BonA.0 = BonA.dat %>% 
  filter(!is.na(inputName),
         inputName %in% lowStressProjs) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 0, lty=2, color='grey') +
  geom_point(alpha=ALPHA,
             color='grey') +
  # geom_smooth(method='lm', se=FALSE) +
  theme(legend.position=legendPos,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(x='',
       y='',
       subtitle='cluster B',
       fill='',
       shape='') +
  lims(y=YLIMS)


# lowPlot = add_custom_go_overlay(lowstress, sub, more.sub)
BonA.1 = add_custom_go_overlay(BonA.0, sub, more.sub) +
  scale_shape_manual(values=c(21, 22, 24, 25)) +
  scale_fill_manual(values=c('white', 'firebrick', 'olivedrab1', 'goldenrod'))

# lowPlotColored



#### FINAL PLOT
ylab = ggdraw() + draw_label('bioproject delta rank', angle=90)
xlab = ggdraw() + draw_label('all cluster A delta rank')
plts = plot_grid(highPlot+theme(legend.position='none'),
                 lowPlotColored+theme(legend.position='none'))
top=plot_grid(ylab, plts, nrow=1, rel_widths=c(0.03, 1))
topx = plot_grid(top, xlab, nrow=2,rel_heights = c(1, 0.05))
plot_grid(topx, high.legend, nrow=2, rel_heights=c(1, 0.05))



#or plot stacked
plts = plot_grid(highPlot+theme(legend.position='none'),
                 lowPlotColored+theme(legend.position='none'),
                 nrow=2)
top=plot_grid(ylab, plts, nrow=1, rel_widths=c(0.03, 1))
topx = plot_grid(top, xlab, nrow=2,rel_heights = c(1, 0.05))
plot_grid(topx, high.legend, nrow=2, rel_heights=c(1, 0.05))


#New Final plot
pltList = list(AonA.1, BonB.1, BonA.1)
plot_grid(plotlist=pltList, nrow=3)


pltList2 = lapply(pltList, function(x) return(x+lims(x=XLIMS)))
plot_grid(plotlist=pltList, nrow=3)

# BUILD WITH BIOPROJECTS LABELED ------------------------------------------
#plot for all
lowstress = ipdat %>% 
  filter(!is.na(inputName),
         inputName %in% lowStressProjs) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 0, lty=2, color='grey') +
  geom_point(alpha=ALPHA,
             color='grey') +
  # geom_smooth(method='lm', se=FALSE) +
  theme(legend.position=legendPos,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(x='',
       y='',
       subtitle='low-stress',
       fill='',
       shape='') +
  lims(y=YLIMS)
add_custom_go_overlay_bioproject(lowstress, sub, more.sub) + theme(legend.position='right')


#subset for a project
SELECT_PROJECT='H_matz_heatTolLat_PRJNA279192'
low.ipdat = ipdat %>% 
  filter(inputName == SELECT_PROJECT)


upcoords = build_annotation_coords(upSelect,
                                   low.ipdat,
                                   sigOnly=FALSE,
                                   xlabel_left_bound=10,
                                   xlabel_right_bound=1500,
                                   ylabel_bottom_bound=YLIMS[1],
                                   ylabel_top_bound = -1300)
downcoords = build_annotation_coords(downSelect,
                                     low.ipdat,
                                     sigOnly=FALSE,
                                     xlabel_left_bound=-1800,
                                     xlabel_right_bound=-800,
                                     ylabel_bottom_bound=1500,
                                     ylabel_top_bound = 2000)

xcoords = rbind(upcoords[[1]], downcoords[[1]])
sub = rbind(upcoords[[2]], downcoords[[2]])
more.sub = rbind(upcoords[[3]], downcoords[[3]])

lowstress = ipdat %>% 
  filter(!is.na(inputName),
         inputName == SELECT_PROJECT) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 0, lty=2, color='grey') +
  geom_point(alpha=ALPHA,
             color='grey') +
  # geom_smooth(method='lm', se=FALSE) +
  theme(legend.position=legendPos,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(x='',
       y='',
       subtitle='low-stress',
       fill='',
       shape='') +
  lims(y=YLIMS)
add_custom_go_overlay_bioproject(lowstress, sub, more.sub) + theme(legend.position='right')





# ORGANIZE EXPRESSION DATA  -----------------------------------------------------

#upload the individual deseq results
ll=load('deseqResults/all_by_bioproject.Rdata')
ll
head(idat)

#upload the 'all', 'high', and 'low' results
alldat = read_deseq_res('deseqResults/stress_deseqResults.Rdata', 'all')
highdat = read_deseq_res('deseqResults/corStress_deseqResults.Rdata', 'high')
lowdat = read_deseq_res('deseqResults/lowStress_deseqResults.Rdata', 'low')

#merge these together
cdfList = list(idat, alldat, highdat, lowdat)
cdat = cdfList %>% 
  purrr::reduce(full_join, by = c('gene')) %>% 
  as_tibble()

#gather into long format
lcs = cdat %>% 
  select(gene, colnames(cdat)[grep('_lfc', colnames(cdat))])
ldat = lcs %>% 
  gather(key='project', value='lfc', 2:ncol(lcs))
projects = unique(ldat$project)



# LOOK AT EXPRESSION RESULTS BY FUNCITON ----------------------------------

sgenes = ribosome.genes
sgenes = red.genes
sgenes = white.genes
sgenes = green.genes
sgenes = antiox.genes
sgenes = ros.genes
sgenes = hsp.genes     #up only in high stress
sgenes = folding.genes #doesn't show anything
sgenes = prolif.genes


sgenes = get_go_genes(lgo, 'GO:0044183') #protein folding chaperone


dat = ldat
sub = dat %>% 
  filter(gene %in% sgenes,
         !is.na(lfc)) 



pltList = list()
for (p in projects){
  plt = dat %>% 
    filter(project == p,
           gene %in% sgenes,
           !is.na(lfc))  %>% 
    ggplot(aes(x=project, y=lfc)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty=2)
  pltList[[p]]=plt
}

plot_grid(plotlist = pltList,
          nrow=5)















# attmepted delta rank correlations for one-tailed, didn't work well --------

inputFiles = c('corStress_down_For_MWU.csv',
               'corStress_up_For_MWU.csv',
               'lowStress_down_For_MWU.csv',
               'lowStress_up_For_MWU.csv',
               'redMembership_ForMWU.csv',
               'stressDAPC_ForMWU.csv')


#combinations of bioprojects by stress type
ll=load('./go_mwu/groupedInputs.Rdata')
groupedNames = c('bleached', 'heatNOBLEACH', 'immune', 'salinityNOBLEACH', 'ph')
inputFiles = append(paste(groupedNames, 'down_For_MWU.csv', sep='_'),
                    paste(groupedNames, 'up_For_MWU.csv', sep='_'))

#individual bioprojects
ll=load('./go_mwu/individualInputs.Rdata')
inputFiles = append(paste(individualNames, 'down_For_MWU.csv', sep='_'),
                    paste(individualNames, 'up_For_MWU.csv', sep='_'))


ll=load('metadata/corStressProjs.Rdata')
ll
sum(corStressProjs %in% individualNames)==length(corStressProjs)
sum(lowStressProjs %in% individualNames)==length(lowStressProjs)


goDivision = 'BP'
combind_up_and_down_go_res = function(inputName){
  upName = paste(paste('./go_mwu/MWU', goDivision, sep = "_"), paste(inputName,'up_For_MWU.csv',sep='_'), sep = "_")
  downName = paste(paste('./go_mwu/MWU', goDivision, sep = "_"), paste(inputName, 'down_For_MWU.csv',sep='_'), sep = "_")
  upRes = read.table(upName, header = TRUE, stringsAsFactors=FALSE) %>% 
    mutate(enrichment = -log(pval, 10))
  downRes = read.table(downName, header = TRUE, stringsAsFactors=FALSE) %>% 
    mutate(enrichment = log(pval, 10))
  res = rbind(upRes, downRes)
  res$inputName = inputName
  return(res)
}

plotenrich = function(cdat){
  cdat %>% 
    ggplot(aes(x=enrichment, y=delta.rank)) +
    geom_point() 
}





cdatList = lapply(individualNames, function(x) combind_up_and_down_go_res(x))
names(cdatList) = individualNames
cdat = cdatList %>% 
  purrr::reduce(rbind)


xInputName = 'corStress'
xdat = combind_up_and_down_go_res(xInputName) %>% 
  select(enrichment, name)

pdat = xdat %>% 
  left_join(cdat, by = 'name')



head(pdat)

pdat %>% 
  mutate(sGroup = if_else(inputName %in% corStressProjs,
                          'high-stress',
                          'low-stress')) %>% 
  ggplot(aes(x=enrichment.x, y=enrichment.y, color=sGroup)) +
  geom_point() +
  geom_smooth(method='lm') 




# PLOT HIGH-STRESS BY ITSELF ----------------------------------------------

mnames = read_tsv('metadata/detailed_tables/stressNamesModified.txt') %>% 
  dplyr::rename('inputName'=id) %>% 
  filter(!grepl("^All", treat2)) %>% 
  mutate(treat2=factor(treat2, levels=c('bleached', 'heat', 'hyposalinity', 'immune', 'pH')))

high.ipdat = ipdat %>% 
  filter(inputName %in% corStressProjs)


upcoords = build_annotation_coords(upSelect,
                                   high.ipdat,
                                   sigOnly=TRUE,
                                   xlabel_left_bound=10,
                                   xlabel_right_bound=1500,
                                   ylabel_bottom_bound=-1800,
                                   ylabel_top_bound = -500)
downcoords = build_annotation_coords(downSelect,
                                     high.ipdat,
                                     sigOnly=TRUE,
                                     xlabel_left_bound=-1800,
                                     xlabel_right_bound=-800,
                                     ylabel_bottom_bound=500,
                                     ylabel_top_bound = 1800) 

xcoords = rbind(upcoords[[1]], downcoords[[1]])
sub = rbind(upcoords[[2]], downcoords[[2]]) %>% 
  filter(summary!='ribosomes')
more.sub = rbind(upcoords[[3]], downcoords[[3]]) %>% 
  filter(summary!='ribosomes')

highPlot = add_custom_go_overlay(highstress, sub, more.sub) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  scale_fill_manual(values=c('white', 'firebrick', 'dodgerblue', 'olivedrab1', 'goldenrod')) +
  labs(x='All together',
       y='Individual studies',
       fill='',
       shape='',
       subtitle='Functional enrichment') +
  theme(legend.position='bottom') 

high.legend <- cowplot::get_legend(highPlot) %>% 
  ggdraw()

