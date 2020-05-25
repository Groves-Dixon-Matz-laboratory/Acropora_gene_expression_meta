#plot_go_correlations.R
#plots scatterplots of delta ranks for groups of projects and individual projects
#purpose is to show how consistent GO enrichment is accross projects
#for the main plot, the X axis will show enrichment for a large group of projects' samples together
#the Y axis will show data from individual projects. Tight correlation shows how consistent the
#enrichment patterns are.

#SETUP
rm(list=ls())
library(tidyverse)
library(DESeq2)
source('figurePlotting/rna_functions.R')


# PLOT CORRELATION OF TWO-TAILED DELTA RANKS FROM GO MWU FOR COMBINATIONS ------------------

#set alpha level for all plots
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
  


# PLOT

grouped = pdat %>% 
  filter(!is.na(inputName)) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point(alpha=ALPHA) +
  geom_smooth(method='lm') +
  labs(x='All type A delta rank',
       y='Specific stress type delta rank',
       color='stress type')
  



# REPEAT DELTA RANK CORRELATION WITH INDIVIDUAL PROJECTS ------------------

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

xInputName = 'corStress'
xdat = read_in_two_tailed_go_results(xInputName) %>% 
  select(delta.rank, p.adj, term, name)



#merge with the long-formatted GO enrichment results from each individual project
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


# RED MODULE --------------------------------------------------------------

ll=load('wgcna/moduleAssignment.Rdata')
ll
module.genes = rownames(geneModuleMembership)
red.genes = module.genes[moduleColors=='red']
length(red.genes)


# WHITE MODULE ------------------------------------------------------------

white.genes = module.genes[moduleColors=='white']
length(white.genes)
  

# GREEN MODULE ------------------------------------------------------------

green.genes = module.genes[moduleColors=='green4']
length(green.genes)



# REPLOT DELTA RANK CORS WITH OVERLAY -------------------------------------
#Here we plot the scatterplots from above, but with key GO terms overlaid on top
#These key GO terms are cherry picked based on previous publications and the 
#enrichment results from the 'all type A' enrichment (x axis read in as corStress above)

#LOAD MODIFIED NAMES FOR LABELING
#altered names from the arbitrary labels used before 
mnames = read_tsv('metadata/detailed_tables/stressNamesModified.txt') %>% 
  mutate('inputName'=my_title) %>% 
  filter(!grepl("^All", treat2)) %>% 
  mutate(treat2=factor(treat2, levels=c('bleached', 'heat', 'hyposalinity', 'immune', 'pH')))


#LOAD CUSTOM GO GROUPINGS FROM creat_custom_go_groups.R
ll=load('figurePlotting/selectGoGroups.Rdata')
ll

#choose coordinates for teh up and down sets of GO terms
upcoords = build_annotation_coords(upSelect,
                                   pdat,
                                   sigOnly=FALSE,
                                   xlabel_left_bound=300,
                                   xlabel_right_bound=1700,
                                   ylabel_bottom_bound=-1800,
                                   ylabel_top_bound = -500)
downcoords = build_annotation_coords(downSelect,
                                     pdat,
                                     sigOnly=FALSE,
                                     xlabel_left_bound=-1800,
                                     xlabel_right_bound=-800,
                                     ylabel_bottom_bound=500,
                                     ylabel_top_bound = 1800)

xcoords = rbind(upcoords[[1]], downcoords[[1]])
sub = rbind(upcoords[[2]], downcoords[[2]])
more.sub = rbind(upcoords[[3]], downcoords[[3]])



#re-build original
ALPHA=0.1
grouped = pdat %>% 
  filter(!is.na(inputName)) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 0, lty=2, color='grey') +
  geom_point(alpha=ALPHA) +
  # geom_smooth(method='lm', se=FALSE) +
  labs(x='All high-stress delta rank',
       y='Specific stress type delta rank',
       color='stress type') +
  theme(legend.position='right')


#add new points and annotations for select processes
textSize=6
overlayPointSize = 5
add_custom_go_overlay(grouped, sub, more.sub)
STICKADD=200


# REPLOT DELTA RANK CORS WITH OVERLAY FOR INDIVIDUAL PROJECTS ------------------------------

YLIMS = c(-2600, 2000)


corStressProjsSepBEWWs = corStressProjs[corStressProjs != 'j1_thisStudy_PRJNA559404'] %>% 
  append(c('bewwCold', 'bewwHeat', 'bewwHyposalinity', 'bewwMulti'))

high.ipdat = ipdat %>% 
  filter(inputName %in% corStressProjs)


upcoords = build_annotation_coords(upSelect,
                                   high.ipdat,
                                   sigOnly=FALSE,
                                   xlabel_left_bound=400,
                                   xlabel_right_bound=1550,
                                   ylabel_bottom_bound=YLIMS[1],
                                   ylabel_top_bound = -600,
                                   stickAdd=STICKADD,
                                   pointTop=FALSE,
                                   center=FALSE)
downcoords = build_annotation_coords(downSelect,
                                     high.ipdat,
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


highstress = high.ipdat %>% 
  filter(!is.na(inputName)) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 0, lty=2, color='grey') +
  geom_point(alpha=ALPHA,
             color='grey') +
  # geom_smooth(method='lm', se=FALSE) +
  lims(y=YLIMS)

highPlot = add_custom_go_overlay(highstress, sub, more.sub, hjust=0) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  scale_fill_manual(values=c('white', 'firebrick', 'dodgerblue', 'olivedrab1', 'goldenrod')) +
  labs(x='',
       y='',
       subtitle='type A',
       fill='',
       shape='') +
  theme(legend.position=legendPos,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 
# highPlot

high.legend <- cowplot::get_legend(highPlot) %>% 
  ggdraw()


#REPEAT FOR LOW STRESS

low.ipdat = ipdat %>% 
  filter(inputName %in% lowStressProjs)


upcoords = build_annotation_coords(upSelect,
                                   low.ipdat,
                                   sigOnly=FALSE,
                                   xlabel_left_bound=100,
                                   xlabel_right_bound=1600,
                                   ylabel_bottom_bound=YLIMS[1],
                                   ylabel_top_bound = -1500,
                                   stickAdd=STICKADD,
                                   pointTop=FALSE,
                                   center=FALSE)
downcoords = build_annotation_coords(downSelect,
                                     low.ipdat,
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
       subtitle='type B',
       fill='',
       shape='') +
  lims(y=YLIMS)


lowPlot = add_custom_go_overlay(lowstress, sub, more.sub)
lowPlotColored = add_custom_go_overlay(lowstress, sub, more.sub) +
  scale_shape_manual(values=c(21, 22, 24, 25)) +
  scale_fill_manual(values=c('white', 'firebrick', 'olivedrab1', 'goldenrod'))

# lowPlotColored



#### FINAL PLOT
ylab = ggdraw() + draw_label('independent bioproject delta rank', angle=90)
xlab = ggdraw() + draw_label('full type A dataset delta rank')
plts = plot_grid(highPlot+theme(legend.position='none'),
                 lowPlotColored+theme(legend.position='none'))
top=plot_grid(ylab, plts, nrow=1, rel_widths=c(0.03, 1))
topx = plot_grid(top, xlab, nrow=2,rel_heights = c(1, 0.05))
# plot_grid(topx, high.legend, nrow=2, rel_heights=c(1, 0.05))



#or plot stacked
plts = plot_grid(highPlot+theme(legend.position='none'),
                 lowPlotColored+theme(legend.position='none'),
                 nrow=2)
top=plot_grid(ylab, plts, nrow=1, rel_widths=c(0.03, 1))
topx = plot_grid(top, xlab, nrow=2,rel_heights = c(1, 0.05))
plot_grid(topx, high.legend, nrow=2, rel_heights=c(1, 0.05))




#=================================================================
# REPLOT WITH TYPE-B STRESS AS X AXIS -------------------


# READ IN DATA TO BE X-AXIS 
#now read enrichment results for all type-B samples at once (lowStress input for GO_MWU)
xInputName = 'lowStress'
xdat = read_in_two_tailed_go_results(xInputName) %>% 
  select(delta.rank, p.adj, term, name)



# MERGE 
#merge up with the individual project enrichment results
#note replacing object ipdat
ipdat = xdat %>% 
  left_join(icdat, by = 'name')


# PLOT BASIC SCATTERPLOTS 
#color coded by biopject
legendPos = 'right'
highstress = ipdat %>% 
  filter(!is.na(inputName),
         inputName %in% corStressProjs) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point(alpha=ALPHA) +
  geom_smooth(method='lm', se=FALSE) +
  theme(legend.position=legendPos) #now shows negative associations of typeA projects with all of type-B at once

lowstress = ipdat %>% 
  filter(!is.na(inputName),
         inputName %in% lowStressProjs) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point(alpha=ALPHA) +
  geom_smooth(method='lm', se=FALSE) +
  theme(legend.position=legendPos) #now shows positive associations of individual type B with all type B



# REPLOT DELTA RANK CORS WITH OVERLAY FOR INDIVIDUAL PROJECTS ------------------------------
#THIS TIME WITH TYPE-B AS X AXIS


YLIMS = c(-2600, 2200)


corStressProjsSepBEWWs = corStressProjs[corStressProjs != 'j1_thisStudy_PRJNA559404'] %>% 
  append(c('bewwCold', 'bewwHeat', 'bewwHyposalinity', 'bewwMulti'))

high.ipdat = ipdat %>% 
  filter(inputName %in% corStressProjs)


upcoords = build_annotation_coords(upSelect,
                                   high.ipdat,
                                   sigOnly=FALSE,
                                   xlabel_left_bound=-2500,
                                   xlabel_right_bound=250,
                                   ylabel_bottom_bound=YLIMS[1],
                                   ylabel_top_bound = -1000,
                                   stickAdd=STICKADD,
                                   pointTop=FALSE,
                                   center=FALSE)
downcoords = build_annotation_coords(downSelect,
                                     high.ipdat,
                                     sigOnly=FALSE,
                                     xlabel_left_bound=-250,
                                     xlabel_right_bound=1500,
                                     ylabel_bottom_bound=1800,
                                     ylabel_top_bound = 2200,
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


highstress = high.ipdat %>% 
  filter(!is.na(inputName)) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 0, lty=2, color='grey') +
  geom_point(alpha=ALPHA,
             color='grey') +
  # geom_smooth(method='lm', se=FALSE) +
  lims(y=YLIMS)

highPlot = add_custom_go_overlay(highstress, sub, more.sub, hjust=0) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  scale_fill_manual(values=c('white', 'firebrick', 'dodgerblue', 'olivedrab1', 'goldenrod')) +
  labs(x='',
       y='',
       subtitle='type A',
       fill='',
       shape='') +
  theme(legend.position=legendPos,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 
# highPlot

high.legend <- cowplot::get_legend(highPlot) %>% 
  ggdraw()


#REPEAT FOR LOW STRESS

low.ipdat = ipdat %>% 
  filter(inputName %in% lowStressProjs)


upcoords = build_annotation_coords(upSelect,
                                   low.ipdat,
                                   sigOnly=FALSE,
                                   xlabel_left_bound=-2100,
                                   xlabel_right_bound=400,
                                   ylabel_bottom_bound=500,
                                   ylabel_top_bound = 2200,
                                   stickAdd=STICKADD,
                                   pointTop=TRUE,
                                   center=TRUE)
downcoords = build_annotation_coords(downSelect,
                                     low.ipdat,
                                     sigOnly=FALSE,
                                     xlabel_left_bound=-300,
                                     xlabel_right_bound=1500,
                                     ylabel_bottom_bound=-1990,
                                     ylabel_top_bound = -2000,
                                     stickAdd=STICKADD,
                                     pointTop=FALSE,
                                     center=TRUE)

xcoords = rbind(upcoords[[1]], downcoords[[1]])
sub = rbind(upcoords[[2]], downcoords[[2]])
more.sub = rbind(upcoords[[3]], downcoords[[3]])

#plot
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
       subtitle='type B',
       fill='',
       shape='') +
  lims(y=YLIMS)


lowPlot = add_custom_go_overlay(lowstress, sub, more.sub)
lowPlotColored = add_custom_go_overlay(lowstress, sub, more.sub) +
  scale_shape_manual(values=c(21, 22, 24, 25)) +
  scale_fill_manual(values=c('white', 'firebrick', 'olivedrab1', 'goldenrod'))

# lowPlotColored



#### FINAL PLOT
ylab = ggdraw() + draw_label('independent bioproject delta rank', angle=90)
xlab = ggdraw() + draw_label('full type B dataset delta rank')
plts = plot_grid(highPlot+theme(legend.position='none'),
                 lowPlotColored+theme(legend.position='none'))
top=plot_grid(ylab, plts, nrow=1, rel_widths=c(0.03, 1))
topx = plot_grid(top, xlab, nrow=2,rel_heights = c(1, 0.05))
# plot_grid(topx, high.legend, nrow=2, rel_heights=c(1, 0.05))


#or plot stacked
plts = plot_grid(highPlot+theme(legend.position='none'),
                 lowPlotColored+theme(legend.position='none'),
                 nrow=2)
top=plot_grid(ylab, plts, nrow=1, rel_widths=c(0.03, 1))
topx = plot_grid(top, xlab, nrow=2,rel_heights = c(1, 0.05))
plot_grid(topx, high.legend, nrow=2, rel_heights=c(1, 0.05))





