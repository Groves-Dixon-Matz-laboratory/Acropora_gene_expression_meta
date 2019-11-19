#frontloading_comparison.R

rm(list=ls())
library(tidyverse)
library(cowplot)
library(ggplot2)
library(DESeq2)
library(adegenet)
source('./figurePlotting/rna_functions.R')



# upload coldata for stressed samples -------------------------------------


coldata = read_csv('metadata/subset_tables/allStress_Coldata.csv') %>% 
  mutate(bleached = if_else(is.na(bleached),
                            'no',
                            bleached))
coldata$my_letter = sapply(coldata$my_title, function(x) strsplit(x, '_')[[1]][1])


# SPLIT INTO CLUSTER A AND CLUSTER B --------------------------------------

ll=load('metadata/corStressProjs.Rdata')
ll

#isolate the runs from each
acoldata = coldata %>% 
  filter(my_title %in% corStressProjs)

bcoldata = coldata %>% 
  filter(my_title %in% lowStressProjs)

dim(acoldata)
dim(bcoldata)
nrow(acoldata) / (nrow(acoldata)+nrow(bcoldata))

#isolate group vectors
agroups = paste('A', acoldata$stress, sep='_')
bgroups = paste('B', bcoldata$stress, sep='_')


#isolate rld data for each group
ll=load('largeIgnored/fullDataset_project_controlled.Rdata')
ll
rld.df = t(datExpr)

adat = rld.df[,acoldata$Run]
dim(adat)
bdat = rld.df[,bcoldata$Run]
dim(bdat)


# RUN DAPC ----------------------------------------------------------------

#function to run dapc
#rld.data = data frame of variance stabilized counts
#group.vector = vector of two grouping strings matching the colnames of rld.data
run_2group_dapc = function(rld.data, group.vector){
  print('running clustering...')
  clus=find.clusters(t(rld.data),max.n.clus=15, n.clust=2, n.pca=10)
  clus$grp=group.vector  #set the transplantation site as the groups
  cdf = data.frame(Run=colnames(rld.data),
                   group = clus$grp)
  print('running dapc...')
  dp=dapc(t(rld.data),clus$grp, n.da=1, perc.pca=80)
}


assemble_dp_coords = function(dp.object){
  res = tibble(Run = rownames(dp.object$ind.coord),
                   LD1 = dp.object$ind.coord[,'LD1'],
                   group = as.character(dp.object$grp))
}

assemble_pred_coords = function(pred.object, group.vector){
  res = tibble(Run = rownames(pred.object$ind.scores),
                   LD1 = pred.object$ind.scores[,'LD1'],
                   group = group.vector)
}

adp = run_2group_dapc(adat, agroups)
bdp = run_2group_dapc(bdat, bgroups)

#extract dataframes of the coordinates
adat2 = assemble_dp_coords(adp)
bdat2 = assemble_dp_coords(bdp)

#look at the initial distributions (both ways to double-check functions)
ALPHA=0.75
scatter(adp,bg="white",scree.da=FALSE,legend=TRUE,solid=ALPHA)
adat2 %>% 
  ggplot(aes(x=LD1,fill=group)) +
  geom_density(alpha=ALPHA)
scatter(bdp,bg="white",scree.da=FALSE,legend=TRUE,solid=ALPHA)
bdat2 %>% 
  ggplot(aes(x=LD1,fill=group)) +
  geom_density(alpha=ALPHA)


# APPLY DAPC TO EACH OTHER --------------------------------------------

#run predictions
pred.AonB<-predict.dapc(adp,newdata=(t(bdat)))
pred.BonA<-predict.dapc(bdp,newdata=(t(adat)))

#extract the dataframes
apb = assemble_pred_coords(pred.AonB, bgroups)
bpa = assemble_pred_coords(pred.BonA, agroups)


#look at distributions
COLORS = c('blue', 'firebrick', 'dodgerblue', 'orange')

#for the group A dapc
rbind(adat2, apb) %>% 
  ggplot(aes(x=LD1,fill=group)) +
  geom_density(alpha=ALPHA) +
  scale_fill_manual(values=COLORS)


#for the group B dapc
rbind(bdat2, bpa) %>% 
  ggplot(aes(x=LD1,fill=group)) +
  geom_density(alpha=ALPHA) +
  scale_fill_manual(values=COLORS)



# USE FOR FRONTLOADING PREDICTIONS ----------------------------------------

ww = c('DC', 'CD', 'CA', 'CB', 'DA', 'DB')
cc = c('AB', 'BA', 'AC', 'AD', 'BC', 'BD')
ww = c('DC', 'CD')
cc = c('AB', 'BA')

#first append the larvae reads from Dixon et al. 2015
larvae.coldata = read_csv('metadata/ALL_Coldata.csv') %>% 
  filter(my_title=='H_matz_heatTolLat_PRJNA279192',
         my_stage=='larval') %>% 
  mutate(my_letter='H',
         stress='larvae',
         parentage=substr(treatDescription, start=1, stop=2),
         lat = if_else(parentage %in% ww,
                       'ww',
                       'unassigned'),
         lat = if_else(parentage %in% cc,
                       'cc',
                       lat)) %>% 
  filter(lat %in% c('ww', 'cc'))



#ADD THE SURIVIVAL DATA FROM DIXON ET AL. 2015
surv = read_csv('metadata/project_specific_tables/dixon15_larvalHeatSurvival.csv') %>% 
  mutate(treatDescription = paste(family, rep, sep=''),
         prop=count/20) %>% 
  left_join(larvae.coldata, by = 'treatDescription')


#set up color coding for families
fams = unique(survProp$family)
cc = 'dodgerblue'
ch = 'black'
hc = 'firebrick'
hh = 'red'
colors = c(cc, ch, ch, cc, ch, ch, hc, hc, hh, hh)

#check the original surivival plot
survProp %>% 
  ggplot(aes(x=time, y=mnProp, color=family)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values=colors)

#make a modified coldata for B
bcoldata.w.larvae = rbind(bcoldata,
                          larvae.coldata[,colnames(bcoldata)])

#make one for A
acoldata.w.larvae = rbind(acoldata,
                          larvae.coldata[,colnames(acoldata)])

#modified counts
bdat.wl = rld.df[,bcoldata.w.larvae$Run]
#modified group vector
bgroups.wl = bcoldata.w.larvae$stress

#repeat the predictions 
pred.AonBwl<-predict.dapc(adp,newdata=(t(bdat.wl)))
pred.BonBwl<-predict.dapc(bdp,newdata=(t(bdat.wl)))

#extract the dataframes
apbwl = assemble_pred_coords(pred.AonBwl, bgroups.wl)
bpbwl = assemble_pred_coords(pred.BonBwl, bgroups.wl)


# PLOT DAPC PREDICTIONS ---------------------------------

#FOR CLUSTER A'S DAPC

#general
rbind(adat2, apbwl) %>% 
  ggplot(aes(x=LD1,fill=group)) +
  geom_density(alpha=ALPHA) +
  labs(x='Cluster A DAPC') +
  scale_fill_manual(values=c('blue', 'red', 'dodgerblue', 'green', 'orange'))


#merge surival data with DAPC A
sdapcA = surv %>% 
  left_join(apbwl, by = 'Run') %>% 
  filter(time > 34)


#DAPC FROM CLUSTER A AND LARVAL SURIVIVAL

#for all points
colors=c(cc, cc, hh, hh)
sdapcA %>% 
  filter(!is.na(LD1)) %>% 
  ggplot(aes(x=LD1, y=prop)) +
  geom_point(aes(color=family)) +
  geom_smooth(method='lm') +
  scale_color_manual(values=colors) +
  labs(x='Cluster A DAPC')
lm1=lm(sdapcA$prop ~ sdapcA$LD1)
summary(lm1)

#for means
sdapcAmns =  sdapcA %>% 
  filter(!is.na(LD1)) %>% 
  group_by(Run, family) %>% 
  summarize(mnProp=mean(prop),
            LD1=mean(LD1))
sdapcAmns %>% 
  ggplot(aes(x=LD1, y=mnProp)) +
  geom_point(aes(color=family)) +
  geom_smooth(method='lm') +
  scale_color_manual(values=colors) +
  labs(x='Cluster A DAPC')
lm1=lm(sdapcA$prop ~ sdapcA$LD1)
summary(lm1)



#DAPC CLUSTER B
sdapcB = surv %>% 
  left_join(bpbwl, by = 'Run') %>% 
  filter(time > 28)
lm1=lm(sdapcB$prop ~ sdapcB$LD1)
summary(lm1)

#plot
sdapcB %>% 
  filter(!is.na(LD1)) %>% 
  ggplot(aes(x=LD1, y=prop)) +
  geom_point(aes(color=family)) +
  geom_smooth(method='lm') +
  scale_color_manual(values=colors) +
  labs(x='cluster B DAPC',
       y='survival at 37 hr')

#density
sdapcB %>% 
  filter(!is.na(LD1)) %>% 
  ggplot(aes(x=LD1, fill=family)) +
  geom_density() +
  scale_fill_manual(values=colors) +
  labs(x='Cluster B DAPC')


#with parentage
aWithPar = rbind(adat2, apbwl) %>% 
  left_join(larvae.coldata[,c('Run', 'lat')], by = 'Run') %>% 
  mutate(group2 = if_else(group == 'larvae',
                          lat,
                          group))
aWithPar %>% 
  ggplot(aes(x=LD1,fill=group2)) +
  geom_density(alpha=ALPHA) +
  labs(x='Cluster A DAPC')

aWithPar %>% 
  filter(group2 %in% c('cc', 'ww')) %>% 
  ggplot(aes(x=LD1,fill=group2)) +
  geom_density(alpha=ALPHA) +
  labs(x='Cluster A DAPC')


#for the group B dapc
rbind(bpbwl) %>% 
  ggplot(aes(x=LD1,fill=group)) +
  geom_density(alpha=ALPHA) +
  labs(x='Cluster B DAPC')


#replot with parentage
bpbwl %>% 
  left_join(larvae.coldata[,c('Run', 'lat')], by = 'Run') %>% 
  mutate(group2 = if_else(group == 'larvae',
                          lat,
                          group)) %>% 
  ggplot(aes(x=LD1,fill=group2)) +
  geom_density(alpha=ALPHA) +
  labs(x='Cluster B DAPC')


# PLOT EIGS FOR WGCNA MODULES ---------------------------------------------

ll=load('wgcna/moduleEigengenes.Rdata')
ll
mEigs$Run=rownames(mEigs)
head(mEigs)


#plot correlation with surivial
surv %>% 
  filter(time>35,
         family %in% c('CD', 'DC', 'AB', 'BA')) %>% 
  left_join(mEigs, by = 'Run') %>% 
  ggplot(aes(x=MEred, y=prop)) +
  geom_point(aes(color=family)) +
  geom_smooth(method='lm') +
  scale_color_manual(values=colors) +
  labs(x='red module')



#plot larvae's red levels
acoldata.w.larvae %>% 
  left_join(mEigs, by = 'Run') %>% 
  ggplot(aes(x=MEred,fill=stress)) +
  geom_density(alpha=ALPHA) +
  labs(x='red module eigengene')

#plot their red levels alone
larvae.coldata %>% 
  left_join(mEigs, by = 'Run') %>% 
  ggplot(aes(x=MEred,fill=lat)) +
  geom_density(alpha=ALPHA) +
  labs(x='red module eigengene')


#plot larvae green levels
larvae.coldata %>% 
  left_join(mEigs, by = 'Run') %>% 
  ggplot(aes(x=MEgreen4,fill=lat)) +
  geom_density(alpha=ALPHA) +
  labs(x='green module eigengene')



# PLOT CORE STRESS GENES --------------------------------------------------

cgenes = read_tsv('results_tables/core_stress_genes.tsv') %>% 
  pull(gene)

get_expression_ranks = function(rld.data){
  trsub = t(rld.data)
  rranks = data.frame(row.names=rownames(trsub))
  for (gene in cgenes){
    glevels = trsub[,gene]
    grank = rank(trsub[,gene])
    rranks[,gene]=grank
  }
  mnRanks=apply(rranks, 1, mean)
  res = data.frame(Run=names(mnRanks),
                   mnRank=mnRanks)
  return(res)
}



rsub = rld.df[cgenes,larvae.coldata$Run]
lrank=get_expression_ranks(rsub)
lrank %>% 
  left_join(larvae.coldata, by='Run') %>% 
  ggplot(aes(x=mnRank,fill=lat)) +
  geom_density(alpha=ALPHA) +
  labs(x='core gene ranks')

lrank %>% 
  left_join(surv, by='Run') %>% 
  filter(time>35) %>% 
  ggplot(aes(x=mnRank, y=prop, color=family)) +
  geom_point() +
  labs(x='core gene ranks')


# LOOK AT POOLS -----------------------------------------------------------

#upload pc data
ll=load('figurePlotting/stress_pca_LD.Rdata')
ll


#POOLS IN k1 (Seneca and Palumbi 2015)

#note that 300 is the highly variable pool, 400 is the medium variable pool
k1_pools = coldata %>% 
  filter(my_title=='k1_Palumbi_lab_heat_resilience_PRJNA274410') %>% 
  separate(treatDescription, into=c('geno', 'num', 'poolNum', 'cnt.treat', 'time')) %>% 
  mutate(origin=if_else(poolNum=='300',
                        'HVP',
                        'MVP'),
         treat2 = paste(treat, origin, sep='_'))


#plot with red module
k1_pools %>% 
  left_join(mEigs, by = 'Run') %>% 
  ggplot(aes(x=MEred,fill=treat2)) +
  geom_density(alpha=ALPHA) +
  labs(x='red module eigengene')

#plot with white module
k1_pools %>% 
  left_join(mEigs, by = 'Run') %>% 
  ggplot(aes(x=MEwhite,fill=treat2)) +
  geom_density(alpha=ALPHA) +
  labs(x='white module eigengene')

#plot PC1
k1_pools %>% 
  left_join(pdat, by = 'Run') %>% 
  ggplot(aes(x=PC1,fill=treat2)) +
  geom_density(alpha=ALPHA) +
  labs(x='PC1')


#plot pools against cluster A dapc
adat2 %>% 
  inner_join(k1_pools, by = 'Run') %>% 
  mutate(treat2 = paste(treat, origin, sep='_')) %>% 
  ggplot(aes(x=LD1, fill=treat2)) +
  geom_density(alpha=ALPHA) +
  labs(x='Cluster A DAPC')

#plot core gene ranks
rsub = rld.df[cgenes,k1_pools$Run]
krank=get_expression_ranks(rsub)
krank %>% 
  left_join(k1_pools, by='Run') %>% 
  ggplot(aes(x=mnRank,fill=treat2)) +
  geom_density(alpha=ALPHA) +
  labs(x='core gene ranks')


#POOLS IN L (Barshis et al. 2013)
L_pools = coldata %>% 
  filter(my_title=='L_Barshis_bleachResillience_PRJNA177515') %>% 
  separate(treatDescription, into=c('geno', 'cnt.treat', 'pool')) %>% 
  mutate(origin=if_else(pool=='HV',
                        'HVP',
                        'MVP'),
         treat2 = paste(treat, origin, sep='_'))


#red module
L_pools %>% 
  left_join(mEigs, by = 'Run') %>% 
  ggplot(aes(x=MEred,fill=treat2)) +
  geom_density(alpha=ALPHA) +
  labs(x='red module eigengene')

#white
L_pools %>% 
  left_join(mEigs, by = 'Run') %>% 
  ggplot(aes(x=MEwhite,fill=treat2)) +
  geom_density(alpha=ALPHA) +
  labs(x='white module eigengene')


#pc1
L_pools %>% 
  left_join(pdat, by = 'Run') %>% 
  ggplot(aes(x=PC1,fill=treat2)) +
  geom_density(alpha=ALPHA) +
  labs(x='PC1')



#plot pools against cluster A dapc
adat2 %>% 
  inner_join(L_pools, by = 'Run') %>% 
  ggplot(aes(x=LD1, fill=treat2)) +
  geom_density(alpha=ALPHA) +
  labs(x='Cluster A DAPC')






