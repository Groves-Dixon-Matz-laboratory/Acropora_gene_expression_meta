#plot_individual_projects.R
library(tidyverse)
library(pheatmap)
library(DESeq2)
source('figurePlotting/rna_functions.R')

# READ IN THE DATA FROM DESEQ ---------------------------------------------

#gather the individual files
iFileList = list.files(path = 'deseqResults/individual_projects/', pattern = '*deseqResults.Rdata', full.names=TRUE)
inames = sub('_deseqResults.Rdata', '', list.files(path = 'deseqResults/individual_projects/', pattern = '*deseqResults.Rdata', full.names=FALSE))
names(iFileList)=inames

#gather grouped files
gFileList = list.files(path = 'deseqResults', pattern = '*deseqResults.Rdata', full.names=TRUE)
gnames0 = sub('_deseqResults.Rdata', '', list.files(path = 'deseqResults/', pattern = '*deseqResults.Rdata', full.names=FALSE))
names(gFileList)=gnames0
toRemove = c('bleached', 'cold_NoBEWW', 'filteredImmune', 'heat_NoBEWW', 'salinity_NoBEWW')
gnames = gnames0[!gnames0 %in% toRemove]
gFileList = gFileList[gnames]
length(gFileList)
names(gFileList)


#read in the files for individual studies
idfList = list()
for (i in 1:length(iFileList)){
  x=iFileList[[i]]
  load(x)
  n=inames[i]
  idfList[[n]]=read_deseq_res(x,n)
}

#read in the stress desea results
gdfList = list()
for (i in 1:length(gFileList)){
  x=gFileList[[i]]
  load(x)
  n=gnames[i]
  gdfList[[n]]=read_deseq_res(x,n)
}



#combine them
cdfList = append(idfList, gdfList)
cdfList=idfList


#merge them
cdat = cdfList %>% 
  purrr::reduce(full_join, by = c('gene'))
rownames(cdat)=cdat$gene
ldat = cdat %>% 
  dplyr::select(grep('_lfc', colnames(cdat)))
colnames(ldat) = sub('_lfc', '', colnames(ldat))

#save write out the table
order = colnames(cdat)
order = order[order!='gene']
order = append(c('gene'), order)
idat = cdat[,order]
save(idat, file= 'deseqResults/all_by_bioproject.Rdata')


# ARRANAGE CORRELATION MATRIX ---------------------------------------------

#get correlations and modify names
noAllBeww = ldat %>% 
  dplyr::select(-j1_thisStudy_PRJNA559404)
c=cor(noAllBeww, use="pairwise.complete.obs")
save(c, file='results_tables/bioproject_stress_correlation_matrix.Rdata')
colnames(c) %>% 
  write.table('./metadata/detailed_tables/stressNamesRaw.tsv', quote=FALSE, row.names=FALSE)
#manually fill in nicer names for plotting
modNames = read.table('./metadata/detailed_tables/stressNamesModified.txt',
                      sep='\t',
                      header=TRUE,
                      row.names='my_title',
                      stringsAsFactors=FALSE)
modNames=modNames[colnames(c),]
projName = paste(modNames$Bioproject, modNames$ref)
projTreat = paste(modNames$treat, modNames$Bioproject)
treat = modNames$treat

#----ARRANGE ANNOTATIONS FOR HEATMAP/CLUSTERING

#setup colors
library(RColorBrewer)
library(vegan)
COLOR = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100) # display.brewer.all()
# COLOR = inferno(100)


#setup bleached projects
adat=read_csv('metadata/subset_tables/allStress_Coldata.csv')
bleachedProjs = adat %>% 
  filter(bleached=='yes' & !is.na(bleached)) %>% 
  pull(my_title) %>% 
  unique()
j1Exps = colnames(c)[grep('beww', colnames(c))]
bleachedProjs = append(bleachedProjs, j1Exps)


##### BUILD MATRIX WITH BLEACHING LABELS FOR HEATMAP
labels = colnames(c)
corWbleach = c[,"bewwMulti"]
rcorWbleach = 1-corWbleach
d <- as.dist(1-c)
h=hclust(d, method='average')
ho=reorder(h, rcorWbleach)
hmLabs = h$labels
blch = hmLabs %in% bleachedProjs
m=matrix(nrow=length(blch), ncol=length(blch))
for (i in 1:length(blch)){
  for (j in 1:length(blch)){
    r=blch[i] & blch[j]
    m[i,j]=r
  }
}
m[m==TRUE]<-'b'
m[m==FALSE]<-''

##### MAKE CALLS ON BIOPROJECT GROUPINGS

plot(h)
cut.height = 1
abline(h = cut.height, col = "red", lty = 2);

# Determine cluster under the line
library(WGCNA)
clust = cutreeStatic(h, cutHeight = cut.height, minSize = 4)
table(clust)
corStressProjs = rownames(c)[(clust==1)]

#look at minimum correlations within these groups
cStrong = c[corStressProjs,corStressProjs]
mean(cStrong)
min(cStrong)

#combine experiments from this study into one
corStressProjs = append(c('j1_thisStudy_PRJNA559404'),
                        corStressProjs[!grepl('^beww', corStressProjs)])
lowStressProjs = rownames(c)[(clust==2)]
save(corStressProjs, lowStressProjs, file='metadata/corStressProjs.Rdata')

#output coldata files for each group
cs = adat %>% 
  filter(my_title %in% corStressProjs)
corStressProjs %in% unique(cs$my_title)
cs %>% 
  write_csv('metadata/subset_tables/corStress_Coldata.csv')

ls = adat %>% 
  filter(my_title %in% lowStressProjs)
lowStressProjs %in% unique(ls$my_title)
ls %>% 
  write_csv('metadata/subset_tables/lowStress_Coldata.csv')


#look at minimum correlations within these groups
cStrong = c[rownames(c)[(clust==1)],rownames(c)[(clust==1)]]
strong0 = as.vector(cStrong)
strong = strong0[strong0 < 0.99999999999]
mean(strong)
min(strong)
cWeak = c[lowStressProjs, lowStressProjs]
weak0 = as.vector(cWeak)
weak = weak0[weak0 < 0.99999999999]
mean(weak)
min(weak)


#make annotaiton labels for heatmap (this didn't look as good but kept here just in case)
anots = data.frame(row.names=colnames(c),
                   cluster = if_else(colnames(c) %in% corStressProjs | grepl('^beww', colnames(c)),
                                   'type A',
                                   'type B'))
anots
cols=gg_color_hue(3)
controlColor = cols[1]
clusterAcol = COLOR[70]
clusterBcol = COLOR[30]
my_colour = list(cluster = c(`type A`=clusterAcol, `type B`=clusterBcol))


#format annotations for cluster and bleaching
blch = colnames(c) %in% bleachedProjs
blch[blch==TRUE]<-'True'
blch[blch==FALSE]<-'False'
clusterA = rownames(anots)[anots$cluster=='type A']
clusterB = rownames(anots)[anots$cluster=='type B']
bleachAnnot = data.frame(bleached = blch,
                         row.names=colnames(c),
                         stringsAsFactors=FALSE)
bleachAnnot$bleached[rownames(bleachAnnot)=='F_Uqueensland_ph_PRJNA269992']<-'mild'
bleachAnnot$bleached[bleachAnnot$bleached=='True']<-'yes'
bleachAnnot$bleached[bleachAnnot$bleached=='False']<-'no'
bleachAnnot$cluster<-NULL
my_colour = list(bleached = c(yes='white', no='coral3', mild = 'coral'),
                 cluster = c(`type A`=clusterAcol, `type B`=clusterBcol))

#check
bleachAnnot %>% 
  mutate(proj = colnames(c),
         label=projName)



#plot 
BREAKS=seq(-.1, 0.5, length.out=length(COLOR)+1)
pheatmap(c,
         labels_row = treat,
         labels_col=projName,
         na_col=COLOR[100],
         color=COLOR,
         breaks=BREAKS,
         cluster_rows = h,
         cluster_cols = h,
         fontsize_number = 10,
         annotation_row = bleachAnnot,
         annotation_col = anots,
         annotation_colors = my_colour,
         treeheight_row=25,
         treeheight_col=0)


#save order for kog heatmap
orderedProjects = h$labels[h$order]
save(orderedProjects, bleachedProjs, anots, my_colour, modNames, projName, bleachAnnot, file='kog_mwu/heatmapObjects.Rdata')



# plot heatmap for random forest with A and B clusters --------------------

#------------- continue here on PC
library(caret)
library(glmnet)
ll=load('category_prediction/stress_prediction_wAandB/stratified_allStress_train.csv__stratified_allStress_test__WITH_A_AND_B.csv_predictionResults_wA_and_B.Rdata')
ll

#format a proportion accuracy table
Ns = table(rfres$obs)
cm=confusionMatrix(data=rfres$pred,
                   reference=rfres$obs,
                   positive="stress_A")
tab = cm$table
colnames(tab)
propCm = tab
propCm[,'control'] = propCm[,'control'] / Ns['control']
propCm[,'stress_A'] = propCm[,'stress_A'] / Ns['stress_A']
propCm[,'stress_B'] = propCm[,'stress_B'] / Ns['stress_B']
propCm
propCm2 = round(propCm, digits=2)


library(pheatmap)
library(RColorBrewer)
genLabs = c('control', 'type A', 'type B')
rannots = data.frame(prediction=genLabs,
                     row.names=rownames(propCm),
                     stringsAsFactors=FALSE)
cannots = data.frame(reference=genLabs,
                     row.names=colnames(propCm),
                     stringsAsFactors=FALSE)
my_colour = list(prediction = c(control=controlColor, `type A`=clusterAcol, `type B` = clusterBcol),
                 reference = c(control=controlColor, `type A`=clusterAcol, `type B` = clusterBcol))


COLOR = colorRampPalette(c("grey", "black"))(100)
pheatmap(propCm,
         labels_row = genLabs,
         labels_col=genLabs,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = propCm2,
         fontsize_number = 15,
         fontsize = 12,
         number_color = 'white',
         color=COLOR,
         legend=FALSE,
         annotation_row = rannots,
         annotation_col = cannots,
         annotation_colors = my_colour,
         annotation_legend = FALSE)


# COMPARE STATS FOR TWO GROUPS --------------------------------------------

#First show that log2 fold differences are greater in cluster A
head(ldat)
al = ldat[,corStressProjs] %>% 
  gather(key='proj',
         value='lfc') %>% 
  mutate(cluster='A')
bl = ldat[,lowStressProjs] %>% 
  gather(key='proj',
         value='lfc') %>% 
  mutate(cluster='B')
c=rbind(al, bl)
c %>% 
  ggplot(aes(x=cluster,y=abs(lfc))) +
  geom_boxplot() +
  lims(y=c(0,1.5)) +
  labs(y=bquote("|"*log[2]~' fold change|'))

alfc = abs(al$lfc)
blcf = abs(bl$lfc)
t.test(alfc, blcf)
median(alfc, na.rm=TRUE)
median(blcf, na.rm=TRUE)
median(alfc, na.rm=TRUE) / median(blcf, na.rm=TRUE)


#compare the medians for the individual bioprojects
med_c = c %>% 
  mutate(alfc = abs(lfc)) %>% 
  group_by(proj, cluster) %>% 
  summarize(med_alfc = median(alfc, na.rm=TRUE))

t.test(med_c$med_alfc ~ med_c$cluster) #not significant

med_c %>% 
  mutate(cluster=factor(cluster, levels = c('A', 'B')),
         num = as.numeric(cluster)) %>% 
  ggplot(aes(x=cluster, y=med_alfc)) +
  geom_violin() +
  geom_jitter(aes(x=num, y=med_alfc), size=3, width=0.1) +
  labs(y = bquote('median abs log'[2]~'fold difference'))

med_c %>% 
  arrange(cluster, med_alfc)

#next show temperatures were higher in cluster A
temps = read_csv('./metadata/detailed_tables/detailedHeat.csv')
mnt = temps %>% 
  mutate(cluster = if_else(my_title %in% corStressProjs,
                           'A',
                           'B')) %>% 
  filter(stress=='stressed') %>% 
  group_by(my_title, cluster) %>% 
  summarize(mnTemp=mean(temp),
            mnOver=mean(degOverAmb))
fmns = mnt %>% 
  group_by(cluster) %>% 
  summarize(mnT=mean(mnTemp),
            medT=median(mnTemp),
            mnO=mean(mnOver),
            medO=median(mnOver))
fmns

mnt %>% 
  ggplot(aes(x=mnTemp, fill=cluster)) +
  geom_density(alpha=0.5) +
  geom_vline(xintercept = fmns$medT, aes(color=fmns$cluster))

mnt %>% 
  filter(my_title!='f1_parkinson_hotCold_PRJNA423227') %>% 
  ggplot(aes(x=mnTemp, fill=cluster)) +
  geom_density(alpha=0.5)


# compare high and low stress against pca and dapc -------------------------------------
source('figurePlotting/rna_functions.R')
input = 'largeIgnored/stress_project_controlled.Rdata'
ll=load(input)
ll


#reformat with stress lvl based on the projects identified in heatmap
sdt = adat %>% 
  filter(!my_title=='d1_mohamed_microalga_PRJNA398338') %>% 
  mutate(stressLvl = if_else(my_title %in% corStressProjs & stress=='stressed',
                              'type A',
                              'unassigned'),
         stressLvl = if_else(my_title %in% corStressProjs & stress=='control',
                             'control',
                             stressLvl),
         stressLvl = if_else(my_title %in% lowStressProjs & stress=='control',
                             'control',
                             stressLvl),
         stressLvl = if_else(my_title %in% lowStressProjs & stress=='stressed',
         'type B',
         stressLvl),
         stressLvl = factor(stressLvl, levels=c('control', 'type A', 'type B')))
srld = datExpr[sdt$Run,] %>% 
  t() 

u=sdt %>% 
  filter(stressLvl=='unassigned')
unique(u$my_title)

#look where they fall on stress PCA and DAPC
ALPHA=0.85
#pca
ll=load('figurePlotting/stress_pca_LD.Rdata')
ll
head(pdat)
pcaLvl = sdt %>% 
  select(stressLvl, Run) %>% 
  left_join(pdat, by = 'Run') %>% 
  ggplot(aes(x=PC1, fill=stressLvl)) +
  geom_density(alpha=ALPHA) +
  labs(x='PC1',
       fill='') +
  scale_fill_manual(values=c(controlColor, clusterAcol, clusterBcol)) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


#dapc
dapcLvl = sdt %>% 
  select(stressLvl, Run) %>% 
  left_join(pdat, by = 'Run') %>% 
  ggplot(aes(x=LD1, fill=stressLvl)) +
  geom_density(alpha=ALPHA) +
  labs(fill='') +
  scale_fill_manual(values=c(controlColor, clusterAcol, clusterBcol)) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

#repeat for bleached
pcaLvlb = sdt %>% 
  select(Run) %>% 
  left_join(pdat, by = 'Run') %>% 
  mutate(bleached = ifelse(is.na(bleached),
                            'no',
                            bleached),
         blchLvl = if_else(bleached=='yes' & stress.x=='stressed',
                           'stress bleached',
                           'unassigned'),
         blchLvl = if_else(bleached!='yes' & stress.x=='stressed',
                           'stress not bleached',
                           blchLvl),
         blchLvl = if_else(stress.x=='control',
                           'control',
                           blchLvl)) %>% 
  ggplot(aes(x=PC1*-1, fill=blchLvl)) +
  geom_density(alpha=0.5) +
  labs(x='PC1',
       fill='none') +
  scale_fill_manual(values=c(controlColor, clusterAcol, clusterBcol)) 


#plot together
plot_grid(pcaLvl+theme(legend.position='none'), dapcLvl, rel_widths=c(0.75,1))



#original dapc
head(pdat)
pdat %>% 
  ggplot(aes(x=LD1, fill=stress.y)) +
  geom_density()


#correlate with zoox loss
ll=load('metadata/coldataWithSymbionts.Rdata')
ll
psdat  = pdat %>% 
  left_join(coldata, by = 'Run') 

psdat %>% 
  mutate(stress = if_else(Treatment=='control',
                          'control',
                          'stressed'),
         PC1=PC1*-1) %>% 
  ggplot() +
  geom_smooth(aes(x=PC1, y=zSym), method='lm') +
  geom_point(aes(x=PC1, y=zSym, color=stress)) +
  labs(y='symbiont amount (z-score)')


# look for high and low within projects -----------------------------------

projects = unique(sdt$my_title)
projects
pr="b1_Aguilar_hypoosmotic stress_PRJNA380267"
psdat = sdt %>% 
  select(stressLvl, Run) %>% 
  left_join(pdat, by = 'Run') %>% 
  rename(stress.x='stress') %>% 
  mutate(PC1=PC1*-1)
plotList = list()
for (pr in projects){
  if (pr %in% corStressProjs){
    lvl='High'
  } else {
    lvl='Low'
  }
  plt=psub = psdat %>% 
    filter(my_title==pr) %>% 
    ggplot(aes(x=PC1, fill=stress)) +
    geom_density(alpha=0.5) +
    labs(title=pr, subtitle=lvl)
  plotList[[pr]]=plt
}

plotList[[pr]]
plot_grid(plotlist=plotList, nrow=4)


projects
#b1_Aguilar_hypoosmotic stress_PRJNA380267
pr="b1_Aguilar_hypoosmotic stress_PRJNA380267"
psub = psdat %>% 
  filter(my_title==pr) %>% 
  separate(treatDescription, into=c('time', 'psu', 'geno', 'spp'), sep='_')

pcaPlain=psub %>% 
  ggplot(aes(x=PC1, fill=stress)) +
  geom_density(alpha=0.5) +
  lims(x=c(-15,15)) +
  labs(subtitle='all')


pca1h=psub %>% 
  filter(time=='1h') %>% 
  mutate(treatTime = paste(stress, time, sep='.')) %>% 
  ggplot(aes(x=PC1, fill=psu)) +
  geom_density(alpha=0.6)+
  lims(x=c(-15,15)) +
  labs(subtitle='1h')

pca24h=psub %>% 
  filter(time=='24h') %>% 
  mutate(treatTime = paste(stress, time, sep='.')) %>% 
  ggplot(aes(x=PC1, fill=psu)) +
  geom_density(alpha=0.6)+
  lims(x=c(-15,15)) +
  labs(subtitle='24h')

pca48h=psub %>% 
  filter(time=='48h') %>% 
  mutate(treatTime = paste(stress, time, sep='.')) %>% 
  ggplot(aes(x=PC1, fill=psu)) +
  geom_density(alpha=0.6)+
  lims(x=c(-15,15)) +
  labs(subtitle='48h')

plot_grid(pcaPlain, pca1h, pca24h, pca48h, nrow=4)



# look at k1_palumbi ------------------------------------------------------

head(pdat)
pdat %>% 
  filter(my_title=='k1_Palumbi_lab_heat_resilience_PRJNA274410') %>% 
  mutate(time = if_else(grepl('_5h', treatDescription),
                        '5h',
                        '20hr'),
         treatTime = paste(treat, time, sep='.')) %>% 
  ggplot(aes(x=PC1, fill=treatTime)) +
  geom_density()
  
