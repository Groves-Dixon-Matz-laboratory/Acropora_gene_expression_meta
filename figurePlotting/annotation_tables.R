#annotation_tables.R

#organize large table with all the annotations from the project:
#1. log2 and p-values from each invididual project
#2. log2 and p-values from stress types
  #cluster A (high stress)
  #cluster B (low stress)
  #bleached
  #stres without bleaching
  #heat
  #heat no bleached
  #hyposalinity
  #hyposalinity no bleached (this is just Aguilar et al. 2017)
  #immune
  #pH
#3. WGCNA module membership values
#4. variable importance measures for:
  #stress PC1
  #DAPC
  #lasso regression
  #Random Forest


#other tables to output:
#1. raw counts against reference genome
#2. raw vsd table
#3. counts controlled vsd counts
#4. project controlled vsd counts
#5. zoox counts


library(tidyverse)
source('figurePlotting/rna_functions.R')



# LOG2 AND P-VALUES FOR INDIVIDUAL PROJECTS -------------------------------

#here we want all the deseq2 results for each individual file
mnames = read_tsv('metadata/detailed_tables/stressNamesModified.txt')
ll=load('deseqResults/all_by_bioproject.Rdata')
ll
#swap out for concise names
titles0 = colnames(idat)[colnames(idat)!='gene']
titles1=titles0[grep('_lfc', titles0)]
titles = sub('_lfc', '', titles1)
newNames = swap_proj_for_author_Bioproj(titles)
newColNames = colnames(idat)
for (i in 1:length(titles)){
  t=titles[i]
  n=newNames[i]
  newColNames = sub(t,n,newColNames)
}
data.frame(colnames(idat),newColNames)
colnames(idat)=newColNames
head(idat)

#look at Sacsin
ss = idat %>% 
  filter(gene %in% c('Amillepora11972',
                     'Amillepora11973',
                     'Amillepora11974'))
ss


# LOG2 AND P-VALUES FOR STRESS TYPES --------------------------------------

#in addition to the individual BioProjects, we want all type A and all type B

#get the high-stress and low-stress gruops
clusterPaths = c('deseqResults/corStress_deseqResults.Rdata',
                 'deseqResults/lowStress_deseqResults.Rdata',
                 'deseqResults/stress_deseqResults.Rdata')
clusterNames = c('clusterA',
                 'clusterB',
                 'allStressStudies')
names(clusterPaths)=clusterNames

#read in the stress desea results
cgList = list()
for (i in 1:length(clusterPaths)){
  x=clusterPaths[[i]]
  load(x)
  n=clusterNames[i]
  cgList[[n]]=read_deseq_res(x,n)
}

#merge them
cgdat = cgList %>% 
  purrr::reduce(full_join, by = c('gene'))
orderCols = c('gene', colnames(cgdat)[colnames(cgdat)!='gene'])
cgdat=cgdat[,orderCols]
rownames(cgdat)=cgdat$gene

#merge with results from individual bioprojects
head(idat)
head(cgdat)
deseq_dat = cgdat %>% 
  full_join(idat, by = 'gene') %>% 
  arrange(clusterA_p)
head(deseq_dat)
dim(deseq_dat)

#write out
write_tsv(deseq_dat, path='results_tables/all_deseqResults.tsv')


# WGCNA MODULE MEMBERSHIP -------------------------------------------------

#load the module membership data from wgcan
ll=load('wgcna/moduleAssignment.Rdata')
ll
head(geneModuleMembership)
head(moduleColors)
gmdat = geneModuleMembership
gmdat$gene=rownames(geneModuleMembership)
orderedCols = c('gene', colnames(geneModuleMembership))
gmdat=gmdat[,orderedCols]
gmdat$assignment=moduleColors
head(gmdat)
write_tsv(gmdat, path='results_tables/wgcna_moduleMembership.tsv')


# CONTEXTUAL ANNOTATIONS  -------------------------------------------------------

#FIRST GET CORE SET
adat = read_tsv('metadata/Amillepora_euk.emapper.annotations.tsv') %>% 
  dplyr::rename('gene'=query_name)
#pull candidates from idat
ll=load('metadata/corStressProjs.Rdata')
newCorProj = swap_proj_for_author_Bioproj(corStressProjs)
ildat = idat[,paste(newCorProj, 'lfc', sep='_')]
lfcCut=0.59 #a 1.5 fold change
nup = apply(ildat, 1, function(x) sum(x > lfcCut, na.rm=TRUE))
ndown = apply(ildat, 1, function(x) sum(x < -lfcCut, na.rm=TRUE))
CUT=9
either = nup >= CUT #decided to just go up
table(nup)
table(ndown)
sum(either)
sidat = ildat[either, ]
sidat$gene=rownames(sidat)

#important module membership
smdat = gmdat %>% 
  filter(assignment=='red') %>% 
  select(gene, MMred)

#merge them all up
datList = list(sidat, smdat)
fdat = datList %>% 
  purrr::reduce(full_join, by='gene') %>% 
  filter(gene %in% sidat$gene) %>% 
  left_join(adat, by = 'gene') %>% 
  dplyr::select('gene', 'seed_eggNOG_ortholog', 'predicted_gene_name', 'eggNOG.annot', colnames(sidat))
dim(fdat)
fdat$eggNOG.annot[is.na(fdat$eggNOG.annot)]<-'none'
t=table(fdat$eggNOG.annot)
t
t['none']/sum(t)
sum(t)


# add less stringent set --------------------------------------------------
lfcCut=0.32 #a 1.25 fold change
nup = apply(ildat, 1, function(x) sum(x > lfcCut, na.rm=TRUE))
ndown = apply(ildat, 1, function(x) sum(x < -lfcCut, na.rm=TRUE))
CUT=8
upgenes = names(nup)[nup>=CUT]
downgenes = names(ndown)[ndown>=CUT]
length(upgenes)
length(downgenes)
cgenes = rownames(ildat)
context = data.frame(gene = cgenes,
                     general_stress_stringent = ifelse(cgenes %in% fdat$gene,
                                             'up',
                                             NA),
                     general_stress_relaxed = NA) %>% 
  mutate(general_stress_relaxed = ifelse(cgenes %in% upgenes,
                                         'up',
                                         general_stress_relaxed),
         general_stress_relaxed = ifelse(cgenes %in% downgenes,
                                         'down',
                                         general_stress_relaxed))
table(context$general_stress_stringent)
table(context$general_stress_relaxed)



#GET UNIVERSAL IN BLEACHED
ll=load('kog_mwu/heatmapObjects.Rdata')
bleachedProjs = bleachedProjs[bleachedProjs != 'F_Uqueensland_ph_PRJNA269992']
bleachedProjs = bleachedProjs[bleachedProjs != 'F_Uqueensland_ph_PRJNA269992']
colnames(ildat)
blch.ildat = ildat[,c('thisStudyAll_PRJNA559404_lfc',
                      "Seneca2015_PRJNA274410_lfc",
                      "Barshis2013_PRJNA177515_lfc")]
head(blch.ildat)
bgenes = rownames(blch.ildat)

lfcCut=0.59 #a 1.5 fold change
nup = apply(blch.ildat, 1, function(x) sum(x > lfcCut, na.rm=TRUE))
ndown = apply(blch.ildat, 1, function(x) sum(x < -lfcCut, na.rm=TRUE))
CUT=3
upgenes = names(nup)[nup>=CUT]
downgenes = names(ndown)[ndown>=CUT]
length(upgenes)
length(downgenes)
blch.ildat$bleaching = NA
blch.ildat$bleaching[bgenes %in% upgenes]<-'up'
blch.ildat$bleaching[bgenes %in% downgenes]<-'down'
table(blch.ildat$bleaching)
bres = data.frame(gene = bgenes,
                  bleaching = blch.ildat$bleaching)

#merge with context
context2 = context %>% 
  left_join(bres, by = 'gene')

#add general annotations
context3 = context2 %>% 
  left_join(adat, by = 'gene')
dim(context3)


#write out
context3 %>% 
  write_tsv(path='results_tables/contextual_annots.tsv')
