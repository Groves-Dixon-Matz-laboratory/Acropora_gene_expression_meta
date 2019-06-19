#plot_stress_scatterplots.R

library(DESeq2)
library(tidyverse)
library(cowplot)
source('figurePlotting/rna_functions.R')

# load the deseq data for each stresser --------------------

stress=load_deseq('./deseqResults/stress_deseqResults.Rdata', 'stress')
heat = load_deseq('./deseqResults/heat_deseqResults.Rdata', 'heat')
heatNob = load_deseq('./deseqResults/heatNoBEWW_deseqResults.Rdata', 'heatNoB')
cold = load_deseq('./deseqResults/cold_deseqResults.Rdata', 'cold')
coldNob = load_deseq('./deseqResults/cold_NoBEWW_deseqResults.Rdata', 'coldNoB')
salinity = load_deseq('./deseqResults/salinity_deseqResults.Rdata', 'salinity')
salinityNob = load_deseq('./deseqResults/salinity_NoBEWW_deseqResults.Rdata', 'salinityNoB')
bleach = load_deseq('./deseqResults/bleached_deseqResults.Rdata', 'bleached')
immune = load_deseq('./deseqResults/immune_deseqResults.Rdata', 'immune')
ph = load_deseq('./deseqResults/ph_deseqResults.Rdata', 'ph')



# write out GO MWU input files --------------------------------------------
write_out_go(stress, './go_mwu/stressForGO.csv')
write_out_go(heat, './go_mwu/heatForGO.csv')
write_out_go(cold, './go_mwu/coldForGO.csv')
write_out_go(salinity, './go_mwu/salinityForGO.csv')
write_out_go(immune, './go_mwu/immuneForGO.csv')



# GET CORE STRESS GENES ------------------------------------------------
#genes that are significant in the general stress set, and each specific one (with same log2fc direction)

pvalueType = 'padj'
sigCut = 0.1

#MERGE THEM ALL TOGETHER AS adat
adat = data.frame(stress)
colnames(adat)=paste('allStress', colnames(adat), sep='_')
colnames(adat)[1]='gene'
specifics = list(heat, cold, salinity, immune)
names = c('heat', 'cold', 'salinity', 'immune')

for (i in 1:length(specifics)){
  s=specifics[[i]]
  n=names[i]
  colnames(s)=paste(n, colnames(s), sep='_')
  colnames(s)[1]='gene'
  adat = adat %>% 
    full_join(s, by = 'gene')
}

#THEN SUBSET FOR FOLDCHANGES AND PVALUES
specificCols = colnames(adat)[!grepl('allStress', colnames(adat))]
ldat = adat[,specificCols[grep('log2Fold', specificCols)]]
pdat = adat[,specificCols[grep(pvalueType, specificCols)]]
head(ldat)
head(pdat)


#SET UP DFS FOR WHETHER SIGN ON LOG CHANGE MATCHES GENERAL AND P-VALUE IS SIGNIFICANT
genUp = adat$allStress_log2FoldChange >=0
sUp = ldat>=0
sMatch = data.frame(sUp[,1]==genUp,
                    sUp[,2]==genUp,
                    sUp[,3]==genUp,
                    sUp[,4]==genUp)
colnames(sMatch) = sapply(colnames(ldat), function(x) strsplit(x, '_')[[1]][1])
sSig = pdat < 0.1

#SET UP FINAL DF WITH WHETHER BOTH ARE TRUE
fdat = data.frame(sMatch[,1] & sSig[,1],
                  sMatch[,2] & sSig[,2],
                  sMatch[,3] & sSig[,3],
                  sMatch[,4] & sSig[,4])
colnames(fdat) = colnames(sMatch)

#GET TOTALS
matchSums = apply(fdat, 1, sum)
table(matchSums)
head(adat)
genSig = adat$allStress_padj < sigCut
adat$coreStress = matchSums > 3 & genSig
sum(adat$coreStress, na.rm=TRUE)
coreGenes = adat %>% 
  filter(coreStress) %>% 
  pull(gene)
length(coreGenes)

# GET SPECIFIC STRESS GENES -----------------------------------------------

#genes that are only significant for the particular stress
head(pdat)
head(sSig)

#isolate the genes significant for single dataset
sigSums = apply(sSig, 1, sum)
sigSpecific = sigSums <= 0.5*max(sigSums, na.rm=T)
sigSpecific = sigSums==1

#make stress-specific calls
adat$heatSpecific = sSig[,paste('heat', pvalueType,sep='_')] & sigSpecific 
adat$coldSpecific = sSig[,paste('cold', pvalueType,sep='_')] &  sigSpecific
adat$salinitySpecific = sSig[,paste('salinity', pvalueType,sep='_')] & sigSpecific
adat$immuneSpecific = sSig[,paste('immune', pvalueType,sep='_')] &  sigSpecific

#how many are there?
sum(adat$heatSpecific, na.rm=T)
sum(adat$coldSpecific, na.rm=T)
sum(adat$salinitySpecific, na.rm=T)
sum(adat$immuneSpecific, na.rm=T)
sum(adat$heatSpecific & adat$coldSpecific, na.rm=T)
sum(adat$heatSpecific & adat$salinitySpecific, na.rm=T)
sum(adat$heatSpecific & adat$immuneSpecific, na.rm=T)


save(adat, file='deseqResults/coreStressDf.Rdata')


#grab vectors for each set of specific genes
heatGenes = adat %>% 
  filter(heatSpecific) %>% 
  pull(gene)
coldGenes= adat %>% 
  filter(coldSpecific) %>% 
  pull(gene)
saltGenes = adat %>% 
  filter(salinitySpecific) %>% 
  pull(gene)
immuneGenes = adat %>% 
  filter(immuneSpecific) %>% 
  pull(gene)
sum(heatGenes %in% coreGenes)
sum(coldGenes %in% coreGenes)
sum(saltGenes %in% coreGenes)
sum(immuneGenes %in% coreGenes)



# BUILD VOLCANO PLOTS -----------------------------------------------------


vheat = plot_overlay_volcano(heat, 'heat', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vheatNob = plot_overlay_volcano(heatNob, 'heat', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vcold = plot_overlay_volcano(cold, 'cold', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vcoldNob = plot_overlay_volcano(coldNob, 'cold', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vsalinity = plot_overlay_volcano(salinity, 'salinity', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vsalinityNob = plot_overlay_volcano(salinityNob, 'salinity', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vbleach = plot_overlay_volcano(bleach, 'bleached', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vimmune = plot_overlay_volcano(immune, 'immune challenge', overlayGenes1 = coreGenes, overlayGenes2 =FALSE)
vph = plot_overlay_volcano(ph, 'pH', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)

plot_grid(vheat, vcold, vsalinity, vimmune, vph, vbleach, nrow=6)
plot_grid(vheatNob, vcoldNob, vsalinityNob, vbleach, vimmune, vph, nrow=6)

# Run PCAs ----------------------------------------------------------------

#load data
vsheat = load_vsd('./normCounts/heat_project_controlled.Rdata')
vsheatNob = load_vsd('./normCounts/heat_NOBEWW_project_controlled.Rdata')
vscold = load_vsd('./normCounts/cold_project_controlled.Rdata')
vscoldNob = load_vsd('./normCounts/cold_NOBEWW_vsd.Rdata')
vssalt = load_vsd('./normCounts/salinity_project_controlled_adultOnly.Rdata')
vssaltNob = load_vsd('./normCounts/salinity_NOBEWW_vsd_adultOnly.Rdata')
vsbleach = load_vsd('./normCounts/bleached_vsd.Rdata')
vsimmune = load_vsd('./normCounts/immune_project_controlled.Rdata')
vsph = load_vsd('./normCounts/ph_project_controlled.Rdata')

#set plotting vars
fraction = 1/10
NTOP = round(fraction * nrow(vscold[[1]]), digits=0)


pcheat = plotStressPCA(df = vsheat[[1]], coldat = vsheat[[2]], intgroup = 'treat', ntop=NTOP, main = 'heat', xInvert=-1)
pcheatNob = plotStressPCA(df = vsheatNob[[1]], coldat = vsheatNob[[2]], intgroup = 'treat', ntop=NTOP, main = 'heat', xInvert=-1)
pccold = plotStressPCA(df = vscold[[1]], coldat = vscold[[2]], intgroup = 'treat', ntop=NTOP, main = 'cold', xInvert=1)
pccoldNob = plotStressPCA(df = vscoldNob[[1]], coldat = vscoldNob[[2]], intgroup = 'treat', ntop=NTOP, main = 'cold', xInvert=1, pcs=4)
pcsalt = plotStressPCA(df = vssalt[[1]], coldat = vssalt[[2]], intgroup = 'treat', ntop=NTOP, main = 'hyposalinity')
pcsaltNob = plotStressPCAwShape(df = vssaltNob[[1]], coldat = vssaltNob[[2]], intgroup = 'treat', ntop=NTOP, main = 'hyposalinity')
pcbleach = plotStressPCA(df = vsbleach[[1]], coldat = vsbleach[[2]], intgroup = 'bleached', ntop=NTOP, main = 'bleached')
pcimmune = plotStressPCA(df = vsimmune[[1]], coldat = vsimmune[[2]], intgroup = 'treat', ntop=NTOP, main = 'immune challenge', xInvert=-1)
pcph = plotStressPCA(df = vsph[[1]], coldat = vsph[[2]], intgroup = 'treat', ntop=NTOP, main = 'pH')



# PLOT EACH PARTICULAR AGAINST GENERAL MINUS THAT STRESS -----------------------

mheat = load_deseq('./deseqResults/stressNoHeat_deseqResults.Rdata', 'stress_minus_heat')
mcold = load_deseq('./deseqResults/stressNoCold_deseqResults.Rdata', 'stress_minus_heat')
msalinity = load_deseq('./deseqResults/stressNoSalinity_deseqResults.Rdata', 'stress_minus_salinity')
mbleach = load_deseq('./deseqResults/stressNoBEWW_deseqResults.Rdata', 'stress_minus_beww')
mimmune= load_deseq('./deseqResults/stressNoImmune_deseqResults.Rdata', 'stress_minus_immune')
mph = load_deseq('./deseqResults/stressNoPH_deseqResults.Rdata', 'stress_minus_ph')

XLAB='\n'
YLAB='\n'
hplt = do_scatter(heat, mheat, XLAB, YLAB, 'heat vs general') #all heat against all non-heat stress
cplt = do_scatter(cold, mcold, XLAB, YLAB, 'cold vs general') #all cold against all non-cold stress
splt = do_scatter(salinity, msalinity, XLAB, YLAB, 'salinity vs general') #all salinity against all non-salinity
iplt = do_scatter(immune, mimmune, XLAB, YLAB, 'immune vs general') #all immune against all non-immune
pplt = do_scatter(ph, mph, XLAB, YLAB, 'pH vs general')             #all pH against all non -H
hpltNob = do_scatter(heatNob, mheat, XLAB, YLAB, 'heat vs general')  #heat without BEWW against all non-heat
cpltNob = do_scatter(coldNob, mcold, XLAB, YLAB, 'cold vs general')  #cold without BEWW against all non-cold
spltNob = do_scatter(salinityNob, msalinity, XLAB, YLAB, 'salinity vs general') #salinity without BEWW against all non-salinity
bplt = do_scatter(bleach, mbleach, XLAB, YLAB, 'bleach vs general')          #BEWW samples against all non-BEWW stress


#example of overlay plot
hplt = do_scatter_overlay(heat, mheat, XLAB, YLAB, 'heat', coreGenes)
cplt = do_scatter_overlay(cold, mcold, XLAB, YLAB, 'cold', coreGenes)
splt = do_scatter_overlay(salinity, msalinity, XLAB, YLAB, 'salinity', coreGenes)
iplt = do_scatter_overlay(immune, mimmune, XLAB, YLAB, 'immune challenge', coreGenes)
pplt = do_scatter_overlay(ph, mph, XLAB, YLAB, 'low pH', coreGenes)
bplt = do_scatter_overlay(bleach, mbleach, XLAB, YLAB, 'bleach', coreGenes)


#plot main figure with chosen set
plot_grid(pcheat, hplt,
          pccold, cplt,
          pcsalt, splt,
          pcimmune, iplt,
          pcph, pplt,
          nrow=5, 
          rel_widths = c(1,0.75))

#plot supplemental
plot_grid(pcbleach, bplt,
          pcheatNob, hpltNob,
          pcimmune, iplt,
          pcsaltNob, spltNob,
          pcph, pplt,
          nrow=5, 
          rel_widths = c(1,0.75))






# SCATTERPLOTS AGAINST ADULT V LARVA FOR REFERENCE ------------------------

avl = load_deseq('./deseqResults/adultVlarva_deseqResults.Rdata', 'adult_v_larva')
vdev = plot_volcano(avl, 'adult vs larva')
hplt = do_scatter(heat, avl, NULL, NULL, 'heat vs age')
cplt = do_scatter(cold, avl, NULL, NULL, 'cold vs age')
splt = do_scatter(salinity, avl, NULL, NULL, 'salinity vs age')
iplt = do_scatter(immune, avl, NULL, NULL, 'immune vs age')
pplt = do_scatter(ph, avl, NULL, NULL, 'pH vs age')

plot_grid(vdev,
          hplt,
          cplt,
          splt,
          iplt,
          pplt,
          nrow=6)

# BUILD SCATTERPLOTS OF EACH PARTICULAR AGAINST EACH -------------------------

datList = list(heat, cold, salinity, immune, ph)
names = c('heat', 'cold', 'salinity', 'immune', 'pH')
pltList = list()

for (i in 1:length(datList)){
    for (j in 1:length(datList)){
      df1=datList[[i]]
      df2=datList[[j]]
      name1=names[i]
      name2=names[j]
      compareString = paste(name1, name2, sep=' vs ')
      print(compareString)
      # plt=do_scatter(df1, df2, NULL, NULL, compareString)
      plt=do_scatter_overlay(df1, df2, NULL, NULL, compareString, coreGenes)
      pltList[[compareString]]=plt
  }
}
names(pltList)

#replace self-comparisons with volcano plots
pltList[['heat vs heat']] = vheat
pltList[['cold vs cold']] = vcold
pltList[['salinity vs salinity']] = vsalinity
pltList[['immune vs immune']] = vimmune
pltList[['pH vs pH']] = vph


plot_grid(plotlist=pltList, nrow=5)











# apply DAPCs to each dataset ---------------------------------------------

#load the rld.dfs for each
pathList = list('./heat/data/project_controlled.Rdata',
                './cold/data/project_controlled.Rdata',
                './salinity/data/project_controlled.Rdata',
                './immune/data/project_controlled.Rdata',
                './ph/data/project_controlled.Rdata')
names=c('heat', 'cold', 'salinity','immune', 'pH')
expList = list()
traitList = list()
allGenes = c()

for (i in 1:length(names)){
  path=pathList[[i]]
  name=names[i]
  print('-------')
  print(name)
  print(path)
  load(path)
  expList[[name]]=datExpr
  traitList[[name]]=datTraits
  allGenes = append(allGenes, colnames(datExpr))
}
sharedGenes = names(table(allGenes))[table(allGenes)==length(pathList)]
length(expList)
length(traitList)
length(allGenes)
length(sharedGenes)
save(sharedGenes, file='./sharedGenes.Rdata')

# REDUCE TO SHARED GENES --------------------------------------------------
for (n in names){
  print('--------')
  print(n)
  dat = expList[[n]]
  sub = dat[,sharedGenes]
  expList[[n]]=sub
  print('before:')
  print(ncol(dat))
  print('after:')
  print(ncol(sub))
}


# RUN DAPC ON EACH --------------------------------------------------------

library(adegenet)

# #modify stress trait table to resemble the others
# stressTraits = traitList[['stress']]
# stressTraits$treat=stressTraits$stress
# traitList[['stress']]=stressTraits
# 
#modify the immune dataset to remove NAs (This is temporary fix, have removed these earlier on)
icoldata = traitList[['immune']] %>% 
  filter(!is.na(treat))
iexp = expList[['immune']][icoldata$Run,]
sum(rownames(iexp)==icoldata$Run)==nrow(iexp)
traitList[['immune']]=icoldata
expList[['immune']]=iexp

#loop through datasets and run DAPC
dpList = list()
for (n in names){
  print('---------')
  print(n)
  dat=expList[[n]]
  coldata=traitList[[n]]
  print('Running PCA...')
  pcp=prcomp(dat, retx=TRUE, center=TRUE, scale=TRUE)
  scores=pcp$x
  print('Clustering...')
  clus=find.clusters(dat,max.n.clus=15, n.clust=2, n.pca=10)
  cdf = data.frame(clus$grp)
  clus$grp=coldata$treat
  cdf = data.frame(clus$grp)
  cdf$Run = coldata$Run
  print('Running DAPC...')
  dp=dapc(dat,clus$grp, n.da=1, perc.pca=80)
  scatter(dp,bg="white",scree.da=FALSE,legend=TRUE,solid=.4)
  dpList[[n]]=dp
}

save(dpList, file='./dapcList.Rdata')


# APPLY DAPCs TO EACH OTHER DATASET ---------------------------------------

altList = list()

name1 = 'immune'
name2 = 'heat'
pair = paste(name1,name2,sep='-')

#loop through each pair of datasets and apply their dapcs against each other
for (i in 1:(length(names)-1)){
  name1 = names[i]
  rest=(i+1):length(names)
  for (j in rest){
    name2 = names[j]
    pair = paste(name1,name2,sep='-')
    print('-----------')
    print(pair)
    #set up data
    dp1 = dpList[[name1]]
    dat2 = expList[[name2]]
    coldat2 = traitList[[name2]]
    #apply dapc
    pred.sup<-predict.dapc(dp1,newdata=dat2) 
    #format results
    pred.dp<-dp
    pred.dp$ind.coord<-pred.sup$ind.scores
    pred.dp$posterior<-pred.sup$posterior
    pred.dp$assign<-pred.sup$assign
    pred.dp$grp<-as.factor(coldat2$treat)
    #save in list
    altList[[pair]]=pred.dp
  }
}


plot_predicted_ghost = function(pair){
  pdp = altList[[pair]]
  df = data.frame(pdp$ind.coord)
  df$grp = as.vector(pdp$grp)
  df %>% 
    ggplot(aes(x=LD1,color=grp, fill=grp)) +
    geom_density(alpha=0.5) +
    labs(subtitle=pair)
}

pltList = lapply(names(altList), function(x) plot_predicted_ghost(x))


plot_grid(plotlist = pltList)


# IDENTIFY GENERAL VS SPECIFIC STRESS GENES -------------------------------

#for a general stress gene, it must be significant in all the datasets
#for a stress-specific gene, it must be significant in an individual dataset and NOT significant in the all stress dataset
#assumption here is that the ALL stress dataset captured every single general stress gene as significant,
#but also had false positives for general stress

