#plot_stress_scatterplots.R

# PICK DATASET YOU WANT TO PLOT -------------------------------------------

#standard dataset
source('figurePlotting/load_stress_for_scatterplots.R')

#or

#the correlated stress projecst only set
source('figurePlotting/load_corrStress_for_scatterplots.R')


# LIBS --------------------------------------------------------------------
library(DESeq2)
library(tidyverse)
library(cowplot)
library(caret)
library(gtools)
source('figurePlotting/rna_functions.R')


# write out GO MWU input files --------------------------------------------
sgoOut = paste(datasetPrefix, 'stressForGO.csv', sep='_')
write_out_go(stress, paste('./go_mwu', sgoOut, sep='/'))


# GET CORE STRESS GENES ------------------------------------------------
#genes that are significant in the general stress set, and each specific one (with same log2fc direction)

allFoldChangeCut = 5/4 #all stress types must show this


#convert to log2 scales
allLcut = log(allFoldChangeCut, 2) 
mostLcut = log(mostFoldChangeCut, 2)
sigCut = 0.01

#get base of significant up and down in all set
baseSig = stress %>% 
  filter(padj<sigCut) %>% 
  pull(gene)
length(baseSig)

#ASSEMBLE LOG2 FOLD CHANGES FOR EACH STRESS TYPE
adat = data.frame(stress)
colnames(adat)=paste('allStress', colnames(adat), sep='_')
colnames(adat)[1]='gene'
specifics = list(heat, cold, salinity, immune, ph); names = c('heat', 'cold', 'salinity', 'immune', 'ph') #with bleached distributed
#specifics = list(bleach, heatNob, salinityNob, immune, ph); names = c('bleach', 'heat', 'salinity', 'immune', 'ph') #with bleached segregated

adat = stress %>% 
  select(gene, log2FoldChange) %>% 
  dplyr::rename('allStress'=log2FoldChange)
for (i in 1:length(specifics)){
  s=specifics[[i]]
  n=names[i]
  print(n)
  colnames(s)[1]='gene'
  colnames(s)[2]=n
  adat = adat %>% 
    full_join(s[,c('gene', n)], by = 'gene')
}

#ASSIGN CORE UPREGULATED
ldat = data.frame(adat, row.names='gene')
basePass = rownames(ldat) %in% baseSig
allUpSums = apply(ldat, 1, function(x) return(sum(x > allLcut)))
allDownSums = apply(ldat, 1, function(x) return(sum(x < -allLcut)))
nAll = ncol(ldat)


#core up
coreUp = ldat %>% 
  mutate(gene=rownames(ldat)) %>% 
  filter(allUpSums >= nAll & basePass)
nrow(coreUp)

#core down
coreDown = ldat %>% 
  mutate(gene=rownames(ldat)) %>% 
  filter(allDownSums >= nAll & basePass)
nrow(coreDown)

#final core gene assignment
coreGenes = append(coreUp$gene, coreDown$gene)
length(coreGenes)

# BUILD PCAS ----------------------------------------------------------------

#set plotting vars
fraction = 1/10
NTOP = round(fraction * nrow(vscold[[1]]), digits=0)


pcheat = plotStressPCA(df = vsheat[[1]], coldat = vsheat[[2]], intgroup = 'treat', ntop=NTOP, main = 'heat', xInvert=-1)
pcheatNob = plotStressPCA(df = vsheatNob[[1]], coldat = vsheatNob[[2]], intgroup = 'treat', ntop=NTOP, main = 'heat', xInvert=hNobInvert)
pccold = plotStressPCA(df = vscold[[1]], coldat = vscold[[2]], intgroup = 'treat', ntop=NTOP, main = 'cold', xInvert=1)
if (datasetPrefix=='standard'){
  pccoldNob = plotStressPCA(df = vscoldNob[[1]], coldat = vscoldNob[[2]], intgroup = 'treat', ntop=NTOP, main = 'cold', xInvert=1, pcs=4)
}
pcsalt = plotStressPCA(df = vssalt[[1]], coldat = vssalt[[2]], intgroup = 'treat', ntop=NTOP, main = 'hyposalinity', xInvert=-1)
pcsaltNob = plotStressPCA(df = vssaltNob[[1]], coldat = vssaltNob[[2]], intgroup = 'treat', ntop=NTOP, main = 'hyposalinity')
pcbleach = plotStressPCA(df = vsbleach[[1]], coldat = vsbleach[[2]], intgroup = 'bleached', ntop=NTOP, main = 'bleached')
pcimmune = plotStressPCA(df = vsimmune[[1]], coldat = vsimmune[[2]], intgroup = 'treat', ntop=NTOP, main = 'immune challenge', xInvert=-1)
pcph = plotStressPCA(df = vsph[[1]], coldat = vsph[[2]], intgroup = 'treat', ntop=NTOP, main = 'pH', xInvert=-1)
#handle unique case where points get cut off
if (datasetPrefix=="corrProjects"){
  pcphYlims=c(-12,10)
  pcph = pcph + lims(y=pcphYlims) #handle unique case where points get cut off
} else{
  pcph = plotStressPCA(df = vsph[[1]], coldat = vsph[[2]], intgroup = 'treat', ntop=NTOP, main = 'pH', xInvert=1)
}


# PLOT EACH PARTICULAR AGAINST GENERAL MINUS THAT STRESS -----------------------

#plot overlay plots
XLAB=bquote(log[2]~'diff. general')
YLAB=bquote(log[2]~'diff. specific')
XLAB=''
YLAB=''
hplt = do_scatter_overlay(dfy=heat, dfx=mheat, XLAB, 'heat', ' ', coreGenes)
cplt = do_scatter_overlay(cold, mcold, XLAB, 'cold', ' ', coreGenes)
splt = do_scatter_overlay(salinity, msalinity, XLAB, 'hyposal.', ' ', coreGenes)
iplt = do_scatter_overlay(immune, mimmune, XLAB, 'immune', ' ', coreGenes)
pplt = do_scatter_overlay(ph, mph, 'other stresses', 'pH', ' ', coreGenes, )
bplt = do_scatter_overlay(bleach, mbleach, XLAB, 'bleached', ' ', coreGenes)
hpltNob = do_scatter_overlay(heatNob, mheat, XLAB, 'heat', ' ', coreGenes)
spltNob = do_scatter_overlay(salinityNob, msalinity, XLAB, 'hyposal.', ' ', coreGenes)
if (datasetPrefix=='standard'){
  cpltNob = do_scatter_overlay(coldNob, mcold, XLAB, 'cold', ' ', coreGenes)
} #(doesn't exist for corrStress)


# PLOT THE STRES MINUS PLOTS ---------------------------------------------

#get the data for reference
genStressLab='other'
hdat = get_pred_dat(mh.ppath, hm.ppath, genStressLab, 'heat')
cdat = get_pred_dat(mc.ppath, cm.ppath, genStressLab, 'cold')
sdat = get_pred_dat(ms.ppath, sm.ppath, genStressLab, 'hyposal.')
idat = get_pred_dat(mi.ppath, im.ppath, genStressLab, 'immune')
pdat = get_pred_dat(mp.ppath, pm.ppath, genStressLab, 'pH')
all = rbind(hdat, cdat, sdat, idat, pdat)
all %>% 
  filter(Method=='RF',
         Train!='other')


#for main figure
hbar = plot_rf_acc_for_specifics(mh.ppath, hm.ppath, genStressLab, 'heat',' ') + labs(x=' ', y=' ')
cbar = plot_rf_acc_for_specifics(mc.ppath, cm.ppath, genStressLab, 'cold', ' ') + labs(x=' ', y=' ')
sbar = plot_rf_acc_for_specifics(ms.ppath, sm.ppath, genStressLab, 'hyposal.', ' ') + labs(x=' ', y=' ')
ibar = plot_rf_acc_for_specifics(mi.ppath, im.ppath, genStressLab, 'immune', ' ') + labs(x=' ', y=' ')
pbar = plot_rf_acc_for_specifics(mp.ppath, pm.ppath, genStressLab, 'pH', ' ') + labs(x='training set', y=' ')
#for supplement with separate beww
bbar = plot_rf_acc_for_specifics(mb.ppath, bm.ppath, genStressLab, 'bleached', ' ') + labs(x=' ', y=' ')
hbarNob = plot_rf_acc_for_specifics(mhNob.ppath, hmNob.ppath, genStressLab, 'heat', ' ') + labs(x=' ', y=' ')
sbarNob = plot_rf_acc_for_specifics(msNob.ppath, smNob.ppath, genStressLab, 'hyposal.', ' ') + labs(x=' ', y=' ')




#plot with bleaching distributed accross stressors
plot_grid(pcheat, hplt, hbar,
          pccold, cplt,cbar,
          pcsalt, splt,sbar,
          pcimmune, iplt,ibar,
          pcph+theme(legend.position='bottom',
                     legend.title=element_blank()) +
            scale_color_discrete( labels = c("control", "stressed")), 
          pplt+theme(legend.position='bottom'),
          pbar+theme(legend.position='bottom',
                     legend.title=element_blank()),
          nrow=5, 
          rel_heights=c(1,1,1,1,1.4),
          labels=LETTERS[1:15])

#get mean accuracy



#plot with bleaching segregated
plot_grid(pcbleach, bplt, bbar,
          pcheatNob, hpltNob, hbarNob,
          pcsaltNob, spltNob, sbarNob,
          pcimmune, iplt, ibar,
          pcph+theme(legend.position='bottom',
                     legend.title=element_blank()) +
            scale_color_discrete( labels = c("control", "stressed")), 
          pplt+theme(legend.position='bottom'),
          pbar+theme(legend.position='bottom',
                     legend.title=element_blank()),
          nrow=5, 
          rel_heights=c(1,1,1,1,1.4),
          labels=LETTERS[1:15])


# SCATTERPLOTS AGAINST ADULT V LARVA FOR REFERENCE ------------------------

avl = load_deseq('./deseqResults/adultVlarva_deseqResults.Rdata', 'adult_v_larva')
vdev = plot_overlay_volcano(avl, 'adult vs larva', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
hplt = do_scatter_overlay(heat, avl, NULL,NULL, 'heat vs age', coreGenes)
cplt = do_scatter_overlay(cold, avl, NULL, NULL, 'cold vs age', coreGenes)
splt = do_scatter_overlay(salinity, avl, NULL, NULL, 'salinity vs age', coreGenes)
iplt = do_scatter_overlay(immune, avl, NULL, NULL, 'immune vs age', coreGenes)
pplt = do_scatter_overlay(ph, avl, NULL, NULL, 'pH vs age', coreGenes)

plot_grid(vdev,
          hplt,
          cplt,
          splt,
          iplt,
          pplt,
          nrow=2,
          labels=LETTERS[1:6])


# BUILD VOLCANO PLOTS -----------------------------------------------------


vheat = plot_overlay_volcano(heat, 'heat', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vheatNob = plot_overlay_volcano(heatNob, 'heat', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vcold = plot_overlay_volcano(cold, 'cold', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vcoldNob = plot_overlay_volcano(coldNob, 'cold', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vsalinity = plot_overlay_volcano(salinity, 'hyposalinity', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vsalinityNob = plot_overlay_volcano(salinityNob[salinityNob$gene!='Amillepora11267',], 'hyposalinity', overlayGenes1 = coreGenes, overlayGenes2 = FALSE) #note remove single gene with extreme value to fit
vbleach = plot_overlay_volcano(bleach, 'bleached', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)
vimmune = plot_overlay_volcano(immune, 'immune', overlayGenes1 = coreGenes, overlayGenes2 =FALSE)
vph = plot_overlay_volcano(ph, 'pH', overlayGenes1 = coreGenes, overlayGenes2 = FALSE)

# BUILD SCATTERPLOTS OF EACH PARTICULAR AGAINST EACH -------------------------

#set up list of datasets and name it
datList = list(heat, cold, salinity, immune, ph)
nHeat = 'Ht'
nCold = 'Cold'
nSalinity = 'Hsal'
nImmune = 'Imm'
nph = 'pH'
names = c(nHeat, nCold, nSalinity, nImmune, nph)


#loop through and create scatterplots
pltList = list()
R2s = c()
for (i in 1:length(datList)){
  for (j in 1:length(datList)){
    df1=datList[[i]]
    df2=datList[[j]]
    name1=names[i]
    name2=names[j]
    compareString = paste(name1, name2, sep=' vs ')
    print(compareString)
    plt=do_scatter_overlay(df1, df2, NULL, NULL, ' ', coreGenes, writeR2=FALSE)
    pltList[[compareString]]=plt
    if(name1 != name2){
      r2=getLog2R2(df1, df2)
      R2s=append(R2s, r2)
    }
  }
}
names(pltList)

#replace self-comparisons with volcano plots
pltList[[paste(nHeat, nHeat, sep=' vs ')]] = vheat
pltList[[paste(nCold, nCold, sep=' vs ')]] = vcold
pltList[[paste(nSalinity, nSalinity, sep=' vs ')]] = vsalinity
pltList[[paste(nImmune, nImmune, sep=' vs ')]] = vimmune
pltList[[paste(nph, nph, sep=' vs ')]] = vph



#REPLACE BOTTOM TRIANGLE WITH THE PREDICTION ACCURACIES
#heat-cold
hc_acc = plot_rf_acc_for_specifics(hc.ppath, ch.ppath, nHeat, nCold)
pltList[[paste(nCold, nHeat, sep=' vs ')]] = hc_acc

#heat-salt
hs_acc = plot_rf_acc_for_specifics(hs.ppath, sh.ppath, nHeat, nSalinity)
pltList[[paste(nSalinity, nHeat, sep=' vs ')]] = hs_acc

#cold-salt
cs_acc = plot_rf_acc_for_specifics(cs.ppath, sc.ppath, nCold, nSalinity)
pltList[[paste(nSalinity, nCold, sep=' vs ')]] = cs_acc

#heat-immune
hi_acc = plot_rf_acc_for_specifics(hi.ppath, ih.ppath, nHeat, nImmune)
pltList[[paste(nImmune, nHeat, sep=' vs ')]] = hi_acc

#cold-immune
ci_acc = plot_rf_acc_for_specifics(ci.ppath, ic.ppath, nCold, nImmune)
pltList[[paste(nImmune, nCold, sep=' vs ')]] = ci_acc

#salt-immune
si_acc = plot_rf_acc_for_specifics(si.ppath, is.ppath, nSalinity, nImmune)
pltList[[paste(nImmune, nSalinity, sep=' vs ')]] = si_acc

#heat-ph
hp_acc = plot_rf_acc_for_specifics(hp.ppath, ph.ppath, nHeat, nph)
pltList[[paste(nph, nHeat, sep=' vs ')]] = hp_acc

#cold-ph
cp_acc = plot_rf_acc_for_specifics(cp.ppath, pc.ppath, nCold, nph)
pltList[[paste(nph, nCold, sep=' vs ')]] = cp_acc

#salinity-ph
sp_acc = plot_rf_acc_for_specifics(sp.ppath, ps.ppath, nSalinity, nph)
pltList[[paste(nph, nSalinity, sep=' vs ')]] = sp_acc

#immune-ph
ip_acc = plot_rf_acc_for_specifics(ip.ppath, pi.ppath, nImmune, nph)
pltList[[paste(nph, nImmune, sep=' vs ')]] = ip_acc



#BUILD THE PLOT
top=plot_grid(plotlist = pltList)
mean(R2s)
max(R2s)
min(R2s)

# REPEAT FOR NO BEWW SETS -------------------------

#set up list of datasets and name it
datList = list(bleach, heatNob, salinityNob, immune, ph)
nBleach = 'Blch'
nHeat = 'Ht'
nSalinity = 'Hsal'
nImmune = 'Imm'
nph = 'pH'
names = c(nBleach, nHeat, nSalinity, nImmune, nph)


#loop through and create scatterplots
pltList = list()
R2s = c()
for (i in 1:length(datList)){
  for (j in 1:length(datList)){
    df1=datList[[i]]
    df2=datList[[j]]
    name1=names[i]
    name2=names[j]
    compareString = paste(name1, name2, sep=' vs ')
    print(compareString)
    plt=do_scatter_overlay(df1, df2, NULL, NULL, ' ', coreGenes, writeR2=FALSE)
    pltList[[compareString]]=plt
    if(name1 != name2){
      r2=getLog2R2(df1, df2)
      R2s=append(R2s, r2)
    }
  }
}
names(pltList)

#replace self-comparisons with volcano plots
pltList[[paste(nBleach, nBleach, sep=' vs ')]] = vbleach
pltList[[paste(nHeat, nHeat, sep=' vs ')]] = vheatNob
pltList[[paste(nSalinity, nSalinity, sep=' vs ')]] = vsalinityNob
pltList[[paste(nImmune, nImmune, sep=' vs ')]] = vimmune
pltList[[paste(nph, nph, sep=' vs ')]] = vph


#REPLACE BOTTOM TRIANGLE WITH THE PREDICTION ACCURACIES
#bleach-heat
bh_accNob = plot_rf_acc_for_specifics(bhNob.ppath, hbNob.ppath, nBleach, nHeat)
pltList[[paste(nHeat, nBleach, sep=' vs ')]] = bh_accNob

#bleach-salt
bs_accNob = plot_rf_acc_for_specifics(bsNob.ppath, sbNob.ppath, nBleach, nSalinity)
pltList[[paste(nSalinity, nBleach, sep=' vs ')]] = bs_accNob

#heat-salt
hs_accNob = plot_rf_acc_for_specifics(hsNob.ppath, shNob.ppath, nHeat, nSalinity)
pltList[[paste(nSalinity, nHeat, sep=' vs ')]] = hs_accNob

#bleach-immune
bi_accNob = plot_rf_acc_for_specifics(biNob.ppath, ibNob.ppath, nBleach, nImmune)
pltList[[paste(nImmune, nBleach, sep=' vs ')]] = bi_accNob

#heat-immune
hi_accNob = plot_rf_acc_for_specifics(hiNob.ppath, ihNob.ppath, nHeat, nImmune)
pltList[[paste(nImmune, nHeat, sep=' vs ')]] = hi_accNob

#salt-immune
si_accNob = plot_rf_acc_for_specifics(siNob.ppath, isNob.ppath, nSalinity, nImmune)
pltList[[paste(nImmune, nSalinity, sep=' vs ')]] = si_accNob

#bleach-ph
bp_accNob = plot_rf_acc_for_specifics(bpNob.ppath, pbNob.ppath, nBleach, nph)
pltList[[paste(nph, nBleach, sep=' vs ')]] = bp_accNob

#heat-ph
hp_accNob = plot_rf_acc_for_specifics(hpNob.ppath, phNob.ppath, nHeat, nph)
pltList[[paste(nph, nHeat, sep=' vs ')]] = hp_accNob

#salinity-ph
sp_accNob = plot_rf_acc_for_specifics(spNob.ppath, psNob.ppath, nSalinity, nph)
pltList[[paste(nph, nSalinity, sep=' vs ')]] = sp_accNob

#immune-ph (this one can stay as it was since they didn't have any BEWW)
ip_accNob = plot_rf_acc_for_specifics(ip.ppath, pi.ppath, nImmune, nph)
pltList[[paste(nph, nImmune, sep=' vs ')]] = ip_accNob



#BUILD THE PLOT
plot_grid(plotlist = pltList)
mean(R2s)
max(R2s)
min(R2s)


#old junk below kept just in case



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

