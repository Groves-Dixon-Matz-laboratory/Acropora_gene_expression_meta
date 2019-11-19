#plot_all_stress.R


library(tidyverse)
library(cowplot)
library(ggplot2)
library(DESeq2)
source('./figurePlotting/rna_functions.R')


# load the normalized counts ----------------------------------------------

input = 'largeIgnored/stress_project_controlled.Rdata'
ll=load(input)
ll
rld.df = t(datExpr)


#upload additional coldata with bleaching info
coldata = read_csv('metadata/subset_tables/allStress_Coldata.csv') %>% 
  mutate(bleached = if_else(is.na(bleached),
                            'no',
                            bleached))
coldata$my_letter = sapply(coldata$my_title, function(x) strsplit(x, '_')[[1]][1])
rld.df=rld.df[,coldata$Run]
dim(rld.df)
dim(coldata)

# modify treat names and get summary table -------------------------------------------------------

#modify treatment names
coldata = coldata %>% 
  mutate(Treatment = if_else(treat=='cold_then_heat',
                             'multiple',
                             treat),
         Treatment = if_else(treat=='heat_then_cold',
                             'multiple',
                             Treatment),
         Treatment = if_else(Treatment=='challenge',
                             'immune challenge',
                             Treatment),
         Treatment = if_else(Treatment=='low_pH',
                             'low pH',
                             Treatment),
         Treatment = if_else(Treatment=='low_salinity',
                             'hyposalinity',
                             Treatment),
         projType = if_else(projType=='pH+heat',
                            'multiple',
                            projType))

#make ordered factors for plotting downstream
unique(coldata$Treatment)
treatLevels = c('control', 'heat', 'cold', 'hyposalinity', 'immune challenge', 'low pH', 'multiple')
coldata$Treatment = factor(coldata$Treatment, levels=treatLevels)

#identify studies with too few reps
coldata %>% 
  filter(stress=='stressed') %>% 
  group_by(my_title) %>% 
  summarize(N=n())


#upload additional coldata with bleaching info
coldata2 = read_csv('metadata/subset_tables/allStress_Coldata.csv') %>% 
  mutate(bleached = if_else(is.na(bleached),
                            'no',
                            bleached))

coldata2 %>% 
  group_by(treat) %>% 
  summarize(N=n())


#get summary table
sumTab=coldata %>% 
  group_by(Treatment) %>% 
  summarize(`N projects` = length(unique(my_title)),
            `N samples` = n(),
            `N bleached` = sum(bleached=='yes')) %>% 
  mutate(`N bleached` = ifelse(Treatment=='control',
                                0,
                                `N bleached`))

sumTab
apply(sumTab[2:4], 2, sum)


# plot volcano ------------------------------------------------------------

ll=load('deseqResults/stress_deseqResults.Rdata')
volc=plot_volcano(data.frame(res), TITLE=NULL) 

# plot PCA ----------------------------------------------------------------

fraction = 1/10
NTOP = round(fraction * nrow(rld.df), digits=0)
point.size=1
treat.df = plotStressPCA(df = rld.df, coldat = coldata, intgroup = 'stress', ntop=NTOP, main = NULL, SIZE=point.size, returnData = T) #with modified one
treatPCA = plotStressPCA(df = rld.df, coldat = coldata, intgroup = 'stress', ntop=NTOP, main = NULL, SIZE=point.size, returnData = F) #with modified one


# RUN DAPC --------------------------------------------------------------------

######## PREPARE DATA FOR DAPC ########
library(adegenet)
pcp=prcomp(t(rld.df), retx=TRUE, center=TRUE, scale=TRUE)
scores=pcp$x
screeplot(pcp,bstick=T) # only 1st PC is non-random when using diff exp genes only; higher PCs non-random when using all data...

# adegenet: finding clusters (even though we know what clusters we want) - choose 4 PCs and 2 groups
clus=find.clusters(t(rld.df),max.n.clus=15, n.clust=2, n.pca=10) #[degs10,]

#Use clus$grp to rename
clus$grp=coldata$stress #set the transplantation site as the groups
cdf = data.frame(clus$grp)
head(cdf)
cdf$Run = coldata$Run

#run dapc using the assigned groups
dp=dapc(t(rld.df),clus$grp, n.da=1, perc.pca=80)
dpf = data.frame(dp$ind.coord)
dpf$Run = rownames(dpf)


#output the gene loading values for dapc
varCont = dp$var.contr
varCont %>% 
  data.frame() %>% 
  ggplot(aes(x=LD1)) +
  geom_density()
goOut = data.frame(gene=rownames(varCont),
                   contribution = varCont[,'LD1'])
head(goOut)
goOut %>% 
  write_csv(path='./go_mwu/stressDAPC_ForMWU.csv')


#output the variable importance for PC1 and DAPC as single dataframe
head(varCont)
pca.dapc.varImp = varCont

pc1Load = dp$pca.loadings %>% 
  data.frame() %>% 
  dplyr::select(1)

pca.dapc.varImp = merge(pca.dapc.varImp, pc1Load, by = 0)
colnames(pca.dapc.varImp)[1]='gene'
head(pca.dapc.varImp)
save(pca.dapc.varImp, file='results_tables/pcaDapcLoads.Rdata')


###### CORRELATE DAPC AND PCA
head(coldata)
treat.df$Run = rownames(treat.df)

#load high and low stress from heatmaps in plot_individual_projects.R
ll=load('metadata/corStressProjs.Rdata')
ll

#merge up the PCA and DAPC data
pdat = dpf %>% 
  full_join(treat.df, by = 'Run') %>% 
  full_join(coldata, by = 'Run') %>% 
  mutate(my_title=if_else(grepl("^j1_", my_title),
                          "j1_thisStudy_PRJNA559404",
                          my_title),
         slvl = if_else(my_title %in% corStressProjs & stress.x=='stressed',
                        'high',
                        'control'),
         slvl = if_else(my_title %in% lowStressProjs & stress.x=='stressed',
                        'low',
                        slvl),
         slvl = factor(slvl, levels=c('control', 'low', 'high')))
save(pdat, file='figurePlotting/stress_pca_LD.Rdata')


# plot stress prediction results -------------------------------------------
ll=load('./category_prediction/stratifiedRandom/stratified_allStress_train.csv__stratified_allStress_test.csv_predictionResults.Rdata')
ll
library(caret)
head(cdat)
head(lres)
head(rfres)

#look at category prediction models
ll=load('./category_prediction/stratifiedRandom/stratified_allStress_train.csv__stratified_allStress_test.csv_predictionResults.Rdata')
cdat$gene=rownames(cdat)
var.imp$gene=rownames(var.imp)
predDat = load_pred_stats('./category_prediction/stratifiedRandom/stratified_allStress_train.csv__stratified_allStress_test.csv_predictionResults.Rdata')
predBp = plot_pred_stats(predDat, NULL) + theme(legend.position='top')


# final replots -----------------------------------------------------------

#set up table object
library(gridExtra)
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(sumTab, rows=NULL, theme=ttheme_minimal())

#volcano
volc=plot_volcano(data.frame(res), TITLE=NULL) + theme(legend.position='top') + labs(color='FDR<0.1', subtitle="\n") 

#pca
mcoldata=coldata %>% 
  mutate(stress = if_else(stress=='control','C','S'))
treatPCA = plotStressPCA(df = rld.df, coldat = mcoldata, intgroup = 'stress', ntop=NTOP, main = NULL, SIZE=point.size, returnData = F, xInvert=1) +theme(legend.position='right')

#pca density plots
pdens = treat.df %>% 
  mutate(stress = if_else(stress=='control',
                            'C',
                            'S')) %>% 
  ggplot(aes(x=PC1*-1, fill=stress)) +
  geom_density(alpha=0.5) +
  theme(legend.position='top') +
  labs(fill=NULL, x='PC1')


#density plots for high stress and low stress on PC1
pcDens = pdat %>% 
  ggplot(aes(x=PC1*-1, fill=slvl)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(breaks=c(-30,0,30)) +
  labs(fill='stress\nlevel',
       subtitle='\n',
       x='PC1')

#dapc density
ddens = pdat %>% 
  mutate(stress = if_else(stress.x=='control',
                          'C',
                          'S')) %>% 
  ggplot(aes(x=LD1, fill=stress)) +
  geom_density(alpha=0.5) +
  labs(fill=NULL,
       x='DAPC') +
  theme(legend.position='top')



#correlation PCA and DAPC by treat type
pdCor = pdat %>% 
  mutate(ttype = if_else(treat=='challenge',
                         'immune challenge',
                         treat),
         ttype = if_else(treat %in% c('cold_then_heat', 'heat_then_cold'),
                         'multiple',
                         ttype),
         ttype = if_else(treat == 'low_pH',
                         'low pH', 
                         ttype),
         ttype = if_else(treat=='low_salinity',
                         'hyposalinity',
                         ttype),
         ttype = factor(ttype, levels=c('control',
                                        'heat',
                                        'cold',
                                        'hyposalinity',
                                        'immune challenge',
                                        'low pH',
                                        'multiple')),
         blch = if_else(bleached=='yes',
                        'bleached',
                        'not bleached'),
         blch = if_else(ttype=='control',
                        'not bleached',
                        blch),
         blch=factor(blch, levels=c('not bleached', 'bleached'))) %>% 
  ggplot(aes(x=PC1, y=LD1)) +
  geom_smooth(span=0.3, se=FALSE, lwd=1) +
  geom_point(aes(color=ttype,
                 shape=blch), size=1) +
  labs(x='PC1',
       y='DAPC',
       color='',
       shape='') +
  theme(legend.title=element_blank())


#correlation with label by stress level
pdCorStsLvl = pdat %>% 
  ggplot(aes(x=PC1, y=LD1)) +
  geom_smooth(span=0.3, se=FALSE, lwd=1) +
  geom_point(aes(color=slvl), size=1) +
  labs(x='PC1',
       y='DAPC',
       color='stress\nlevel')



#new one with prediction barplot (don't like this as much)
quartz()
#main version of figure 2
plot_grid(tbl, volc+theme(legend.position='none'), treatPCA,  ddens, pdCor, predBp, nrow=3, rel_widths = c(1, 0.5), labels=LETTERS[1:6], axis='b')


#alternative versions:

#new version of figure 2
row1 = plot_grid(tbl, volc+theme(legend.position='none'), nrow=1, labels=LETTERS[1:2], rel_widths=c(1, 0.5))
row2 = plot_grid(treatPCA+theme(legend.position='none'), ddens, nrow=1, rel_widths=c(1, 0.75), labels=LETTERS[3:4])
row3 = plot_grid(pcDens+theme(legend.position='none'), pdCorStsLvl, nrow=1, labels=LETTERS[5:6], rel_widths=c(0.75, 1))
plot_grid(row1, row2, row3, nrow=3)


#simpler one that doesn't bother with DAPC
row1 = plot_grid(tbl, volc+theme(legend.position='none'), nrow=1, labels=LETTERS[1:2], rel_widths=c(1, 0.5))
row2 = plot_grid(treatPCA+theme(legend.position='top'), pcDens, nrow=1, rel_widths=c(1, 1), labels=LETTERS[3:4])
plot_grid(row1, row2, nrow=2)



#simpler one that doesn't bother with DAPC
row1 = plot_grid(tbl, volc+theme(legend.position='none'), nrow=1, labels=LETTERS[1:2], rel_widths=c(1, 0.5))
row2 = plot_grid(treatPCA+theme(legend.position='top'), pcDens, nrow=1, rel_widths=c(1, 1), labels=LETTERS[3:4])
plot_grid(row1, row2, nrow=2)



# plot stress type densities along PC1 and DA1 ----------------------------

#pca density plots
pdens2 = treat.df %>% 
  left_join(coldata) %>% 
  ggplot(aes(x=PC1*-1, fill=Treatment)) +
  geom_density(alpha=0.75) +
  theme(legend.position='none') +
  labs(fill=NULL, x='PC1')

#dapc density
ddens2 = pdat %>% 
  ggplot(aes(x=LD1, fill=Treatment)) +
  geom_density(alpha=0.75) +
  labs(fill=NULL, x='DAPC')


plot_grid(pdens2, ddens2, rel_widths=c(.75, 1))



# plot non-zero coefficients from lasso regression and variable importance from RF ------------------

#NON-ZERO COEFFICIENTS FROM LASSO LOGISTIC REG.
cdat$gene = rownames(cdat)
cdat %>% 
  mutate(annot=convert_to_annot(gene)) %>% 
  filter(abs(coef)>0) %>% 
  arrange(coef) %>% 
  filter(!grepl('(Intercept)', gene)) %>% 
  mutate(geneA=if_else(annot=='none',
                      gene,
                      annot),
         geneA=factor(geneA, levels=geneA)) %>% 
  ggplot() +
  geom_bar(stat='identity', aes(x=geneA, y=coef)) +
  coord_flip() +
  labs(y='Coefficent', x='Gene')


#double-check this against deseq results
ll=load('deseqResults/stress_deseqResults.Rdata')
ll
res=data.frame(res) %>% 
  mutate(gene=rownames(res))
crdat = cdat %>% 
  filter(gene != '(Intercept)',
         abs(coef)>0) %>% 
  left_join(res, by = 'gene')

#gene naames
ndat = convert_to_annot(crdat$gene)
labs = paste(crdat$gene, ndat, sep='_')


#plot
crdat %>% 
  ggplot(aes(x=log2FoldChange,y=coef)) +
  geom_point() +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  geom_text(aes(label=ndat), hjust=0)


#look at gene names
acdat = cdat %>% 
  mutate(annot=convert_to_annot(gene)) %>% 
  filter(abs(coef)>0,
         !grepl('(Intercept)', gene)) %>% 
  arrange(coef)

sum(acdat$annot=='none') / nrow(acdat)


#VARIABLE IMPORTANCE FOR RF
select = 99.7 #percentile to use as cutoff
select = 97.1
var.imp=var.imp %>% 
  arrange(desc(MeanDecreaseAccuracy)) %>% 
  mutate(rank=seq(nrow(var.imp), 1, -1))
selectStr = paste(as.character(select), '%', sep='')
pctile = quantile(var.imp$MeanDecreaseAccuracy, probs=seq(0,1, by = 0.001))[selectStr]
fullplt=var.imp %>% 
  ggplot(aes(y=rank,x=MeanDecreaseAccuracy)) +
  geom_vline(xintercept = pctile, lty=2) +
  geom_point() +
  labs(y='Gene rank', x='Mean Decrease Accuracy')

sel.imp = var.imp %>% 
  filter(MeanDecreaseAccuracy>=pctile) %>% 
  mutate(labs=convert_to_annot(gene),
         labs=if_else(labs=='none',
                      gene,
                      labs),
         labLengths = sapply(labs, function(x) nchar(x)),
         splitLabs = sapply(labs, function(x) strsplit(x, "(", fixed=TRUE)[[1]][1]),
         labs = if_else(labLengths > 40,
                        paste(substr(splitLabs, start=1, stop=30), '...', sep=''),
                        labs))
#replace weird annotaiton
sel.imp$labs[sel.imp$labs==')-reductase']<-'reductase'

selplt = sel.imp %>% 
  ggplot(aes(y=rank,x=MeanDecreaseAccuracy)) +
  geom_point() +
  labs(x='Mean Decrease Accuracy', y='Gene') +
  scale_y_continuous(breaks = sel.imp$rank, labels=sel.imp$labs)


#draw with inset
inset = fullplt + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ggdraw(selplt) +
  draw_plot(inset, .55, .12, .4, .4)


#### SAVE THE NON-ZERO COEFFS FROM LASSO AND VAR.IMP FROM RANDOM FOREST
save(cdat, var.imp, file='results_tables/classificationImportance.Rdata')


# check dapc against stress intensity -------------------------------------

#load detailed heat trait table
heat = read_csv('./metadata/detailed_tables/detailedHeat.csv') %>% 
  left_join(pdat, by='Run') %>% 
  select(my_title.x, Run, LD1, PC1, temp, degOverAmb, hoursExposed, treat.x)


heat %>% 
  filter(temp != 'not_clear') %>%
  mutate(temp=as.numeric(temp)) %>% 
  ggplot(aes(x=temp, y=LD1, color=my_title.x)) +
  geom_point()


heat %>% 
  ggplot(aes(x=degOverAmb, y=LD1, color=my_title.x)) +
  geom_point()

heat %>% 
  mutate(cumExpose = degOverAmb*hoursExposed) %>% 
  ggplot(aes(x=cumExpose, y=LD1, color=my_title.x)) +
  geom_point()

