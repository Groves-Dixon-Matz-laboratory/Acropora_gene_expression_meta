#kog_MWU.R

source('figurePlotting/rna_functions.R')
library(KOGMWU)


#run single file 
inputFile = 'wgcna/go_mwu_module_style/MMred_goInput.csv'
inputFile = 'go_mwu/clusterAstress_For_MWU.csv'
kogFile = 'kog_mwu/Amil.v2.eggnogWebsite.gene2kog.tsv' 
inputData = read.csv(inputFile)
gene2kog = read.table(kogFile, sep='\t', header = TRUE)

M = kog.mwu(data=read.csv(inputFile), gene2kog=gene2kog, Alternative='t')
M
################# Generate heatmap
# compiling a table of delta-ranks to compare these results:
kognames=M$term
kogtable=data.frame("M" = M$delta.rank)
row.names(kogtable)=kognames

#sort the dataframe
y = rownames(kogtable)
x = kogtable$M
z = data.frame("M" = x, 'names' = y)
z = z[with(z, order(M)),]
kogtable = data.frame('M' = z$M)
row.names(kogtable) = z$names


# making a heatmap. Un-remark and run the next two lines if you don't have packages pheatmap and RColorBrewer installed. 
library(pheatmap)
library(RColorBrewer)

# you may want to adust these colors to make sure 0 is in the yellow; this will depend on your scale of delta-ranks
color = colorRampPalette(rev(c(brewer.pal(n = 7, name = "RdYlBu"),"darkblue","darkblue")))(100)

quartz()
pheatmap(as.matrix(kogtable), clustering_distance_rows = "euclidean", cluster_rows = F, cluster_cols = F, color=color, treeheight_row=15, treeheight_col=20)


#designate kogs to remove based on size
kogsToRemove = c('Function Unknown',
                 'Cell motility',
                 'Nuclear structure')



# RUN FOR MULTIPLE FILES --------------------------------------------------
ll=load('./go_mwu/groupedInputs.Rdata')
ll
ll=load('./go_mwu/individualInputs.Rdata')
ll
groupedNames
individualNames
names = groupedNames
names = individualNames
names = append(individualNames, c('clusterAstress', 'clusterBstress'))
inputs = paste(paste('./go_mwu', names,sep='/'), '_For_MWU.csv', sep='')

# #optionally add the modules
# modNames = c('red_module', 'white_module', 'green_module')
# modPaths = c('wgcna/go_mwu_module_style/MMred_goInput.csv',
#              'wgcna/go_mwu_module_style/MMwhite_goInput.csv',
#              'wgcna/go_mwu_module_style/MMgreen4_goInput.csv')
# names = append(names, modNames)
# inputs = append(inputs, modPaths)


kogResults = list()
for (i in 1:length(inputs)){
  inputFile=inputs[i]
  inputData = read.csv(inputFile)
  name = names[i]
  print(paste(name,'...',sep=''))
  cname= paste(name, 'DR', sep='_')
  M=kog.mwu(data = inputData, gene2kog = gene2kog, Alternative='t')
  M[,cname]=M$delta.rank
  M$term=as.character(M$term)
  kogResults[[cname]]=M
}

kdat = kogResults %>% 
  purrr::reduce(full_join, by='term') %>% 
  select('term', ends_with("_DR"))
colnames(kdat)=sub('_DR', '', colnames(kdat))
colnames(kdat)
rownames(kdat)=kdat$term

#format matrix of delta ranks for building heatmap
kmat = kdat
kmat = kmat[!kmat$term %in% kogsToRemove,]
kmat$term<-NULL
kmat = as.matrix(kmdat)

#order by median
meds = apply(kmat, 1, median)
orderedNames = names(meds)[order(meds)]

#order by A
obyA = kmat %>% 
  data.frame() %>% 
  mutate(kog=rownames(kmat)) %>% 
  arrange(clusterAstress) %>% 
  pull(kog)
orderedNames=obyA


kmat=kmat[orderedNames,]
plot(density(kmat))
BREAKS=seq(-1500, 1500, length.out=length(color)+1)
kmat.ordered = kmat[orderedNames,]
pheatmap(kmat.ordered,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "maximum",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color=color,
         treeheight_row=15,
         treeheight_col=20,
         breaks = BREAKS
         )




#plot to match other heatmap
ll=load('kog_mwu/heatmapObjects.Rdata')
ll
c=kmat[orderedNames,orderedProjects]
bleachAnnot = data.frame(bleached = bleachAnnot[orderedProjects,],
                         cluster = anots[orderedProjects,],
                         row.names=orderedProjects)
head(bleachAnnot)
dim(c)
m=c

#set modified project names
orderedProjects %in% colnames(kmat)
modNames2 = modNames[orderedProjects,]
projName = paste(modNames2$Bioproject, modNames2$ref)


BREAKS=seq(-1200, 1200, length.out=length(color)+1)
pheatmap(c,
         labels_col=projName,
         color=color,
         breaks=BREAKS,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_number = 10,
         annotation_col = bleachAnnot,
         annotation_colors = my_colour,
         treeheight_row=25,
         treeheight_col=0)



#subset
ll=load('metadata/corStressProjs.Rdata')
ll
corStressProjs
keep=corStressProjs
keep = append(keep, names[grep('beww', names)])
ksub=kmat[,keep]
BREAKS=seq(-1000, 1000, length.out=length(color)+1)
meds = apply(ksub, 1, median)
orderedNames = names(meds)[order(meds)]
ksub.ordered = ksub[orderedNames,]
pheatmap(ksub.ordered, 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "maximum",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color=color,
         treeheight_row=15,
         treeheight_col=20,
         breaks=BREAKS
         # labels_row=kdat$term
)



# RUN FOR YEAST DATA ------------------------------------------------------
inputs0 = list.files('kog_mwu/yeast_inputs', '*_ForMWU.csv', full.names = TRUE)
inputs = inputs0[!grepl('MWU_', inputs0) &
                   !grepl('.tmp', inputs0) &
                   !grepl('MF_', inputs0) &
                   !grepl('CC_', inputs0) &
                   !grepl('BP_', inputs0)]
names0 = sub('_ForMWU.csv', '', inputs)
names=sub('kog_mwu/yeast_inputs/', '', names0, fixed=TRUE)
ykogFile = 'kog_mwu/Lachancea_kluyveri_cds.fasta.emapperWebsite.gene2kog.tsv'
ygene2kog = read.table(ykogFile, sep='\t', header = TRUE)

kogResults = list()
for (i in 1:length(inputs)){
  inputFile=inputs[i]
  inputData = read.csv(inputFile)
  name = names[i]
  cname= paste(name, 'DR', sep='_')
  M=kog.mwu(data = inputData, gene2kog = ygene2kog, Alternative='t')
  M[,cname]=M$delta.rank
  M$term=as.character(M$term)
  kogResults[[cname]]=M
}

y.kdat = kogResults %>% 
  purrr::reduce(full_join, by='term') %>% 
  select('term', ends_with("_DR"))


y.kdat %>% 
  select(grep('_DR', colnames(y.kdat)))

colnames(y.kdat)=sub('_DR', '', colnames(y.kdat))
colnames(y.kdat)

y.kmat = y.kdat %>% 
  select(-'term') %>% 
  as.matrix()
rownames(y.kmat)=y.kdat$term
meds = apply(y.kmat, 1, median)
orderedNames = names(meds)[order(meds)]
y.kmat=y.kmat[orderedNames,]

#build heatmap
library(pheatmap)
library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name = "RdYlBu"),"darkblue","darkblue")))(100)
BREAKS=seq(-1500, 1500, length.out=length(color)+1)
pheatmap(y.kmat, clustering_distance_rows = "euclidean",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color=color,
         treeheight_row=15,
         treeheight_col=20,
         breaks=BREAKS
         # labels_row=kdat$term
)



# plot together -----------------------------------------------------------

head(y.kmat)
head(kmat)

mat = merge(kmat, y.kmat, by = 0)
rownames(mat)=mat$Row.names
mat$Row.names<-NULL
dim(mat)
meds = apply(mat, 1, median)
orderedNames = names(meds)[order(meds)]
mat=mat[orderedNames,]
BREAKS=seq(-1000, 1000, length.out=length(color)+1)
pheatmap(mat,
         clustering_distance_rows = "correlation",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color=color,
         treeheight_row=15,
         treeheight_col=20,
         breaks=BREAKS
         # labels_row=kdat$term
)


adat = read_tsv('metadata/Amillepora_euk.emapper.annotations.tsv')

adat %>% 
  filter(`COG cat`=='Y') %>% 
  write_csv(path='~/Desktop/nucStruc.csv')

# plot a selection --------------------------------------------------------

#order by A
obyA = kmat %>% 
  data.frame() %>% 
  mutate(kog=rownames(kmat)) %>% 
  arrange(clusterAstress) %>% 
  pull(kog)
orderedNames=obyA


colnames(mat)
selection = c('All_yeast',
              'H_matz_heatTolLat_PRJNA279192',
              'clusterAstress',
              'clusterBstress',
              'bewwHeat')
selection = c('All_yeast',
              'clusterAstress',
              'clusterBstress')
smat=mat[orderedNames,selection]

BREAKS=seq(-1000, 1000, length.out=length(color)+1)
smat=smat[,c("clusterAstress", "clusterBstress", "All_yeast")]
COLLABS = c('type A', 'type B', 'yeast')
pheatmap(smat,
         clustering_distance_cols = "correlation",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color=color,
         treeheight_row=15,
         treeheight_col=20,
         labels_col = COLLABS,
         breaks=BREAKS,
         fontsize=12
         
)

plotdf = data.frame(smat)
colnames(plotdf) = c('typeA', 'typeB', 'yeast')


ab = plot_scatter(plotdf, xcol='typeA', ycol='typeB',  addStars=TRUE) + labs(x='type A', y='type B')
lmab = lm(typeB~typeA, data=plotdf)
summary(lmab)
ay = plot_scatter(plotdf, xcol='typeA', ycol='yeast') + labs(x='type A', y='yeast')
lmay = lm(typeA~yeast, data=plotdf)
summary(lmay)
by = plot_scatter(plotdf, xcol='typeB', ycol='yeast') + labs(x='type B', y='yeast')
lmby = lm(typeB~yeast, data=plotdf)
summary(lmby)
pltList = list(ab, ay, by)
BREAKS=c(-1500, 0, 1500)
LIMS=c(-2100, 2100)

pltList2 = lapply(pltList, function(x) return(x +
                                                scale_x_continuous(breaks=BREAKS, limits=LIMS) +
                                                scale_y_continuous(breaks = BREAKS, limits=LIMS) +
                                                geom_smooth(method='lm', se=FALSE)))
plot_grid(plotlist=pltList2, nrow=3)




dat = smat
pairs(smat)
dat %>% 
  gather(key='study', value='deltaRank', -All_yeast) %>% 
  ggplot(aes(x=All_yeast, y=deltaRank, group=study, color=study)) +
  geom_point()

pairs(smat)


# RUN FOR INDIVIDUAL PROJECTS ---------------------------------------------

inputs = list.files('go_mwu', '*_IndStudy_ForMWU.csv', full.names = TRUE)
names = sub('_IndStudy_ForMWU.csv', '', list.files('go_mwu', '*_IndStudy_ForMWU.csv'))
kogFile = 'kog_mwu/Amillepora.v2.00.chrs_gene2kogClass.tab' 
kogResults = list()
for (i in 1:length(inputs)){
  inputFile=inputs[i]
  name = names[i]
  cname= paste(name, 'DR', sep='_')
  M=kog.mwu(inname = inputFile, gene2kog = kogFile, Alternative='two.sided')
  M[,cname]=M$delta.rank
  M$term=as.character(M$term)
  kogResults[[cname]]=M
}

kdat = kogResults %>% 
  purrr::reduce(full_join, by='term') %>% 
  select('term', ends_with("_DR"))


kdat %>% 
  select(grep('_DR', colnames(kdat)))

colnames(kdat)=sub('_DR', '', colnames(kdat))
colnames(kdat)

kmat = kdat %>% 
  select(-'term') %>% 
  as.matrix()
rownames(kmat)=kdat$term
meds = apply(kmat, 1, median)
orderedNames = names(meds)[order(meds)]
kmat=kmat[orderedNames,]

#build heatmap
library(pheatmap)
library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name = "RdYlBu"),"darkblue","darkblue")))(100)
pheatmap(kmat, clustering_distance_rows = "euclidean",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color=color,
         treeheight_row=15,
         treeheight_col=20,
         # labels_row=kdat$term
)


#subset
keep=c('adultVlarva',
       'stress',
       'heat',
       'cold',
       'salinity',
       'immune',
       'ph'
)
ksub=kmat[,keep]
pheatmap(ksub, clustering_distance_rows = "euclidean",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color=color,
         treeheight_row=15,
         treeheight_col=20,
         # labels_row=kdat$term
)


# plot together yeast and individual -----------------------------------------------------------

head(y.kmat)
head(kmat)

mat = merge(kmat, y.kmat, by = 0)
rownames(mat)=mat$Row.names
mat$Row.names<-NULL
dim(mat)
meds = apply(mat, 1, median)
orderedNames = names(meds)[order(meds)]
mat=mat[orderedNames,]
pheatmap(mat, clustering_distance_rows = "euclidean",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color=color,
         treeheight_row=15,
         treeheight_col=20,
         # labels_row=kdat$term
)



selectCoral = colnames(kmat)


smat = mat[,selectCoral]
ymat = mat[,colnames(y.kmat)]
ymeds = apply(ymat, 1, median)



smat$yeast = ymeds

dat = smat
dat %>% 
  gather(key='stress', value='deltaRank', selectCoral) %>% 
  filter(stress=='heat_NoBEWW') %>% 
  ggplot(aes(x=yeast, y=deltaRank, group=stress, color=stress)) +
  geom_point()



head(kmat)


lm1=lm(mat$H_matz_heatTolLat_PRJNA279192~mat$All_yeast)
summary(lm1)
p1=mat %>% 
  ggplot(aes(x=H_matz_heatTolLat_PRJNA279192, y=All_yeast)) +
  geom_point() +
  geom_smooth(method='lm')


lm1=lm(mat$j1_BLEACH_EWW~mat$All_yeast)
summary(lm1)
p2=mat %>% 
  ggplot(aes(x=j1_BLEACH_EWW, y=All_yeast)) +
  geom_point()+
  geom_smooth(method='lm')



lm1=lm(mat$j1_BLEACH_EWW~mat$H_matz_heatTolLat_PRJNA279192)
summary(lm1)
p3=mat %>% 
  ggplot(aes(x=j1_BLEACH_EWW, y=H_matz_heatTolLat_PRJNA279192)) +
  geom_point()+
  geom_smooth(method='lm')


plot_grid(p1,p2,p3, nrow=1)
