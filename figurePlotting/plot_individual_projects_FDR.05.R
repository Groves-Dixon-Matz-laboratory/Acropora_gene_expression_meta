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


#read in the files for individual studies
idfList = list()
for (i in 1:length(iFileList)){
  x=iFileList[[i]]
  load(x)
  n=inames[i]
  idfList[[n]]=read_deseq_res(x,n)
}

######### FILER THE INDIVIDUAL PROJECTS BY P-VALUE
CUT=0.05
filter_by_fdr = function(dat){
  p = dat[,2]
  padj = p.adjust(p)
  fdat = dat %>% 
    filter(padj < CUT)
  print('------')
  print(paste('N significant =', sum(padj < CUT, na.rm=TRUE)))
  print(paste("N rows =", nrow(fdat)))
  return(fdat)
}
idfList2 = map(idfList, function(x) filter_by_fdr(x))

# cdfList = append(idfList, gdfList)
cdfList=idfList2


#merge them
cdat = cdfList %>% 
  purrr::reduce(full_join, by = 'gene')
rownames(cdat)=cdat$gene
ldat = cdat %>% 
  dplyr::select(grep('_lfc', colnames(cdat)))
colnames(ldat) = sub('_lfc', '', colnames(ldat))

#save write out the table
order = colnames(cdat)
order = order[order!='gene']
order = append(c('gene'), order)
idat = cdat[,order]


# ARRANAGE CORRELATION MATRIX ---------------------------------------------

#get correlations and modify names (note have to remove the datasets with no significant genes here)
noAllBeww = ldat %>% 
  dplyr::select(-j1_thisStudy_PRJNA559404,
                -A_moya_acid_PRJNA149513,
                -F_Uqueensland_ph_PRJNA269992,
                -J_weiss_immune_PRJNA200542,
                -f1_parkinson_hotCold_PRJNA423227)
c=cor(noAllBeww, use="pairwise.complete.obs")
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

BREAKS=seq(-.1, 0.5, length.out=length(COLOR)+1)
pheatmap(c,
         labels_col=projName,
         labels_row = treat,
         na_col=COLOR[100],
         color=COLOR,
         breaks=BREAKS,
         cluster_rows = h,
         cluster_cols = h)


#get mean cor for group A
bs = c('N_Hainan_U_heat_PRJNA308355', 'H_matz_heatTolLat_PRJNA279192')
a_cors = c[!rownames(c) %in% bs,
           !colnames(c) %in% bs]
a_cors[a_cors==1]<-NA
b_cors = c[rownames(c) %in% bs,
           colnames(c) %in% bs]
b_cors[b_cors==1]<-NA
mean(a_cors, na.rm=TRUE)
mean(b_cors, na.rm=TRUE)
