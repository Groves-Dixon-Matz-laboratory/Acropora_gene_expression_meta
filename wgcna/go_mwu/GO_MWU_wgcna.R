
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################

#swtich to the WGCNA go mwu working dir and assign file names
setwd("./wgcna/go_mwu")
goDatabase="go.obo"
goAnnotations="amil_zachFullerV2_gos.tsv"
source("gomwu.functions.R")

#set variables to run for each module
ll=load('moduleInputFiles.Rdata')
ll
inputFiles
divisions=c('CC')
divisions=c('BP')

c('input1.csv', 'input2.csv')

# run go_mwu for each -----------------------------------------------------


for (goDivision in divisions){
  print('==============')
  print('==============')
  print('==============')
  print(goDivision)
  #set BP smallest to 50 so you don't get too many
  if (goDivision=='BP'){
    SMALLEST=50
  } else {
    SMALLEST=10
  }
  for (input in inputFiles){
    print('--------------')
    print(paste(input, '...', sep='.'))
    # Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
    gomwuStats(input, goDatabase, goAnnotations, goDivision,
               perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
               largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
               smallest=SMALLEST,   # a GO category should contain at least this many genes to be considered
               clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
               Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
               # Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
               #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
    )
  }
}




# record how many significant and save thsoe ---------------------------------------------
library(tidyverse)
divRec = c()
inRec = c()
sigRec = c()
for (goDivision in divisions){
  for (input in inputFiles){
    resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
    go.res=read.table(resName, header = T)
    totSig = sum(go.res$p.adj < 0.1)
    divRec = append(divRec, goDivision)
    inRec = append(inRec, input)
    sigRec = append(sigRec, totSig)
    sig=go.res[go.res$p.adj < 0.1,]
    sigOut=paste( c('./resultsFiles/', goDivision, input, '_sigGos.tsv'), collapse='')
    if (nrow(sig)>0){
      sig %>% 
        write_tsv(path=sigOut)
    }
  }
}
res = tibble('goDivision'=divRec,
                 'input'=inRec,
                 'nSig'=sigRec)
res %>% 
  write_tsv(path='./resultsFiles/gomwu_results_summary.tsv')


# OUTPUT SUPPLEMENTARY TABLES -------------------------------------------
allDf = data.frame()
forSuppTables = inputFiles[inputFiles!='orange_moduleInput.csv']
for (input in forSuppTables){
  inCol = strsplit(input, '_')[[1]][1]
  goDivisions=c('BP', 'MF', 'CC')
  go.res = data.frame()
  for (goDivision in goDivisions){
    print(goDivision)
    resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
    divDf = read.table(resName, header = T) %>% 
      filter(p.adj<0.1)
    divDf$goDivision=goDivision
    go.res=rbind(go.res, divDf)
  }
  go.res$module=inCol
  allDf = rbind(allDf, go.res)
}
allDf %>% 
  write_tsv(path='../../results_tables/module_go_mwu_results.tsv')



# PLOT A SINGLE DATASET ---------------------------------------------------
input='red_moduleInput.csv';COL='red';CUT=0.005;goDivision='BP'
input='red_moduleInput.csv';COL='red';CUT=0.01;goDivision='MF'
input='red_moduleInput.csv';COL='red';CUT=0.05;goDivision='CC'

input='white_moduleInput.csv';COL='black';CUT=0.0000002;goDivision='BP'
input='white_moduleInput.csv';COL='black';CUT=0.05;goDivision='MF'
input='white_moduleInput.csv';COL='black';CUT=0.01;goDivision='CC'

input='green4_moduleInput.csv';COL='green4';CUT=0.05;goDivision='MF'
input='green4_moduleInput.csv';COL='green4';CUT=0.01;goDivision='CC'
input='green4_moduleInput.csv';COL='green4';CUT=0.01;goDivision='BP'

quartz()
gomwuPlot(input,goAnnotations,goDivision,
          absValue=0.00001,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
          level1=CUT, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
          level2=CUT, # FDR cutoff to print in regular (not italic) font.
          level3=CUT, # FDR cutoff to print in large bold font.
          txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
          treeHeight=0.5, # height of the hierarchical clustering tree
          colors=c(COL,COL,COL,COL) # these are default colors, un-remar and change if needed
)


# plot results for each module with any significant enrichment -------------------
sig = res %>% 
  filter(nSig>=2)

for (i in 1:nrow(sig)){
  row=sig[i,]
  goDivision=row['goDivision']
  input = row['input']
  figFileName = paste('./resultsFiles/', sep='', paste(paste(goDivision, input, sep='_'), 'tree.pdf', sep='_'))
  pdf(figFileName)
  gomwuPlot(input,goAnnotations,goDivision,
            absValue=0.00001,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
            level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
            level2=0.05, # FDR cutoff to print in regular (not italic) font.
            level3=0.001, # FDR cutoff to print in large bold font.
            txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
            treeHeight=0.5, # height of the hierarchical clustering tree
            #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
  )
  dev.off()
}
    

#Code below is included in separate script look_at_genes_within_GO_terms.R

#------------ LOOK AT GENES WITHIN GO TERMS ------------#
library(cowplot)
library(tidyverse)
#upload the results for each genes output by the functions above
#re-pick input if you need to
input = 'brown_moduleInput.csv'
goDivision = 'CC'

#upload GO enrichment results
resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
go.res=read.table(resName, header = T)
go.res = go.res %>% 
  arrange(pval)
sig = go.res[go.res$p.adj < 0.1,]
head(sig, n=20)


#upload the gene-GO associations
geneResName=paste(goDivision, input, sep='_')
gene.res=read.table(geneResName, header=T)
head(gene.res)


#search the GO results for a particular term
searchString = 'ribo'
sig[grep(searchString, sig$name),] 


#select the GO term from your results file (sigGo above)
go="GO:0006979" #oxidative stress
go="GO:0000302" #response to reactive oxygen species
go='GO:0004930' #G-protein coupled receptor activity
go='GO:0007409;GO:0048812' #neuron projection morphogenesis
go='GO:0043009;GO:0009792;GO:0009790'

#subset for that GO term
go.genes = gene.res[gene.res$term == go, 'seq']
length(go.genes)

#get gene names
geneSet = as.character(go.genes)

#GATHER GENE NAMES
adat = read.table('../../metadata/Amillepora_euk.emapper.annotations.tsv',
                  sep='\t',
                  header = TRUE)
rownames(adat) = adat$query_name
annots = adat[geneSet,] %>% 
  dplyr::select(query_name,eggNOG.annot) %>% 
  dplyr::rename(gene=query_name)
nrow(annots)
length(geneSet)


#merge with deseq results
ll=load('../../correlated_only/deseqResults/stress_deseqResults.Rdata')
ll
mdat = res %>% 
  data.frame() %>% 
  mutate(gene = rownames(res)) %>% 
  right_join(annots, by = 'gene')
head(mdat)
geneSet = mdat$gene
geneNames=mdat$eggNOG.annot

#check the log2 values match expectation from heatmap
mdat %>% 
  ggplot(aes(y=log2FoldChange)) +
  geom_boxplot()


#------------ BUILD HEATMAP FOR A GIVEN GO TERM ------------#
library(pheatmap)


#you'll need the variance stabilized counts
ll=load('../../largeIgnored/strongStressOnly_project_controlled.Rdata')
ll
rld.df=data.frame(t(datExpr))
head(rld.df)


#subset the variance stabilized counts to include only the GO of interest
g=rld.df[geneSet,]
g[1:10,1:10]
cdat2 = coldata[colnames(g),] %>% 
  mutate(stress=if_else(treat=='control',
                        'control',
                        'stress'))
cdat2=cdat2[with(cdat2,order(cdat2$stress, cdat2$my_title)),]
orderedRuns = cdat2$Run
orderedTreats = cdat2$stress
orderedG = g[,orderedRuns]

pheatmap(test,
         cluster_cols=FALSE,
         cluster_rows=TRUE,
         labels_col=paste(orderedTreats, orderedRuns, sep='.')[1:10])









