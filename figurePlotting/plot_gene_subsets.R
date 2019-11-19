#plot_gene_subsets.R


#------------ LOOK AT GENES WITHIN GO TERMS ------------#
library(cowplot)
library(tidyverse)
setwd('wgcna/go_mwu/')
#upload the results for each genes output by the functions above
#re-pick input if you need to
input = 'brown_moduleInput.csv'
input = 'blue_moduleInput.csv'
goDivision = 'CC'

#upload GO enrichment results
resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
go.res=read.table(resName, header = T, stringsAsFactors=FALSE)
go.res = go.res %>% 
  arrange(pval)
sig = go.res[go.res$p.adj < 0.1,]
head(sig, n=20)


#upload the gene-GO associations
geneResName=paste(goDivision, input, sep='_')
gene.res=read.table(geneResName, header=T)
head(gene.res)


#search the GO results for a particular term
searchString = 'ribosom'
search=go.res[grep(searchString, go.res$name),] 
search

#set gos to the select terms
go=search$term
go

#subset for that GO term
goSub = gene.res[gene.res$term %in% go,] %>% 
  dplyr::rename('gene'=seq)



#GATHER GENE NAMES
adat = read_tsv('../../metadata/Amillepora_euk.emapper.annotations.tsv') %>% 
  dplyr::rename('gene'=query_name)
head(adat)
annots = goSub %>% 
  left_join(adat, by = 'gene')


#merge with deseq results
ll=load('../../deseqResults/corStress_deseqResults.Rdata')
ll
resdf = res %>% 
  data.frame() %>% 
  mutate(gene = rownames(res)) 
mdat = annots %>% 
  left_join(resdf, by = 'gene') %>% 
  mutate(riboType = if_else(grepl('ribosomal protein', eggNOG.annot),
                            'ribosome',
                            'other'),
         riboType = if_else(grepl('mitochondrial', eggNOG.annot),
                            'mitochondrial',
                            riboType))
head(mdat)
geneSet = mdat$gene
geneNames=mdat$eggNOG.annot

#check the log2 values match expectation from heatmap
mdat %>% 
  ggplot(aes(y=log2FoldChange)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2)


#check for different types of ribosomal genes
mdat %>% 
  ggplot(aes(x=riboType, y=log2FoldChange)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2)


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
cdat2 = datTraits[colnames(g),]
cdat2=cdat2[with(cdat2,order(cdat2$stress, cdat2$my_title)),]
orderedRuns = cdat2$Run
orderedTreats = cdat2$stress
orderedG = g[,orderedRuns]
catLabs = paste(orderedTreats, orderedRuns, sep='.')

pheatmap(orderedG,
         cluster_cols=FALSE,
         cluster_rows=TRUE,
         labels_col=orderedTreats)






# COMPARE GO AND KOG ANNOTATIONS ------------------------------------------
#goal here is to figure out why the kog annotations don't show 
#what we expect for ribosomal stuff
setwd('~/gitreps/Acropora_gene_expression_meta/')
kog = read_tsv('kog_mwu/Amillepora.v2.00.chrs_gene2kogClass.tab', col_names=c('gene', 'kog')) %>% 
  left_join(adat, by = 'gene') %>% 
  select(gene, kog, GO_terms, KEGG_KOs, eggNOG.annot)

selectKog = "Translation, ribosomal structure and biogenesis"
skog = kog %>% 
  filter(kog==selectKog) %>% 
  left_join(resdf, by = 'gene')


#look at gene names for ribosomes
allNames = skog %>% 
  pull(eggNOG.annot) %>% 
  unique()


#set some rules for true ribosomal genes
skog = skog %>% 
  mutate(riboType = if_else(grepl('ribosomal protein', eggNOG.annot),
                            'ribosome',
                            'other'),
         riboType = if_else(grepl('Ribosomal protein', eggNOG.annot),
                            'ribosome',
                            riboType),
         riboType = if_else(grepl('mitochondrial', eggNOG.annot),
                            'mitochondrial',
                            riboType),
         inGo = gene %in% mdat$gene)



#write out to look at it
skog %>% 
  filter(riboType=='other') %>% 
  write_csv(path='~/Desktop/other.csv')


#pull out true ribosomes
skog %>% 
  ggplot(aes(x=riboType, y = log2FoldChange)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2)


#pull out true ribosomes
skog %>% 
  ggplot(aes(x=inGo, y = log2FoldChange)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2)



skog %>% 
  filter(log2FoldChange > 0) %>% 
  write_csv(path='~/Desktop/up.csv')
