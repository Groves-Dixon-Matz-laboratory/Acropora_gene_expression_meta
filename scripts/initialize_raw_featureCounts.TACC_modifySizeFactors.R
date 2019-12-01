#!/usr/bin/env Rscript
#intialize_raw_featureCounts_TACC.R
#clean up and format raw counts
#input is expected from feature counts or merge_feature_counts.R 
#can be run with subsets of data on both sides, will take overlap
#and output missing Runs from either the coldata or counts table

library(optparse)



option_list = list(
  
  make_option(c("--counts"), type="character", default=NULL, 
              help="Counts file name"),

  make_option(c("--coldat"), type="character", default=NULL, 
              help="Coldata file name"),

  make_option(c("--mnCountCut"), type="integer", default=5, 
              help="Mean read count. Genes with < this value will be removed."),

  make_option(c("--treeCut"), type="integer", default=50000, 
              help="Tree cut level"),

  make_option(c("--callOutliers"), type="logical", default=TRUE, 
              help="Logical for whether to do outlier analysis"),

  make_option(c("--o"), type="character",  
              help="output prefix")


)

print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
countsIn = 'all_featureCounts_geneCounts_laneDupsRemd.tsv'
coldataIn = 'ALL_Coldata.csv'
mnCountCut = 3
treeCut = 5000
doOutliers = FALSE
outPrefix = 'ALL'



suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(WGCNA))
suppressMessages(library(pheatmap))
suppressMessages(library(tidyverse))


#upload col data
print("Reading in coldata file...")
coldata0 = read.csv(coldataIn, stringsAsFactors=FALSE)
rownames(coldata0)=coldata0$Run
print("Head:")
print(head(coldata0))


#upload read counts
print(paste("Reading in counts file"))
counts = read.table(countsIn, sep="\t", header = T, row.names='Geneid', stringsAsFactors=FALSE) %>%
  select(-Chr, -Start, -End, -Strand, -Length)
print(paste(ncol(counts), 'total samples'))



#SET UP COLDATA
#first revise the miss-sexed individuals
traitRuns = coldata0$Run
countRuns = colnames(counts)
inCountNotTrait = countRuns[!countRuns %in% traitRuns]
inTraitNotCount = traitRuns[!traitRuns %in% countRuns]
overlap = countRuns[countRuns %in% traitRuns]

#match up the runs from coldata and counts table
print("Matching up the coldata and counts data...")
print(paste("Total runs in counts =", length(countRuns)))
print(paste("Total runs in trait data =", length(traitRuns)))
print(paste("Total overlap=", length(countRuns[countRuns %in% traitRuns])))
coldata=coldata0[overlap,]
counts=counts[,overlap]

#print out info on matching up Runs
countNotTrait = paste(outPrefix, 'in_count_not_trait.tsv', sep='_')
traitNotCount = paste(outPrefix, 'in_trait_not_count.tsv', sep='_')

if (length(inCountNotTrait) >0){
  print('writing runs found in count file but not in traits file to in_count_not_trait.tsv')
  data.frame(inCountNotTrait) %>%
    write_tsv(path=countNotTrait)
}
if (length(inTraitNotCount) >0){
  print('writing trait data for runs not found in count file to in_trait_not_count.tsv')
  coldata0 %>%
    filter(Run %in% inTraitNotCount) %>%
    write_tsv(path=traitNotCount)
}
if (length(inCountNotTrait)==0 & length(inTraitNotCount)==0){
  print('counts runs and trait table runs matched up.')
}
if (nrow(coldata) != length(overlap)){
  exit('Error, matching up files broke somehow.')
}



#get count totals
tots = apply(counts, 2, sum)
# pdf(file="total_read_counts_hist.pdf")
# hist(tots, main='total read counts')
# dev.off()
# write.table(data.frame(tots), file="counted_on_genes.tsv", sep="\t", quote=F) 
mean(tots)/1e6
median(tots)/1e6


#remove genes with low coverage
print("Removing genes with low coverage...")
cc=counts
means=apply(cc,1,mean)
print(paste('Cutoff =', mnCountCut))
print('Passing/Failing counts:')
print(table(means>mnCountCut))
print(paste('Removed percentage =', round(sum(means<mnCountCut)/length(means), digits=3)*100 ))
counts=cc[means>mnCountCut,]


###########
#upload the module membership data
mm=read.table('wgcna_moduleMembership.tsv', header = TRUE)
head(mm)
green = mm[mm$assignment=='green4' & mm$MMgreen4<0.4,]
dim(green)
genesForsizeFactor = green$gene
length(genesForsizeFactor)
countsForSizeFactor = counts[genesForsizeFactor,]
dim(countsForSizeFactor)


#------- GET RAW VARIANCE STABILIZED COUNTS ------------#

#First estimate size factors
print("Setting data for DESeq")
ddsHTSeqGreen<-DESeqDataSetFromMatrix(countsForSizeFactor,
	colData = coldata,
	design = formula(~1))
ddsGreen = estimateSizeFactors(ddsHTSeqGreen)
greenSizeFactors = sizeFactors(ddsGreen)
head(greenSizeFactors)


#now get variance stabilized for entire set using those size factors
dds<-DESeqDataSetFromMatrix(counts,
                                 colData = coldata,
                                 design = formula(~1))
sizeFactors(dds)<-greenSizeFactors
dds <- estimateDispersions(dds)
vsd <- varianceStabilizingTransformation(dds)
rld.df=assay(vsd)
save(vsd, file='greenVsd.Rdata')

#=====================================================================================
#
#  Code chunk 2
# transpose the dataset you have samples as rows and genes as columns
#=====================================================================================

datExpr0 = as.data.frame(t(rld.df));

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

#check that the dataset doesn't have geneswith too many missing values
#these would likely represent lowly expressed genes and under sequenced samples
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK



#=====================================================================================
#
#  Code chunk 4

#=====================================================================================
#removing genes that were flagged with too many missing values
#note how many genes we have right now
before = ncol(datExpr0)
print(before)


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
rld.df=t(datExpr0)
vsd=vsd[rownames(rld.df),]
counts=counts[rownames(rld.df),]
dim(rld.df)
dim(vsd)
nrow(datExpr0)
after = ncol(datExpr0)
print(paste(before - after, "Genes With Too Many Missing Values Were Removed"))

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

###build sample heatmaps 
#print("Building heatmap...")
#pdf(file='raw_heatmap.pdf')
#pheatmap(cor(rld.df), labels_row=coldata$treat, labels_col=coldata$my_title)
#dev.off()

# #now cluster samples based on gene expression to identify outliers
# print("Building Sample Clustering Tree...")
# sampleTree = hclust(dist(datExpr0), method = "average");
# # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# # The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(12,9)
# #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
# pdf(file='sample_clustering_tree.pdf')
# par(cex = 0.6);
# par(mar = c(0,5,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# dev.off()

#=====================================================================================
#
#  Code chunk 6
# 
#=====================================================================================


if (doOutliers){
  print('Doing outlier analyses...')


  #Remove outliers by setting a branch cut threshold
  # Plot a line to show the cut
  print(paste(c("Replotting cut tree with cutoff at ", treeCut, '...'), collapse=''))
  cut.height = treeCut
  pdf('cut_clustering_tree.pdf')
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  abline(h = cut.height, col = "red", lty = 2);
  dev.off()
  # Determine cluster under the line
  clust = cutreeStatic(sampleTree, cutHeight = cut.height, minSize = 4)
  table(clust)
  # clust 1 contains the samples we want to keep.
  keepSamples = (clust==1)
  keepSampleNames = rownames(datExpr0)[keepSamples]
  outlierNames = rownames(datExpr0)[clust==0]
  datExpr = datExpr0[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr) #number of samples left after outlier removal
  print(paste(length(outlierNames), "samples were flagged as outliers and removed:"))
  outlierNames
  print(paste(nSamples, "samples were kept"))


  #replot heatmap without outlier
  rld.df = rld.df[, !colnames(rld.df) %in% outlierNames]
  counts = counts[, !colnames(counts) %in% outlierNames]

  #re-build sample heatmaps 
  print("re-Building heatmap after outlier removal...")
  pdf(file='outlierRemoved_heatmap.pdf')
  pheatmap(cor(rld.df))
  dev.off()
  counts=counts[,!colnames(counts) %in% outlierNames]
  coldata=coldata[!coldata$sample %in% outlierNames,]
} else {
  datExpr=datExpr0
  print('Skipping outlier analyses.')
}



#save the outlier names so you can optionally remove them in other analyses
# save(outlierNames, file = 'datasets/outliers.Rdata')
print('writing out results')
save(counts, coldata, file=paste(outPrefix, "deseqBaselineInput.Rdata", sep='_'))
save(vsd, rld.df, coldata, file=paste(outPrefix, "vsd.Rdata", sep='_'))

#also save input in format for WGCNA
datTraits=coldata
save(datExpr, datTraits, file=paste(outPrefix, 'wgcna_input_green4SizeFactors.Rdata', sep='_'))

