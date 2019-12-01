#wgcna4_module-correlations.R
#This script calculates and plots relationships between 
#module eigengenes and sample traits. Inputs come from:
#get_variance_stabilized_counts.R and wgcna3b_step-wise_network_construction.R
#code is adapted from examples given here: https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/


#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
rm(list=ls())

# Load the WGCNA package
library(WGCNA)
library(tidyverse)
source('wgcna/scripts/wgcna_functions.R')
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


# Load the expression and trait data saved in the first part
lnames = load(file = "largeIgnored/project_controlled_varCut0.75_expCut0.75.Rdata");   
#The variable lnames contains the names of loaded variables.
lnames #(when you load a .Rdata object, if can have multiple variables assigned, the list of them is saved in this lnames variable, in this case the two objects you loaded are "datExpr" and "datTraits")


# Load network data saved in the second part.
lnames = load(file = "wgcna/inputs/wgcna3b_manual_sft12_minModSize30_cutDpth0.3_signed.Rdata")
# lnames = load(file = "wgcna/temp/wgcna3b_manual_sft12_minModSize30_cutDpth0.3_signed.Rdata")
lnames


#match up the rownames in the expression 
rownames(datExpr)[1:10]     #this tells us the rownames as given to WGCNA
rownames(datTraits)[1:10]   #this is for the entire set of datTraits and includes samples that were not put into WGCNA
rownames(MEs) = rownames(datExpr)
datTraits = datTraits[rownames(datTraits) %in% rownames(datExpr),]  #reduce datTraits to those we have WGCNA results for
#if everything is matched up this should say TRUE
sum(rownames(datTraits) == rownames(datExpr)) == length(rownames(datExpr))
head(datTraits)
dim(datTraits)


#load extra column data
coldata = read.csv('metadata/ALL_Coldata.csv') #replaced with R oject made with symbiont_count_analysis.R
rownames(coldata)=coldata$Run
# ll=load('metadata/coldataWithSymbionts.Rdata')
# ll
scdata = read.csv('metadata/subset_tables/allStress_Coldata.csv')
stressRuns = scdata$Run


#check order of sample IDs in coldata
c2 = coldata[rownames(datExpr),]
sum(c2$Run==rownames(datExpr))==nrow(datExpr)
c2$bleached[is.na(c2$bleached)]<-'no'
c2$blch=c2$bleached


#load high vs low stress projects
ll=load('metadata/corStressProjs.Rdata')
ll
corStressProjs
lowStressProjs


datTraits2=c2 %>% 
  mutate(stress = if_else(Run %in% scdata$Run & treat!='control',
                          'stressed',
                          'control'),
         stress = ifelse(!Run %in% scdata$Run,
                         NA,
                         stress),
         highStress = ifelse(my_title %in% corStressProjs,
                             1,
                             NA),
         highStress = if_else(my_title %in% corStressProjs & stress=='control',
                              0,
                              highStress),
         lowStressVsAll = if_else(my_title %in% lowStressProjs & stress=='stressed',
                                  1,
                                  0),
         lowStress = ifelse(my_title %in% lowStressProjs,
                            1,
                            NA),
         lowStress = if_else(my_title %in% lowStressProjs & stress=='control',
                             0,
                             lowStress),
         #bleached
         bleached = ifelse(stress=='stressed' & blch=='yes' & (my_title %in% corStressProjs),
                           1,
                           NA),
         bleached = if_else(stress=='control' & blch=='yes' & (my_title %in% corStressProjs),
                            0,
                            bleached),
         #heat
         heat = if_else(stress=='stressed',
                        1,
                        0),
         heat = ifelse(projType=='temp' & !grepl('^Ch', treatDescription),
                       heat,
                       NA),
         #heat no beww
         heat = ifelse(blch=='yes',
                       NA,
                       heat),
         #high heat
         high.heat = ifelse(my_title %in% corStressProjs,
                            heat,
                            NA),
         low.heat = ifelse(my_title %in% lowStressProjs,
                           heat,
                           NA),
         high.heatNoBleach = ifelse(blch=='yes',
                                    NA,
                                    high.heat),
         #low salinity
         hyposalinity = if_else(stress=='stressed',
                                1,
                                0),
         hyposalinity = ifelse(projType == 'salinity',
                               hyposalinity,
                               NA),
         #low salinity no beww
         hyposal.noBEWW = ifelse(my_title=='j1_thisStudy_PRJNA559404',
                                 NA,
                                 hyposalinity),
         #immune
         immune = if_else(stress=='stressed',
                          1,
                          0),
         immune = ifelse(projType == 'immune',
                         immune,
                         NA),
         high.immune = ifelse(my_title %in% corStressProjs,
                              immune,
                              NA),
         low.immune = ifelse(my_title %in% lowStressProjs,
                             immune,
                             NA),
         #low pH
         pH = if_else(stress=='stressed',
                      1,
                      0),
         pH = ifelse(projType == 'pH',
                     pH,
                     NA),
         high.pH = ifelse(my_title %in% corStressProjs,
                          pH,
                          NA),
         low.pH = ifelse(my_title %in% lowStressProjs,
                         pH,
                         NA),
         stress = if_else(stress=='control',
                          0,
                          1)
  )

#choose the set of sample traits to plot and reformat
selection = c('Run', 'bleached', 'high.heatNoBleach', 'hyposal.noBEWW', 'high.immune', 'high.pH', 'low.heat', 'low.immune', 'low.pH')
datTraits2 = datTraits2%>% 
  dplyr::select(selection)
rownames(datTraits2) = datTraits2$Run
datTraits2$Run<-NULL
apply(datTraits2, 2, function(x) sum(x, na.rm=TRUE))
head(datTraits2)

#change names for plotting
datTraits3 = datTraits2 %>% 
  dplyr::rename('heat'=high.heatNoBleach,
                'hyposalinity'=hyposal.noBEWW,
                'immune'=high.immune,
                'pH'=high.pH,
                ' heat'=low.heat,
                ' immune'=low.immune,
                ' pH'=low.pH)

#check column classes
datTraits3 %>% 
  as_tibble()

# datTraits3 %>% 
#   write_csv(path='~/Desktop/traits.csv')

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate module eigengenes with color labels
#This will provide the first principal component for expression
#behavior for the genes in each module. On that principal component,
#each sample will have a loading value. These values can be corelated
#with our known sample traits to get an idea what biological mechanism
#the co-reguated genes in a given module might be responding to
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs1 = MEs0[rownames(datTraits2),]


########################################
############ SUBSET BY SIZE ############
########################################

#plot all modules
minModSizeToPlot = 0;outputMM=TRUE

# plot only large modules
minModSizeToPlot = 600;outputMM=FALSE

########################################
########################################
########################################


module.sizes = table(moduleColors)
passing=names(module.sizes)[module.sizes>minModSizeToPlot]
MEs = orderMEs(MEs1[,paste('ME', passing, sep='')])


########################################
############ COLOR SWAPPING ############
########################################

# #replot the dendrogram
# plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors),
#                     c("Dynamic Tree Cut", "Merged dynamic"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)


#choose colors to swap out
swap = data.frame(oldCols=c('brown',
                            'yellow',
                            'blue',
                            'magenta'),
                  newCols=c('red',
                            'white',
                            'green4',
                            'orange'))

#make swaps
modCols =  swap_cols(swap$oldCols, swap$newCols, moduleColors, colnames(MEs), dynamicColors)
moduleColors=modCols[[1]]
colnames(MEs) = modCols[[2]]
dynamicColors = modCols[[3]]


# #replot the dendrogram after color changes
# plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors),
#                     c("Dynamic Tree Cut", "Merged dynamic"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)

########################################
########################################

#use the cor() function to get the correlations between the module eigengenes and the trait data
moduleTraitCor = cor(MEs, datTraits3, use = "p");
#get p values as well
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#write out module loadings
mEigs=MEs
rownames(mEigs) = rownames(datTraits3)
save(mEigs, file='wgcna/moduleEigengenes.Rdata')



#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

sizeGrWindow(10,6)

#now replot the heatmap
sizeGrWindow(10,6)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
COLORS = greenWhiteRed(50)
blueWhiteRed(50)


# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
COLORS = greenWhiteRed(50)
blueWhiteRed(50)

#set up additional formatting variables
rows = rownames(moduleTraitCor)
sub.colors = substr(rownames(moduleTraitCor), 3, 50) #trim away the ME from the
module.sizes = paste(sub.colors, table(moduleColors)[sub.colors], sep = "\n")

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits3),
               yLabels = rownames(moduleTraitCor),
               ySymbols = module.sizes,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



#=====================================================================================
#REPLOT THE MODULE EIGENGENE CLUSTERING TREE
#This shows you how the modules are similar to one another.
#Pushing the merging theshold (argument 3 for wgcna3b_step-wise_network_construction.R) will join modules together
#This is like moving a horizontal line up this tree figure, if the line is above a node, modules below that node will be joined into one
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)

#plot them
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=0.35, lty=2, col='red')

#=====================================================================================


#=====================================================================================
#
#  Code chunk 4 Gather module membership data for the genes in each module
#
#=====================================================================================
#"module membership" is a measure of how strongly individual genes correlate with the 
#module eigengene. The gene that matches best with the module eigengene can be thought
#of as a hub gene, (ie it's variation in expression across the samples is most exemplary 
#of the module)

#SELECT A TRAIT AND P-VALUE CUTOFF
TRAIT="bleached"
TRAIT='stress'
PCUT=1

#SUBSET FOR MODULES SIGNIFICANT FOR THAT TRAIT
traitPvalues = moduleTraitPvalue[,TRAIT]
keep = traitPvalues < PCUT
subCor = moduleTraitCor[keep,]
subP = moduleTraitPvalue[keep,]
rows = rownames(moduleTraitCor)[keep]
sub.colors = substr(rownames(subCor), 3, 50) #trim away the ME from the
module.sizes = paste(sub.colors, table(moduleColors)[sub.colors], sep = "\n")
length(sub.colors)



# # Define dataframe trait.df containing the a trait of interest from datTraits3
trait.df = as.data.frame(datTraits3[,TRAIT], row.names=rownames(datTraits3));
names(trait.df) = TRAIT
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr[rownames(MEs),], MEs, use = "p"));
head(geneModuleMembership)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr[rownames(trait.df),], trait.df, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(trait.df), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.df), sep="");
modules = sub.colors
plot.cols = modules 

#save the module assignment data for other plotting
if (outputMM){
  save(geneModuleMembership, moduleColors, file='wgcna/moduleAssignment.Rdata')
}

redGo = data.frame(gene=rownames(geneModuleMembership),
                   redMembership = geneModuleMembership[,'MMred'])
head(redGo)
redGo %>% 
  write_csv(path='./go_mwu/redMembership_ForMWU.csv')
  


ll=load('wgcna/moduleAssignment.Rdata')
ll
head(geneModuleMembership)


#=====================================================================================
#
#  Code chunk 5 plot scatterplots of module membership and trait correlations
#
#=====================================================================================   
#The code below loops through each of the modules that were significant for the TRAIT of 
#interest that you chose above. Each point in the scatterpot is a gene in the module. 
#For each one of the it plots the the genes' correlation
#with the trait of interest against the genes' module memberships,
#(the correlation between the gene's variation accross samples and the module eigengene)
#A tight correlation here is suggestive that the correlated variation in gene expression 
#captured by the module is truely associated with the trait of interest.
#A tight correlation basically means, the better a gene fits into this module, the more strongly
#it correlates with the trait.



# quartz()
length(modules)
par(mfrow=c(2,2))
ggplotList = list()
for (m in modules){
	column = match(m, modNames);
	moduleGenes = moduleColors==m;
	
	# sizeGrWindow(7, 7);
	# par(mfrow = c(1,1));
	verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	                   abs(geneTraitSignificance[moduleGenes, 1]),
	                   xlab = paste("Module Membership in", m, "module"),
	                   ylab = paste("Correlation with", TRAIT),
	                   main = paste("Module membership vs. gene significance\n"),
	                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'black', bg = m, pch = 21, cex = 1.5)
	ggplotList[[m]]=ggVerboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	                                     abs(geneTraitSignificance[moduleGenes, 1]),
	                                     xlab = paste("membership in", m, "module"),
	                                     ylab = paste("correlation with", TRAIT),
	                                     main = paste(""),
	                                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'black', bg = m, pch = 21, cex = 1.5)
}
# ggplotList[['red']]
plot_grid(plotlist = ggplotList, nrow=1)


##################################################
################# GO MWU OUTPUTS #################
##################################################

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
modules = sub('MM', '', colnames(geneModuleMembership))
genes = rownames(geneModuleMembership)
inputFiles=c()

for (m in modules){
  moduleGenes = genes[moduleColors==m]
  godf = data.frame(gene=genes,
                    inMod=if_else(genes %in% moduleGenes,
                                  1,
                                  0))
  outname = paste(c('./wgcna/go_mwu/', m, '_moduleInput.csv'), collapse='')
  print(paste('module = ', m))
  print(paste('total genes =', length(moduleGenes)))
  write.csv(godf, file=outname, row.names=F, quote=F)
  inputFiles = append(inputFiles, paste(m, '_moduleInput.csv', sep=''))
}

save(inputFiles, file='./wgcna/go_mwu/moduleInputFiles.Rdata')











###################################################
#### LOOK AT RESULTS FOR A A PARTICULAR MODULE ####
###################################################

source("scripts/zebrafish_RNAseq_functions.R")
library("biomaRt")
embl = useMart("ensembl", dataset="drerio_gene_ensembl")
# listAttributes(embl)


##### BUILD BOXPLOTS #####
#gather genes for modules of interest
genes = rownames(geneModuleMembership)
modules = c('mediumpurple4', 'darkolivegreen4')
moduleGenes = c()
mod = c()
for (m in modules){
	print("---------------------")
	print(m)
	toAdd=genes[moduleColors==m]
	mlabs=rep(m, length(toAdd))
	print(paste(length(toAdd), "genes added"))
	print(toAdd)
	print(mlabs)
	moduleGenes = append(moduleGenes, toAdd)
	mod=append(mod, mlabs)
}
mres=data.frame(moduleGenes, mod)
rownames(mres) = mres[,1]
head(mres)

#merge with DESeq results
lnames=load("results/ethanol_full_LRT_results.Rdata")
head(res.eth)
r=data.frame(res.eth[rownames(res.eth) %in% moduleGenes,])
mdat=merge(mres,r,by=0)
mdat
nrow(mdat) == nrow(mres)

# b=boxplot(mdat$log2FoldChange~mdat$mod, axes=F, notch=T, outline=F, ylab="Log2 Fold Difference", ylim=c(-2,2.1))
mdat$xs=match(mdat$mod, rev(modules))
plot(mdat$log2FoldChange~jitter(mdat$xs, factor=1.2), bg=mdat$mod, pch=21, cex=1.2, ylim=c(-2,2.1), xlim=c(.5,2.5), axes=F, ylab="Log2 Fold Difference", xlab='')
axis(2)
axis(1, at=1:length(modules), label=rev(modules))
b=boxplot(mdat$log2FoldChange~mdat$mod, axes=F, outline=F, add=T, lwd=1.5, xlab="")
abline(h=0, lty=2, col='black', lwd=1.5)


#do ggplot boxplots
library(ggplot2)
p <- ggplot(mdat, aes(x=mod, y= log2FoldChange)) + geom_boxplot() + coord_cartesian(ylim = c(-2, 2)) + theme_bw() + xlab("Module") + ylab("Log2 Fold Difference")
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, fill=rev(mod)) + geom_hline(yintercept=0, linetype='dashed', color='black', size=1)




#pick the module to look at
m='mediumpurple4'
m='darkolivegreen4'



#gather gene names for this module
genes = rownames(geneModuleMembership)
moduleGenes = genes[moduleColors==m]
moddat = data.frame(geneModuleMembership[moduleGenes, paste("MM",m,sep='')]);rownames(moddat) = moduleGenes;colnames(moddat)=c(m)
head(moddat)
dat=merge_gene_names(moddat, sort.column=m)
head(dat)
length(moduleGenes)
dim(dat)
dat


#write out the gene names
fileName=paste(m, 'gene_data.tsv',sep="_")
write.table(dat, file=paste('results', fileName, sep="/"), sep="\t", quote=F, row.names=F)


#LOOK AT DISTRIBUTION ACCROSS CHROMOSOME 4
datc4=na.omit(dat[dat$chromosome_name==4,])
print("Percentage of module on Chromosome 4:")
nrow(datc4)/nrow(dat)*100
chrom4_length = 76625712 #taken from here http://www.ensembl.org/Danio_rerio/Location/Chromosome?r=4:59726329-59826329
plot(density(datc4$start_position/1e6, bw=4), xlim=c(0,chrom4_length/1e6+1), main="Module Gene Locations", xlab='Chromosome 4 position (Mb)', axes=F)
axis(1, at = seq(10, 70, by = 10))
axis(2)


#BUILD HEATMAPS TO SHOW EXPRESSION DIFFERENCES
lnames=load("datasets/raw_rld.Rdata")
head(rld.df)

#select for this module
x=rld.df[rownames(rld.df) %in% moduleGenes, colnames(rld.df) %in% rownames(datTraits3)]
dim(x)
length(moduleGenes)
head(x)
sum(colnames(x) == rownames(datTraits3)) == ncol(x)


#sort the dataframe by ethanol treatment, and time
datTraits3$ordTime = as.numeric(as.character(datTraits3$time))
y=datTraits3[with(datTraits3, order(ordTime, ethanol)), ]
y=datTraits3[with(datTraits3, order(ethanol, ordTime)), ]
# y=y[y$time==2,]
x=x[,rownames(y)]
colnames(x)


#get z-scores for genes
head(x)
gmeans = apply(x, 1, mean)
gsds = apply(x, 1, sd)

z=(x-gmeans)/gsds
head(z)

library(pheatmap)
library(plotrix)
labs = get_gene_names(rownames(z))
labels = labs$description
for (i in 1:length(labels)){
	if(labels[i] == ""){
		labels[i]<-labs$external_gene_name[i]
	}
}

pheatmap(z,cluster_cols=F, treeheight_row=0, border_color=NA,clustering_distance_rows="correlation", labels_row=labels)

#BUILD BARPLOTS
lnames=load("results/ethanol_full_LRT_results.Rdata")
head(res.eth)
r=res.eth[rownames(res.eth) %in% moduleGenes,]
mn=mean(r$log2FoldChange)
se = std.error(r$log2FoldChange)
sd = sd(r$log2FoldChange)
plotCI(x=0,y=mn,uiw=se, axes=F, xlab=m, ylab='Mean log2 Fold Difference');axis(2)
YLIM=c(-.4, .4)
YLIM=c(-1.3,1.3)
plotCI(x=0,y=mn,uiw=se, ylim=YLIM, axes=F, xlab=m, ylab='Mean log2 Fold Difference', lwd=1.5)
axis(2)
abline(h=0, lty=2)


#store for plotting both
keepmn=mn
keepse=se
keepcol=m
YLIM=c(-.4, .4)
plotCI(x=0,y=keepmn,uiw=keepse, ylim=YLIM, axes=F, xlab=keepcol, ylab='Mean log2 Fold Difference', lwd=1.5)
axis(2)
abline(h=0, lty=2)



#plot both on same axis
mns = c(mn, keepmn)
ses = c(se, keepse)
cols = c(m, keepcol)
plotCI(x=1:2,y=mns,uiw=ses, axes=F, col=cols, pch=26, lwd=2, xlab='', ylab='Mean log2 Fold Difference', ylim=c(-1.3,1.3), xlim=c(0,3));axis(2)
abline(h=0, lty=2)
text(cex=1, x=1:2, y=-1.5, cols, xpd=TRUE, srt=45)
abline(h=0, lty=2)


#OUTPUT FOR GO ENRICHMENT
up = c('mediumpurple4')
down = c('darkolivegreen4')

#output for upregulated modules
mgenes=c()
for (m in up){
	mgenes=append(mgenes, genes[moduleColors==m])
}
mgenes=unique(mgenes)
o1=data.frame(genes, 0)
colnames(o1) = c('gene', 'stat')
o1$stat[o1$gene %in% mgenes]<-1
head(o1)
sum(as.numeric(o1$stat))
dim(o1)
write.csv(o1, file='~/gitreps/zebrafish_early_ethanol_RNASeq/wgcna/goMWU/ethanol_up_mod_genes.csv', row.names=F, quote=F)


#output for downregulated modules
#output for upregulated modules
mgenes=c()
for (m in down){
	mgenes=append(mgenes, genes[moduleColors==m])
}
mgenes=unique(mgenes)
o1=data.frame(genes, 0)
colnames(o1) = c('gene', 'stat')
o1$stat[o1$gene %in% mgenes]<-1
head(o1)
sum(as.numeric(o1$stat))
dim(o1)
write.csv(o1, file='~/gitreps/zebrafish_early_ethanol_RNASeq/wgcna/goMWU/ethanol_down_mod_genes.csv', row.names=F, quote=F)


