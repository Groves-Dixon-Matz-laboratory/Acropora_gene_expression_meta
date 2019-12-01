#!/usr/bin/env Rscript
# wgcna3b_step-wise_network_construction.R

#example Usage:
#wgcna3b_step-wise_network_construction.R 18 10 0 wgcna-01_output.RData signed



library(optparse)
options(stringsAsFactors = FALSE) #important for WGCNA

option_list = list(
  
  make_option(c("--softPower"), type="integer",
              help="Soft power threshold to use. (Get this with wgcna2_get_soft_threshold.R)"),

  make_option(c("--minSize"), type="integer", default=30, 
              help="Minimum module size"),

  make_option(c("--mergeCutoff"), type="double", default=0, 
              help="Module Merging Threshold"),

  make_option(c("--input"), type="character", default='wgcna_input.RData', 
              help="Network type to use for getting soft threshold"),

  make_option(c("--networkType"), type="character", default='signed', 
              help="Network type to use (signed or unsigned)"),

  make_option(c("--nCores"), type="integer", default=1, 
              help="Number of cores to use")
)


print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
softPower = opt$softPower
minModuleSize = opt$minSize
MEDissThres = opt$mergeCutoff
input = opt$input
NETWORK.TYPE = opt$networkType
nCores = as.numeric(opt$nCores)


#print out arguments to user
print("----------------------------")
print("----------------------------")
print(paste("Using Soft Threshold Power =", softPower))
print(paste("Using Minimum Module Size =", minModuleSize))
print(paste("Using Module Merging Threshold =", MEDissThres))
print(paste("Using Network type =", NETWORK.TYPE))
print(paste("Number of cores to use =", nCores))

#assemble the output file name specific to the chosen parameters
soft = paste("sft",softPower,sep = "")
minModSize = paste("minModSize",minModuleSize,sep = "")
cutDpth = paste("cutDpth", MEDissThres, sep = "")
stats = paste(paste(soft, minModSize, sep = "_"), cutDpth, sep = "_")
outName = paste( paste(paste("wgcna3b_manual", stats, sep = "_"), NETWORK.TYPE, sep = "_"), "Rdata", sep = ".")
print(paste("Saving results as", outName))


suppressMessages(library(WGCNA))
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
if (nCores>1){
  allowWGCNAThreads(nThreads=nCores)
}
# Load the data saved in the first part
lnames = load(file = input);


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


#calculate adjacency using the soft power threshold from wgcna2 figure
adjacency = adjacency(datExpr, power = softPower, type=NETWORK.TYPE);

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType=NETWORK.TYPE);
dissTOM = 1-TOM

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================



# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)


#plot the gene tree
pdf(file="geneTree.pdf", width=40, height=20)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()



#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


#Note minModuleSize was given as argument when running script

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = TRUE, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)

#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
                    

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================
                  
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors, softPower = softPower)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)

#plot them
pdf(file="Clustering_module_eigengenes.pdf")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================

#tutorial
# Plot the cut line into the dendrogram
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


sizeGrWindow(12, 9)
pdf(file = "Plots_geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================

#tutorial
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;


# Save module colors and labels for use in subsequent parts
sampleNames = rownames(datExpr)
save(sampleNames,  MEs, moduleLabels, moduleColors, geneTree, dynamicColors, file = outName)



