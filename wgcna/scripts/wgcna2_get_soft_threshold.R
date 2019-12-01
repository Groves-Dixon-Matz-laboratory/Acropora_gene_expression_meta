#!/usr/bin/env Rscript
#wgcna2_get_soft_threshold.R

#this script is for deciding what soft theshold to use.
#it is based off the first half of 'FemaleLiver-02-networkConstr-auto.R'
#in the WGCNA tutorial set.


#the process is too hard for desktop computers, so this
#this is intended to be run on TACC.

#We want the scale-free topology fit index to reach at least
# 0.8 at a power less than 15 for unsigned networks
# or less than 30 for a signed network (WGCNA FAQs). 


#upload libarary
library(WGCNA)
library(optparse)
options(stringsAsFactors = FALSE) #important for WGCNA

option_list = list(
  
  make_option(c("--input"), type="character", default='wgcna_input.RData', 
              help="R object with datExpr and datTraits variables"),
  make_option(c("--networkType"), type="character", default='signed', 
              help="Network type to use for getting soft threshold")
)

print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input = opt$input
NETWORK.TYPE = opt$networkType
print(paste("Loading File", input))
print(paste("Getting Soft Threshold for network type =", NETWORK.TYPE))

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
#enableWGCNAThreads()



# Load the data saved in the first part
lnames = load(file = input);
#The variable lnames contains the names of loaded variables.
lnames



# Choose a set of soft-thresholding powers

powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = NETWORK.TYPE, verbose = 5)



# Plot the results:

pdf("soft_threshold_plot1.pdf")
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
# abline(h=0.80,col="red")
#dev.off()
#pdf("soft_threshold_plot2.pdf")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()





