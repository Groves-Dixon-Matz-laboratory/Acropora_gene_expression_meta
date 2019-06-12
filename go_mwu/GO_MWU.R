
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.

setwd("./go_mwu/")

inputFiles = c('stressForGO.csv',
               'heatForGO.csv',
               'coldForGO.csv',
               'salinityForGO.csv',
               'immuneForGO.csv')


# Edit these to match your data file names: 
goAnnotations="amil_zachFullerV2_gos.tsv" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
source("gomwu.functions.R")

#select GO category
goDivision="CC"
goDivision="MF"
goDivision="BP"


for (input in inputFiles){
  print(paste(c('Running ', input, '...'), collapse=''))
  # Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
  gomwuStats(input, goDatabase, goAnnotations, goDivision,
  	perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
  	largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
  	smallest=5,   # a GO category should contain at least this many genes to be considered
  	clusterCutHeight=0.0, # threshold for merging similar (gene-sharing) terms. See README for details.
    # Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
  #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
  #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
  )
  # do not continue if the printout shows that no GO terms pass 10% FDR.
}



# IDENTIFY CORE STRESS GO TERMS -------------------------------------------

#funtion to read in results from gomwuStats()
library(tidyverse)
read_go_results = function(inputName){
  resName = paste(paste('MWU', goDivision, sep = "_"), inputName, sep = "_")
  res=read.table(resName, header = T, stringsAsFactors=F)
  res=res[order(res$pval),]
  genName = sub('ForGO.csv', '', inputName)
  res[,genName] = res$p.adj<0.1
  res %>% 
    select(c('term', genName))
}


#select GO category
goDivision="CC"
goDivision="MF"
goDivision="BP"



#These are the ones that show up for all five datasets
gdat <- inputFiles %>%
  map(function(x) read_go_results(x)) %>% 
  purrr::reduce(full_join, by = c('term'))
head(gdat)
nSig = apply(gdat[,2:ncol(gdat)], 1, sum)
table(nSig)
passingTerms= gdat %>% 
  filter(nSig==5) %>% 
  pull(term)




# BUILD MODIFIED RESULTS FILES FOR PLOTTING THE CONCENSUS -----------------

#first upload results from the all stress set
allStressInput = 'stressForGO.csv'
resName = paste(paste('MWU', goDivision, sep = "_"), allStressInput, sep = "_")
res=read.table(resName, header = T)
res=res[order(res$pval),]

#modify it to be significant only for the passing terms
modInput = 'coreForGO.csv'
modResName = paste(paste('MWU', goDivision, sep = "_"), modInput, sep = "_")
modRes = res %>% 
  mutate(p.adj = if_else(term %in% passingTerms,
                         p.adj,
                         1.0))
write.table(modRes, file=modResName, row.names=F)

#also copy the version of the secondary results file to give it new modified name
#(note reading this in and writing it back out somehow screwed things up)
resName2 = paste(goDivision, allStressInput, sep='_')
modResName2 = paste(goDivision, modInput, sep='_')
copyCommand = paste(c('cp', resName2, modResName2), collapse=' ')
system(copyCommand)


# PLOT RESULTS FOR PASSING TERMS ONLY -------------------------------------

# Plotting results
quartz()
results=gomwuPlot(modInput,goAnnotations,goDivision,
#	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
	absValue=1,
	level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
	level2=0.05, # FDR cutoff to print in regular (not italic) font.
	level3=0.001, # FDR cutoff to print in large bold font.
	txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.5, # height of the hierarchical clustering tree
#	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results

