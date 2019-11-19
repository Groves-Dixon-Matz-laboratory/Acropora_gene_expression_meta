

##The unexpected result of a subset of stress projects showing opposite response as the rest made me think there may be a bug
#This was to double-check that this was not the case
#everything appears to have checked out, but keeping to script to remind myself I checked.

library(tidyverse)
ll=load('largeIgnored/allStress_deseqBaselineInput.Rdata')
ll
head(coldata)

runIndividual=TRUE
variable = 'stress'
treatString='stressed'
controlString = 'control'
# loop through projects individually --------------------------------------
if (runIndividual==TRUE){
  print('Running individual projects...')
  formula = paste('~', variable,sep='')
  projects = unique(coldata$my_title)
  for (proj in projects){
    print('----------')
    print(paste(c('Running project', proj,'...'), collapse=''))
    subColdata = coldata %>% 
      filter(my_title==proj)
    print(paste('total samples in subcoldata = ', nrow(subColdata)))
    subCounts = counts %>% 
      select(subColdata$Run)
    print(paste('total samples in subcounts = ', ncol(subCounts)))
    if (nrow(subColdata)<50){
      dds<-DESeqDataSetFromMatrix(subCounts,
                                  colData = subColdata,
                                  design = formula(formula))
      dds <- DESeq(dds)
      res = results(dds, contrast = c(variable, treatString, controlString))
      outName = paste(proj,'deseqResults.Rdata', sep='_')
      save(res, file=outName)
    }
  }
}



# check an individual one -------------------------------------------------


counts = read.table('largeIgnored/all_featureCounts_geneCounts_laneDupsRemd.tsv', sep='\t', header = TRUE)
rownames(counts)=counts$Geneid
coldata=read_csv('metadata/ALL_Coldata.csv')
unique(coldata$my_title)
sproj = 'H_matz_heatTolLat_PRJNA279192'

subColdata = coldata %>% 
  filter(my_title==sproj,
         my_stage=='adult')
subCounts = counts[,subColdata$Run]
head(subCounts)
dim(subColdata)
dim(subCounts)


ddsTEST<-DESeqDataSetFromMatrix(subCounts,
                            colData = subColdata,
                            design = ~stress)
ddsTEST <- DESeq(ddsTEST)
resTEST = results(ddsTEST, contrast = c('stress', 'stressed', 'control'))
head(resTEST)


ll=load('~/gitreps/Acropora_gene_expression_meta/deseqResults/individual_projects/H_matz_heatTolLat_PRJNA279192_deseqResults.Rdata')
ll

m=merge(data.frame(resTEST), data.frame(res), by = 0)
m %>% 
  ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_point()


