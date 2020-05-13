#full_dataset_pca.R


library(tidyverse)
library(cowplot)
library(ggplot2)
library(DESeq2)
source('figurePlotting/rna_functions.R')




# upload the count controlled vst counts ----------------------------------
# input = 'largeIgnored/fullDataset_vsd.Rdata' #for reference
input = 'largeIgnored/fullDataset_count_controlled.Rdata'
ll=load(input)
ll
rld.df = t(datExpr)
coldata=datTraits

#modify coldata for plotting
pullProj = function(x){
  y=strsplit(x, '_')[[1]]
  return(y[length(y)])
}
coldata$project = sapply(coldata$my_title, function(x) pullProj(x))
coldata$project[coldata$project=='EWW']<-'This study'
unique(coldata$project)
coldata$my_letter = substr(coldata$my_title, start=1, stop=1)
sampleCounts = coldata %>% 
  group_by(project, my_letter) %>% 
  summarize(N=n())
sampleCounts$forLeg = paste(sampleCounts$project, paste(sampleCounts$N, ')',sep=''), sep=' (N=')
head(sampleCounts)
coldata$forLegend = 'notAssigned'
for (i in 1:nrow(sampleCounts)){
  p=sampleCounts$project[i]
  fl = sampleCounts$forLeg[i]
  coldata = coldata %>% 
    mutate(forLegend = if_else(project==p,
                               fl,
                               forLegend))
}
coldata$my_stage[coldata$my_stage=='larval']<-'larva'
coldata$my_stage = factor(coldata$my_stage, levels=c('gamete', 'embryo', 'larva', 'recruit', 'adult'))
head(coldata)

# build PCA ---------------------------------------------------------------

fraction = 1/10
NTOP = round(fraction * nrow(rld.df), digits=0)


####### uncomment IF YOU WANT TO SEE CORRELATION WITH READ COUNTS
# pca.proj = plotProjectPCA(df = rld.df, coldat = coldata, intgroup = 'forLegend', ntop=NTOP, main = "Project", SIZE=1, returnData = T)
# lmdat = pca.proj %>% 
#   mutate(Run = rownames(pca.proj)) %>% 
#   left_join(coldata, by = 'Run')
# lmdat %>% 
#   ggplot(aes(x=log(rawCounts, 10), y=PC1)) +
#   geom_point()
# lm1=lm(log(lmdat$rawCounts, 10)~lmdat$PC1)
# summary(lm1)
##########################################

#get the project pca and legend separately
pca.proj = plotProjectPCA(df = rld.df, coldat = coldata, intgroup = 'forLegend', ntop=NTOP, main = "Project", SIZE=1, returnData = F) 
legend <- cowplot::get_legend(pca.proj + theme(legend.text=element_text(size=5)))
pca.proj.noleg = pca.proj + theme(legend.position='none')

#just for kicks, what does it look like if you remove effect of developmental stage?
contStage <- removeBatchEffect(rld.df, batch = coldata$my_stage)
pca.contStage = plotProjectPCA(df = contStage, coldat = coldata, intgroup = 'forLegend', ntop=NTOP, main = "Project", SIZE=1, returnData = F) 


#get the stage pca and density
pca.stage = plotStressPCA(df = rld.df, coldat = coldata, intgroup = 'my_stage', ntop=NTOP, main = "Developmental stage", SIZE=1, returnData = F, xInvert=-1)
stage.df = plotStressPCA(df = rld.df, coldat = coldata, intgroup = 'my_stage', ntop=NTOP, main = "Developmental stage", SIZE=1, returnData = T, xInvert=-1)

g <- ggplot_build(pca.stage)
colors = g$data[[1]] %>% 
  select(colour, group) %>% 
  filter(group %in% 3:5) %>% 
  group_by(group) %>% 
  unique() %>% 
  arrange(group) %>% 
  pull(colour)
  

stageDens = stage.df %>% 
  filter(!my_stage %in% c('embryo', 'gamete')) %>% 
  ggplot(aes(x=PC1, fill=my_stage)) +
  geom_density(alpha=0.8) +
  labs(fill=NULL, subtitle='\n', title='\n') + 
  scale_fill_manual(values=colors) +
  theme(legend.position = 'none')

#stick the point legend onto the density plot
dev.legend <- cowplot::get_legend(pca.stage + 
                                    theme(legend.position='right') +
                                    guides(colour = guide_legend(override.aes = list(size=5))))
stageDensL = plot_grid(stageDens, dev.legend, rel_widths = c(2.5, 1), nrow=1)
legend <- cowplot::get_legend(pca.proj + theme(legend.text=element_text(size=9)))
plot_grid(pca.proj.noleg, legend, pca.stage, stageDensL, nrow=2, axis='b', rel_widths=c(1, 1, 1, 0.75), labels=c('A', '', '\nB', '\nC'))



