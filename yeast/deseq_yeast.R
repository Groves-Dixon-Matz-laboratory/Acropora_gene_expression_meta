#!/usr/bin/env Rscript

#load the data
library(DESeq2)
library(tidyverse)

#format counts
counts = read.table('yeast/yeast_cds_counts.tsv')
samples=colnames(counts)
length(samples)

#build coldata
cdat0 = read_tsv('metadata/project_specific_tables/PRJEB10946_yeast_sraRunTable.txt') %>% 
  mutate(Treat=sub('Lachancea kluyveri (CBS3082_a), ', '', Title, fixed=TRUE),
         Treat = if_else(Treat=='glucose 2%',
                         'control',
                         Treat)) %>% 
  select(Run, Treat)
t=sapply(cdat0$Treat, function(x) strsplit(x, ' ')[[1]][1])
coldata = cdat0 %>% 
  mutate(Treat=t,
         Treat = if_else(Treat=='low',
                         'lowTemp',
                         Treat),
         Treat = if_else(Treat=='high',
                         'highTemp',
                         Treat),
         Stress=if_else(Treat=='control',
                        'control',
                        'stress')) %>% 
  data.frame()
rownames(coldata)=coldata$Run
coldata=coldata[samples,]
sum(coldata$Run==colnames(counts))==ncol(counts)


# run for full dataset ----------------------------------------------------

dds<-DESeqDataSetFromMatrix(counts,
                            colData = coldata,
                            design = ~Treat)
dds <- DESeq(dds)

#GATHER THE LOG2 FC FOR EACH STRESS

#get the stressors
stressors = coldata %>% 
  filter(Treat != 'control') %>% 
  pull(Treat) %>% 
  unique()

#function to get results with contrast for each stressor
get_res = function(x){
  print(paste(x,'...',sep=''))
  results(dds, contrast = c('Treat', s, 'control'))
}

#run through list
resList=list()
for (s in stressors){
  print(paste(s,'...',sep=''))
  resList[[s]]=results(dds, contrast = c('Treat', s, 'control'))
}


# run final time for all stress vs control --------------------------------
ddss<-DESeqDataSetFromMatrix(counts,
                            colData = coldata,
                            design = ~Stress)
ddss <- DESeq(ddss)
ress=results(ddss, contrast=c('Stress', 'stress', 'control'))
resList[['All stress']]=ress
stressors = append(stressors, 'All stress')

save(stressors, resList, file='yeast/deseqResults.Rdata')

# plot volcanos -----------------------------------------------------------
library(cowplot)

#funcituon to build volcano
svolcano = function(s){
  r=resList[[s]]
  data.frame(r) %>% 
    mutate(sig=padj<0.1) %>% 
    ggplot(aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
    geom_point(alpha=0.5) +
    scale_color_manual(values=c('black', 'red')) +
    labs(subtitle=s) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position='none') +
    lims(x=c(-5,5), y=c(0,60))
    # labs(x=bquote(log[2]~'Fold Difference'), y=bquote('-log'[10]~'p'))
}

#run for each stressor
vlist = lapply(stressors, function(s) svolcano(s))
plts = plot_grid(plotlist=vlist)


#set up labels and legend
leg <- cowplot::get_legend(vlist[[1]]+theme(legend.position='bottom'))
legplt = plot_grid(leg)


#set up labels
ylab = ggdraw() + draw_label(bquote('-log'[10]~'p-value'), angle=90)
xlab = plot_grid(ggdraw() + draw_label(bquote(log[2]~'Fold Difference')))

#plot
top1 = plot_grid(ylab, plts, nrow=1, rel_widths=c(0.03, 1))
top2 = plot_grid(top1, xlab, nrow=2, rel_heights=c(1, 0.05))
final= plot_grid(top2, legplt, nrow=2, rel_heights=c(1, 0.05))
final




# get stress-stress correlations ------------------------------------------
library(pheatmap)
#assemble into single dataframe
extract_logfcs = function(s){
  r=resList[[s]] %>% 
    data.frame() %>% 
    select(log2FoldChange)
  colnames(r)=s
  r$gene=rownames(r)
  return(r)
}
extract_ps = function(s){
  r=resList[[s]] %>% 
    data.frame() %>% 
    select(pvalue)
  colnames(r)=s
  r$gene=rownames(r)
  return(r)
}

lfcList = lapply(stressors, function(s) extract_logfcs(s))
pList = lapply(stressors, function(s) extract_ps(s))


#organize log2 fold changes
ldat = lfcList %>% 
  purrr::reduce(full_join, by='gene')
colnames(ldat)
lcs = ldat %>% 
  select(-`All stress`, -gene) %>% 
  as_tibble()
head(lcs)

c=cor(lcs, use="pairwise.complete.obs")
c
pheatmap(c)
plot(density(c))

#look at some scatterplots
head(lcs)

xcol='glycerol'
ycol='highTemp'

pltList=list()
for (xcol in stressors[1:3]){
  for (ycol in stressors[4:6]){
    pair=paste(xcol,ycol,sep='-')
    plt=lcs %>% 
      ggplot(aes_string(x=xcol, y=ycol)) +
      geom_point(alpha=0.5) +
      geom_smooth(method='lm')
    pltList[[pair]]=plt
  }
}
plot_grid(plotlist=pltList)

lcs %>% 
  ggplot(aes_string(x=xcol, y=ycol)) +
  geom_point(alpha=0.5)

lcs %>% 
  gather(key='stressType', value='lfc') %>% 
  ggplot(aes) +
    geom_line() +
    facet_wrap(facets = vars(genus))


#organize pvalues
pdat = pList %>% 
  purrr::reduce(full_join, by='gene')
colnames(pdat)
ps = pdat %>% 
  select(-`All stress`, -gene)
head(ps)
sigCount = apply(ps, 1, function(x) sum(x<0.1))
table(sigCount)



# get variance stabilized -------------------------------------------------

rld = vst(dds)
rld.df = data.frame(assay(rld))
colnames(rld.df)==coldata$Run
pheatmap(cor(rld.df), labels_row=coldata$Treat)


# output for mwu ----------------------------------------------------------

source('figurePlotting/rna_functions.R')

for (s in stressors){
  print(paste(s,'...',sep=''))
  df=data.frame(resList[[s]])
  df$gene=rownames(df)
  outName=paste(s, 'ForMWU.csv', sep='_')
  write_out_go(df, paste('~/Desktop/kog_mwu', outName, sep='/'))
}







