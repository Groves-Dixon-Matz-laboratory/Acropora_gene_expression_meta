library(gtools)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())


#plot a scatterplot
library(gtools)
plot_scatter = function(df, xcol, ycol, addStars = FALSE){
  bad = c(NA, NaN, Inf, -Inf)
  badx = df[,xcol] %in% bad
  bady = df[,ycol] %in% bad
  badr = badx | bady
  subdf = df[,c(ycol, xcol)] %>% 
    filter(!badr)
  lm1 = lm(subdf[,ycol] ~ subdf[,xcol])
  r2 = round(summary(lm1)$r.squared, digits=2)
  pvalue = summary(lm1)$coefficients[2,4]
  df %>% 
    ggplot(aes_string(x=xcol, y=ycol)) +
    geom_point(alpha=ALPHA) +
    labs(subtitle=bquote(R^2*' = '*.(r2)))
  if (addStars==TRUE){
    print('Adding stars for significance')
    starString=stars.pval(pvalue)[1]
    print(stars.pval(pvalue))
    if (starString %in% c('.', ' ')){
      starString=''
    }
    plt=df %>% 
      ggplot(aes_string(x=xcol, y=ycol)) +
      geom_point(alpha=ALPHA) +
      labs(subtitle=bquote(R^2*' = '*.(r2)*.(starString)))
  } else{
    plt=df %>% 
      ggplot(aes_string(x=xcol, y=ycol)) +
      geom_point(alpha=ALPHA) +
      labs(subtitle=bquote(R^2*' = '*.(r2)))
  }
  return(plt)
  
}


#ggplot version of pairs()
my_pairs = function(df){
  print(colnames(df))
  pltList = list()
  for (xcol in colnames(df)){
    for (ycol in colnames(df)){
      pair=paste(xcol, ycol, sep='-')
      print(pair)
      if (xcol==ycol){
        plt= ggdraw() + draw_label(xcol)
      } else {
        plt = plot_scatter(df, xcol, ycol)
      }
      pltList[[pair]]=plt
    }
  }
  return(pltList)
}


#get ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#function to swap project names

swap_proj_names = function(df, modNamesPath='./metadata/detailed_tables/stressNamesModified.txt'){
  mnames = read_tsv(modNamesPath)
  mdf = df %>% 
    left_join(mnames, by = 'my_title')
  newNames = paste(mdf$Bioproject, mdf$ref)
  return(newNames)
}


swap_proj_for_author_Bioproj= function(my_title_vector, modNamesPath='./metadata/detailed_tables/stressNamesModified.txt'){
  mnames = read.table(modNamesPath, sep='\t', header = TRUE)
  rownames(mnames)=mnames$my_title
  ordered=mnames[my_title_vector,]
  newNames = paste(ordered$author, ordered$Bioproject, sep='_')
  return(newNames)
}


#function to read in the res files
read_deseq_res = function(x, n){
  print(paste(n, '...', sep=''))
  load(x)
  d=data.frame(res[,c('log2FoldChange', 'pvalue')])
  colnames(d) = c('lfc', 'p')
  colnames(d) = paste(n, colnames(d), sep='_')
  d$gene = rownames(d)
  return(d)
}


#get go genes
get_go_genes = function(lgo, goterm){
  lgo %>% 
    filter(go==goterm) %>% 
    pull(gene) %>% 
    unique()
}


#function to read in GO results made in standard two-tailed way for delta-rank correlations
read_in_two_tailed_go_results = function(inputName){
  infile = paste(paste('./go_mwu/MWU', goDivision, sep = "_"), paste(inputName,'For_MWU.csv',sep='_'), sep = "_")
  res = read.table(infile, header = TRUE, stringsAsFactors=FALSE)
  res$inputName = inputName
  return(res)
}



#Funciton to convert a vector of gene names to one of annotations based on Amillepora_euk.emapper.annotations.tsv
convert_to_annot = function(geneVector, annotDfFile='results_tables/Amillepora_euk.emapper.annotations.tsv', geneCol='query_name', annotToSwapCol='eggNOG_annot'){
  #THIS DIDNT SEEM TO BE WORKING RIGHT
}


#write out mwu input from deseq res 
write_out_go = function(deseqRes, outPath){
  df=data.frame(deseqRes)
  df$gene=rownames(df)
  go = df %>% 
    mutate(upregulated = log2FoldChange>0,
           logp= if_else(upregulated,
                         -log(pvalue, 10),
                         -log(pvalue, 10)*-1)) %>%
    select(gene,logp)
  out=paste(outPath, 'For_MWU.csv', sep='_')
  go %>% 
    write_csv(out)
}



#write out mwu input for up and downreulgated genes from deseq res objects
write_out_up_and_down_go = function(deseqRes, outPrefix){
  print('---')
  print('writing up and down mwu inputs')
  df=data.frame(deseqRes)
  df$gene=rownames(df)
  df$logp = -log(df$pvalue, 10)
  upOut = paste(outPrefix, 'up_For_MWU.csv', sep='_')
  downOut = paste(outPrefix, 'down_For_MWU.csv', sep='_')
  print(paste('Writing upregulated genes for MWU to', upOut))
  print(paste('Writing downregulated genes for MWU to', downOut))
  up = df %>% 
    filter(log2FoldChange>0)
  down = df %>% 
    filter(log2FoldChange<0)
  up %>% 
    select(gene, logp) %>% 
    write_csv(path=upOut)
  down %>% 
    select(gene, logp) %>% 
    write_csv(path=downOut)
}


plot_volcano = function(df, TITLE){
  sdf = df %>% 
    arrange(pvalue) %>% 
    mutate(sig = !is.na(padj) & padj < 0.1) %>% 
    as_tibble()
  sdf
  
  g=ggplot(data=sdf, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
    geom_point(alpha=0.1) +
    scale_color_manual(values=c('black', 'red')) + 
    labs(title=TITLE, x=bquote(log[2]~difference), y=bquote("-"*log[10]~"p")) +
    theme(legend.position='none') +
    lims(x=c(-3,3))
  return(g)
}

plot_multicolor_volcano = function(df, TITLE, overlayGenes1=FALSE, overlayGenes2=FALSE){
  sdf = df %>% 
    arrange(pvalue) %>% 
    mutate(sig = !is.na(padj) & padj < 0.1,
           color = if_else(sig,
                           'red',
                           'black'),
           color = if_else(gene %in% overlayGenes1,
                           'blue',
                           color),
           color = if_else(gene %in% overlayGenes2,
                           'darkgreen',
                           color)) %>% 
    as_tibble()
  sdf
  
  g=ggplot(data=sdf, aes(x=log2FoldChange, y=-log(pvalue, 10))) +
    geom_point(color=sdf$color) +
    labs(title=TITLE, x=NULL, y=NULL) +
    theme(legend.position='none') +
    lims(x=c(-3,3))
  return(g)
}



plot_overlay_volcano = function(df, TITLE, overlayGenes1=FALSE, overlayGenes2=FALSE){
  sdf = df %>% 
    arrange(pvalue) %>% 
    mutate(sig = !is.na(padj) & padj < 0.1,
           color = if_else(sig,
                           'red',
                           'black'),
           color = if_else(gene %in% overlayGenes1,
                           'blue',
                           color),
           color = if_else(gene %in% overlayGenes2,
                           'darkgreen',
                           color)) %>% 
    as_tibble()
  sdf
  
  g=ggplot(data=sdf, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
    geom_point(alpha=0.5) +
    scale_color_manual(values=c(backgroundPtCol, 'red')) + 
    labs(title=TITLE, x=NULL, y=NULL) +
    theme(legend.position='none') +
    scale_x_continuous(limits=SCATTER_LIMX, breaks=X_BREAKS)
  if (overlayGenes1[1] != FALSE){
    bdf = sdf %>% 
      filter(gene %in% overlayGenes1)
    g=g+geom_point(data=bdf, aes(x=log2FoldChange, y=-log(pvalue, 10)), color=foregroundPtCol, alpha=0.2)
  }
  if (overlayGenes2[1] != FALSE){
    bdf = sdf %>% 
      filter(gene %in% overlayGenes2)
    g=g+geom_point(data=bdf, aes(x=log2FoldChange, y=-log(pvalue, 10)), color='darkgreen')
  }
  return(g)
}





plotStressPCA <- function (df, coldat, intgroup = "condition", ntop = 25000, returnData = F, pcs = 2, pc1 = 1, pc2 = 2, main = "\n", SIZE = 2, legendTitle=NULL, xInvert=1) 
{
  rv <- rowVars(df)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  intgroup.df <- as.data.frame(coldat[, intgroup, 
                                      drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    coldat[[intgroup]]
  }
  d <- data.frame(pca$x[,1:pcs], group = group, 
                  intgroup.df, name = colnames(df))
  
  d$my_letter = coldat$my_letter
  attr(d, "percentVar") <- percentVar[1:2]
  #Invert X if needed
  d[,paste('PC', pc1, sep = '')]=d[,paste('PC', pc1, sep = '')]*xInvert
  g = ggplot(data = d, 
             aes_string(x = paste('PC', pc1, sep = ''),
                        y = paste('PC', pc2, sep = ''), color = "group")) + 
    geom_point(size = SIZE) +
    xlab(paste0(paste0(paste0("PC", pc1), ": "), 
                round(percentVar[pc1] * 100), "%")) + 
    ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 100), "%")) + 
    coord_fixed() 
  g = g + labs(subtitle=main, color=legendTitle)
  g = g + theme(legend.position='none',
                axis.ticks.x=element_blank(),
                axis.ticks.y=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank())
  
  #add labels
  # g=g+geom_text(aes(label=my_letter),hjust=1, vjust=0, show.legend = FALSE)
  
  # g = g + theme_bw()
  print(g)
  if (returnData == T){
    return(d)
  }
  else{
    return(g)
  }
}

plotStressPCAwShape <- function (df, coldat, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n", SIZE = 2, SHAPE='my_stage', xInvert=1) 
{
  rv <- rowVars(df)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  intgroup.df <- as.data.frame(coldat[, intgroup, 
                                      drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    coldat[[intgroup]]
  }
  d <- data.frame(pca$x[,1:pcs], group = group, 
                  intgroup.df, name = colnames(df))
  
  d$my_letter = coldat$my_letter
  d[,SHAPE] = coldat[,SHAPE]
  
  attr(d, "percentVar") <- percentVar[1:2]
  #Invert X if needed
  d[,paste('PC', pc1, sep = '')]=d[,paste('PC', pc1, sep = '')]*xInvert
  g = ggplot(data = d, 
             aes_string(x = paste('PC', pc1, sep = ''),
                        y = paste('PC', pc2, sep = ''), 
                        color = "group",
                        shape=SHAPE)) + 
    geom_point(size = SIZE) + 
    xlab(paste0(paste0(paste0("PC", pc1), ": "), 
                round(percentVar[pc1] * 100), "%")) + 
    ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 100), "%")) + 
    coord_fixed()
  g = g + labs(subtitle=main)
  
  #add labels
  # g=g+geom_text(aes(label=my_letter),hjust=0, vjust=0, show.legend = FALSE)
  
  # g = g + theme_bw()
  print(g)
  if (returnData == T){
    return(d)
  }
  else{
    return(g)
  }
}



#modified version of the function provided in the DESeq package that can be run from a dataframe
#instead of DESeqTransform object output from
plotProjectPCA <- function (df, coldat, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n", SIZE = 5) 
{
  rv <- rowVars(df)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  intgroup.df <- as.data.frame(coldat[, intgroup, 
                                      drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    coldat[[intgroup]]
  }
  d <- data.frame(pca$x[,1:pcs], group = group, 
                  intgroup.df, name = colnames(df))
  
  d$my_letter = coldat$my_letter
  d$PC1=d$PC1*-1
  
  attr(d, "percentVar") <- percentVar[1:2]
  g = ggplot(data = d, 
             aes_string(x = paste('PC', pc1, sep = ''),
                        y = paste('PC', pc2, sep = ''), color = "group")) + 
    geom_point(size=SIZE) +
    xlab(paste0(paste0(paste0("PC", pc1), ": "), 
                round(percentVar[pc1] * 100), "%")) + 
    ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 100), "%")) + 
    coord_fixed()
  g = g + labs(subtitle=main, color=NULL) 
  g = g + theme(axis.ticks.x=element_blank(),
                axis.ticks.y=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank())
  
  #add labels
  # g=g+geom_text(aes(label=my_letter),hjust=0, vjust=0, show.legend = FALSE)
  
  # g = g + theme_bw()
  print(g)
  if (returnData == T){
    return(d)
  }
  else{
    return(g)
  }
}


#Function to load deseq data in plot_deseq_scatterplots.R
load_deseq = function(path, tag){
  ll=load(path)
  df=as_tibble(data.frame(res)) %>% 
    mutate(gene=rownames(res)) %>% 
    select(gene,log2FoldChange,pvalue,padj)
  return(df)
}


#Function to load variance stabilized counts from deseq
load_vsd = function(path){
  ll=load(path)
  datTraits$my_letter = substr(datTraits$my_title, start=1, stop=1)
  return(list(t(datExpr), datTraits))
}



#SCATTERPLOT VARIABLES
SCATTER_LIMX = c(-4.5, 4.5)
SCATTER_LIMY = c(-4.5, 4.5)
X_BREAKS = c(-3, 0, 3)
Y_BREAKS = c(-3, 0, 3)
ALPHA = 0.75
backgroundPtCol = 'black'
foregroundPtCol = 'dodgerblue'
# #ggplot default
# backgroundPtCol = '#F8766D'
# foregroundPtCol = '#00BFC4'
# #blue-gold
# backgroundPtCol = '#143D59'
# foregroundPtCol = '#F4B41A'

#Function to build scatterplots of log2 fold changes for 2 dataframes
do_scatter= function(df1, df2, xlab, ylab, title){
  mdat=df1 %>% 
    left_join(df2, by='gene')
  lm1=lm(mdat$log2FoldChange.x~mdat$log2FoldChange.y)
  r2=round(summary(lm1)$r.squared, digits=2)
  mdat %>% 
    ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
    geom_point(alpha=ALPHA, color=backgroundPtCol) +
    # geom_smooth(method='lm') +
    labs(x=xlab,y=ylab,title=title,subtitle=bquote(R^2*"="*.(r2))) +
    lims(x=SCATTER_LIMX, y=SCATTER_LIMY)
}


#Function to build scatterplots of log2 fold changes for 2 dataframes
do_scatter_overlay= function(dfy, dfx, xlab, ylab, title, overlay, r2aCoords=c(2, -4), r2bCoords=c(-2, 4), writeR2=TRUE){
  mdat=dfx %>% 
    left_join(dfy, by='gene')
  lm1=lm(mdat$log2FoldChange.x~mdat$log2FoldChange.y)
  r2a=round(summary(lm1)$r.squared, digits=2)
  sub=mdat %>% 
    filter(gene %in% overlay)
  lm2=lm(sub$log2FoldChange.x~sub$log2FoldChange.y)
  r2b=round(summary(lm2)$r.squared, digits=2)
  r2aString = bquote(R^2*"="*.(r2a))
  mdat$select = mdat$gene %in% overlay
  PLOT = mdat %>% 
    mutate(select = if_else(gene %in% overlay,
           'core stress',
           'all genes')) %>% 
    ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y, color=select)) +
    geom_point(alpha=ALPHA) +
    # geom_smooth(method='lm', lwd=0.5, se=FALSE) +
    labs(x=xlab,
         y=ylab,
         subtitle=title) +
    scale_x_continuous(breaks=X_BREAKS, limits=SCATTER_LIMX) +
    scale_y_continuous(breaks=Y_BREAKS, limits=SCATTER_LIMY) +
    scale_color_manual(values=c(backgroundPtCol, foregroundPtCol)) +
    theme(legend.title=element_blank(),
          legend.position='none') 
    
    #optionally overlay the blue points for visibility (It's too much though)
    # geom_point(data=sub,
    #            aes(x=log2FoldChange.x, y=log2FoldChange.y),
    #            color=foregroundPtCol,
    #            alpha=ALPHA/2)
    
    if (writeR2==TRUE){
      PLOT=PLOT+annotate("text", x = r2aCoords[1], y = r2aCoords[2],
               label = paste('italic(R)^2 == ', r2a), parse=TRUE, color=backgroundPtCol) +
      annotate("text", x = r2bCoords[1], y = r2bCoords[2],
               label = paste('italic(R)^2 == ', r2b), parse=TRUE, color=foregroundPtCol)
    } else {
      PLOT=PLOT+annotate("text", x = r2aCoords[1], y = r2aCoords[2],
             label = paste(r2a), parse=TRUE, color=backgroundPtCol) +
      annotate("text", x = r2bCoords[1], y = r2bCoords[2],
             label = paste(r2b), parse=TRUE, color=foregroundPtCol)
    }
  return(PLOT)
}

getLog2R2 = function(dfy, dfx){
  mdat=dfx %>% 
    left_join(dfy, by='gene')
  lm1=lm(mdat$log2FoldChange.x~mdat$log2FoldChange.y)
  r2a=round(summary(lm1)$r.squared, digits=2)
  return(r2a)
}

#Function to convert a contingency table from confusion matrix into plotable accuracy percentages
get_acc_pct = function(df, method){
  cm=confusionMatrix(data=df$pred,
                     reference=df$obs,
                     positive='stressed')
  tab=cm$table
  pcts = tab %>% 
    sweep(2,colSums(tab),`/`) %>% 
    data.frame() %>% 
    mutate(Agree=if_else(Prediction==Reference,
                         'Agree',
                         'Disagree'),
           Method = method,
           Reference = if_else(Reference=='stressed',
                                'S',
                                'C'))
  return(list(pcts, cm))
}


load_pred_stats = function(fileName){
  load(fileName)
  lcmdat = get_acc_pct(lres, 'LR')
  rfcmdat = get_acc_pct(rfres, 'RF')
  lpcts = lcmdat[[1]]
  lcm = lcmdat[[2]]
  rfpcts = rfcmdat[[1]]
  rfcm = rfcmdat[[2]]
  return(list('pcts'=list(lpcts, rfpcts), 'cms'=list(lcm, rfcm)))
}


plot_pred_stats = function(datList, subtitle){
  dat = datList[["pcts"]] %>% 
    purrr::reduce(rbind)
  cmdat = datList[["cms"]]
  lpval = cmdat[[1]]$overall['AccuracyPValue']
  rfpval = cmdat[[2]]$overall['AccuracyPValue']
  lpsymb = as.character(stars.pval(lpval)[1])
  rfpsymb = as.character(stars.pval(rfpval)[1])
  symdf = tibble(x=c(1,2), y=c(1.01,1.01), lab=c(lpsymb, rfpsymb))
  PSIZE=5
  pnudge=0
  #change capitalization to match figs
  dat$Agree[dat$Agree=='Agree']<-'agree'
  dat$Agree[dat$Agree=='Disagree']<-'disagree'
  dat$Agree = factor(dat$Agree, levels=c('disagree', 'agree'))
  dat %>% 
    group_by(Method, Agree) %>% 
    summarize(pct = sum(Freq)/2) %>% 
    ggplot() +
    geom_bar(aes(x=Method, y=pct, fill=Agree), stat='identity') + 
    geom_text(data=symdf, aes(x=x, y=y, label=lab), size=PSIZE, nudge_y=pnudge) +
    labs(y='Freq', subtitle=subtitle) +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position='none') +
    scale_y_continuous(breaks=c(0,0.5, 1), limits=c(0, 1.1))
}



plot_2bar_pred_stats = function(datList){
  lpcts=datList[[1]]
  rfpcts=datList[[2]]
  marginL = -0.1
  lplt = lpcts %>% 
    ggplot(aes(x=Reference, y=Freq, fill=Agree)) +
    geom_bar(stat='identity') +
    theme(legend.position='none',
          plot.margin=margin(l=marginL,unit="cm")) +
    labs(x='LR')
  rfplt = rfpcts %>% 
    ggplot(aes(x=Reference, y=Freq, fill=Agree)) +
    geom_bar(stat='identity') +
    theme(legend.title=element_blank(),
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin=margin(l=marginL,unit="cm")) +
    labs(x='RF')
  predBp = plot_grid(lplt,rfplt, nrow=1, rel_widths=c(.75, 1))
  return(predBp)
}




single_label = function(LAB, SIZE, ANGLE){
  ggplot() +
    geom_text(aes(label=LAB, x=0,y=0), size=SIZE, angle=ANGLE) + 
    theme(axis.line.x=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin=margin(l=-0.1,unit="cm")) +
    lims(x=c(-.1,.1), y=c(-.1, .1))
}



#used in plot_rf_acc_for_specifics to plot accuracy barplots for specifics multipanel
extract_stats = function(cms){
  l=cms[[1]]$overall
  r=cms[[2]]$overall
  la = l['Accuracy']
  lp = l['AccuracyPValue']
  ra = r['Accuracy']
  rp = r['AccuracyPValue']
  return(list('la'=la,
              'lp'=lp,
              'ra'=ra,
              'rp'=rp))
}

plot_rf_acc_for_specifics = function(predDat1, predDat2, name1, name2, SUBTITLE=NULL){
  cms1=load_pred_stats(predDat1)[[2]]
  cms2=load_pred_stats(predDat2)[[2]]
  stats1 = extract_stats(cms1)
  stats2 = extract_stats(cms2)
  pvals = c(stats1$lp, stats1$rp, stats2$lp, stats2$rp)
  psymbols = sapply(pvals, function(x) as.character(stars.pval(x)))
  psymbols[psymbols=='.']<- ''
  dat = data.frame('props'=c(stats1$la,
                             stats1$ra,
                             stats2$la,
                             stats2$ra),
                   'Method'=c('LR', 'RF', 'LR', 'RF'),
                   'Train'=c(name1, name1, name2, name2),
                   pvalues = pvals,
                   psymbols = psymbols,
                   xsym = c(1,1,2,2),
                   ysym = rep(1.01, 4)
  ) %>% 
    mutate(agree = props,
           disagree = 1-agree)
  print(dat)
  dat %>% 
    gather(agree, Acc, c('agree', 'disagree')) %>% 
    filter(Method=='RF') %>% 
    mutate(Train=factor(Train, levels=c(name1, name2)),
           agree=factor(agree, levels=c('disagree', 'agree'))) %>% 
    ggplot() +
    geom_bar(aes(x=Train, y=Acc, fill=agree), stat='identity') +
    geom_text(aes(x=xsym, y=ysym, label=psymbols), size=5, nudge_y=0)+
    theme(legend.position='none') +
    labs(x=NULL, y=NULL, subtitle=SUBTITLE) +
    scale_y_continuous(breaks=c(0,0.5, 1), limits=c(0, 1.1)) 
    
}



get_pred_dat = function(predDat1, predDat2, name1, name2){
  cms1=load_pred_stats(predDat1)[[2]]
  cms2=load_pred_stats(predDat2)[[2]]
  stats1 = extract_stats(cms1)
  stats2 = extract_stats(cms2)
  pvals = c(stats1$lp, stats1$rp, stats2$lp, stats2$rp)
  psymbols = sapply(pvals, function(x) as.character(stars.pval(x)))
  psymbols[psymbols=='.']<- '+'
  dat = data.frame('props'=c(stats1$la,
                             stats1$ra,
                             stats2$la,
                             stats2$ra),
                   'Method'=c('LR', 'RF', 'LR', 'RF'),
                   'Train'=c(name1, name1, name2, name2),
                   pvalues = pvals,
                   psymbols = psymbols,
                   xsym = c(1,1,2,2),
                   ysym = rep(1.01, 4))
  return(dat)
}

#funciton to add Library Layout data from original sra table
add_single_or_pe = function(df){
  adat = read_csv('all_acropora_sra_runInfoTable.csv') %>% 
    select(Run, LibraryLayout)
  df %>% 
    left_join(adat, by = 'Run')
}







#function for organizing annotation coordinates for GO correlation plot
#(used in core_stress_go_correlations.R)
build_annotation_coords = function(selection,
                                   goData,
                                   sigOnly,
                                   xlabel_left_bound,
                                   xlabel_right_bound,
                                   ylabel_bottom_bound,
                                   ylabel_top_bound,
                                   stickAdd=0,
                                   pointTop=FALSE,
                                   center=FALSE){
  if (sigOnly){
    goData = goData %>% 
      filter(p.adj.y < 0.1)
  }
  xcoords = selection %>% 
    left_join(goData, by='name') %>% 
    filter(!is.na(inputName)) %>% 
    group_by(summary) %>% 
    summarize(mnDelta.x=mean(delta.rank.x, na.rm=TRUE))
  sub = selection %>% 
    left_join(goData, by='name') %>% 
    filter(!is.na(inputName)) %>% 
    group_by(inputName, summary) %>% 
    summarize(mnDelta.y=mean(delta.rank.y, na.rm=TRUE)) %>% 
    left_join(xcoords, by = 'summary') %>% 
    left_join(mnames, by = 'inputName') %>% 
    mutate(treatShape=as.numeric(treat2) + 20)
  more.sub = sub %>% 
    group_by(summary) %>% 
    summarize(text.x=mean(mnDelta.x),
              text.y=mean(mnDelta.y),
              segBottom = min(mnDelta.y-stickAdd),
              segTop = max(mnDelta.y)+stickAdd)
  if (pointTop){
    more.sub$pointerYend=more.sub$segTop
  } else {
    more.sub$pointerYend=more.sub$segBottom
  }
  if (center){
    more.sub$hjust=0.5
  } else {
    more.sub$hjust=0
  }
  xpositions = seq(xlabel_left_bound, xlabel_right_bound, length.out = nrow(more.sub))
  ypositions = seq(ylabel_bottom_bound, ylabel_top_bound, length.out = nrow(more.sub))
  more.sub = more.sub %>% 
    arrange(text.x) %>% 
    mutate(text.y.leveled = ypositions,
           text.x.shifted = xpositions)
  return(list(xcoords,
              sub,
              more.sub))
}


add_custom_go_overlay = function(original, sub, more.sub, hjust=0){
  original +
    #add the pointers from labels to center of each custom go group
    annotate("segment",
             x = more.sub$text.x,
             xend = more.sub$text.x.shifted,
             y = more.sub$pointerYend,
             yend = more.sub$text.y.leveled,
             colour = "grey40",
             lwd=0.2) +
    #add vertical lines connecting summarized functions
    annotate("segment",
             x = more.sub$text.x,
             xend = more.sub$text.x,
             y = more.sub$segBottom,
             yend = more.sub$segTop,
             colour = "black",
             lwd=0.5) +
    #add point indicating mean for each summary grouping
    annotate("point",
             x = more.sub$text.x,
             y = more.sub$text.y,
             colour = "black",
             pch='-',
             size=12) +
    #add points for projects or stress types
    geom_point(data=sub,
               aes(x=mnDelta.x,
                   y=mnDelta.y,
                   fill=treat2,
                   shape=treat2),
               size=overlayPointSize,
               color='black') +
    #add labels
    annotate("text",
             x=more.sub$text.x.shifted,
             y=more.sub$text.y.leveled,
             label=more.sub$summary,
             hjust = more.sub$hjust)
}





add_custom_go_overlay_bioproject = function(original, sub, more.sub, text.size=5){
  original +
    #add the pointers from labels to center of each custom go group
    annotate("segment",
             x = more.sub$text.x,
             xend = more.sub$text.x.shifted,
             y = more.sub$text.y,
             yend = more.sub$text.y.leveled,
             colour = "grey",
             lwd=0.3) +
    #add vertical lines connecting summarized functions
    annotate("segment",
             x = more.sub$text.x,
             xend = more.sub$text.x,
             y = more.sub$segBottom,
             yend = more.sub$segTop,
             colour = "black",
             lwd=0.5) +
    #add point indicating mean for each summary grouping
    annotate("point",
             x = more.sub$text.x,
             y = more.sub$text.y,
             colour = "black",
             pch='-',
             size=12) +
    #add points for projects or stress types
    geom_point(data=sub,
               aes(x=mnDelta.x,
                   y=mnDelta.y,
                   fill=inputName,
                   shape=inputName),
               size=overlayPointSize,
               color='black') +
    #add labels
    annotate("text",
             x=more.sub$text.x.shifted,
             y=more.sub$text.y.leveled,
             label=more.sub$summary,
             hjust = 0,
             size=text.size)
}



#same as r2 version but for pearson correlation 
plot_scatter_pearsonCor_annotated = function(dat, xcol, ycol, xlab, ylab, ALPHA=0.1){
  dat=data.frame(dat)
  bad = c(NA, NaN, Inf, -Inf)
  x=dat[,xcol]
  y=dat[,ycol]
  rem = (x %in% bad) | (y %in% bad)
  totBad = sum(rem)
  if (totBad>0){
    print(paste('removing', totBad, 'rows with NA/NaN/Inf'))
  }
  dat=dat[!rem,]
  lm1=lm(dat[,ycol]~dat[,xcol])
  r2 = round(summary(lm1)$r.squared, digits=2)
  pearsonCor=cor(x=dat[,xcol],
                 y=dat[,ycol])
  r=round(pearsonCor, digits=2)
  print(summary(lm1))
  plt = dat %>% 
    ggplot(aes_string(x=xcol, y=ycol)) +
    geom_point(alpha = ALPHA) +
    labs(x=xlab,
         y=ylab)
  pbuild = ggplot_build(plt)
  yrange = pbuild$layout$panel_params[[1]]$y.range
  xrange = pbuild$layout$panel_params[[1]]$x.range
  plt +
    annotate("text", x = xrange[1], y = yrange[2],
             label = paste('italic(r) ==', r), parse=TRUE, color='black',
             hjust=0)
}
