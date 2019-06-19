
write_out_go = function(df, outPath){
  go = df %>% 
    mutate(upregulated = log2FoldChange>0,
           logp= if_else(upregulated,
                         -log(pvalue, 10),
                         -log(pvalue, 10)*-1)) %>%
    select(gene,logp)
  go %>% 
    write_csv(outPath)
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
    scale_color_manual(values=c('black', 'red')) + 
    labs(title=TITLE, x=NULL, y=NULL) +
    theme(legend.position='none') +
    lims(x=c(-3,3))
  if (overlayGenes1[1] != FALSE){
    bdf = sdf %>% 
      filter(gene %in% overlayGenes1)
    g=g+geom_point(data=bdf, aes(x=log2FoldChange, y=-log(pvalue, 10)), color='blue', alpha=0.2)
  }
  if (overlayGenes2[1] != FALSE){
    bdf = sdf %>% 
      filter(gene %in% overlayGenes2)
    g=g+geom_point(data=bdf, aes(x=log2FoldChange, y=-log(pvalue, 10)), color='darkgreen')
  }
  return(g)
}





plotStressPCA <- function (df, coldat, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n", SIZE = 2, legendTitle=NULL, xInvert=1) 
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
  
  #add labels
  g=g+geom_text(aes(label=my_letter),hjust=1, vjust=0, show.legend = FALSE)
  
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



#Function to build scatterplots of log2 fold changes for 2 dataframes
do_scatter= function(df1, df2, xlab, ylab, title){
  mdat=df1 %>% 
    left_join(df2, by='gene')
  lm1=lm(mdat$log2FoldChange.x~mdat$log2FoldChange.y)
  r2=round(summary(lm1)$r.squared, digits=2)
  mdat %>% 
    ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
    geom_point(alpha=0.5) +
    geom_smooth(method='lm') +
    labs(x=xlab,y=ylab,title=title,subtitle=bquote(R^2*"="*.(r2))) +
    lims(x=c(-3,3), y=c(-3,3))
}


#Function to build scatterplots of log2 fold changes for 2 dataframes
do_scatter_overlay= function(df1, df2, xlab, ylab, title, overlay, r2aCoords=c(2.5, -1), r2bCoords=c(2.5, 2.5)){
  mdat=df1 %>% 
    left_join(df2, by='gene')
  lm1=lm(mdat$log2FoldChange.x~mdat$log2FoldChange.y)
  r2a=round(summary(lm1)$r.squared, digits=2)
  sub=mdat %>% 
    filter(gene %in% overlay)
  lm2=lm(sub$log2FoldChange.x~sub$log2FoldChange.y)
  r2b=round(summary(lm2)$r.squared, digits=2)
  r2aString = bquote(R^2*"="*.(r2a))
  mdat %>% 
    ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
    geom_point(alpha=0.5) +
    geom_smooth(method='lm', color='black', lwd=0.5, se=F) +
    labs(x=xlab,
         y=ylab,
         subtitle=title
         ) +
    lims(x=c(-3,3), y=c(-3,3)) +
    geom_point(data=sub, aes(x=log2FoldChange.x, y=log2FoldChange.y), color='blue', alpha=0.1) +
    geom_smooth(data=sub, aes(x=log2FoldChange.x, y=log2FoldChange.y), method='lm', color='blue', lwd=0.5, se=F) +
    theme(plot.title = element_text(size=10, hjust=0),
          plot.subtitle = element_text(color = "black", hjust=0, size=10)) +
    annotate("text", x = r2aCoords[1], y = r2aCoords[2],
             label = paste('italic(R) ^ 2 ==', r2a), parse=TRUE, color='black') +
    annotate("text", x = r2bCoords[1], y = r2bCoords[2],
             label = paste('italic(R) ^ 2 ==', r2b), parse=TRUE, color='blue')
    
}










