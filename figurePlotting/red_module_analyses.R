#red_module_analyses.R
#Here we do some further exploring of the red module
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
rm(list = ls())


# load red module ---------------------------------------------------------

#load the genes' module memberships
memb_dat = read_tsv('results_tables/wgcna_moduleMembership.tsv') %>% 
  dplyr::select(gene, assignment, MMred)

red_genes = memb_dat %>% 
  filter(assignment=='red') %>% 
  pull(gene)
length(red_genes)
  

# red module and log2 for full dataset ------------------------------------

#load the log2 differences from full dataset analysis
ll=load('deseqResults/stress_deseqResults.Rdata')
ll
rdat = data.frame(res) %>% 
  rownames_to_column('gene') %>% 
  mutate(red = gene %in% red_genes) %>% 
  left_join(memb_dat, by = 'gene') %>% 
  as_tibble
rdat

#show log2 fold differences by red or not red
basic_box = rdat %>% 
  ggplot(aes(x=red, y=log2FoldChange, fill=red)) +
  geom_boxplot() +
  scale_fill_manual(values=c('grey', 'red')) +
  labs(y = bquote('stress log'[2]~'fold difference'),
       fill='red module') +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
#plot(basic_box)

#function to plot correlations with red module membership
do_red_mm_scatterplot = function(dat, col1, col2){
  x=dat %>% 
    pull(col1)
  y=dat %>% 
    pull(col2)
  lm1 = lm(y~x)
  r2 = round(summary(lm1)$r.square, digits=2)
  r2_str = paste("R^2==", r2, sep='')
  
  plt=dat %>% 
    ggplot(aes_string(x=col1, y=col2)) +
    geom_point(pch=21, fill='red', color='black')
  pbuild = ggplot_build(plt)
  yrange = pbuild$layout$panel_params[[1]]$y.range
  xrange = pbuild$layout$panel_params[[1]]$x.range
  plt +
    annotate("text", x = xrange[1], y = yrange[2],
             label = r2_str, parse=TRUE, color='black',
             hjust=0)
}

#correlate log2 and module membership
red_log2_plt = rdat %>% 
  filter(log2FoldChange < 3) %>% 
  do_red_mm_scatterplot('MMred', 'log2FoldChange') +
  labs(y = bquote('stress log'[2]~'fold difference'),
       x = 'red membership')


# build log2 scatterplot for all bioprojects ------------------------------

#load log2 values from each project
ll=load('deseqResults/all_by_bioproject.Rdata')
idat = idat %>% 
  dplyr::select(!ends_with('_p')) %>% 
  pivot_longer(-gene,
               names_to = 'project',
               values_to = 'lfc') %>% 
  mutate(my_title = sub('_lfc', '', project)) %>% 
  dplyr::select(-project)
head(idat)


idat %>% 
  group_by(my_title) %>% 
  summarize(nonNA = sum(!is.na(lfc)))


#merge with membership
irdat = idat %>% 
  left_join(memb_dat, by = 'gene')

r2s = c()
iplts = list()
for (p in unique(irdat$my_title)){
  sub = irdat %>% 
    filter(my_title==p)
  lm1 = lm(sub$lfc ~ sub$MMred)
  r2 = summary(lm1)$r.square
  r2s = append(r2s, r2)
  plt = sub %>% 
    ggplot(aes(x=MMred, y=lfc)) +
    geom_point(pch=19, fill='red') +
    labs(subtitle = p)
  iplts[[p]]=plt
}
names(r2s) = unique(irdat$my_title)
#plot_grid(plotlist = iplts)



# boxplot for individual projects -----------------------------------------

#get the prettier names
name_df = read.table('./metadata/detailed_tables/stressNamesModified.txt',
                      sep='\t',
                      header=TRUE,
                      row.names='my_title',
                      stringsAsFactors=FALSE) %>% 
  rownames_to_column('my_title') %>% 
  unite(fname, Bioproject, ref, treat, sep=' ')
mod_names = name_df %>% 
  pull(fname)
names(mod_names) = name_df$my_title


#make module assignments
rbdat = irdat %>% 
  filter(my_title != "j1_thisStudy_PRJNA559404") %>% 
  mutate(redAssigned = if_else(assignment=='red',
                               'red',
                               'not_red'),
         redAssigned = if_else(is.na(redAssigned),
                               'not_red',
                               redAssigned),
         fname = mod_names[my_title])

#load the type A stress projects
ll=load('metadata/corStressProjs.Rdata')
corStressProjs

#get the median log2 fold change for red module accross projects
rmeds = rbdat %>% 
  group_by(my_title, fname, redAssigned) %>% 
  summarize(medRed = mean(lfc, na.rm=TRUE),
            nonNA = sum(!is.na(lfc))) %>% 
  filter(redAssigned=='red') %>% 
  mutate(response = if_else(my_title %in% corStressProjs,
                        'type A',
                        'type B'),
         response = if_else(grepl('beww', my_title),
                        'type A',
                        response),
         response = factor(response, levels = c('type A', 'type B'))) %>% 
  arrange(desc(medRed))

#get their order
med_order = rmeds %>% 
  pull(fname)


#red only
rbdat$fname = factor(rbdat$fname, levels = med_order)
rmeds$fname = factor(rmeds$fname, levels = med_order)
violin_plt = rbdat %>% 
  filter(redAssigned=='red') %>% 
  ggplot(aes(x=fname, y = lfc)) +
  geom_hline(yintercept = 0, lty=2) +
  geom_violin(fill='red') +
  geom_point(data=rmeds, aes(x=as.numeric(fname),
                             y=medRed,
                             shape = response),
             size=3,
             fill='black') +
  scale_shape_manual(values = c(21,25)) +
  lims(y=c(-3,3))  +
  labs(x = 'BioProject',
       y = bquote('log'[2]~'fold change stress treatment'),
       shape = '') +
  theme(legend.position = 'top',
        legend.justification = 0,
        axis.text.y = element_text(size=11)) +
  coord_flip() 
  
  

# red module gene contributions ------------------------------------------------------

#load the gene's loading values for dapc
ord_dat = read_tsv('results_tables/gene_contibutions.tsv') %>% 
  left_join(memb_dat, by = 'gene') %>% 
  mutate(PCA.pa.1 = PCA.pa.1*-1) %>% 
  as_tibble
ord_dat

#correlate pc1 loading values and module membership
red_pc1_plt = do_red_mm_scatterplot(ord_dat, 'MMred', 'PCA.pa.1') +
  labs(y = 'gene loading PC1',
       x = 'red module membership')

#correlate DAPC gene loading values and module membership
red_dapc_plt = do_red_mm_scatterplot(ord_dat, 'MMred', 'LD1') +
  labs(y = 'gene loading DAPC',
       x = 'red membership')

#non-zero lasso coefficients
red_coef_plt = ord_dat %>% 
  filter(abs(coef) > 0) %>% 
  dplyr::select(coef, MMred) %>% 
  na.omit() %>% 
  do_red_mm_scatterplot('MMred', 'coef')
  
#MeanDecreaseAccuracy from RF
red_rfImp_plt = ord_dat %>% 
  do_red_mm_scatterplot('MMred', 'MeanDecreaseAccuracy')




# get red module mean expression levels -----------------------------------

#load traits
coldata = read_csv('metadata/subset_tables/allStress_Coldata.csv') %>% 
  dplyr::select(Run, my_title, stress, treatDescription)

#swap out individual experiments for this study
beww_exps = coldata %>% 
  separate(treatDescription, into = c('exp')) %>% 
  pull(exp)
coldata = coldata %>% 
  mutate(my_title = if_else(my_title=='F_Uqueensland_ph_PRJNA269992',
                            'F_Uqueensland_ph_PRJNA269992',
                            my_title),
         my_title = if_else(my_title=="j1_thisStudy_PRJNA559404",
                            beww_exps,
                            my_title),
         my_title = if_else(my_title %in% c('Ht1', 'Ht2'),
                            'bewwHeat',
                            my_title),
         my_title = if_else(my_title == 'Ls',
                            'bewwHyposalinity',
                            my_title),
         my_title = if_else(my_title == 'Ch',
                            'bewwCold',
                            my_title),
         my_title = if_else(my_title %in% c('Tw', 'CtH', 'HtC'),
                            'bewwMulti',
                            my_title))
unique(coldata$my_title)

#load the normalized expression for all genes used in WGCNA
ll=load('largeIgnored/project_controlled_varCut0.75_expCut0.75.Rdata')
ll
rld.df = data.frame(t(datExpr)) %>% 
  rownames_to_column('gene') %>% 
  left_join(memb_dat, by = 'gene') %>% 
  as_tibble

#get the mean expression level for the red module
red_rld = rld.df %>% 
  filter(assignment=='red') %>% 
  as_tibble

#get the mean red expression
red_mns = red_rld %>% 
  dplyr::select(-gene,
                -assignment,
                -MMred) %>% 
  apply(2,mean) %>% 
  data.frame() %>% 
  set_names(c('mnRed')) %>% 
  rownames_to_column('Run') %>% 
  left_join(coldata, by = 'Run') %>% 
  as_tibble


# correlate mean red expression and sample loadings for orination axes -------

#load the pca and dapc loading values
ll=load('results_tables/pca_dapc_sample_loads.Rdata')
ll
pdat = pdat %>% 
  left_join(red_mns, by = 'Run') %>% 
  as_tibble

#plot pC1 vs mnRED
red_samplePC1_plt = do_red_mm_scatterplot(pdat, 'mnRed', 'PC1') +
  labs(y = 'sample PC1 score',
       x = 'mean expression level red module genes')
#red_samplePC1_plt


#plot pC1 vs mnRED
red_dapc_plt = do_red_mm_scatterplot(pdat, 'mnRed', 'LD1') +
  labs(y = 'sample DAPC score',
       x = 'mean expression level red module genes')



# correlation within type A -----------------------------------------------

#build a correlation matrix of mean red expression level among 
long_rld = red_rld %>% 
  dplyr::select(-gene,
                -assignment,
                -MMred) %>% 
  pivot_longer(everything(),
               names_to = 'Run',
               values_to = 'expression') %>% 
  left_join(coldata, by = 'Run')

#get mean red module expression level for stressed samples for each Bioproject
p_mns = long_rld %>% 
  filter(!is.na(my_title),
         stress=='stressed') %>% 
  group_by(my_title) %>% 
  summarize(mnRed = mean(expression))

#get mean correlation with all other bioprojects for each Bioproject
ll=load('results_tables/bioproject_stress_correlation_matrix.Rdata')
ll
ll=load('metadata/corStressProjs.Rdata')
corStressProjs
typeA = colnames(c)[colnames(c) %in% corStressProjs | grepl('beww', colnames(c))]
longc = data.frame(c) %>% 
  rownames_to_column('my_title.x') %>% 
  pivot_longer(-my_title.x,
               names_to = 'my_title.y',
               values_to = 'cor') %>% 
  mutate(type.x = if_else(my_title.x %in% typeA,
                          'A',
                          'B'),
         type.y = if_else(my_title.y %in% typeA,
                          'A',
                          'B'))


#get mean correlation with type A
mn_corA = longc %>% 
  filter(type.y =='A',
         my_title.x != my_title.y) %>% 
  group_by(my_title.x) %>% 
  summarize(mn_cor = mean(cor)) %>% 
  dplyr::rename(my_title = my_title.x) %>% 
  full_join(p_mns, by = 'my_title')

#plot pC1 vs mnRED
red_corr_plt = do_red_mm_scatterplot(mn_corA, 'mnRed', 'mn_cor') +
  labs(y = 'mean correlation Type A',
       x = 'mean red level')
#red_corr_plt



# delta ranks for red module and type A stress ----------------------------

#we have module membership for each gene
memb_dat

#load the results for the full stress dataset
typeA_delta = read.table('go_mwu/MWU_BP_corStress_For_MWU.csv',
                         header = TRUE) %>% 
  as_tibble()

#we have gene-go associations
gene_go = read_tsv('go_mwu/amil_zachFullerV2_gos.tsv') %>% 
  set_names(c('gene', 'go')) %>% 
  filter(!is.na(go))

# #make a long format version
# long_go = data.frame()
# for (i in 1:nrow(gene_go)){
#   rsub = gene_go[i,]
#   g = rsub['gene']
#   gos = rsub['go']
#   res=data.frame('gene' = g,
#                  'go' = unlist(str_split(gos, ';')))
#   long_go = rbind(long_go, res)
# }
#or just load it
ll=load('metadata/long_go.Rdata')
ll

#merge with module membership
head(long_go)
mm_go = long_go %>% 
  inner_join(memb_dat, by = 'gene') %>% 
  as_tibble


#make long format of delta ranks
long_delta = data.frame()
adf = data.frame(typeA_delta)
adf$term = as.character(adf$term)
for (i in 1:nrow(typeA_delta)){
  rsub = adf[i,c('delta.rank', 'term', 'name')]
  dr = rsub['delta.rank']
  name = rsub['name']
  gos = as.character(rsub['term'])
  split_gos = unlist(str_split(gos, ';'))
  res=data.frame('name' = name,
                 'go' = split_gos,
                 'delta.rank' = dr)
  long_delta = rbind(long_delta, res)
}

#merge that with the module membership
mm_delta = long_delta %>% 
  left_join(mm_go, by = 'go') %>% 
  as_tibble

#now average the moduel membership by go term
mn_mm = mm_delta %>% 
  group_by(name, delta.rank) %>% 
  summarize(mnRed = mean(MMred))

#plot initial
red_delta_plt = mn_mm %>% 
  ggplot(aes_string(x='mnRed', y='delta.rank')) +
  geom_point(pch=21, fill='red', color='black')

#load custom groupings from plot_go_correlations.R
ll=load('figurePlotting/selectGoGroups.Rdata')
ll
upSelect$dir = 'up'
downSelect$dir = 'down'
add_horiz = 0.2
add_vert = 500
custom = rbind(upSelect, downSelect)
overlay_df = mn_mm %>% 
  inner_join(custom, by = 'name') %>% 
  group_by(summary) %>% 
  summarize(mnDelta = mean(delta.rank),
            mnRed = mean(mnRed)) %>% 
  arrange(mnRed)

overlay_df = overlay_df %>% 
  mutate(side = if_else(mnRed > 0.13,
                        'right',
                        'left'),
         nudge_sign = if_else(side == 'left',
                              -1,
                              1),
         nudge_horz = add_horiz*nudge_sign,
         hjust = if_else(side=='left',
                         1,
                         0),
         nudge_vert = add_vert,
         nudge_vert = if_else(summary=='ROS',
                              add_vert - 100,
                              nudge_vert),
         nudge_vert = if_else(summary=='NF-KB',
                              add_vert - 300,
                              nudge_vert),
         nudge_vert = if_else(summary=='immune activition',
                              add_vert - 200,
                              nudge_vert),
         nudge_vert = if_else(summary=='membrane vesicle',
                              add_vert - 50,
                              nudge_vert))

red_delta_plt2 = red_delta_plt +
  geom_point(data = overlay_df,
             aes(x=mnRed,
                 y=mnDelta), size = 2, color = 'black', pch=19) +
  geom_text(data = overlay_df,
            aes(x=mnRed + overlay_df$nudge_horz,
                y=mnDelta + overlay_df$nudge_vert,
                label = summary), hjust = overlay_df$hjust, size=3, fontface='bold') +
  annotate("segment",
           x = overlay_df$mnRed + overlay_df$nudge_horz,
           xend = overlay_df$mnRed,
           y = overlay_df$mnDelta + overlay_df$nudge_vert,
           yend = overlay_df$mnDelta,
           colour = "black",
           lwd=0.5) +
  coord_cartesian(x=c(-0.8, 0.8)) +
  labs(y = 'all type A delta rank',
       x = "GO term's mean red membership")



#add R2 annoation
pbuild = ggplot_build(red_delta_plt2)
yrange = pbuild$layout$panel_params[[1]]$y.range
xrange = pbuild$layout$panel_params[[1]]$x.range
lm1 = lm(mn_mm$delta.rank~mn_mm$mnRed)
r2 = round(summary(lm1)$r.square, digits=2)
r2_str = paste("R^2==", r2, sep='')

red_delta_plt3 = red_delta_plt2 +
  annotate("text", x = xrange[1]+0.05, y = yrange[2],
           label = r2_str, parse=TRUE, color='black',
           hjust=0)



# replot red corr plot with shapes matching violin ------------------------

red_corr_plt = do_red_mm_scatterplot(mn_corA, 'mnRed', 'mn_cor') +
  labs(y = 'mean correlation Type A',
       x = 'mean red level')

red_corr_plt0 = mn_corA %>% 
  mutate(response = if_else(my_title %in% corStressProjs,
                            'type A',
                            'type B'),
         response = if_else(grepl('beww', my_title),
                            'type A',
                            response)) %>% 
  ggplot(aes_string(x='mnRed', y='mn_cor')) +
  geom_point(aes(shape = response),
             fill='red', color='black', size=5) +
  scale_shape_manual(values=c(21,25))
#add R2 annoation
pbuild = ggplot_build(red_corr_plt0)
yrange = pbuild$layout$panel_params[[1]]$y.range
xrange = pbuild$layout$panel_params[[1]]$x.range
lm1 = lm(mn_corA$mn_cor~mn_corA$mnRed)
r2 = round(summary(lm1)$r.square, digits=2)
r2_str = paste("R^2==", r2, sep='')
red_corr_plt = red_corr_plt0 +
  annotate("text", x = xrange[1]+0.05, y = yrange[2],
           label = r2_str, parse=TRUE, color='black',
           hjust=0) +
  theme(legend.position = 'none') +
  labs(y=c('mean correlation type A'))

# assemble final figure ---------------------------------------------------

log2_axis_lab = bquote(log[2]~"fold change stress")

small = plot_grid(red_log2_plt + labs(y=log2_axis_lab,
                                      x="gene red membership"),
                  red_corr_plt + labs(x="BioProject mean red level"),
                  nrow=1,
                  labels = c('B', 'C'))
right = plot_grid(small, 
                  plot_grid(red_delta_plt3 + 
                              labs(x='GO term mean red membership'),
                            labels=c('D')),
                  nrow=2)

plot_grid(plot_grid(violin_plt, labels=c('A')), right,
          nrow=1, rel_widths = c(1, 0.9))



# wihtout log2 ------------------------------------------------------------

thousandths = function(x){
  return(x/1e3)
}
right = plot_grid(red_corr_plt + labs(x="BioProject mean red level"), 
                  red_delta_plt3 + labs(x='GO term mean red membership',
                                        y ='type A enrichment') +
                    scale_y_continuous(labels = thousandths),
                  labels = c('B', 'C'),
                  nrow=2,
                  align='h',
                  axis = 'l')

plot_grid(plot_grid(violin_plt, labels=c('A')),
          right,
          nrow=1,
          rel_widths = c(1, 0.75))

