#plot_all_stress.R


library(tidyverse)
library(cowplot)
library(ggplot2)
library(DESeq2)
source('./figurePlotting/rna_functions.R')


# load the normalized counts ----------------------------------------------

input = 'largeIgnored/stress_project_controlled.Rdata'
ll=load(input)
ll
rld.df = t(datExpr)
coldata=datTraits
coldata$my_letter = substr(coldata$my_title, start=1, stop=1)


# modify treat names and get summary table -------------------------------------------------------

#modify treatment names
coldata = coldata %>% 
  mutate(Treatment = if_else(treat=='cold_then_heat',
                             'multiple',
                             treat),
         Treatment = if_else(treat=='heat_then_cold',
                             'multiple',
                             Treatment),
         Treatment = if_else(Treatment=='challenge',
                             'immune challenge',
                             Treatment),
         Treatment = if_else(Treatment=='low_pH',
                             'low pH',
                             Treatment),
         Treatment = if_else(Treatment=='low_salinity',
                             'hyposalinity',
                             Treatment))

#make ordered factors for plotting downstream
unique(coldata$Treatment)
treatLevels = c('control', 'heat', 'cold', 'hyposalinity', 'immune challenge', 'low pH', 'multiple')
coldata$Treatment = factor(coldata$Treatment, levels=treatLevels)


#get summary table
sumTab=coldata %>% 
  group_by(Treatment) %>% 
  summarize(`N projects` = length(unique(my_title)),
            `N samples` = n())

sumTab


# plot volcano ------------------------------------------------------------

ll=load('deseqResults/stress_deseqResults.Rdata')
volc=plot_volcano(data.frame(res), TITLE=NULL) 

# plot PCA ----------------------------------------------------------------

fraction = 1/10
NTOP = round(fraction * nrow(rld.df), digits=0)
point.size=1
treat.df = plotStressPCA(df = rld.df, coldat = coldata, intgroup = 'stress', ntop=NTOP, main = NULL, SIZE=point.size, returnData = T) #with modified one
treatPCA = plotStressPCA(df = rld.df, coldat = coldata, intgroup = 'stress', ntop=NTOP, main = NULL, SIZE=point.size, returnData = F) #with modified one


# RUN DAPC --------------------------------------------------------------------

######## PREPARE DATA FOR DAPC ########
library(adegenet)
pcp=prcomp(t(rld.df), retx=TRUE, center=TRUE, scale=TRUE)
scores=pcp$x
screeplot(pcp,bstick=T) # only 1st PC is non-random when using diff exp genes only; higher PCs non-random when using all data...

# adegenet: finding clusters (even though we know what clusters we want) - choose 4 PCs and 2 groups
clus=find.clusters(t(rld.df),max.n.clus=15, n.clust=2, n.pca=10) #[degs10,]

#Use clus$grp to rename
clus$grp=coldata$stress #set the transplantation site as the groups
cdf = data.frame(clus$grp)
head(cdf)
cdf$Run = coldata$Run

#run dapc using the assigned groups
dp=dapc(t(rld.df),clus$grp, n.da=1, perc.pca=80)
dpf = data.frame(dp$ind.coord)
dpf$Run = rownames(dpf)

###### CORRELATE DAPC AND PCA
head(coldata)
treat.df$Run = rownames(treat.df)

#merge up the PCA and DAPC data
pdat = dpf %>% 
  full_join(treat.df, by = 'Run') %>% 
  full_join(coldata, by = 'Run')


# plot stress prediction results -------------------------------------------
ll=load('./category_prediction/stratified_allStress_train.csv__stratified_allStress_test.csv_predictions_predictionResults.Rdata')
library(caret)
head(cdat)
head(lres)
head(rfres)

lpcts = get_acc_pct(lres, 'logistic')
rfpcts = get_acc_pct(rfres, 'RF')


marginL = -0.1
lplt = lpcts %>% 
  ggplot(aes(x=Reference, y=Freq, fill=Agree)) +
  geom_bar(stat='identity') +
  theme(legend.position='none',
        plot.margin=margin(l=marginL,unit="cm")) +
  labs(x='Logistic reg.')
rfplt = rfpcts %>% 
  ggplot(aes(x=Reference, y=Freq, fill=Agree)) +
  geom_bar(stat='identity') +
  theme(legend.title=element_blank(),
        axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin=margin(l=marginL,unit="cm")) +
  labs(x='Rand. forest')
predBp = plot_grid(lplt,rfplt, nrow=1, rel_widths=c(.7, 1))

# final replots -----------------------------------------------------------

#set up table object
library(gridExtra)
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(sumTab, rows=NULL, theme=ttheme_minimal())

#volcano
volc=plot_volcano(data.frame(res), TITLE=NULL) 

#pca
mcoldata=coldata %>% 
  mutate(stress = if_else(stress=='control','C','S'))
treatPCA = plotStressPCA(df = rld.df, coldat = mcoldata, intgroup = 'stress', ntop=NTOP, main = NULL, SIZE=point.size, returnData = F) #with modified one


#pca density plots
pdens = treat.df %>% 
  mutate(stress = if_else(stress=='control',
                            'C',
                            'S')) %>% 
  ggplot(aes(x=PC1*-1, fill=stress)) +
  geom_density(alpha=0.5) +
  theme(legend.position='top') +
  labs(fill=NULL, x='PC1')

#dapc density
ddens = pdat %>% 
  ggplot(aes(x=LD1, fill=stress.x)) +
  geom_density(alpha=0.5) +
  labs(fill=NULL, x='DAPC') +
  theme(legend.position='none')

#correlation
pdCor = pdat %>% 
  ggplot(aes(x=PC1*-1, y=LD1)) +
  geom_smooth(span=0.3, se=FALSE, lwd=1) +
  geom_point(aes(color=Treatment), size=1) +
  labs(x='PC1', y='DAPC', color=NULL)

plot_grid(tbl, volc, treatPCA,  pdens, pdCor, ddens, nrow=3, rel_widths = c(1, 0.5), labels=LETTERS[1:6])
plot_grid(volc, tbl,  pdens, treatPCA, ddens, pdCor, nrow=3, rel_widths = c(0.5, 1), labels=LETTERS[1:6])

#new one with prediction barplot (don't like this as much)
plot_grid(tbl, volc, treatPCA,  ddens, pdCor, predBp, nrow=3, rel_widths = c(1, 0.5), labels=LETTERS[1:6])


# check dapc against stress intensity -------------------------------------

#load detailed heat trait table
heat = read_csv('./metadata/detailed_tables/detailedHeat.csv') %>% 
  left_join(pdat, by='Run') %>% 
  select(my_title.x, Run, LD1, PC1, temp, degOverAmb, hoursExposed, treat.x)


heat %>% 
  filter(temp != 'not_clear') %>%
  mutate(temp=as.numeric(temp)) %>% 
  ggplot(aes(x=temp, y=LD1, color=my_title.x)) +
  geom_point()


heat %>% 
  ggplot(aes(x=degOverAmb, y=LD1, color=my_title.x)) +
  geom_point()

heat %>% 
  mutate(cumExpose = degOverAmb*hoursExposed) %>% 
  ggplot(aes(x=cumExpose, y=LD1, color=my_title.x)) +
  geom_point()

