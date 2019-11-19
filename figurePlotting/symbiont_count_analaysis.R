#symbiont_count_analysis.R
library(tidyverse)
library(cowplot)

# PREPARE DATA ------------------------------------------------------------

#upload sample traits
tdat = read_csv('metadata/ALL_Coldata.csv')

#need to sum lane duplicates again
ddat = read_csv('metadata/lane_duplicate_runs.csv', col_names=c('Run', 'fileName')) %>% 
  mutate(sample=sapply(fileName, function(x) strsplit(x, '_')[[1]][1]))


#upload raw symbiont read counts
sdat0 = read_tsv('metadata/allSymbTranCounts.tsv') %>% 
  rename(Run = sample)


# sdat0 = read_tsv('metadata/allSymbGenomeCounts.tsv') %>% 
#   rename(Run = sample) %>% 
#   filter(Run!='SRR6396378')


#subset the runs that are lane duplicates and sum them
sumedDups = sdat0 %>% 
  filter(Run %in% ddat$Run) %>% 
  merge(ddat, by = 'Run') %>% 
  group_by(sample) %>% 
  summarize(Run=Run,
            all=sum(all),
            nonSym=sum(nonSym),
            cladeA=sum(cladeA),
            cladeB=sum(cladeB),
            cladeC=sum(cladeC),
            cladeD=sum(cladeD)) %>% 
  select(-sample)

#assmelbe back together, get proportions, merge with other sample traits
sdat1 = sdat0 %>% 
  filter(!Run %in% ddat$Run) %>% 
  rbind(sumedDups) %>% 
  mutate(symb = cladeA+cladeB+cladeC+cladeD,
         fracA=cladeA/symb,
         fracB=cladeB/symb,
         fracC=cladeC/symb,
         fracD=cladeD/symb,
         propSym=symb/all,
         propHost=nonSym/all) %>% 
  select(-nonSym,cladeA,cladeB,cladeC,cladeD) %>%
  full_join(tdat, by = 'Run') %>% 
  filter(!is.na(my_stage)) %>% 
  mutate(Organism = if_else(Organism=='Acropora hyacinthus complex sp. E JTL-2012',
                            'Acropora hyacinthus',
                            Organism),
         Organism = sub('Acropora', 'A.', Organism))


# add symbiont amount z-scores within projects ----------------------------

#z-score funciton
z_score = function(vector){
  mn=mean(vector)
  z= (vector - mean(vector)) / sd(vector)
  return(z)
}

#get them
projects = unique(sdat1$my_title)
sdat = data.frame()
for (pj in projects){
  psub = sdat1 %>% 
    filter(my_title==pj) %>% 
    mutate(zSym = z_score(propSym))
  sdat=rbind(sdat,psub)
}
dim(sdat1)
dim(sdat)

#save
coldata=data.frame(sdat, row.names=sdat$Run)
save(coldata, file='metadata/coldataWithSymbionts.Rdata')


# OVERALL STATS -----------------------------------------------------------

#hist
sdat %>% 
  ggplot(aes(x=propSym)) +
  geom_histogram()

#by developmental stage
sdat %>% 
  mutate(my_stage = factor(my_stage, levels=c('gamete', 'embryo', 'larval', 'recruit', 'adult'))) %>% 
  ggplot(aes(x=my_stage, y=propSym)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2) +
  labs(x='developmental stage', y='proportion symbiont reads')


#adult and larvae
sdat %>% 
  filter(my_stage %in% c('adult', 'larval')) %>% 
  ggplot(aes(x=my_stage, y=propSym*100)) +
  geom_boxplot() +
  labs(x='Developmental stage', y='symbiodinium reads (%)') +
  geom_hline(yintercept = 0, lty=2)


#----- CLADE PROPORTIONS BY SPECIES
adat = sdat %>% 
  filter(my_stage %in% c('adult'))

#clade A
order=adat %>% 
  filter(my_stage %in% c('adult')) %>% 
  group_by(Organism) %>% 
  summarize(med=median(fracA)) %>% 
  arrange(desc(med)) %>% 
  pull(Organism)

sdat %>% 
  mutate(Organism = factor(Organism, levels=order)) %>% 
  ggplot(aes(x=Organism, y=fracA*100)) +
  geom_boxplot() +
  labs(x='Species', y='Clade A (%)')

#clade D
order=adat %>% 
  filter(my_stage %in% c('adult')) %>% 
  group_by(Organism) %>% 
  summarize(med=median(fracD)) %>% 
  arrange(desc(med)) %>% 
  pull(Organism)

sdat %>% 
  mutate(Organism = factor(Organism, levels=order)) %>% 
  ggplot(aes(x=Organism, y=fracD*100)) +
  geom_boxplot() +
  labs(x='Species', y='Clade D (%)')


#barplot of all at once
adults = sdat %>% 
  filter(my_stage=='adult')
tprops = adults %>% 
  group_by(Organism) %>% 
  summarize(cA=sum(cladeA),
            cB=sum(cladeB),
            cC=sum(cladeC),
            cD=sum(cladeD)) %>% 
  mutate(tot=cA+cB+cC+cD,
         A=cA/tot,
         B=cB/tot,
         C=cC/tot,
         D=cD/tot) %>% 
  select(Organism,A,B,C,D) %>% 
  gather(key=clade, value='fraction', A,B,C,D) %>% 
  group_by(Organism, clade) %>% 
  summarize(mn=mean(fraction, na.rm=TRUE)) %>% 
  mutate(species=Organism)

order = tprops %>% 
  filter(clade=='C') %>% 
  arrange(desc(mn)) %>% 
  pull(species)

last3 = c('A. aspera', 'A. cervicornis', 'A. palmata')
order2 = c(order[!order %in% last3], last3)

tprops %>% 
  mutate(species=factor(species, levels=order2)) %>% 
  ggplot(aes(x=species, y=mn, fill=clade)) +
  geom_bar(stat='identity') +
  labs(y='proportion', x='species') +
  theme(axis.text.x=element_text(face="italic"))


# LOOK AT BEWW SYMBIONTS --------------------------------------------------

bew = sdat %>% 
  filter(grepl("^s", Run))
bew %>% 
  mutate(bleached = if_else(treat=='control',
                            'control',
                            'bleached')) %>% 
  ggplot(aes(x=fracD, fill=bleached)) +
  geom_density(alpha=0.4) +
  labs(x='Proportion clade D')




# GET BLEACHING ESTIMATES -------------------------------------------------
#goal here is to go through each stress project and 
#get the ratio of symbionts in stress over controls based on read counts



######### optionally add in zooxType_PL.tsv
# propSym=sdat$propSym
# 
# pldat0 = read_tsv('metadata/zooxType_PL.tsv') %>% 
#   mutate(Run=sub('.sam', '', sample))
# 
# pldat = sdat %>% 
#   left_join(pldat0, by = 'Run') %>% 
#   select(Run, fracZ) 
# head(pldat)
# sum(pldat$Run==sdat$Run)==nrow(sdat)
# sdat$propSym = pldat$fracZ


#read in the trait data from wrangle_sample_data.R
stdat = read_csv('metadata/subset_tables/stressColdata.csv')

ssdat = sdat %>% 
  filter(Run %in% stdat$Run) %>% 
  mutate(stress = if_else(treat=='control',
                          'control',
                          'stress'))

#set up modified list like in wrangle_sample_data.R
datList = list()
projects = unique(ssdat$my_title)
for (i in 1:length(projects)){
  p=projects[i]
  datList[[p]]=ssdat %>% filter(my_title==p)
}
length(datList)

#set up empty to fill with modified dataframes
assembleList=list()


# j1_thisStudy_PRJNA559404 -------------------------------------------------
projects
selectProject = "j1_thisStudy_PRJNA559404"
pdat = datList[[selectProject]]

mod0=pdat %>% 
  separate(treatDescription, into=c('exp', 'geno', 'CorT', 'tp', 'rep')) %>% 
  unite('scGroup', c('exp', 'geno',  'tp', 'rep'), sep='.') %>% 
  mutate(scGroup=sub('CtH', 'Tw', scGroup),
         scGroup=sub('HtC', 'Tw', scGroup),
         scGroup=sub('Ht2', 'Tw', scGroup))
stressors = mod0 %>% 
  filter(stress=='stress') %>% 
  pull(treat) %>% 
  unique()

ldat = data.frame()
for (sts in stressors){
  ssub = mod0 %>% 
    filter(treat==sts)
  csub = mod0 %>% 
    filter(treat=='control' & scGroup %in% ssub$scGroup)
  if (sts=='heat'){
    csub=csub %>% 
      mutate(projType='temp')
  }
  r=rbind(ssub, csub) %>% 
    select(scGroup, projType, stress, propSym) %>% 
    spread(key=stress, value=propSym) %>% 
    mutate(sym.lfc=log(stress/control, 2))
  ldat = rbind(ldat, r)
}

ldat %>% 
  mutate(exp=sapply(scGroup, function(x) strsplit(x, '.', fixed=TRUE)[[1]][1]),
         expType=paste(exp, projType, sep='_')) %>% 
  ggplot(aes(x=expType, y=sym.lfc)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2) +
  labs(x='BEWW Experiment', y=bquote(log[2]~FC~bleach~vs~control))

#This was so strange that I didn't finish running this for other projects





# z-scores for symbiont amounts -------------------------------------------
z_score = function(vector){
  mn=mean(vector)
  z= (vector - mean(vector)) / sd(vector)
  return(z)
}

res=data.frame()
pltList=list()
for (sp in projects){
  print(sp)
  sub=ssdat %>% 
    filter(my_title==sp) %>% 
    mutate(propD=cladeD/all)
  sub$zSym=z_score(sub$symb / sub$all)
  res=rbind(res,sub)
  plt=sub %>% 
    ggplot(aes(x=stress, y=zSym)) +
    geom_boxplot() +
    labs(subtitle=sp) +
    lims(y=c(-2,2))
  # plt=sub %>% 
  #   ggplot(aes(x=propSym, fill=stress)) +
  #   geom_density() +
  #   labs(subtitle=sp)
  pltList[[sp]]=plt
}
head(res)

plot_grid(plotlist=pltList, nrow=4)





