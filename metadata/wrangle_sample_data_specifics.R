#wrangle_sample_data_specifics.R
#add more stats to some of the trait tables
library(tidyverse)
source('figurePlotting/rna_functions.R')
rm(list=ls())
# make heat more specific -------------------------------------------------
#include boolean corrStress for whether the study is among the correlated stress BioProjects (from project correlation heatmap)

heat0 = read_csv('./metadata/subset_tables/heat_Coldata.csv')
heatProjects = unique(heat0$my_title)
modHeats = list()

#deal with F_Uqueensland_ph_PRJNA269992 (Kaniewska et al. 2015)
proj='F_Uqueensland_ph_PRJNA269992'
sub0=heat0 %>% 
  filter(my_title==proj)
sub1 = sub0 %>% 
  separate(treatDescription, c('desc', 'temp', 'pco2'), sep='_') %>% 
  mutate(temp=as.numeric(sub('C', '', temp)))
sub = sub0 %>% 
  mutate(temp = sub1$temp,
         degOverAmb = temp - 24,
         hoursExposed =  24*7*5,#5 weeks
         corrStress = FALSE)
modHeats[[proj]]=sub


#deal with E_aspera_PRJNA266455
proj='E_aspera_PRJNA266455'
sub0=heat0 %>% 
  filter(my_title==proj)
exposure = sapply(sub0$treatDescription, function(x) strsplit(x, ' ')[[1]][2]) %>% 
  as.numeric()
sub = sub0 %>% 
  mutate(temp = if_else(treat=='heat',
                        30,
                        23.5),
         degOverAmb = if_else(treat=='heat',
                              6,
                              0),
         hoursExposed = exposure,
         corrStress = TRUE)
modHeats[[proj]]=sub

#deal with "H_matz_heatTolLat_PRJNA279192"   
proj='H_matz_heatTolLat_PRJNA279192'
sub0=heat0 %>% 
  filter(my_title==proj)
sub = sub0 %>% 
  mutate(temp = if_else(treat=='heat',
                        31.5,
                        28),
         degOverAmb = if_else(treat=='heat',
                              3.5,
                              0),
         hoursExposed = 72,
         corrStress = FALSE)
modHeats[[proj]]=sub

#deal with L_Barshis_bleachResillience_PRJNA177515
#information taken from methods in 'Genomic basis for coral resilience to climate change'
#temperature indicate means, because they fluctuated
heatProjects
proj='L_Barshis_bleachResillience_PRJNA177515'
sub0=heat0 %>% 
  filter(my_title==proj)
sub = sub0 %>% 
  mutate(temp = if_else(treat=='heat',
                        31.9,
                        29.2),
         degOverAmb = if_else(treat=='heat',
                              2.7,
                              0),
         hoursExposed = 72,
         corrStress = TRUE)
modHeats[[proj]]=sub

#deal with N_Hainan_U_heat_PRJNA308355
#couldn't find any info for this one
heatProjects
proj='N_Hainan_U_heat_PRJNA308355'
sub0=heat0 %>% 
  filter(my_title==proj)
sub = sub0 %>% 
  mutate(temp = if_else(treat=='heat',
                        32,
                        26),
         degOverAmb = if_else(treat=='heat',
                              6,
                              0),
         hoursExposed = 12,
         corrStress=FALSE)
modHeats[[proj]]=sub

#deal with W_gajigan_thermalMicroRNA_PRJNA298496
heatProjects
proj='W_gajigan_thermalMicroRNA_PRJNA298496'
sub0=heat0 %>% 
  filter(my_title==proj)
sub = sub0 %>% 
  mutate(temp = as.numeric(sub('C', '', treatDescription)),
         degOverAmb = if_else(treat=='heat',
                              5,
                              0),
         hoursExposed = 4,
         corrStress=TRUE)
modHeats[[proj]]=sub

#deal with f1_parkinson_hotCold_PRJNA423227
heatProjects
proj='f1_parkinson_hotCold_PRJNA423227'
sub0=heat0 %>% 
  filter(my_title==proj)
sub = sub0 %>% 
  mutate(temp = if_else(treat=='heat',
                        35,
                        28),
         degOverAmb = if_else(treat=='heat',
                              7,
                              0),
         hoursExposed = 1,
         corrStress=FALSE)
modHeats[[proj]]=sub

#deal with j1_BLEACH_EWW
heatProjects
proj='j1_thisStudy_PRJNA559404'
sub0=heat0 %>% 
  filter(my_title==proj)
sub = sub0 %>% 
  mutate(temp = if_else(treat=='heat',
                        36,
                        28),
         degOverAmb = if_else(treat=='heat',
                              8,
                              0),
         hoursExposed = 3,
         corrStress=TRUE)
modHeats[[proj]]=sub

#deal with k1_Palumbi_lab_heat_resilience_PRJNA274410
heatProjects
proj='k1_Palumbi_lab_heat_resilience_PRJNA274410'
sub0=heat0 %>% 
  filter(my_title==proj)
sub = sub0 %>% 
  mutate(temp = if_else(treat=='heat',
                        35,
                        29),
         degOverAmb = if_else(treat=='heat',
                              6,
                              0),
         hoursExposed = 3,
         corrStress=TRUE)
modHeats[[proj]]=sub

newHeat = modHeats %>% 
  purrr::reduce(rbind) %>% 
  as_tibble()

newHeat %>% 
  write_csv(path='./metadata/detailed_tables/detailedHeat.csv')


# look at immune more specific ---------------------------------------------
#note the main difference for these seems to be HOW the treatment was done
#the two that didn't correlate seems to just not really have stressed the coral
immune0 = read_csv('./metadata/subset_tables/immuneColdata.csv')

idat = immune0 %>% 
  group_by(my_title) %>% 
  summarize(N=n())


# look at pH studies ------------------------------------------------------

ph = read_csv('./metadata/subset_tables/phColdata.csv')

ph %>% 
  group_by(my_title) %>% 
  summarize(N=n())




# make stratified random samples for full stress --------------------------
#idea here is to build a training set and testing set for stress_prediction.R with
#random but equal representation from each type of stress
#note didn't end up using this
stress = read_csv('./metadata/subset_tables/allStress_Coldata.csv')
stress$proj_stress = paste(stress$my_title, stress$stress, sep='_')


probDev = 0.6
probVal = 0.4
proj_stress = unique(stress$proj_stress)
#assemble 
sample.ind = lapply(proj_stress, function(x) sample(2,
                                                    nrow(stress[stress$proj_stress==x,]),
                                                    replace=T,
                                                    prob = c(probDev, probVal)
)
) %>% 
  unlist()

trainSet = stress[sample.ind==1,]
testSet = stress[sample.ind==2,]
nrow(trainSet)/nrow(stress)
nrow(testSet)/nrow(stress)
table(trainSet$stress)
table(testSet$stress)
trainSet %>% 
  group_by(proj_stress) %>% 
  summarize(N=n())
testSet %>% 
  group_by(proj_stress) %>% 
  summarize(N=n())

trainSet %>% 
  write_csv(path='./metadata/subset_tables/stratified_allStress_train.csv')
testSet %>% 
  write_csv(path='./metadata/subset_tables/stratified_allStress_test.csv')
