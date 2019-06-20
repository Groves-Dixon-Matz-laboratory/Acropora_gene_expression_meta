#wrangle_sample_data_specifics.R
#add more stats to some of the trait tables
library(tidyverse)

# make stratified random samples for full stress --------------------------
#idea here is to build a training set and testing set for stress_prediction.R with
#random but equal representation from each type of stress
stress = read_csv('./metadata/subset_tables/stressColdata.csv')
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

# make heat more specific -------------------------------------------------

heat0 = read_csv('./metadata/subset_tables/heatColdata.csv')
heatProjects = unique(heat0$my_title)
modHeats = list()

#deal with E_aspera_PRJNA266455
proj='E_aspera_PRJNA266455'
sub0=heat0 %>% 
  filter(my_title==proj)
exposure = sapply(sub0$treatDescription, function(x) strsplit(x, ' ')[[1]][2]) %>% 
  as.numeric()
sub = sub0 %>% 
  mutate(temp = 'not_clear',
         degOverAmb = if_else(treat=='heat',
                              6,
                              0),
         hoursExposed = exposure)
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
         hoursExposed = 72)
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
         hoursExposed = 72)
modHeats[[proj]]=sub

#deal with N_Hainan_U_heat_PRJNA308355
#couldn't find any info for this one
heatProjects
proj='N_Hainan_U_heat_PRJNA308355'
sub0=heat0 %>% 
  filter(my_title==proj)
sub = sub0 %>% 
  mutate(temp = NA,
         degOverAmb = NA,
         hoursExposed = NA)
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
         hoursExposed = 4)
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
         hoursExposed = 1)
modHeats[[proj]]=sub

#deal with j1_BLEACH_EWW
heatProjects
proj='j1_BLEACH_EWW'
sub0=heat0 %>% 
  filter(my_title==proj)
sub = sub0 %>% 
  mutate(temp = if_else(treat=='heat',
                        36,
                        28),
         degOverAmb = if_else(treat=='heat',
                              8,
                              0),
         hoursExposed = 3)
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
         hoursExposed = 3)
modHeats[[proj]]=sub

newHeat = modHeats %>% 
  purrr::reduce(rbind) %>% 
  as_tibble()

newHeat %>% 
  write_csv(path='./metadata/detailed_tables/detailedHeat.csv')




# make salinity more specific ---------------------------------------------









