#wrangle_sample_data.R

library(tidyverse)
setwd('./metadata')
dat0 = read.csv('my_all_acropora_sra_runInfoTable.csv', stringsAsFactors=FALSE)


# ADD DATA FROM BLEACHING EVERY WHICH WAY PROJECT -------------------------

#upload the tag-seq runs
runs = read_tsv('./project_specific_tables/bleachEWW_runs.txt', col_names='run') %>% 
  separate(run, into=c('s', 'S'), remove=F) %>% 
  mutate(num=as.numeric(sub('s', '', s)))
runs

#upload the trait data for the runs
traits0 = read_csv('./project_specific_tables/bleachEWW_traits.csv') %>% 
  mutate(timepoint = paste('tp', timepoint, sep=''),
         treatment = paste(experiment, treat, sep='_')) %>% 
  unite(col='sampleID', exp,geno,treat,timepoint,rep, sep='.', remove=F)

#merge
traits = runs %>% 
  left_join(traits0, by = 'num')
nrow(traits)
nrow(runs)


#make equivalent dataframe
bdat = data.frame(
  'my_title' = rep('j1_thisStudy_PRJNA559404', nrow(traits)),
  'BioSample' = traits$sampleID,
  'Run' = traits$run,
  'Library_Name' = traits$treatment,
  'Sample_Name' = traits$sampleID,
  'isolate' = traits$sampleID,
  "Organism" = rep('Acropora millepora', nrow(traits)),
  'my_treatment' = traits$treat,
  'my_stage' = rep('adult', nrow(traits))
)

#are the dataframes equivalent?
sum(colnames(bdat) == colnames(dat0))==ncol(dat0)

#assemble
dat = dat0 %>% 
  rbind(bdat)
nrow(bdat)
nrow(dat0)
nrow(dat)



# SET UP LIST OF DATAFRAMES TO ADD PROJECT-SPECIFIC INFO ------------------

datList = list()
projects = unique(dat$my_title)
for (i in 1:length(projects)){
  p=projects[i]
  datList[[p]]=dat %>% filter(my_title==p)
}
length(datList)

#set up empty to fill with modified dataframes
assembleList=list()




### Go through and fill in treatment info for each project
#assign:
#projType = type of project, (eg. pH, temperature, light), this will be same for all samples in each project
#treatDescription = longer description of treatment for refernce, will vary by treatment group
#treat            = coded treatment group for sample

# A_moya_acid_PRJNA149513 -------------------------------------------------
selectProject = 'A_moya_acid_PRJNA149513'
pdat = datList[[selectProject]]
mod = pdat %>% 
  mutate(projType = 'pH',
         treatDescription = sapply(pdat$Library_Name, function(x) strsplit(x, ": ", fixed=TRUE)[[1]][2]),
    treat = if_else(grepl('Control', treatDescription),
                         'control',
                         'low_pH' ),
         special = 'none'
         ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# B_whiteBand_PRJNA222758 -------------------------------------------------
selectProject = 'B_whiteBand_PRJNA222758'
pdat = datList[[selectProject]]
mod = pdat %>% 
  mutate(projType = 'immune',
         treatDescription = if_else(BioSample=='SAMN02380460',
                                    'white band',
                                    'healthy' ),
         treat = if_else(BioSample=='SAMN02380460',
                         'challenge',
                         'control' ),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# C_palumbi_pools_PRJNA242821 ----------------------------------------------
#note that pool 300 is the Highly Variable Pool and pool 400 is the Medium Variable Pool
selectProject = 'C_palumbi_pools_PRJNA242821'
pdat = datList[[selectProject]]
mod = pdat %>% 
  mutate(projType = 'temperature_resilience',
         treatDescription = if_else(isolate=='o300_t300',
                                    'HVP_to_HVP',
                                    'not_assigned'),
         treatDescription = if_else(isolate=='o300_t400',
                                    'HVP_to_MVP',
                                    treatDescription),
         treatDescription = if_else(isolate=='o400_t300',
                                    'MVP_to_HVP',
                                    treatDescription),
         treatDescription = if_else(isolate=='o400_t400',
                                    'MVP_to_MVP',
                                    treatDescription),
         treat = treatDescription,
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# E_aspera_PRJNA266455 ----------------------------------------------------
selectProject = 'E_aspera_PRJNA266455'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'temp',
         treatDescription = my_treatment,
         treat = if_else(grepl('(C', Sample_Name, fixed=T),
                         'control',
                         'heat'),
         special = 'none') %>% 
  filter(!grepl('N', Sample_Name)) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# F_Uqueensland_ph_PRJNA269992 --------------------------------------------
selectProject = 'F_Uqueensland_ph_PRJNA269992'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'pH+heat',
         treatDescription = my_treatment,
         treat = if_else(grepl('RCP', treatDescription),
                         'multiple',
                         'control' ),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# G_barilanU_moonlight_PRJNA276913 ----------------------------------------

selectProject = 'G_barilanU_moonlight_PRJNA276913'
pdat = datList[[selectProject]]
data.frame(pdat)
unique(pdat$Sample_Name)
mod = pdat %>% 
  mutate(projType = 'light',
         treatDescription = Library_Name,
         treat = 'not_done_yet',
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod



# H_matz_heatTolLat_PRJNA279192 -------------------------------------------
#note only adults were exposed to heat stress
#larval heat stress data came from an older dataset based on SOLID reads, so not used in this study
#for the adults, rep numbers 1-3 were controls and 4-6 got heat stress
selectProject = 'H_matz_heatTolLat_PRJNA279192'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = if_else(my_stage=='larval',
                            'temperature_resilience',
                            'temp'),
         repNum = '1',
         repNum = if_else(my_stage=='adult',
                          substr(Sample_Name, 2,2),
                          substr(Sample_Name, 3,3)),
         repNum = as.numeric(repNum),
                          
         treatDescription = if_else(repNum>3,
                                    '3 days at 31.5C',
                                    '3 days at 28C'),
         treatDescription = if_else(my_stage=='larval',
                                    Library_Name,
                                    treatDescription),
         treat = if_else(grepl('31.5C', treatDescription),
                         'heat',
                         'control'),
         treat = if_else(my_stage=='larval',
                         'none',
                         treat),
         special = 'larvae_not_treated'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# I_bertucci_dayNight_PRJNA288809 -----------------------------------------

selectProject = 'I_bertucci_dayNight_PRJNA288809'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'circadian',
         treatDescription = my_treatment,
         treat = if_else(grepl('day', my_treatment),
                         'day',
                         'night'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# J_weiss_immune_PRJNA200542 ----------------------------------------------
selectProject = 'J_weiss_immune_PRJNA200542'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'immune',
         treatDescription = if_else(grepl('mdp', my_treatment),
                                    'bacterial mimic',
                                    'control'),
         treatDescription = if_else(grepl('pic', my_treatment),
                                    'viral mimic',
                                    treatDescription),
         treat = if_else(grepl('mimic', treatDescription),
                         'challenge',
                         'control'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# K_strader_redDiapause_PRJNA292574 ---------------------------------------
selectProject = 'K_strader_redDiapause_PRJNA292574'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'light',
         treatDescription = Library_Name,
         treat = 'not done yet',
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# L_Barshis_bleachResillience_PRJNA177515 ---------------------------------

selectProject = 'L_Barshis_bleachResillience_PRJNA177515'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'temp',
         treatDescription = Sample_Name,
         treat = if_else(grepl('_h_', Sample_Name),
                                    'heat',
                                    'control'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# M_moya_juvenilePh_PRJNA260269 -------------------------------------------

selectProject = 'M_moya_juvenilePh_PRJNA260269'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'pH',
         treatDescription = my_treatment,
         treat = if_else(grepl('380', my_treatment),
                         'control',
                         'low_pH'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# N_Hainan_U_heat_PRJNA308355 ---------------------------------------------
selectProject = 'N_Hainan_U_heat_PRJNA308355'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'temp',
         treatDescription = Sample_Name,
         treat = if_else(grepl('Heat', Sample_Name),
                         'heat',
                         'control'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# O_Mohamed_zoox_PRJNA309168 ----------------------------------------------
selectProject = 'O_Mohamed_zoox_PRJNA309168'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'zoox',
         treatDescription = my_treatment,
         treat = if_else(grepl('infected', my_treatment),
                         'infected',
                         'control'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# P_Vize_timeKeeping_PRJNA293100 ------------------------------------------
#Sample_Names:
#1=tank/field
#2=rep
#3=night/day
#4/5=moon type
#F3NFM
selectProject = 'P_Vize_timeKeeping_PRJNA293100'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'circadian.moon',
         location = if_else(substr(Sample_Name, 1,1)=='F',
                            'field',
                            'tank'),
         timeOfDay=if_else(substr(Sample_Name, 3,3)=='N',
                           'night',
                           'day'),
         moon=substr(Sample_Name, 4,5),
         treatDescription = paste(paste(location,timeOfDay,sep='_'),moon,sep='_'),
         treat = timeOfDay,
         special = 'none'
  ) %>% 
  dplyr::select(-location, -timeOfDay, -moon) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# Q_bar_ilanU_lightSpawn_PRJNA316795 --------------------------------------

selectProject = 'Q_bar_ilanU_lightSpawn_PRJNA316795'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'moon.spawn',
         treatDescription = my_treatment,
         treat = 'not done yet',
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# R_oist_devoStage_PRJDB3244 ----------------------------------------------

selectProject = 'R_oist_devoStage_PRJDB3244'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'development',
         treatDescription = Library_Name,
         treat = if_else(grepl('Adult', Library_Name),
                         'adult',
                         'blastular'),
         treat = if_else(grepl('Gastrula', Library_Name),
                         'gastrula',
                         treat),
         treat = if_else(grepl('Planula', Library_Name),
                         'planula',
                         treat),
         treat = if_else(grepl('Sphere', Library_Name),
                         'sphere',
                         treat),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# S_rachel_immune_PRJNA319662 ---------------------------------------------

selectProject = 'S_rachel_immune_PRJNA319662'
pdat = datList[[selectProject]]
pdat$sam = sub('Sample','',sapply(pdat$Sample_Name, function(x) strsplit(x, '_')[[1]][1]))
head(data.frame(pdat))
annots = read_csv('./project_specific_tables/NewGood_aimsrnaseqsamples.csv', comment='#')
dupSum = read_tsv('./project_specific_tables/laneDupSummary.tsv') #output from sum_lane_dups.R after running feature counts (see data processing walkthough)
sum(dupSum$nlaneDups)==nrow(pdat)
toKeep = dupSum %>% 
  pull(keptName)

mod = pdat %>% 
  left_join(annots, by = 'sam') %>% 
  mutate(projType = 'immune',
         treatDescription = sam,
         treat = if_else(bac=='D',
                             'challenge',
                             'control'),
         treat = if_else(bac=='O',
                             'challenge',
                             treat),
         special = 'laneDups removed with sum_lane_dups.R'
  ) %>% 
  filter(Run %in% toKeep,
         !is.na(bac)) %>% 
  select(colnames(assembleList[[1]])) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# T_yasuoka_brachyury_PRJDB4579 -------------------------------------------
selectProject = 'T_yasuoka_brachyury_PRJDB4579'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'development',
         treatDescription = Library_Name,
         treat = 'not done yet',
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# U_standford_natEnvChange_PRJNA338455 ------------------------------------
selectProject = 'U_standford_natEnvChange_PRJNA338455'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'natural_variation',
         treatDescription = isolate,
         treat = 'none',
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# V_kenkel_co2seep_PRJNA362652 --------------------------------------------
#Libary_Name gives sampling location
#either from a co2 seep or control location
#First letter (I or D) indicates location: I=Upa-Upasina; D=Dobu
#number indicates genotype
selectProject = 'V_kenkel_co2seep_PRJNA362652'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'pH_resilience',
         loc = if_else(substr(Library_Name, 1,1)=='I',
                       'Upa-Upasina',
                       'Dobu'),
         seep = if_else(substr(Library_Name, 2,2)=='B',
                        'seep',
                        'control'),
         geno = substr(Library_Name, 3,5),
         treatDescription = paste(paste(loc,seep,sep='_'), geno,sep='_'),
         treat = if_else(seep=='seep',
                         'low_pH',
                         'control'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# W_gajigan_thermalMicroRNA_PRJNA298496 -----------------------------------
selectProject = 'W_gajigan_thermalMicroRNA_PRJNA298496'
pdat = datList[[selectProject]]
pdat$degC = sapply(pdat$Sample_Name, function(x) strsplit(x, ' ')[[1]][3])
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'temp',
         treatDescription = degC,
         treat = if_else(degC=='34C',
                         'heat',
                         'control'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# X_strader_larvalCompetance_PRJNA379147 ----------------------------------
selectProject = 'X_strader_larvalCompetance_PRJNA379147'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'development',
         treatDescription = 'not done yet',
         treat = 'not done yet',
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# Y_rose_polygenicThermal_PRJNA379450 -------------------------------------

selectProject = 'Y_rose_polygenicThermal_PRJNA379450'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'cryptic_species',
         treatDescription = 'not done yet',
         treat = 'not done yet',
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# Z_takahashiKariyazono_fluorescentCN_PRJDB4562 ---------------------------

selectProject = 'Z_takahashiKariyazono_fluorescentCN_PRJDB4562'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'species',
         spp=if_else(grepl('_D', Library_Name),
                     'digitifera',
                     'tenuis'),
         treatDescription = paste(spp, my_stage,sep='_'),
         treat = treatDescription,
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# a1_Yaun_CO2_earlyDevelop_PRJNA380146 ------------------------------------

selectProject = 'a1_Yaun_CO2_earlyDevelop_PRJNA380146'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'species',
         spp=if_else(grepl('_D', Library_Name),
                     'digitifera',
                     'tenuis'),
         treatDescription = paste(spp, my_stage,sep='_'),
         treat = treatDescription,
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# b1_Aguilar_hypoosmotic stress_PRJNA380267 --------------------------------
selectProject = 'b1_Aguilar_hypoosmotic_stress_PRJNA380267'
ttable = read_tsv('./project_specific_tables/b1_treat_table.tsv')
pdat = datList[[selectProject]] %>% 
  left_join(ttable, by = 'BioSample')
mod = pdat %>% 
  mutate(projType = 'salinity',
         treatDescription = name,
         treat = if_else(grepl('35psu', salt),
                               'control',
                               'low_salinity'),
         special = 'none') %>% 
  select(-num, -name, -time, -salt) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# c1_uga_diseaseTolerance_PRJNA386795 -------------------------------------

selectProject = 'c1_uga_diseaseTolerance_PRJNA386795'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'disease_tolerance',
         treatDescription = 'emailed John Wares',
         treat = 'not done yet',
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# d1_mohamed_microalga_PRJNA398338 ----------------------------------------

selectProject = 'd1_mohamed_microalga_PRJNA398338'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'maybe_immune',
         treatDescription = my_treatment,
         treat = if_else(grepl('Chromera_infected', my_treatment),
                         'control',
                         'challenge'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# e1_bayPalumbi_adaptWarm_PRJNA411943 -------------------------------------

selectProject = 'e1_bayPalumbi_adaptWarm_PRJNA411943'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'genotyping',
         treatDescription = my_treatment,
         treat = 'not done yet',
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# f1_parkinson_hotCold_PRJNA423227 ----------------------------------------

#note, the cold treated samples in these surprisingly look less stressed than the ambient
#not sure what's going on there, but double-checked the treatment calls against those
#on SRA database for the samples on 6-17-19

selectProject = 'f1_parkinson_hotCold_PRJNA423227'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'temp',
         treatDescription = my_treatment,
         treat = if_else(grepl('_35C_', my_treatment),
                         'heat',
                         'control'),
         treat = if_else(grepl('_10C_', my_treatment),
                         'cold',
                         treat),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# g1_unknown_earlySymbiosis -----------------------------------------------

selectProject = 'g1_unknown_earlySymbiosis_PRJDB4715'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'zoox_or_immune',
         treatDescription = my_treatment,
         treat = if_else(grepl('clade C symbiotic', my_treatment),
                         'cladeC',
                         'control'),
         treat = if_else(grepl('clade D symbiotic', my_treatment),
                         'cladeD',
                         treat),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# h1_rocker_waterQuality_PRJNA476311 --------------------------------------
#From NCBI sample identifiers: tag-based RNAseq; PI=Pelorus Island, GB=Geoffrey Bay; 1st number indicates genotype; 2nd number indicates deployment locale (1.=home, .2=away)
#From publication: 
#  Palm Island (PI) is considered 'very good' water quality
#  Geoffrey Bay (GB) is considered 'moderate' water quality
selectProject = 'h1_rocker_waterQuality_PRJNA476311'
pdat = datList[[selectProject]]
data.frame(pdat)

pdat %>% 
  filter(grepl('GB_24', isolate)) %>% 
  data.frame()

mod = pdat %>% 
  separate(Library_Name, into=c('origin', 'colony'), sep='_') %>% 
  separate(colony, into=c('genotypeNumber', 'transplantNumber')) %>%
  mutate(projType = 'transplant_or_temp',
         transplantedTo = if_else(transplantNumber==1,
                                  'home',
                                  'away'),
         treatDescription = paste(paste(origin,genotypeNumber, sep='_'),sep='_'),
         treat = if_else(origin=='GB',
                         'moderate',
                         'veryGood'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod

# i1_TakahashiKariyazono_fluorescentPolymorph -----------------------------

selectProject = 'i1_TakahashiKariyazono_fluorescentPolymorph_PRJDB6468'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'fluorescence',
         treatDescription = 'not done yet',
         treat = 'not done yet',
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# j1_thisStudy_PRJNA559404 -----------------------------------------------------------

selectProject = 'j1_thisStudy_PRJNA559404'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  #remove mismatched samples (This was based on PCA data where they are clearly mixed up)
  filter(!Run %in% c('s37_S74', 's38_S75', 's39_S76', 's40_S77', 's41_S78')) %>% 
  mutate(projType = if_else(grepl('low_salinity', Library_Name),
                            'salinity',
                            'temp'),
         projType = if_else(grepl('three_way', Library_Name) | grepl('hot_to_cold', Library_Name) | grepl('cold_to_hot', Library_Name),
                            'multiple',
                            projType),
         treatDescription = Sample_Name,
         treat = 'control',
         treat = if_else(grepl('chill', Library_Name) & my_treatment=='T',
                         'cold',
                         treat),
         treat = if_else(grepl('heat', Library_Name) & my_treatment=='T',
                         'heat',
                         treat),
         treat = if_else(grepl('hot_to_cold', Library_Name) & my_treatment=='T',
                         'heat_then_cold',
                         treat),
         treat = if_else(grepl('cold_to_hot', Library_Name) & my_treatment=='T',
                         'cold_then_heat',
                         treat),
         treat = if_else(projType=='salinity' & my_treatment=='T',
                         'low_salinity',
                         treat),
         special = 'none'
  ) %>% 
  data.frame()
assembleList[[selectProject]] = mod
dim(mod)

#---OUTPUT SPECIAL TABLE FOR SRA UPLOAD
sra0 = pdat %>% 
  filter(!Run %in% c('s37_S74', 's38_S75', 's39_S76', 's40_S77', 's41_S78'))
head(sra0)
exp = sub('_C', '', sub('_T', '', sra0$Library_Name))
breed = sapply(sra0$isolate, function(x) strsplit(x, '.', fixed=TRUE)[[1]][2])
treat = sapply(sra0$sample_name, function(x) strsplit(x, '.', fixed=TRUE)[[1]][3])
rep = sapply(sra0$isolate, function(x) strsplit(x, '.', fixed=TRUE)[[1]][5])
treatRep = paste(treat, rep, sep='')
iso = paste(paste(exp, breed, sep='_'), treatRep, sep='_')


#make unique iso by replacing sample name periods with underscores
iso = sub('.', '_', sra0$BioSample, fixed=TRUE)
iso = sub('.', '_', iso, fixed=TRUE)
iso = sub('.', '_', iso, fixed=TRUE)
iso = sub('.', '_', iso, fixed=TRUE)

sra = data.frame(sample_name = sra0$Sample_Name,
                 sample_title = sra0$Run,
                 organism=sra0$Organism,
                 isolate = iso,
                 breed = breed,
                 isolation_source = if_else(breed=='L1',
                                            'Little Pioneer Bay near Orpheus Island Research Station',
                                            'North East Orpheus near Orpheus Island Research Station'),
                 collection_date = '11/25/18',
                 geo_loc_name = 'Australia',
                 tissue='nubbin',
                 age='adult', stringsAsFactors=FALSE)

#build description strings
experiment = paste('experiment', exp, sep='=')
geno = paste('colony', breed, sep='=')
replicate = paste('replicate', sapply(sra$sample_name, function(x) strsplit(x, '.', fixed=TRUE)[[1]][5]), sep='=')
t = if_else(treat=='C',
                'control',
                'treated')
treatment= paste('treatment', t, sep='=')
timepoint = sub('tp', 'timepoint=', sapply(sra$sample_name, function(x) strsplit(x, '.', fixed=TRUE)[[1]][4]))

strdf = data.frame(experiment,
                   geno,
                   treatment,
                   replicate,
                   timepoint,
                   stringsAsFactors=FALSE)

#concatenate them
d = strdf %>% 
  unite('d', colnames(strdf), sep=';') %>% 
  pull(d)

#write out
sra$description = d
sra$method = paste('Method=Tagseq', d, sep=';')
sra$fileName = paste(sra$sample_title, 'fastq.gz', sep='.')
head(sra)
sra %>% 
  write_csv('project_specific_tables/sraUploadData.csv')

#This table has the information to paste into the SRA upload tables for both bioSamples and Runs
#Project was uploaded with a unique bioSample for each Run
#File upload was lane-concateanted fastq files for each sample



# k1_Palumbi_lab_heat_resilience_PRJNA274410 ------------------------------

selectProject = 'k1_Palumbi_lab_heat_resilience_PRJNA274410'
pdat = datList[[selectProject]]
data.frame(pdat)
mod = pdat %>% 
  mutate(projType = 'temp',
         treatDescription = Sample_Name,
         treat = if_else(grepl('_c_', treatDescription),
                         'control',
                         'heat'),
         special = 'none'
  ) %>% 
  data.frame()
mod
assembleList[[selectProject]] = mod


# ASSEMBLE TOGETHER -------------------------------------------------------
names(datList)
length(assembleList)
length(datList)
names(datList)[!names(datList) %in% names(assembleList)]

#funtion to select for particular columns in each 
keepCols = c('my_title', 'Run', 'my_stage', 'Organism', 'projType', 'treat', 'treatDescription', 'special')
selectCols = function(x){
  x %>% 
    select(keepCols)
}

#apply to each dataframe so they can be rbound
revList = lapply(assembleList, function(x) selectCols(x))
names(revList)

#rbind into single df
rdat = revList %>% 
  purrr::reduce(rbind) %>% 
  as_tibble()


projects = unique(rdat$my_title)

# GET SUMMARY INFO --------------------------------------------------------

#how many total runs
nrow(rdat)

#check projects
projects = unique(rdat$my_title)
projects
length(projects)


#how many projects of each type?
rdat %>% 
  group_by(projType) %>% 
  summarize(length(unique(my_title)))
sum(table(rdat$Run)==1)==nrow(rdat) #NOTE IF THIS ISN'T TRUE, A PROJECT PROBABLY GOT DROPPED, NEED TO RUN THROUGH MORE SLOWLY


# COMBINE WITH READ COUNT DATA --------------------------------------------
cdat = read_tsv('../pipeline_counts/all_pipeline_counts.txt', col_names=c('Run', 'value', 'stat'))
cdat$statRun = paste(cdat$stat, cdat$Run,sep='_')
dim(cdat)
cdat = cdat[!duplicated(cdat$statRun),]
dim(cdat)
cdat = cdat %>% 
  filter(stat %in% c("rawCounts", "geneCounted")) %>% 
  select(-statRun) %>% 
  spread(stat, value)

frdat0 = rdat %>% 
  left_join(cdat, by = 'Run')

nrow(frdat0)
nrow(rdat)



# COMBINE WITH BLEACHING DATA ---------------------------------------------

#upload bleaching info (this table was built by hand by looking at the publications)
#note that right now it includes only the general bleached or not bleached
#note that bleached call is assigned to controls and stressed corals so they can be isolated easily by experiment
blchDat = read_csv('./detailed_tables/bleaching_info.csv') %>% 
  select(Run, bleachGeneral)
nrow(blchDat)
unique(blchDat$bleachGeneral)

#recode messy calls into yes, no and NA
NAcalls = c('not stated')
notCalls = c('none', 'not yet')
bleachedCalls = c('mild', 'yes')
blchDat = blchDat %>% 
  mutate(bleached = ifelse(bleachGeneral %in% NAcalls,
                           NA,
                           'notAssigned'),
         bleached = if_else(bleachGeneral %in% notCalls,
                            'no',
                            bleached),
         bleached = if_else(bleachGeneral %in% bleachedCalls,
                            'yes',
                            bleached)) %>% 
  select(Run,bleached)
blchDat
unique(blchDat$bleached)


#add stress column
stressTypes = c('pH', 'temp', 'immune', 'salinity', 'multiple')

frdat = frdat0 %>% 
  left_join(blchDat,
            by='Run') %>% 
  mutate(stress = if_else(projType %in% stressTypes,
                         'unassigned',
                         'not_a_stress_experiment'),
         stress = if_else(projType %in% stressTypes & treat=='control',
                          'control',
                          stress),
         stress = if_else(projType %in% stressTypes & treat!='control',
                          'stressed',
                          stress))
unique(frdat$stress)
dim(frdat0)
dim(frdat)





# WRITE OUT FULL TABLE ----------------------------------------------------
source('../figurePlotting/rna_functions.R')
frdat %>% 
  add_single_or_pe() %>% 
  write_csv(path='./ALL_Coldata.csv')


frdat %>% 
  add_single_or_pe() %>% 
  pull(LibraryLayout) %>% 
  table()

# WRITE OUT FORMATTED TABLE FOR PUBLICATION -------------------------------

#gather bioProjects
titleSplits = sapply(frdat$my_title, function(x) strsplit(x, '_'))
bpList = list()
for (ts in titleSplits){
  bp=ts[[length(ts)]][1]
  bpList = append(bpList, bp)
}
bps=unlist(bpList)



#get reference IDs
refs = read_csv('detailed_tables/BioProject_Publication_Table.csv')

format = frdat %>% 
  mutate(BioProject=bps) %>% 
  left_join(refs, by = 'BioProject') %>% 
  rename(Developmental_stage=my_stage,
         Reference_year = year,
         Reference_DOI = DOI,
         Project_type = projType,
         Treatment = treat,
         Treat_description = treatDescription,
         Raw_counts = rawCounts,
         Gene_counted = geneCounted)

#write out
format[,c('Run',
          'BioProject',
          'Reference',
          'Reference_year',
          'Reference_DOI',
          'Organism',
          'Developmental_stage',
          'Project_type',
          'Treatment',
          'Treat_description',
          'Raw_counts',
          'Gene_counted')] %>% 
  write_csv('detailed_tables/my_sample_trait_table.csv')

#save



# ORGANIZE STRESS SET -----------------------------------------------

stressTypes = c('pH', 'temp', 'immune', 'salinity', 'multiple', 'pH+heat')
stress = frdat %>% 
  filter(projType %in% stressTypes,
         !is.na(treat)) %>% 
  mutate(stress = if_else(treat=='control',
                         'control',
                         'stressed'))
#check counts
table(stress$my_title)

#look at bleaching calls
stress %>% 
  group_by(my_title, bleached) %>% 
  summarize(N=n()) 

countTable = stress %>% 
  group_by(treat) %>% 
  summarize(`N Projects`=length(unique(my_title)),
            `N Samples`=n()) 
countTable


stress %>% 
  write_csv(path='./subset_tables/allStress_Coldata.csv')



# UPLOAD THE 'CORRELATED STRESS PROJECTS' ---------------------------------
#note these projects are identified with the script plot_individual_projects.R
#to run that script, you'll need to have run DESeq on all stress projects in the INITIALIZE COUNTS/CONTROL FOR COVARIATES/DESEQ section in walkthrough
ll=load('corStressProjs.Rdata')
ll
corStressProjs
lowStressProjs


#check new sample total table
stress %>% 
  filter(my_title %in% corStressProjs) %>% 
  group_by(treat) %>% 
  summarize(`N Projects`=length(unique(my_title)),
            `N Samples`=n()) 

#write out general stress for correlated projects
corAllStress = stress %>% 
  filter(my_title %in% corStressProjs)

corAllStress %>% 
  write_csv(path='./correlated_stress_project_tables/allStress_Coldata.csv')


# heat --------------------------------------------------------------------

#temp calls because come share controls
#make heat group as temps minus any cold treated dudes
#also remove BEWW samples that are not 'Ht1' or 'Ht2' in destriptionn
#also grab the three-way 'Tw' controls to go with Ht2 from this study

#set up the samples from j1 to remove:
j1keepers = frdat %>% 
  filter(my_title=='j1_thisStudy_PRJNA559404',
         grepl('Ht', treatDescription) | grepl('Tw', treatDescription),
         !grepl('HtC.', treatDescription)) %>% 
  pull(Run)
j1Remove = frdat %>% 
  filter(my_title=='j1_thisStudy_PRJNA559404',
         !Run %in% j1keepers) %>% 
  pull(Run)
  
  
  
#isolate heat set
heat = frdat %>% 
  mutate(projType = if_else(my_title=='j1_thisStudy_PRJNA559404' & grepl('^Tw', treatDescription),
                            'temp',
                            projType)) %>% 
  filter(projType %in% c('temp', 'pH+heat'),
         treat !='cold',
         !(my_title =='j1_thisStudy_PRJNA559404' & Run %in% j1Remove))
table(heat$my_title)
unique(heat$treat)



#standard table
heat %>% 
  write_csv(path='./subset_tables/heat_Coldata.csv')


#output the stress_noHeat set
stress %>% 
  filter(!Run %in% heat$Run,
         !projType=='multiple') %>% 
  write_csv('./subset_tables/stressNOHEAT_Coldata.csv')



#output version without bleaching samples
#This is to double-check that correlations between cold and salinity are not driven just by that
#note bleaching has to have been explicitly mentioned, otherwise no bleaching is assumed
heat %>% 
  filter(bleached == 'no' | is.na(bleached)) %>% 
  write_csv(path='./subset_tables/heatNOBLEACH_Coldata.csv')



# heat for correlated projects only ---------------------------------------


#write out for correlated projects only
corHeat = heat %>% 
  filter(my_title %in% corStressProjs) 
corHeat %>% 
  write_csv(path='./correlated_stress_project_tables/heat_Coldata.csv')

#write out the stress no heat set
corStresNoHeat = stress %>% 
  filter(my_title %in% corStressProjs,
         !Run %in% heat$Run,
         !projType=='multiple')
corStresNoHeat %>% 
  write_csv('./correlated_stress_project_tables/stressNOHEAT_Coldata.csv')



#and for correlated projects only without bleaching
corHeatNoBleach = heat %>% 
  filter(my_title %in% corStressProjs,
         bleached == 'no' | is.na(bleached))
corHeatNoBleach %>% 
  write_csv(path='./correlated_stress_project_tables/heatNOBLEACH_Coldata.csv')


# ORAGNIZE BLEACHED -------------------------------------------------------

bleached = frdat %>% 
  filter(!is.na(bleached) & bleached=='yes') %>% 
  mutate(bleached = if_else(treat=='control',
                            'control',
                            'bleached'))
unique(bleached$bleached)
bleached %>% 
  write_csv(path='./subset_tables/bleached_Coldata.csv')


#make a stress no-bleached set
blchRuns = bleached %>% 
  pull(Run)
stress %>% 
  filter(!Run %in% blchRuns) %>% 
  write_csv('./subset_tables/stressNOBLEACH_Coldata.csv')


# organize bleached for correlated ----------------------------------------
#This turns out to be exactly the same, but outputting for clarity

bleached %>% 
  filter(my_title %in% corStressProjs) %>% 
  write_csv(path='./correlated_stress_project_tables/bleached_Coldata.csv')


#make a stress no-bleached set
stress %>% 
  filter(my_title %in% corStressProjs,
         !Run %in% blchRuns) %>% 
  write_csv('./correlated_stress_project_tables/stressNOBLEACH_Coldata.csv')


# ORGANIZE pH SET ---------------------------------------------------------

#subset
stressTypes = c('pH')
ph = frdat %>% 
  filter(projType %in% stressTypes)

#check counts
table(ph$my_title)
ph %>% 
  group_by(projType) %>% 
  summarize(N=n())

#write out
ph %>% 
  write_csv(path='./subset_tables/ph_Coldata.csv')


#output the stress no pH set
stress %>% 
  filter(!Run %in% ph$Run) %>% 
  write_csv('./subset_tables/stressNOPH_Coldata.csv')


# ORGANIZE pH FOR CORRELATED PROJECTS -------------------------------------

#write out
ph %>% 
  filter(my_title %in% corStressProjs) %>% 
  write_csv(path='./correlated_stress_project_tables/ph_Coldata.csv')


#output the stress no pH set
stress %>% 
  filter(!Run %in% ph$Run,
         my_title %in% corStressProjs) %>% 
  write_csv('./correlated_stress_project_tables/stressNOPH_Coldata.csv')



# ORGANIZE DISEASE SET --------------------------------------------------------

#subset
stressTypes = c('immune')
immune = frdat %>% 
  filter(projType %in% stressTypes)

#check counts
table(immune$my_title)
immune %>% 
  group_by(my_title, treat) %>% 
  summarize(N=n())

#write out
immune %>% 
  write_csv(path='./subset_tables/immune_Coldata.csv')


#output the stress no immune
stress %>% 
  filter(!Run %in% immune$Run) %>% 
  write_csv('./subset_tables/stressNOIMMUNE_Coldata.csv')



# ORGANIZE DISEASE FOR CORRELATED PROJECTS --------------------------------

#write out
immune %>% 
  filter(my_title %in% corStressProjs) %>% 
  write_csv(path='./correlated_stress_project_tables/immune_Coldata.csv')


#output the stress_noHeat set
stress %>% 
  filter(!Run %in% immune$Run,
         my_title %in% corStressProjs) %>% 
  write_csv('./correlated_stress_project_tables/stressNOIMMUNE_Coldata.csv')



# ORGANIZE SALINITY -------------------------------------------------------

#subset
stressTypes = c('salinity')
salt = frdat %>% 
  filter(projType %in% stressTypes)

#check counts
table(salt$my_title)

#write out
salt %>% 
  write_csv(path='./subset_tables/salinity_Coldata.csv')


#write out a no BEWW version
salt %>% 
  filter(!Run %in% blchRuns) %>% 
  write_csv(path='./subset_tables/salinityNOBLEACH_Coldata.csv')
  

#output the stress no Salinity set
stress %>% 
  filter(!Run %in% salt$Run,
         !projType == 'multiple') %>% 
  write_csv('./subset_tables/stressNOSALINITY_Coldata.csv')



# ORGANIZE SALINITY FOR CORRELATED PROJECTS -------------------------------

#write out
salt %>% 
  filter(my_title %in% corStressProjs) %>% 
  write_csv(path='./correlated_stress_project_tables/salinity_Coldata.csv')


#write out a no bleaching version
salt %>% 
  filter(!Run %in% blchRuns,
         my_title %in% corStressProjs) %>% 
  write_csv(path='./correlated_stress_project_tables/salinityNOBLEACH_Coldata.csv')


#output the stress_noSalinity set
stress %>% 
  filter(!Run %in% salt$Run,
         !projType == 'multiple',
         my_title %in% corStressProjs) %>% 
  write_csv('./correlated_stress_project_tables/stressNOSALINITY_Coldata.csv')



# ORGANIZE COLD -----------------------------------------------------------

#first look at which projects have cold
projWithCold = frdat %>% 
  group_by(my_title, treat) %>% 
  summarize(N=n()) %>%
  filter(treat=='cold') %>% 
  pull(my_title)

#subset
cold = frdat %>% 
  filter(my_title %in% projWithCold,
         treat %in% c('cold', 'control'),
         !(my_title=='j1_thisStudy_PRJNA559404' & !grepl('Ch.', treatDescription)))



#check counts
cold %>% 
  group_by(my_title, treat) %>% 
  summarize(N=n())

#write out
cold %>% 
  write_csv(path='./subset_tables/cold_Coldata.csv')

#make a no bleaching version
cold %>% 
  filter(!Run %in% blchRuns) %>% 
  write_csv(path='./subset_tables/coldNOBLEACH_Coldata.csv')


#output the stress no Cold set
stress %>% 
  filter(!Run %in% cold$Run,
         !projType == 'multiple') %>% 
  write_csv('./subset_tables/stressNOCOLD_Coldata.csv')



# ORGANIZE COLD FOR CORRELATED ONLY ---------------------------------------

#write out
cold %>% 
  filter(my_title %in% corStressProjs) %>% 
  write_csv(path='./correlated_stress_project_tables/cold_Coldata.csv')

#no need to make a no BEWW version because only other one is not among correlated projects


#output the stress_noCold set
stress %>% 
  filter(!Run %in% cold$Run,
         !projType == 'multiple',
         my_title %in% corStressProjs) %>% 
  write_csv('./correlated_stress_project_tables/stressNOCOLD_Coldata.csv')



#set back to default working dir
setwd('../')

