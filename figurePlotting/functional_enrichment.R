#funcitonal_enrichment_correlations.R



inputFiles = c('corStress_down_For_MWU.csv',
               'corStress_up_For_MWU.csv',
               'lowStress_down_For_MWU.csv',
               'lowStress_up_For_MWU.csv',
               'redMembership_ForMWU.csv',
               'stressDAPC_ForMWU.csv')


#combinations of bioprojects by stress type
ll=load('./go_mwu/groupedInputs.Rdata')
groupedNames = c('bleached', 'heatNOBLEACH', 'immune', 'salinityNOBLEACH', 'ph')
inputFiles = append(paste(groupedNames, 'down_For_MWU.csv', sep='_'),
                    paste(groupedNames, 'up_For_MWU.csv', sep='_'))

#individual bioprojects
ll=load('./go_mwu/individualInputs.Rdata')
inputFiles = append(paste(individualNames, 'down_For_MWU.csv', sep='_'),
                    paste(individualNames, 'up_For_MWU.csv', sep='_'))


ll=load('metadata/corStressProjs.Rdata')
ll
sum(corStressProjs %in% individualNames)==length(corStressProjs)
sum(lowStressProjs %in% individualNames)==length(lowStressProjs)


goDivision = 'BP'
combind_up_and_down_go_res = function(inputName){
  upName = paste(paste('./go_mwu/MWU', goDivision, sep = "_"), paste(inputName,'up_For_MWU.csv',sep='_'), sep = "_")
  downName = paste(paste('./go_mwu/MWU', goDivision, sep = "_"), paste(inputName, 'down_For_MWU.csv',sep='_'), sep = "_")
  upRes = read.table(upName, header = TRUE, stringsAsFactors=FALSE) %>% 
    mutate(enrichment = -log(pval, 10))
  downRes = read.table(downName, header = TRUE, stringsAsFactors=FALSE) %>% 
    mutate(enrichment = log(pval, 10))
  res = rbind(upRes, downRes)
  res$inputName = inputName
  return(res)
}

plotenrich = function(cdat){
  cdat %>% 
    ggplot(aes(x=enrichment, y=delta.rank)) +
    geom_point() 
}





cdatList = lapply(individualNames, function(x) combind_up_and_down_go_res(x))
names(cdatList) = individualNames
cdat = cdatList %>% 
  purrr::reduce(rbind)


xInputName = 'corStress'
xdat = combind_up_and_down_go_res(xInputName) %>% 
  select(enrichment, name)

pdat = xdat %>% 
  left_join(cdat, by = 'name')



head(pdat)

pdat %>% 
  mutate(sGroup = if_else(inputName %in% corStressProjs,
                          'high-stress',
                          'low-stress')) %>% 
  ggplot(aes(x=enrichment.x, y=enrichment.y, color=sGroup)) +
  geom_point() +
  geom_smooth(method='lm') 



# PLOT CORRELATION OF TWO-TAILED DELTA RANKS FROM GO MWU FOR COMBINATIONS ------------------

source('figurePlotting/rna_functions.R')



########### CHOOSE DATASET YOU WANT TO PLOT 

#COMBINATIONS OF PROJECTS
#combinations of bioprojects by stress type
ll=load('./go_mwu/groupedInputs.Rdata')
groupedNames = c('bleached', 'heatNOBLEACH', 'immune', 'salinityNOBLEACH', 'ph')
inputNames = groupedNames


# READ IN AND BIND INTO SINGLE DATAFRAME 

cdatList = lapply(inputNames, function(x) read_in_two_tailed_go_results(x))
cdat = cdatList %>% 
  purrr::reduce(rbind)



# READ IN DATA TO BE X-AXIS 


xInputName = 'corStress'
xdat = read_in_two_tailed(xInputName) %>% 
  select(delta.rank, p.adj, name)


# MERGE 

pdat = xdat %>% 
  left_join(cdat, by = 'name')


# PLOT

pdat %>% 
  filter(!is.na(inputName)) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point() +
  geom_smooth(method='lm')



# REPEAT DELTA RANK CORRELATION WITH INDIVIDUAL PROJECTS ------------------

#INDIVIDUAL PROJECTS
ll=load('./go_mwu/individualInputs.Rdata')
inputNames = individualNames
ll=load('metadata/corStressProjs.Rdata')
ll

# READ IN AND BIND INTO SINGLE DATAFRAME 

icdatList = lapply(inputNames, function(x) read_in_two_tailed_go_results(x))
icdat = cdatList %>% 
  purrr::reduce(rbind)



# READ IN DATA TO BE X-AXIS 

xInputName = 'corStress'
xdat = read_in_two_tailed(xInputName) %>% 
  select(delta.rank, p.adj, name)

# MERGE 

ipdat = xdat %>% 
  left_join(cdat, by = 'name')

# PLOT

highstress = ipdat %>% 
  filter(!is.na(inputName),
         inputName %in% corStressProjs) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point() +
  geom_smooth(method='lm')

lowstress = ipdat %>% 
  filter(!is.na(inputName),
         inputName %in% lowStressProjs) %>% 
  ggplot(aes(x=delta.rank.x,
             y=delta.rank.y,
             color=inputName)) +
  geom_point() +
  geom_smooth(method='lm')


plot_grid(highstress, lowstress)
