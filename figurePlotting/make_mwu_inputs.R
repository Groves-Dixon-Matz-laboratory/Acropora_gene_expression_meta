#make_mwu_inputs.R
#make input files for GO MWU from the deseq results files output for each of the stress subsets (see INITIALIZE COUNTS/CONTROL FOR COVARIATES/DESEQ in walkthough)
library(DESeq2)
source('figurePlotting/rna_functions.R')



# HIGH AND LOW STRESS BY DIRECTION ----------------------------------------


#First write out the overall high-stress and low-stress go inputs by direction
#write out for corStress
ll=load('deseqResults/corStress_deseqResults.Rdata')
ll
write_out_up_and_down_go(res, './go_mwu/clusterAstress')
write_out_go(res, './go_mwu/clusterAstress')

#write out for lowStress
ll=load('deseqResults/lowStress_deseqResults.Rdata')
ll
write_out_up_and_down_go(res, './go_mwu/clusterBstress')
write_out_go(res, './go_mwu/clusterBstress')


# COMBINATIONS OF PROJECTS -------------------------------

#get the names and paths for the output files
groupedNames = sub('_deseqResults.Rdata', '', list.files('correlated_only/deseqResults', "*deseqResults.Rdata"))
groupedPaths = list.files('correlated_only/deseqResults', "*deseqResults.Rdata", full.names = TRUE)


#look through and make go_mwu inputs for them
for (i in 1:length(groupedPaths)){
  path = groupedPaths[i]
  name = groupedNames[i]
  ll=load(path)
  ll=load(path)
  print(ll)
  write_out_up_and_down_go(res, paste('./go_mwu', name, sep='/'))
  write_out_go(res, paste('./go_mwu', name, sep='/'))
}

groupedNames = append(groupedNames, c('clusterAstress', 'clusterBstress'))
save(groupedNames, file='./go_mwu/groupedInputs.Rdata')


# INDIVIDUAL PROJECTS -----------------------------------------------------

individualNames = sub('_deseqResults.Rdata', '', list.files('deseqResults/individual_projects', "*deseqResults.Rdata"))
individualPaths = list.files('deseqResults/individual_projects', "*deseqResults.Rdata", full.names = TRUE)

for (i in 1:length(individualPaths)){
  path = individualPaths[i]
  name = individualNames[i]
  ll=load(path)
  print(ll)
  write_out_up_and_down_go(res, paste('./go_mwu', name, sep='/'))
  write_out_go(res, paste('./go_mwu', name, sep='/'))
}

save(individualNames, file='./go_mwu/individualInputs.Rdata')
