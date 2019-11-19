#load_stress_for_scatterplots.R

source('figurePlotting/rna_functions.R')

datasetPrefix = 'standard'


# LOAD DESEQ RESULTS --------------------
#for each stress and the nobeww groups
stress=load_deseq('./deseqResults/stress_deseqResults.Rdata', 'stress')
heat = load_deseq('./deseqResults/heat_deseqResults.Rdata', 'heat')
heatNob = load_deseq('./deseqResults/heat_NoBEWW_deseqResults.Rdata', 'heatNoB')
cold = load_deseq('./deseqResults/cold_deseqResults.Rdata', 'cold')
coldNob = load_deseq('./deseqResults/cold_NoBEWW_deseqResults.Rdata', 'coldNoB')
salinity = load_deseq('./deseqResults/salinity_deseqResults.Rdata', 'salinity')
salinityNob = load_deseq('./deseqResults/salinity_NoBEWW_deseqResults.Rdata', 'salinityNoB')
bleach = load_deseq('./deseqResults/bleached_deseqResults.Rdata', 'bleached')
immune = load_deseq('./deseqResults/immune_deseqResults.Rdata', 'immune')
ph = load_deseq('./deseqResults/ph_deseqResults.Rdata', 'ph')

#for the 'minus stresses'
mheat = load_deseq('./deseqResults/stress_noHeat_deseqResults.Rdata', 'stress_minus_heat')
mcold = load_deseq('./deseqResults/stress_noCold_deseqResults.Rdata', 'stress_minus_heat')
msalinity = load_deseq('./deseqResults/stress_noSalinity_deseqResults.Rdata', 'stress_minus_salinity')
mbleach = load_deseq('./deseqResults/stress_noBEWW_deseqResults.Rdata', 'stress_minus_beww')
mimmune= load_deseq('./deseqResults/stress_noImmune_deseqResults.Rdata', 'stress_minus_immune')
mph = load_deseq('./deseqResults/stress_noph_deseqResults.Rdata', 'stress_minus_ph')


# LOAD VSD DATA ----------------------------------------------------------------

#load data
vsheat = load_vsd('./normCounts/heat_project_controlled.Rdata')
vsheatNob = load_vsd('./normCounts/heat_NOBEWW_project_controlled.Rdata');hNobInvert=-1
vscold = load_vsd('./normCounts/cold_project_controlled.Rdata')
vscoldNob = load_vsd('./normCounts/cold_NOBEWW_wgcna_input.Rdata')
vssalt = load_vsd('./normCounts/salinity_project_controlled.Rdata')
vssaltNob = load_vsd('./normCounts/cold_NOBEWW_wgcna_input.Rdata')
vsbleach = load_vsd('./normCounts/bleached_wgcna_input.Rdata')
vsimmune = load_vsd('./normCounts/immune_project_controlled.Rdata')
vsph = load_vsd('./normCounts/ph_project_controlled.Rdata');pcphYlims=NULL



# LOAD PREDICITON DATA ----------------------------------------------------
#stress minus pairings
mh.ppath = './category_prediction/stressMinuses/stress_noHeat_Coldata.csv__heatColdata.csv_predictions_predictionResults.Rdata'
hm.ppath = './category_prediction/stressMinuses/heatColdata.csv__stress_noHeat_Coldata.csv_predictions_predictionResults.Rdata'
mc.ppath = './category_prediction/stressMinuses/stress_noCold_Coldata.csv__coldColdata.csv_predictions_predictionResults.Rdata'
cm.ppath = './category_prediction/stressMinuses/coldColdata.csv__stress_noCold_Coldata.csv_predictions_predictionResults.Rdata'
ms.ppath = './category_prediction/stressMinuses/stress_noSalinity_Coldata.csv__salinityColdata.csv_predictions_predictionResults.Rdata'
sm.ppath = './category_prediction/stressMinuses/salinityColdata.csv__stress_noSalinity_Coldata.csv_predictions_predictionResults.Rdata'
mi.ppath = './category_prediction/stressMinuses/stress_noImmune_Coldata.csv__immuneColdata.csv_predictions_predictionResults.Rdata'
im.ppath = './category_prediction/stressMinuses/immuneColdata.csv__stress_noImmune_Coldata.csv_predictions_predictionResults.Rdata'
mp.ppath = './category_prediction/stressMinuses/stress_noph_Coldata.csv__phColdata.csv_predictions_predictionResults.Rdata'
pm.ppath = './category_prediction/stressMinuses/phColdata.csv__stress_noph_Coldata.csv_predictions_predictionResults.Rdata'
mb.ppath = './category_prediction/stressMinuses/stress_noBEWW_Coldata.csv__BEWW_Coldata.csv_predictions_predictionResults.Rdata'
bm.ppath = './category_prediction/stressMinuses/BEWW_Coldata.csv__stress_noBEWW_Coldata.csv_predictions_predictionResults.Rdata'
#stress minus nobeww
mhNob.ppath = './category_prediction/stressMinuses/stress_noHeat_Coldata.csv__heat_NOBEWW_Coldata.csv_predictions_predictionResults.Rdata'
hmNob.ppath = './category_prediction/stressMinuses/heat_NOBEWW_Coldata.csv__stress_noHeat_Coldata.csv_predictions_predictionResults.Rdata'
msNob.ppath = './category_prediction/stressMinuses/stress_noSalinity_Coldata.csv__salinity_NOBEWW_Coldata.csv_predictions_predictionResults.Rdata'
smNob.ppath = './category_prediction/stressMinuses/salinity_NOBEWW_Coldata.csv__stress_noSalinity_Coldata.csv_predictions_predictionResults.Rdata'
#individual stress pairings
hc.ppath = './category_prediction/specificStresses/heatColdata.csv__coldColdata.csv_predictions_predictionResults.Rdata'
ch.ppath = './category_prediction/specificStresses/heatColdata.csv__coldColdata.csv_predictions_predictionResults.Rdata'
hs.ppath = './category_prediction/specificStresses/heatColdata.csv__salinityColdata.csv_predictions_predictionResults.Rdata'
sh.ppath = './category_prediction/specificStresses/salinityColdata.csv__heatColdata.csv_predictions_predictionResults.Rdata'
cs.ppath = './category_prediction/specificStresses/coldColdata.csv__salinityColdata.csv_predictions_predictionResults.Rdata'
sc.ppath = './category_prediction/specificStresses/salinityColdata.csv__coldColdata.csv_predictions_predictionResults.Rdata'
hi.ppath = './category_prediction/specificStresses/heatColdata.csv__immuneColdata.csv_predictions_predictionResults.Rdata'
ih.ppath = './category_prediction/specificStresses/immuneColdata.csv__heatColdata.csv_predictions_predictionResults.Rdata'
ci.ppath = './category_prediction/specificStresses/coldColdata.csv__immuneColdata.csv_predictions_predictionResults.Rdata'
ic.ppath = './category_prediction/specificStresses/immuneColdata.csv__coldColdata.csv_predictions_predictionResults.Rdata'
si.ppath = './category_prediction/specificStresses/salinityColdata.csv__immuneColdata.csv_predictions_predictionResults.Rdata'
is.ppath = './category_prediction/specificStresses/immuneColdata.csv__salinityColdata.csv_predictions_predictionResults.Rdata'
hp.ppath = './category_prediction/specificStresses/heatColdata.csv__phColdata.csv_predictions_predictionResults.Rdata'
ph.ppath = './category_prediction/specificStresses/phColdata.csv__heatColdata.csv_predictions_predictionResults.Rdata'
cp.ppath = './category_prediction/specificStresses/coldColdata.csv__phColdata.csv_predictions_predictionResults.Rdata'
pc.ppath = './category_prediction/specificStresses/phColdata.csv__coldColdata.csv_predictions_predictionResults.Rdata'
sp.ppath = './category_prediction/specificStresses/salinityColdata.csv__phColdata.csv_predictions_predictionResults.Rdata'
ps.ppath = './category_prediction/specificStresses/phColdata.csv__salinityColdata.csv_predictions_predictionResults.Rdata'
ip.ppath = './category_prediction/specificStresses/immuneColdata.csv__phColdata.csv_predictions_predictionResults.Rdata'
pi.ppath = './category_prediction/specificStresses/phColdata.csv__immuneColdata.csv_predictions_predictionResults.Rdata'
#individual stress pairings nobeww
bhNob.ppath = './category_prediction/specificStresses/BEWW_Coldata.csv__heat_NOBEWW_Coldata.csv_predictions_predictionResults.Rdata'
hpNob.ppath = './category_prediction/specificStresses/heat_NOBEWW_Coldata.csv__BEWW_Coldata.csv_predictions_predictionResults.Rdata'
bsNob.ppath = './category_prediction/specificStresses/BEWW_Coldata.csv__salinity_NOBEWW_Coldata.csv_predictions_predictionResults.Rdata'
bsNob.ppath = './category_prediction/specificStresses/salinity_NOBEWW_Coldata.csv__BEWW_Coldata.csv_predictions_predictionResults.Rdata'
hsNob.ppath = './category_prediction/specificStresses/heat_NOBEWW_Coldata.csv__salinity_NOBEWW_Coldata.csv_predictions_predictionResults.Rdata'
shNob.ppath = './category_prediction/specificStresses/salinity_NOBEWW_Coldata.csv__heat_NOBEWW_Coldata.csv_predictions_predictionResults.Rdata'
biNob.ppath = './category_prediction/specificStresses/BEWW_Coldata.csv__immuneColdata.csv_predictions_predictionResults.Rdata'
ibNob.ppath = './category_prediction/specificStresses/immuneColdata.csv__BEWW_Coldata.csv_predictions_predictionResults.Rdata'
hiNob.ppath = './category_prediction/specificStresses/heat_NOBEWW_Coldata.csv__immuneColdata.csv_predictions_predictionResults.Rdata'
ihNob.ppath = './category_prediction/specificStresses/immuneColdata.csv__heat_NOBEWW_Coldata.csv_predictions_predictionResults.Rdata'
siNob.ppath = './category_prediction/specificStresses/salinity_NOBEWW_Coldata.csv__immuneColdata.csv_predictions_predictionResults.Rdata'
isNob.ppath = './category_prediction/specificStresses/immuneColdata.csv__salinity_NOBEWW_Coldata.csv_predictions_predictionResults.Rdata'
bpNob.ppath = './category_prediction/specificStresses/BEWW_Coldata.csv__phColdata.csv_predictions_predictionResults.Rdata'
pbNob.ppath = './category_prediction/specificStresses/phColdata.csv__BEWW_Coldata.csv_predictions_predictionResults.Rdata'
hpNob.ppath = './category_prediction/specificStresses/heat_NOBEWW_Coldata.csv__phColdata.csv_predictions_predictionResults.Rdata'
phNob.ppath = './category_prediction/specificStresses/phColdata.csv__heat_NOBEWW_Coldata.csv_predictions_predictionResults.Rdata'
spNob.ppath = './category_prediction/specificStresses/salinity_NOBEWW_Coldata.csv__phColdata.csv_predictions_predictionResults.Rdata'
psNob.ppath = './category_prediction/specificStresses/phColdata.csv__salinity_NOBEWW_Coldata.csv_predictions_predictionResults.Rdata'

