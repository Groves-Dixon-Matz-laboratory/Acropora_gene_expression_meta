#load_corrStress_for_scatterplots.R

source('figurePlotting/rna_functions.R')

datasetPrefix = 'corrProjects'


# LOAD DESEQ RESULTS --------------------
#for each stress and the nobeww groups
stress=load_deseq('./correlated_only/deseqResults/allStress_deseqResults.Rdata', 'stress')
heat = load_deseq('./correlated_only/deseqResults/heat_deseqResults.Rdata', 'heat')
heatNob = load_deseq('./correlated_only/deseqResults/heatNOBLEACH_deseqResults.Rdata', 'heatNoB')
cold = load_deseq('./correlated_only/deseqResults/cold_deseqResults.Rdata', 'cold')
# coldNob = load_deseq('./correlated_only/deseqResults/cold_NOBEWW_deseqResults.Rdata', 'coldNoB') #this one doesn't exist for corrStress
salinity = load_deseq('./correlated_only/deseqResults/salinity_deseqResults.Rdata', 'salinity')
salinityNob = load_deseq('./correlated_only/deseqResults/salinityNOBLEACH_deseqResults.Rdata', 'salinityNoB')
bleach = load_deseq('./correlated_only/deseqResults/bleached_deseqResults.Rdata', 'bleached')
immune = load_deseq('./correlated_only/deseqResults/immune_deseqResults.Rdata', 'immune')
ph = load_deseq('./correlated_only/deseqResults/ph_deseqResults.Rdata', 'ph')

#for the 'minus stresses'
mheat = load_deseq('./correlated_only/deseqResults/stressNOHEAT_deseqResults.Rdata', 'stress_minus_heat')
mcold = load_deseq('./correlated_only/deseqResults/stressNOCOLD_deseqResults.Rdata', 'stress_minus_heat')
msalinity = load_deseq('./correlated_only/deseqResults/stressNOSALINITY_deseqResults.Rdata', 'stress_minus_salinity')
mbleach = load_deseq('./correlated_only/deseqResults/stressNOBLEACH_deseqResults.Rdata', 'stress_minus_bleach')
mimmune= load_deseq('./correlated_only/deseqResults/stressNOIMMUNE_deseqResults.Rdata', 'stress_minus_immune')
mph = load_deseq('./correlated_only/deseqResults/stressNOPH_deseqResults.Rdata', 'stress_minus_ph')


# LOAD VSD DATA ----------------------------------------------------------------

#load data
vsheat = load_vsd('./correlated_only/normCounts/heat_project_controlled.Rdata')
vsheatNob = load_vsd('./correlated_only/normCounts/heatNOBLEACH_project_controlled.Rdata');hNobInvert=1
vscold = load_vsd('./correlated_only/normCounts/cold_wgcna_input.Rdata')
# vscoldNob = load_vsd('./correlated_only/normCounts/cold_NOBEWW_vsd.Rdata') #this one doesn't exit for corrStress
vssalt = load_vsd('./correlated_only/normCounts/salinity_project_controlled.Rdata')
vssaltNob = load_vsd('./correlated_only/normCounts/salinityNOBLEACH_wgcna_input.Rdata')
vsbleach = load_vsd('./correlated_only/normCounts/bleached_project_controlled.Rdata')
vsimmune = load_vsd('./correlated_only/normCounts/immune_project_controlled.Rdata')
vsph = load_vsd('./correlated_only/normCounts/ph_wgcna_input.Rdata');pcphYlims=c(-8, 8)



# LOAD PREDICITON DATA ----------------------------------------------------
#stress minus pairings
mh.ppath = './correlated_only/category_prediction/stressMinuses/stressNOHEAT_Coldata.csv__heat_Coldata.csv_predictions_predictionResults.Rdata'
hm.ppath = './correlated_only/category_prediction/stressMinuses/heat_Coldata.csv__stressNOHEAT_Coldata.csv_predictions_predictionResults.Rdata'
mc.ppath = './correlated_only/category_prediction/stressMinuses/stressNOCOLD_Coldata.csv__cold_Coldata.csv_predictions_predictionResults.Rdata'
cm.ppath = './correlated_only/category_prediction/stressMinuses/cold_Coldata.csv__stressNOCOLD_Coldata.csv_predictions_predictionResults.Rdata'
ms.ppath = './correlated_only/category_prediction/stressMinuses/stressNOSALINITY_Coldata.csv__salinity_Coldata.csv_predictions_predictionResults.Rdata'
sm.ppath = './correlated_only/category_prediction/stressMinuses/salinity_Coldata.csv__stressNOSALINITY_Coldata.csv_predictions_predictionResults.Rdata'
mi.ppath = './correlated_only/category_prediction/stressMinuses/stressNOIMMUNE_Coldata.csv__immune_Coldata.csv_predictions_predictionResults.Rdata'
im.ppath = './correlated_only/category_prediction/stressMinuses/immune_Coldata.csv__stressNOIMMUNE_Coldata.csv_predictions_predictionResults.Rdata'
mp.ppath = './correlated_only/category_prediction/stressMinuses/stressNOPH_Coldata.csv__ph_Coldata.csv_predictions_predictionResults.Rdata'
pm.ppath = './correlated_only/category_prediction/stressMinuses/ph_Coldata.csv__stressNOPH_Coldata.csv_predictions_predictionResults.Rdata'
mb.ppath = './correlated_only/category_prediction/stressMinuses/stressNOBLEACH_Coldata.csv__bleached_Coldata.csv_predictions_predictionResults.Rdata'
bm.ppath = './correlated_only/category_prediction/stressMinuses/bleached_Coldata.csv__stressNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'

#stress minus nobeww
mhNob.ppath = './correlated_only/category_prediction/stressMinuses/stressNOHEAT_Coldata.csv__heatNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'
hmNob.ppath = './correlated_only/category_prediction/stressMinuses/heatNOBLEACH_Coldata.csv__stressNOHEAT_Coldata.csv_predictions_predictionResults.Rdata'
msNob.ppath = './correlated_only/category_prediction/stressMinuses/stressNOSALINITY_Coldata.csv__salinityNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'
smNob.ppath = './correlated_only/category_prediction/stressMinuses/salinityNOBLEACH_Coldata.csv__stressNOSALINITY_Coldata.csv_predictions_predictionResults.Rdata'
#individual stress pairings
hc.ppath = './correlated_only/category_prediction/specificStresses/heat_Coldata.csv__cold_Coldata.csv_predictions_predictionResults.Rdata'
ch.ppath = './correlated_only/category_prediction/specificStresses/heat_Coldata.csv__cold_Coldata.csv_predictions_predictionResults.Rdata'
hs.ppath = './correlated_only/category_prediction/specificStresses/heat_Coldata.csv__salinity_Coldata.csv_predictions_predictionResults.Rdata'
sh.ppath = './correlated_only/category_prediction/specificStresses/salinity_Coldata.csv__heat_Coldata.csv_predictions_predictionResults.Rdata'
cs.ppath = './correlated_only/category_prediction/specificStresses/cold_Coldata.csv__salinity_Coldata.csv_predictions_predictionResults.Rdata'
sc.ppath = './correlated_only/category_prediction/specificStresses/salinity_Coldata.csv__cold_Coldata.csv_predictions_predictionResults.Rdata'
hi.ppath = './correlated_only/category_prediction/specificStresses/heat_Coldata.csv__immune_Coldata.csv_predictions_predictionResults.Rdata'
ih.ppath = './correlated_only/category_prediction/specificStresses/immune_Coldata.csv__heat_Coldata.csv_predictions_predictionResults.Rdata'
ci.ppath = './correlated_only/category_prediction/specificStresses/cold_Coldata.csv__immune_Coldata.csv_predictions_predictionResults.Rdata'
ic.ppath = './correlated_only/category_prediction/specificStresses/immune_Coldata.csv__cold_Coldata.csv_predictions_predictionResults.Rdata'
si.ppath = './correlated_only/category_prediction/specificStresses/salinity_Coldata.csv__immune_Coldata.csv_predictions_predictionResults.Rdata'
is.ppath = './correlated_only/category_prediction/specificStresses/immune_Coldata.csv__salinity_Coldata.csv_predictions_predictionResults.Rdata'
hp.ppath = './correlated_only/category_prediction/specificStresses/heat_Coldata.csv__ph_Coldata.csv_predictions_predictionResults.Rdata'
ph.ppath = './correlated_only/category_prediction/specificStresses/ph_Coldata.csv__heat_Coldata.csv_predictions_predictionResults.Rdata'
cp.ppath = './correlated_only/category_prediction/specificStresses/cold_Coldata.csv__ph_Coldata.csv_predictions_predictionResults.Rdata'
pc.ppath = './correlated_only/category_prediction/specificStresses/ph_Coldata.csv__cold_Coldata.csv_predictions_predictionResults.Rdata'
sp.ppath = './correlated_only/category_prediction/specificStresses/salinity_Coldata.csv__ph_Coldata.csv_predictions_predictionResults.Rdata'
ps.ppath = './correlated_only/category_prediction/specificStresses/ph_Coldata.csv__salinity_Coldata.csv_predictions_predictionResults.Rdata'
ip.ppath = './correlated_only/category_prediction/specificStresses/immune_Coldata.csv__ph_Coldata.csv_predictions_predictionResults.Rdata'
pi.ppath = './correlated_only/category_prediction/specificStresses/ph_Coldata.csv__immune_Coldata.csv_predictions_predictionResults.Rdata'
#individual stress pairings nobeww
bhNob.ppath = './correlated_only/category_prediction/specificStresses/bleached_Coldata.csv__heatNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'
hbNob.ppath = './correlated_only/category_prediction/specificStresses/heatNOBLEACH_Coldata.csv__bleached_Coldata.csv_predictions_predictionResults.Rdata'
bsNob.ppath = './correlated_only/category_prediction/specificStresses/bleached_Coldata.csv__salinityNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'
sbNob.ppath = './correlated_only/category_prediction/specificStresses/salinityNOBLEACH_Coldata.csv__bleached_Coldata.csv_predictions_predictionResults.Rdata'
hsNob.ppath = './correlated_only/category_prediction/specificStresses/heatNOBLEACH_Coldata.csv__salinityNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'
shNob.ppath = './correlated_only/category_prediction/specificStresses/salinityNOBLEACH_Coldata.csv__heatNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'
biNob.ppath = './correlated_only/category_prediction/specificStresses/bleached_Coldata.csv__immune_Coldata.csv_predictions_predictionResults.Rdata'
ibNob.ppath = './correlated_only/category_prediction/specificStresses/immune_Coldata.csv__bleached_Coldata.csv_predictions_predictionResults.Rdata'
hiNob.ppath = './correlated_only/category_prediction/specificStresses/heatNOBLEACH_Coldata.csv__immune_Coldata.csv_predictions_predictionResults.Rdata'
ihNob.ppath = './correlated_only/category_prediction/specificStresses/immune_Coldata.csv__heatNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'
siNob.ppath = './correlated_only/category_prediction/specificStresses/salinityNOBLEACH_Coldata.csv__immune_Coldata.csv_predictions_predictionResults.Rdata'
isNob.ppath = './correlated_only/category_prediction/specificStresses/immune_Coldata.csv__salinityNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'
bpNob.ppath = './correlated_only/category_prediction/specificStresses/bleached_Coldata.csv__ph_Coldata.csv_predictions_predictionResults.Rdata'
pbNob.ppath = './correlated_only/category_prediction/specificStresses/ph_Coldata.csv__bleached_Coldata.csv_predictions_predictionResults.Rdata'
hpNob.ppath = './correlated_only/category_prediction/specificStresses/heatNOBLEACH_Coldata.csv__ph_Coldata.csv_predictions_predictionResults.Rdata'
phNob.ppath = './correlated_only/category_prediction/specificStresses/ph_Coldata.csv__heatNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'
spNob.ppath = './correlated_only/category_prediction/specificStresses/salinityNOBLEACH_Coldata.csv__ph_Coldata.csv_predictions_predictionResults.Rdata'
psNob.ppath = './correlated_only/category_prediction/specificStresses/ph_Coldata.csv__salinityNOBLEACH_Coldata.csv_predictions_predictionResults.Rdata'



