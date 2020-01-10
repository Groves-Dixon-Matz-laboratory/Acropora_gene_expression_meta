#stress_prediction_wA_and_B.R

rm(list=ls())
library(tidyverse)
library(caret)
library(glmnet)
library(randomForest)


# SETUP AND RUN RANDOM FOREST ---------------------------------------------



input = 'largeIgnored/stress_project_controlled.Rdata'
trainingFile = 'metadata/subset_tables/stratified_allStress_train.csv'
testFile = 'metadata/subset_tables/stratified_allStress_test.csv'
outcomeVar = 'stress'
posOutcome = 'stressed'
negOutcome = 'control'
outputPrefix = 'category_prediction/stress_prediction_wAandB/stratified_allStress_train.csv__stratified_allStress_test__WITH_A_AND_B.csv'



# LOAD DATA ---------------------------------------------------------------

#load stress categories
ll=load('metadata/corStressProjs.Rdata')
ll


#overall input
print('loading input data...')
ll=load(input)
ll

#make fixes to titles
stressProjs = append(corStressProjs, lowStressProjs)
unique(datTraits$my_title[!datTraits$my_title %in% stressProjs])
datTraits$my_title[datTraits$my_title=="j1_BLEACH_EWW"]<-"j1_thisStudy_PRJNA559404"
datTraits$my_title[datTraits$my_title=="b1_Aguilar_hypoosmotic stress_PRJNA380267"]<-"b1_Aguilar_hypoosmotic_stress_PRJNA380267"

d1Samples = datTraits$Run[datTraits$my_title=="d1_mohamed_microalga_PRJNA398338"]
datTraits=datTraits[datTraits$my_title!="d1_mohamed_microalga_PRJNA398338",]
datExpr = data.frame(datExpr)
datExpr = datExpr[!rownames(datExpr) %in% d1Samples,]


#build 3-way outcome vector
outcome = if_else(datTraits$treat=='control',
                         negOutcome,
                         posOutcome)
outcome = if_else(datTraits$treat!='control' & datTraits$my_title %in% corStressProjs,
                  'stress_A',
                  outcome)
outcome = if_else(datTraits$treat!='control' & datTraits$my_title %in% lowStressProjs,
                  'stress_B',
                  outcome)
unique(outcome)


#check alignment
print('Samples line up in input data?')
print(sum(rownames(datTraits)==rownames(datExpr)))==nrow(datExpr)
datExpr$outcome = outcome


#load training and test sets
trainSamples = read_csv(trainingFile) %>% 
  filter(!Run %in% d1Samples) %>% 
  pull(Run)
testSamples = read_csv(testFile) %>% 
  filter(!Run %in% d1Samples) %>% 
  pull(Run)
print('All training samples in data?')
print(sum(trainSamples %in% rownames(datTraits))==length(trainSamples))
print('All test samples in data?')
print(sum(testSamples %in% rownames(datTraits))==length(testSamples))
print('Number of overlapping samples between traingin and test set:')
print(sum(testSamples %in% trainSamples)) #should be zero


#assign training data and outcome variables
#*note this is done to avoid stack overflow protection error that occur if you use a formula with too many predictors
trainDat = as.matrix(datExpr[trainSamples, !colnames(datExpr)=='outcome'])
trainOutcome = datExpr[trainSamples, colnames(datExpr)=='outcome']
trainOutcome = factor(trainOutcome, levels=c('control', 'stress_A', 'stress_B'))

#assign test data and coutcome
testDat = as.matrix(datExpr[testSamples, !colnames(datExpr)=='outcome'])
testOutcome = datExpr[testSamples, colnames(datExpr)=='outcome']
testOutcome = factor(testOutcome, levels=c('control', 'stress_A', 'stress_B'))

#check agreement
print('All training samples accounted for?')
nrow(trainDat)==length(trainSamples)
print('All test samples accounted for?')
nrow(testDat)==length(testSamples)





# RUN RANDOM FOREST -------------------------------------------------------

#SET UP ALGORITHM
print('Running random forest...')
rf <- randomForest(x = trainDat,
                   y = trainOutcome,
                   ntree=500,
                   importance=T,
                   mtry=1000)

print(rf)
plot(rf)

#LOOK AT VARIABLE IMPORTANCE
varImpPlot(rf,
           sort = T,
           main="Variable Importance",
           n.var=50)
names(rf)
var.imp <- data.frame(rf$importance)
var.imp = var.imp[order(var.imp$MeanDecreaseAccuracy, decreasing=TRUE),]
head(var.imp)


#Now measure model accuracy for full RF model

#predict the response variable using model
predicted.response <- predict(rf, trainDat)
actual.response = trainOutcome
rfres0 = data.frame(pred = predicted.response, obs = actual.response)

#compare with actual from original dataset
print('Random forest predictions for training dataset:')
print(confusionMatrix(data=rfres0$pred,
                      reference=rfres0$obs,
                      positive="stressed"))


#now predict from validation set
predicted.response <- predict(rf, testDat)
actual.response = testOutcome
rfres = data.frame(pred = predicted.response, obs = actual.response)
mean(rfres$pred==rfres$obs)

#how does it do?
print(confusionMatrix(data=rfres$pred,
                      reference=rfres$obs,
                      positive="stress_A"))


#save the results
save(var.imp, rfres, file='category_prediction/stress_prediction_wAandB/stratified_allStress_train.csv__stratified_allStress_test__WITH_A_AND_B.csv_predictionResults_wA_and_B.Rdata')

#send the resulting object to PC for further plotting


#------------- continue here on PC
library(tidyverse)
library(caret)
library(glmnet)
library(randomForest)
ll=load('category_prediction/stress_prediction_wAandB/stratified_allStress_train.csv__stratified_allStress_test__WITH_A_AND_B.csv_predictionResults_wA_and_B.Rdata')
ll

#format a proportion accuracy table
Ns = table(rfres$obs)
cm=confusionMatrix(data=rfres$pred,
                   reference=rfres$obs,
                   positive="stress_A")
tab = cm$table
colnames(tab)
propCm = tab
propCm[,'control'] = propCm[,'control'] / Ns['control']
propCm[,'stress_A'] = propCm[,'stress_A'] / Ns['stress_A']
propCm[,'stress_B'] = propCm[,'stress_B'] / Ns['stress_B']
propCm
propCm2 = round(propCm, digits=2)


library(pheatmap)
library(RColorBrewer)
labs = c('control', 'stress\ncluster A', 'stress\ncluster B')
COLOR = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100) # display.brewer.all()

#PLOT CONFUSION MATRIX AS HEATMAP
#note this is replotted in plot_individual_projects.R with
#more formatting along with other planels for figure 3
pheatmap(propCm,
         labels_row = labs,
         labels_col=labs,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = propCm2,
         fontsize_number = 15,
         color=COLOR)







