


library(tidyverse)
library(caret)
library(glmnet)
library(randomForest)


input = 'stress_project_controlled.Rdata'
trainingFile = 'stress_noBEWW_Coldata.csv'
testFile = 'bleached_Coldata.csv'
outcomeVar = 'stress'
posOutcome = 'stressed'
negOutcome = 'control'
outputPrefix = 'stressNoBEWW_BEWW_predictions'



# LOAD DATA ---------------------------------------------------------------

#overall input
print('loading input data...')
ll=load(input)
datExpr = data.frame(datExpr)
outcome = factor(if_else(datTraits$treat=='control', negOutcome, posOutcome), levels = c(posOutcome, negOutcome))
print('Samples line up in input data?')
print(sum(rownames(datTraits)==rownames(datExpr)))==nrow(datExpr)
datExpr$outcome = outcome


#load training and test sets
trainSamples = read_csv(trainingFile) %>% 
  pull(Run)
testSamples = read_csv(testFile) %>% 
  pull(Run)
print('All training samples in data?')
print(sum(trainSamples %in% rownames(datTraits))==length(trainSamples))
print('All test samples in data?')
print(sum(testSamples %in% rownames(datTraits))==length(testSamples))
print('Number of overlapping samples between traingin and test set:')
print(sum(testSamples %in% trainSamples))


#assign training data and outcome variables
#*note this is done to avoid stack overflow protection error that occur if you use a formula with too many predictors
trainDat = as.matrix(datExpr[trainSamples, !colnames(datExpr)=='outcome'])
trainOutcome = datExpr[trainSamples, colnames(datExpr)=='outcome']

#assign test data and coutcome
testDat = as.matrix(datExpr[testSamples, !colnames(datExpr)=='outcome'])
testOutcome = datExpr[testSamples, colnames(datExpr)=='outcome']

#check agreement
print('All training samples accounted for?')
nrow(trainDat)==length(trainSamples)
print('All test samples accounted for?')
nrow(testDat)==length(testSamples)



# RUN LOGISTIC REGRESSION -------------------------------------------------
set.seed(12345)

#set up input variables from training data
x <- trainDat 
y <- if_else(trainOutcome==posOutcome, 1, 0)

#get optimal lambda value
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
pdf(paste(outputPrefix, 'lassoPlot.pdf', sep='_'))
plot(cv.lasso)
dev.off()

# Fit the final model on the training data
model <- glmnet(x, y, alpha = 1, family = "binomial",
                lambda = cv.lasso$lambda.1se)

# Display regression coefficients
cdat = data.frame(coef(model)[,1])
colnames(cdat) = 'coef'
print('Total predictors:')
print(nrow(cdat))
print('Total non-zero coefficients:')
print(sum(cdat[,1] != 0))

# Make predictions on the test data
x.test <- testDat
probabilities <- model %>% predict(newx = x.test)
predicted.classes <- factor(ifelse(probabilities > 0.5, posOutcome, negOutcome), levels=c(posOutcome, negOutcome))

# Model accuracy basic
observed.classes <- factor(testOutcome, levels=c(posOutcome, negOutcome))
lres = data.frame(pred = predicted.classes, obs = observed.classes)
mean(lres$pred == lres$obs)
#Model accuracy with caret
print(confusionMatrix(data=lres$pred,
                      reference=lres$obs,
                      positive=posOutcome))

# RUN RANDOM FOREST -------------------------------------------------------

#SET UP ALGORITHM
print('Running random forest...')
rf <- randomForest(x = trainDat,
                   y = trainOutcome,
                   ntree=500,
                   importance=T,
                   mtry=1000)

print(rf)
pdf(paste(outputPrefix, 'randomForestPlot.pdf'))
plot(rf)
dev.off()

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
                      positive="stressed"))


#save the results
save(cdat, lres, var.imp, rfres, file=paste(outputPrefix, 'predictionResults.Rdata', sep='_'))