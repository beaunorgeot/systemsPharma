setwd("/Users/beaunorgeot/systemsPharma")

drugs = read.table("drug_data.txt", sep="\t",header = T)
library(dplyr)
drugs = cbind(cLine = rownames(drugs), drugs)
rownames(drugs) = NULL
#remove drugs we won't touch
drugs = drugs %>% select(-c(drug3,drug4,drug7,drug6))


knownDrugs = drugs %>% select(-contains("drug"))
ourDrugs = drugs %>% select(drug1,drug2,drug5)

completeDrug1 = drugs %>% filter(!is.na(drug1)) %>% select(-c(drug2,drug5))
completeDrug2 = drugs %>% filter(!is.na(drug2)) %>% select(-c(drug1,drug5))
completeDrug5 = drugs %>% filter(!is.na(drug5)) %>% select(-c(drug2,drug1))

cnv = read.table("cnv_data.txt", sep="\t",header = T)
cnv = cbind(cLine = rownames(cnv), cnv)
rownames(cnv) = NULL
cnv1 = cnv %>% filter(cnv$cLine %in% completeDrug1$cLine)

geneEx = read.table("exp_data.txt", sep="\t",header = T)
geneX = cbind(cLine = rownames(geneEx), geneEx)
rownames(geneX) = NULL
geneX1 = geneX %>% filter(geneX$cLine %in% completeDrug1$cLine)

# next 
# 1. mutate all column names in cnv adding "c_" to each
names(cnv) = paste("c",names(cnv),sep = "_")
names(geneX) = paste("g",names(geneX),sep = "_")
names(cnv)[1] = "cLine"
names(geneX)[1] = "cLine"
allFeatures = inner_join(cnv,geneX)
allFeatures = inner_join(allFeatures,drugs)
save(allFeatures, file = "allFeatures.RData")

sum(is.na(allFeatures)) # ~ 4% of the data is missing

# are the unknown drugs correlated to any of the known drugs?
only_drugs = drugs %>% select(-cLine)
corrs = cor(only_drugs, use = "pairwise.complete.obs")
sigs = which(corrs > .5 & corrs != 1,arr.ind = T)


####################### Models ####################################
# Beau: 
library(plyr)
library(caret)
library(dplyr)
library(doMC)
registerDoMC(cores=40)
#running on beefy machine
load("allFeatures.RData")
set.seed(666)
inTrain1 = allFeatures %>% select(-c(drug2,drug3,drug5)) %>% filter(!is.na(drug1))

inTraining1 = createDataPartition(inTrain1$drug1, p = 0.8, list = F)
trainData = inTrain1[inTraining1,] %>% select(-cLine)
testData = inTrain1[-inTraining1,] %>% select(-cLine)
trainX = trainData %>% select(-drug1)
trainY = trainData$drug1

myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
modGLMNET1 <- train(trainX, trainY,method = "glmnet", preProcess=c("medianImpute"), trControl = myControl1) 
save(modGLMNET1, file = "modGLMNET1.RData")

myControl1rf <- trainControl(method = "repeatedcv", number = 10, repeats = 5, verboseIter = T, returnResamp = "all")
Forest1 = train(trainX,trainY, method = "rf", metric = "RMSE",ntree = 1029, preProcess=c("medianImpute"), trControl = myControl1rf, importance = T)
save(Forest1, file = "Forest1.RData")

#annoylingly, median imputation isn't working for the GBM in preProcessing, doing it manually..
trainData = trainData %>% mutate_each(funs(ifelse(is.na(.),median(.,na.rm = TRUE),.)))
tnGrid = expand.grid(interaction.depth=c(1,2), n.trees=100, shrinkage=c(0.001, 0.01), n.minobsinnode=10)
gbm1 = train(as.matrix(trainX),trainY, method = "gbm", metric = "RMSE", trControl = myControl1,verbose = F,tuneGrid = tnGrid)
save(gbm1, file = "gbm1.RData")
#below works
gbm1 = train(as.matrix(trainX),trainY, method = "gbm", metric = "RMSE", trControl = myControl1,verbose = F,tuneGrid = jake)
jake = expand.grid(interaction.depth=c(1,2), n.trees=seq(1,15), shrinkage=c(0.001, 0.01), n.minobsinnode=3)

#########giving up ######### trying a different package

#XGBoost supports only numeric matrix data. Converting all training, testing and outcome data to matrix.
testData = inTrain1[-inTraining1,] %>% select(-cLine)
testData = as.matrix(testData)
trainY = trainData$drug1
boost = train(as.matrix(trainX),trainY, method = "xgbLinear", metric = "RMSE", trControl = myControl1)
save(boost, file = "boost.RData")


########
#impute the missing values for the test data
testData_im = testData %>% mutate_each(funs(ifelse(is.na(.),median(.,na.rm = TRUE),.)))


load("modGLMNET1.RData")
#predict.glmnet requires that you remove the reponse column from testing data, it fails with 
#  Error in as.matrix(cbind2(1, newx) %*% nbeta) .... Cholmod error 'X and/or Y have wrong dimensions' at file
# if the column that you want to predict already has values in it
testy = testData_im %>% select(-drug1)
glmnet1_predictions = predict(modGLMNET1,testy)
test_truth = data.frame(obs = testData_im$drug1, pred = glmnet1_predictions)
glmnet1_summary = defaultSummary(test_truth) #RMSE = .33, Rsquared = .84
save(glmnet1_summary, file = "glmnet1_summary.RData")

load("Forest1.RData")
forest1_predictions = predict(Forest1,testData_im)
test_truth = data.frame(obs = testData_im$drug1, pred = forest1_predictions)
forest1_summary = defaultSummary(test_truth)
save(forest1_summary, file = "forest1_summary.RData") #RMSE = .47, R2 = .74

########### final Models and predictions ######################
library(plyr)
library(caret)
library(dplyr)
library(doMC)
registerDoMC(cores=40)
#running on beefy machine
load("allFeatures.RData")
set.seed(666)
inTrain1 = allFeatures %>% select(-c(cLine,drug2,drug3,drug5)) %>% filter(!is.na(drug1))
trainX = inTrain1 %>% select(-drug1)
trainY = inTrain1$drug1

myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
glmnet_drug1 <- train(trainX, trainY,method = "glmnet", preProcess=c("medianImpute"), trControl = myControl1) 
save(glmnet_drug1, file = "glmnet_drug1.RData") #R2: .72

test1 = allFeatures %>% select(-c(cLine,drug2,drug3,drug5)) %>% filter(is.na(drug1)) %>% mutate_each(funs(ifelse(is.na(.),median(.,na.rm = TRUE),.))) %>% select(-drug1)
glmnet1_predictions = predict(glmnet_drug1,test1)
missingDrug1 = allFeatures %>% select(cLine,drug1) %>% filter(is.na(drug1)) %>% select(cLine)
drug1prediction = cbind(missingDrug1,glmnet1_predictions) %>% arrange(glmnet1_predictions)
save(drug1prediction, file = "drug1prediction.RData")

#### 2 ####
inTrain2 = allFeatures %>% select(-c(cLine,drug1,drug3,drug5)) %>% filter(!is.na(drug2))
trainX = inTrain2 %>% select(-drug2)
trainY = inTrain2$drug2

glmnet_drug2 <- train(trainX, trainY,method = "glmnet", preProcess=c("medianImpute"), trControl = myControl1) # RMSE 0.6836211  R2: 0.3924393
save(glmnet_drug2, file = "glmnet_drug2.RData")
## glmnet didn't do nearly as well on drug 2, how does an rf do in sample? ## about the same
rf_drug2 <- train(trainX, trainY,method = "rf", preProcess=c("medianImpute"), trControl = myControl1rf) # RMSE:6952442   R2:0.40
save(rf_drug2, file = "rf_drug2.RData")
###

test2 = allFeatures %>% select(-c(cLine,drug1,drug3,drug5)) %>% filter(is.na(drug2)) %>% mutate_each(funs(ifelse(is.na(.),median(.,na.rm = TRUE),.))) %>% select(-drug2)
glmnet2_predictions = predict(glmnet_drug2,test2)
missingDrug2 = allFeatures %>% select(cLine,drug2) %>% filter(is.na(drug2)) %>% select(cLine)
drug2prediction = cbind(missingDrug2,glmnet2_predictions) %>% arrange(glmnet2_predictions)
save(drug2prediction, file = "drug2prediction.RData")
#### 5 ####
inTrain5 = allFeatures %>% select(-c(cLine,drug1,drug2,drug3)) %>% filter(!is.na(drug5))
trainX = inTrain5 %>% select(-drug5)
trainY = inTrain5$drug5

glmnet_drug5 <- train(trainX, trainY,method = "glmnet", preProcess=c("medianImpute"), trControl = myControl1) # RMSE:0.8136636  R2:0.4132502
save(glmnet_drug5, file = "glmnet_drug5.RData")

test5 = allFeatures %>% select(-c(cLine,drug1,drug2,drug3)) %>% filter(is.na(drug5)) %>% mutate_each(funs(ifelse(is.na(.),median(.,na.rm = TRUE),.))) %>% select(-drug5)
glmnet5_predictions = predict(glmnet_drug5,test5)
missingDrug5 = allFeatures %>% select(cLine,drug5) %>% filter(is.na(drug5)) %>% select(cLine)
drug5prediction = cbind(missingDrug5,glmnet5_predictions) %>% arrange(glmnet5_predictions)
save(drug5prediction, file = "drug5prediction.RData")

############ visualizations ###########
plot(modGLMNET1) # Mixing % vs RMSE
plot(modGLMNET1$finalModel) #coeffecients


############# variable selection for Erlotinib #####################
# There are a couple of ways to go about this
# 1. use RandomForest, build about 20 different models (each w/a different seed), and use a voting system to select the top
# vars that show up the most frequently. Could also do this only 3 times to save time, etc
# 2. Build a few RF's and a few glmnet's (each w/different seeds). Choose vars that only show up as important vars w/both model types, on the logic
# that anything that is predictive in both an RF and a linear model is super important. This can be a good method to get a small number of vars to 
# explore for biological understanding, but maybe isn't a great method for doing future predictions

library(plyr)
library(caret)
library(dplyr)
library(doMC)
registerDoMC(cores=40)
#running on beefy machine
load("allFeatures.RData")
set.seed(111)
inTrainFS = allFeatures %>% select(-contains("drug")) %>% filter(!is.na(Erlotinib)) %>% select(-cLine)

trainX = inTrainFS %>% select(-Erlotinib)
trainY = inTrainFS$Erlotinib

myControl1rf <- trainControl(method = "repeatedcv", number = 10, repeats = 5, verboseIter = T, returnResamp = "all")
Forest_fs = train(trainX,trainY, method = "rf", metric = "RMSE",ntree = 1029, preProcess=c("medianImpute"), trControl = myControl1rf, importance = T)
save(Forest_fs, file = "Forest_fs.RData")
# extract a list of the vars by importance
ImpMeasure<-data.frame(varImp(Forest_fs)$importance)
ImpMeasure$Vars<-row.names(ImpMeasure)
#take the top n-variables
n = 20
forest_vars = ImpMeasure[order(ImpMeasure$Overall,decreasing = T),][1:n,]

# now train a glmnet on erot
glmnet_fs <- train(trainX, trainY,method = "glmnet", preProcess=c("medianImpute"), trControl = myControl1)
save(glmnet_fs, file = "glmnet_fs.RData")
# extract a list of the vars by importance
ImpMeasure<-data.frame(varImp(glmnet_fs)$importance)
ImpMeasure$Vars<-row.names(ImpMeasure)
#take the top n-variables
n = 20
glmnet_vars = ImpMeasure[order(ImpMeasure$Overall,decreasing = T),][1:n,]

#which vars appear in the top 20 for glmnet and forest?
super_vars = intersect(forest_vars,glmnet_vars) #there's no overlap :()