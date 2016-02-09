setwd("/Users/beaunorgeot/systemsPharma")

drugs = read.table("drug_data.txt", sep="\t",header = T)
library(dplyr)
drugs = cbind(cLine = rownames(drugs), drugs)
rownames(drugs) = NULL
#remove drugs we won't touch
drugs = drugs %>% select(-c(drug4,drug7,drug6))


knownDrugs = drugs %>% select(-contains("drug"))
ourDrugs = drugs %>% select(drug1,drug2,drug5)

completeDrug1 = drugs %>% filter(!is.na(drug1)) %>% select(-c(drug2,drug5))
completeDrug2 = drugs %>% filter(!is.na(drug2)) %>% select(-c(drug1,drug5))
completeDrug5 = drugs %>% filter(!is.na(drug1)) %>% select(-c(drug3,drug1))

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
# Beau: start a RF preprocesing w/KNN. How to choose value of K?
library(caret)
library(doMC)
registerDoMC(cores=20)
#running on beefy machine
load("allFeatures.RData")
set.seed(666)
inTrain1 = allFeatures %>% select(-c(drug2,drug3,drug5)) %>% filter(!is.na(drug1))

inTraining1 = createDataPartition(inTrain1$drug1, p = 0.8, list = F)
trainData = inTrain1[inTraining1,] %>% select(-cLine)
testData = inTrain1[-inTraining1,] %>% select(-cLine)
trainX = trainData %>% select(-drug1)
trainY = trainData$drug1

myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)#, verboseIter = T, returnResamp = "all")
modGLMNET1 <- train(trainX, trainY,method = "glmnet", preProcess=c("knnImpute"), trControl = myControl1) 
save(modGLMNET1, ffile = "modGLMNET1.RData")

#Forest = train(trainXAA_all,trainYAA_all, method = "rf", metric = "Kappa",ntree = 1029, trControl = myControl, importance = T)
#save(aa_all_Forest, file = "aa_all_Forest.RData") #Acc = .75  Kapp = .50


runForest = train(trainXp,trainYp, method = "rf", metric = "Kappa",ntree = 1029, trControl = myControl1, importance = T)