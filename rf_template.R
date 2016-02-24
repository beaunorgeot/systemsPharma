
library(caret)
library(dplyr)
bob_1 = train(am ~.,data = mtcars,method = "rf", metric = "Kappa",tuneLength =5,ntree = 50, proximity = T, importance = T, type = 1)
# get the proximity scores for the samples
bob_prox = bob$finalModel$proximity
#look for subclasses within the samples. first split samples by class
# not sure how to do this. Perhaps using RF to cluster? see:http://www.bios.unc.edu/~dzeng/BIOS740/randomforest.pdf

varImp(bob)
bob_2 = train(am ~.,data = mtcars,method = "rf", metric = "Kappa",tuneLength =5,ntree = 50, proximity = T, importance = T, type = 2)
varImp(bob2)
# CHOOSING DIFFERENT TYPE OF IMPORTANCE CHANGES THE ORDER AND VALUE OF THE IMPORTANT FEATURES!!
# type1 = mean decrease in accuracy; type2 = mean decrease in node impurity.
set.seed(3)
bob_3 = train(am ~.,data = mtcars,method = "rf", metric = "Kappa",tuneLength =5,ntree = 50, proximity = T, importance = T, type = 1)

# extract a list of the vars by importance
impVar_bob_1<-data.frame(varImp(bob_1)$importance)
impVar_bob_1$Vars_bob<-row.names(impVar_bob_1)
colnames(impVar_bob_1) = paste(colnames(impVar_bob_1),"1", sep = "_")

impVar_bob_2<-data.frame(varImp(bob_2)$importance)
impVar_bob_2$Vars_bob<-row.names(impVar_bob_2)
colnames(impVar_bob_2) = paste(colnames(impVar_bob_2),"2", sep = "_")

impVar_bob_3<-data.frame(varImp(bob_3)$importance)
impVar_bob_3$Vars_bob<-row.names(impVar_bob_3)
colnames(impVar_bob_3) = paste(colnames(impVar_bob_3),"3", sep = "_")
#join and average the results
bobs_impVar = cbind(impVar_bob_1,impVar_bob_2,impVar_bob_3)
n = 5
bobs_avgTop_vars = bobs_impVar %>% mutate(Vars = Vars_bob_1,avg_imp = (auto_1 + auto_2 + auto_3)/3) %>% arrange(desc(avg_imp)) %>% top_n(n = n,wt=avg_imp) %>% select(Vars,avg_imp)


########################### #take the top n-variables from a single run ################
ImpMeasure<-data.frame(varImp(yourModel)$importance)
ImpMeasure$Vars<-row.names(ImpMeasure)
n = 20
Top_vars = ImpMeasure[order(ImpMeasure$Overall, decreasing = T),][1:n,]

#which vars appear in the top 20 for glmnet and forest?
super_vars = intersect(forest_vars,glmnet_vars)

#set 20 different random seeds, use voting to select the top N features. Use only those features in the final model. 
# less features, less overfitting. 

##randoms ##
featurePlot(train, outcome.org, "strip")
# correlation plot: pretty cool, see https://rpubs.com/flyingdisc/practical-machine-learning-xgboost
corrplot.mixed(cor(train), lower="circle", upper="color", 
               tl.pos="lt", diag="n", order="hclust", hclust.method="complete")
