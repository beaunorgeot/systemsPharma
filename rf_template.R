
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

#### train up a bunch of models with a different random seed ##########
model_list = list()
for (i in 1:3){
  set.seed(i)
  nam = paste('myCars_rf',i,sep="_")
  assign(nam,train(am ~. , data = mtcars, method = "rf", metric = "Kappa",
                   tuneLength = 5, ntree = 5, proximity = T, importance = T))
  model_list[i] = nam
}

# extract a list of the vars by importance as a seperate df for each model
for (i in 1:length(model_list)){
  vars = paste('myCars_rf_vars',i,sep="_")
  assign(vars,data.frame(varImp(get(model_list[[i]]))$importance))
  print(class(get(vars)))
  tempTab = get(vars)
  tempTab[paste("myVars", i, sep="_")] <- row.names(tempTab)
  colnames(tempTab) = paste(colnames(tempTab),i,sep = "_")
  assign(vars,tempTab)
}
#bind the important vars into 1 df
all_vars = myCars_rf_vars_1
for (i in 2:length(model_list)){
  vars = paste('myCars_rf_vars',i,sep="_")
  all_vars = cbind(all_vars, get(vars))
}

n = 5
#you will still need to change the response, in the first select from "auto" to whatever you want
my_avgTop_vars = all_vars %>% select(Vars = myVars_1_1,contains("auto")) %>% mutate(avg_imp = rowSums(.[, -1])/length(model_list)) %>% arrange(desc(avg_imp)) %>% top_n(n = n,wt=avg_imp) %>% select(Vars,avg_imp)


########## below is the line by line that I used to build the programmatic approach above ###############
impVar_bob_1<-data.frame(varImp(bob_1)$importance)
impVar_bob_1$Vars_bob<-row.names(impVar_bob_1)
colnames(impVar_bob_1) = paste(colnames(impVar_bob_1),"1", sep = "_")

impVar_bob_2<-data.frame(varImp(bob_2)$importance)
impVar_bob_2$Vars_bob<-row.names(impVar_bob_2)
colnames(impVar_bob_2) = paste(colnames(impVar_bob_2),"2", sep = "_")

impVar_bob_3<-data.frame(varImp(bob_3)$importance)
impVar_bob_3$Vars_bob<-row.names(impVar_bob_3)
colnames(impVar_bob_3) = paste(colnames(impVar_bob_3),"3", sep = "_")
##################
cars_vars = cbind(myCars_rf_vars_1,myCars_rf_vars_2,myCars_rf_vars_3)
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

################################################
##plotting ##
featurePlot(train, outcome.org, "strip")
# correlation plot: pretty cool, see https://rpubs.com/flyingdisc/practical-machine-learning-xgboost
corrplot.mixed(cor(train), lower="circle", upper="color", 
               tl.pos="lt", diag="n", order="hclust", hclust.method="complete")
