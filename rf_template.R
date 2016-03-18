
library(caret)
library(dplyr)
library(doMC)
registerDoMC(cores = 54)

#set the NAs when reading in a file
final_df = read.csv("df_beau.csv",na.strings=c("",".","NA"))
# read files faster with readr
knock_down_path = "https://github.com/dhimmel/lincs/raw/abcb12f942f93e3ee839e5e3593f930df2c56845/data/consensi/signif/dysreg-knockdown.tsv"
kd_genes = readr::read_tsv(knock_down_path)
kd_genes$direction = as.factor(kd_genes$direction)

library(tidyr) # for separate()
stud_adjust_responders = separate(stud_adjust_responders, into = c("a1","a2","a3") ,col = study_id, sep= "-", remove = T)
# spread() example
answer_df = kd_genes %>%
  group_by(symbol, direction) %>%
  dplyr::summarise(count=n()) %>%
  tidyr::spread(key = direction, value = count, fill = 0) %>%
  # spread: take the direction column, and create a new column for each factor that was in direction. populate each new
  # column with the corresponding "value"s from the column titled "count"
  arrange(desc(up)) %>% top_n(n = 10, wt = desc(up)) # remember that top_n() doesn't do scientific notation. 
# If your values are sci-note then stop after the arrange() call and do
# answer_df[1:10,]

#search within columnnames to find possible matches(without actually selecting them)
colnames(clinical)[grep("bmi",colnames(clinical))]

# how many fast/slow responders received 1 antibiotic vs many antis
# how many fast responders got mutliple antis, how many got only 1 anti? Same for slow responders
someshit = final_df %>% 
  group_by(CAT_response,antitb) %>%
  summarise(n = n())

set.seed(42)
inTrain <- createDataPartition(inTraining$Survived, p = 0.8, list = F)
trainSet <- inTraining[inTrain,]
training = trainSet %>% select(-c(yourResponse))
train_response = trainSet$yourResponse
testSet <- inTraining[-inTrain,]

set.seed(666)
rocControl = trainControl(method = "repeatedcv", number = 10, repeats = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
final_all_rf = train(training, train_response, method= "rf", metric = "ROC",tuneLength =4,ntree = 1500, importance = T, trControl = myControl_f)

kappaControl <- trainControl(method = "repeatedcv", number = 10, repeats = 5 )
bob_1 = train(am ~.,data = mtcars,method = "rf", metric = "Kappa",tuneLength =5,ntree = 500, proximity = T, importance = T, type = 1,trControl = kappaControl)

############# Recursive Feature Elimination ###################
predictors_df <- training # Make sure you're only doing this on your training data!!!!
labels <- train_response
set.seed(1)
rfe(
  predictors_df[complete.cases(predictors_df),],
  labels[complete.cases(predictors_df)],
  #sizes = seq(2, 10, by = 2),
  sizes = seq(2, length(names(training)), by = 4), # start with all the vars(length(names(training))), go all the way down to 2 vars, subtracting 4 each time
  rfeControl = rfeControl(functions = rfFuncs, method = 'cv', number = 5) #rfFuncs means use a randomForest
)


###############################
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

########################### #take the top n-variables from a single run ################
ImpMeasure<-data.frame(varImp(yourModel)$importance)
ImpMeasure$Vars<-row.names(ImpMeasure)
n = 20
Top_vars = ImpMeasure[order(ImpMeasure$Overall, decreasing = T),][1:n,]

#which vars appear in the top 20 for glmnet and forest?
super_vars = intersect(forest_vars,glmnet_vars)

#set 20 different random seeds, use voting to select the top N features. Use only those features in the final model. 
# less features, less overfitting. 

#### train up a bunch of models with a different random seed ##########
mtcars
model_list = list()
for (i in 1:4){
  set.seed(i)
  nam = paste('rf',i,sep="_")
  assign(nam,train(am ~. , data = mtcars, method = "rf", metric = "Kappa",
                   tuneLength = 5, ntree = 5, proximity = T, importance = T))
  model_list[i] = nam
}
# extract a list of the vars by importance as a seperate df for each model
for (i in 1:length(model_list)){
  vars = paste('important_rf_vars',i,sep="_")
  assign(vars,data.frame(varImp(get(model_list[[i]]))$importance))
  print(class(get(vars)))
  tempTab = get(vars)
  tempTab[paste("myVars", i, sep="_")] <- row.names(tempTab)
  colnames(tempTab) = paste(colnames(tempTab),i,sep = "_")
  assign(vars,tempTab)
}

#bind the important vars into 1 df
all_vars = important_rf_vars_1
for (i in 2:length(model_list)){
  vars = paste('important_rf_vars',i,sep="_")
  all_vars = cbind(all_vars, get(vars))
}

# either extract the top n important vars, or remove the top_n() step to simply return a sorted list of vars with their relative importance
# YOU NEED TO CHANGE THE RESPONSE select(contain("your actual variable")),, in the first select from "auto" to whatever you want
# look at all_vars to choose a name
n = 5
my_avgTop_vars = all_vars %>% select(Vars = myVars_1_1,contains("X1")) %>% mutate(avg_imp = rowSums(.[, -1])/length(model_list)) %>% arrange(desc(avg_imp)) %>% select(Vars,avg_imp) %>% top_n(n = n,wt=avg_imp)
#now use only those vars to train a model
training = mtcars[,my_avgTop_vars$Vars] #subset the original df on the columns you want

######################### Choose the number of vars by trying different numbers of vars and assessing out of sample results####
# the idea here is first use my variable importance permutation method above to produce a list of variables ranked by importance
# then take the top 5, 10, 15, etc and build models with them, testing the models out of sample to search for where out of sample
# performance either plateus or begins do decline (overfitting)
confus_acc = list()# assign lists to hold confusionMatrix accurcacy and kappa
confus_kapp = list()# here
for (i in seq(10,200,10)){
  print(paste("starting forest",i/10,sep = " "))
  pvalues$SNP = as.character(pvalues$SNP)
  snps = pvalues %>% arrange(desc(P)) # arrange your vars in order from largest to smallest
  snps = snps[1:i,] # take the top_i variables
  train_response = genotypes_train$bdrGp
  nam = genotypes_train[,snps$SNP] #subset your training data to only include the columns that are the top_i variables
  nam$Pre.FEV1.perc.pred = genotypes_train$Pre.FEV1.perc.pred # I already knew this was my favorite predictor, so I wanted to added it back in, even though it wasn't a ranked variable
  truth = genotypes_test$bdrGp
  test = genotypes_test[,snps$SNP] #must also subset the testSet to only include the top i variables
  test$Pre.FEV1.perc.pred = genotypes_test$Pre.FEV1.perc.pred
  forest = paste('pvalue_rf',i,sep=".") # create a string that will be assigned to forest
  #assign the name pvalue_rf.1 to the train object. this is the programatic way to do: pvalue_rf.1 = train()
  assign(forest,train(nam,train_response, method = "rf", metric = "Kappa", tuneLength = 3, ntree = 1500, importance = T, trControl = myControl))
  pred = paste('pred_rf',i,sep=".") # create a string that will be assigned to the predictions
  assign(pred,predict(get(forest),test)) # do predictions and assign them to the string. forest == "just a string". You must get(forest) to get the train object
  confus = paste('confus_rf',i, sep = ".") # repeat the prediction string assignment for the confusionMatrix
  assign(confus,confusionMatrix(get(pred),truth))
  accKapp = get(confus) # get the confusionMatrix and assign it to a variable
  confus_acc[[i/10]] = accKapp$overall[[1]] # add the ith accuracy value to the list of confusion matrix accuracy values
  confus_kapp[[i/10]] = accKapp$overall[[2]]# same for kappa
  save(list=c(forest),file = paste("pvalue_rf",i,"RData", sep = ".")) # again, forest,pred,confus are just strings. However list=forest actually saves the forest, but the save name needs to be generated, by some R voodoo
  save(list=c(pred), file = paste("pred_rf",i,"RData",sep="."))
  save(list=c(confus), file = paste("confus_rf",i,"RData",sep="."))
}
# the full script this came from was: ~burchardRotation/80_20/pr_pvalue/pr_by_pvalue.R
# to see how to do this with glmnet and dummies vars: ~burchardRotation/80_20/pr_pvalue/glmnet_pr_pvalue.R
save(confus_acc, file = "confus_acc.RData") #here
save(confus_kapp, file = "confus_kapp.RData")# here
###### out of sample plotting ##################
library(tidyr)
accANDkapp = data.frame(numVars = seq(10,200,10),Accuracy = unlist(confus_acc),Kappa = unlist(confus_kapp))
a_long = accANDkapp %>% gather(statistic, value, -numVars)
# This says "Gather all variables except numVars, calling the new key column 'statistic' and the new value column 'value'."
library(ggplot2)
pr_oob_acc_kapp_plot = ggplot() + geom_line(data = a_long, mapping = aes(x = numVars, y = value,colour = statistic)) + geom_point(data = a_long,aes(x = numVars,y = value,colour = statistic))
#options(bitmapType = "cairo") # need this to save in unix env (cesar) #otherwise get X11 related error. Not necessary when in Rstudio
ggsave(pr_oob_acc_kapp_plot,file = "pr_oob_acc_kapp_plot.png", width = 5, height = 5)
###############
############### inSample plotting ##########
resamps <- resamples(list(V_10 = pvalue_rf.10,V_20 = pvalue_rf.20), V_30 = pvalue_rf.30,V_40 = pvalue_rf.40,V_50 = pvalue_rf.50,
                     V_60 = pvalue_rf.60,V_70 = pvalue_rf.70,V_80 = pvalue_rf.80,V_90 = pvalue_rf.90,V_100 = pvalue_rf.100,
                     V_110 = pvalue_rf.110,V_120 = pvalue_rf.120,V_130 = pvalue_rf.130,V_140 = pvalue_rf.140,V_150 = pvalue_rf.150,
                     V_160 = pvalue_rf.160,V_170 = pvalue_rf.170,V_180 = pvalue_rf.180,V_190 = pvalue_rf.190,V_200 = pvalue_rf.200)
bwplot(resamps, layout = c(2, 1))
trellis.par.set(caretTheme())
dotplot(resamps, metric = "Accuracy")

#can get the mean oob error (mean of all trees)
mean(pvalue_rf.10$finalModel$err.rate[,1])
# get the confusion matrix for the final model. This provides misclaffication rates
pvalue_rf.10$finalModel$confusion
# see the misclassification rates for each class in that model
pvalue_rf.10$finalModel$confusion[,3]




##plotting ##
featurePlot(train, outcome.org, "strip")
# correlation plot: pretty cool, see https://rpubs.com/flyingdisc/practical-machine-learning-xgboost
corrplot.mixed(cor(train), lower="circle", upper="color", 
               tl.pos="lt", diag="n", order="hclust", hclust.method="complete")

####################################
# ggploting

library(ggplot2)
#Box plot
b = ggplot() + geom_boxplot(data = c, aes(x=gear, y = mpg, color = gear))
# facet box plot, 1 planel for each gear
bFacet1 <-ggplot() + facet_wrap(~gear) + geom_boxplot(data=c, mapping=aes(x = gear, y= mpg, color = gear))
# facet box plot, 1 planel for each number of cylinders, also fill in the boxes. Notices this makes the median bar disappear :(
bFacet2 <-ggplot() + facet_wrap(~cyl) + geom_boxplot(data=c, mapping=aes(x = gear, y= mpg, color = gear, fill = gear))
# facet by 2 variables, cyl and am(transmission type)
c$am = factor(c$am, levels = c(0,1), labels = c("manualTrans", "autoTrans"))
bFacet3 = ggplot() + facet_grid(cyl ~ am) + geom_boxplot(data = c, aes(x = gear, y = mpg, color = gear))
#flip the x & y coordinates
bFlip = ggplot() + geom_boxplot(data = c, aes(x=gear, y = mpg, color = gear)) + coord_flip()

####factors####
c = mtcars
#mpg by number of gears
c$gear = factor(c$gear)
c$cyl = factor(c$cyl, levels = c(4,6,8),labels=c("4Cyl","6Cyl","8Cyl"))

### optimizing training time #############################################
#The overall complexity of RF is something like ntree⋅mtry⋅(# objects)log(# objects)ntree⋅mtry⋅(# objects)loga(# objects); if you want to speed your computations up, you can try the following:
  
  #1.Use randomForest instead of party, or, even better, ranger or Rborist (although both are not yet battle-tested).
  #2.Don't use formula, i.e. call randomForest(predictors,decision) instead of randomForest(decision~.,data=input).
#3. Use do.trace argument to see the OOB error in real-time; this way you may detect that you can lower ntree.
#4. About factors; RF (and all tree methods) try to find an optimal subset of levels thus scanning 2(# of levels-1)2(# of levels-1) possibilities; to this end it is rather naive this factor can give you so much information -- not to mention that randomForest won't eat factors with more than 32 levels. Maybe you can simply treat it as an ordered one (and thus equivalent to a normal, numeric variable for RF) or cluster it in some groups, splitting this one attribute into several?


######## DUMMY VARS ############
glmnet_raw_predictors = clinical_train %>% select(-c(personID,age,bdrGp))
glmnet_response = clinical_train$bdrGp
#cast into dummies
dummies = dummyVars("~.", data=glmnet_raw_predictors)
# turn dummies into new df
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
#train
myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(666)
glmnet_clinical <- train(dummy_glmnet_train, glmnet_response, method = "glmnet", family = "binomial", trControl = myControl1)
#test it out of sample
testy = clinical_test %>% select(-c(bdrGp,age,personID))
# make dummies dummy
test_dummies = dummyVars("~.", data=testy)
testy = as.data.frame(predict(test_dummies,newdata = testy))
#make predictions
glmnet_clinical_predictions = predict(glmnet_clinical,testy)
glmnet_clinical_confus = confusionMatrix(glmnet_clinical_predictions,clinical_test$bdrGp) 
