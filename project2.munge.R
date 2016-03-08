
setwd("/Users/beaunorgeot/systemsPharma")
dataStud = read.csv("data_project2_pspg245B.csv"); dim(dataStud) # 866 x 71
pkData = read.csv("pk.data.rpt.csv"); dim(pkData) #383 x 7
summary(dataStud) # lots of non data at the bottom
dataStud = dataStud[1:816,] # row 462 col 30 is NA. dataStud$csrclass
table(pkData$STUD) # the 2 studies are 29, 299. 299 must = 29X, the escalting dose model
table(pkData$DOSE) # 5 different doses(450:60,600:211,900:78,1200:30,1500:4)
table(pkData$DOSK) # 3 different doses/kg (10:236, 15:74, 20:73)

stud_adjust = dataStud
stud_adjust$mtb_l_on_days_sc = as.character(stud_adjust$mtb_l_on_days_sc)
sum(stud_adjust$mtb_l_on_days_sc == ".") #47  something's going on here
stud_adjust_responders = stud_adjust %>% filter(!(stud_adjust$mtb_l_on_days_sc == "."))
stud_adjust_responders$mtb_l_on_days_sc = as.numeric(stud_adjust_responders$mtb_l_on_days_sc)

library(tidyr)
stud_adjust_responders = separate(stud_adjust_responders, into = c("a1","a2","a3") ,col = study_id, sep= "-", remove = T) %>% mutate(good_id = paste(a1,a2,a3, sep = "")) %>% select(-c(a1,a2,a3))
pkData = pkData %>% mutate(good_id = as.character(pkData$OID))
#set a threshold for high/low response
fast_folks = stud_adjust_responders %>% filter(mtb_l_on_days_sc <= 70) %>% summarise(n = n()) # 554
stud_adjust_responders= stud_adjust_responders %>% mutate(CAT_response = ifelse(mtb_l_on_days_sc <= 70,"fast","slow")) 
stud_adjust_responders$CAT_response = as.factor(stud_adjust_responders$CAT_response)
#where do decide response threshold
#mtb_l_on_days_sc :: days from first study dose to stable conversion in liquid medium
summary(dataStud$mtb_l_on_days_sc)
par(mfrow=c(2,1))
hist(as.numeric(dataStud$mtb_l_on_days_sc, breaks = 20))
hist(as.numeric(dataStud$mtb_l_on_days_sc, breaks = 10))
#hist(as.numeric(dataStud$mtb_s_on_days_sc))
barchart(dataStud$mtb_s_on_days_sc)
library(ggplot2)
ggplot() + geom_histogram(data=dataStud, mapping=aes(x = mtb_l_on_days_sc, binwidth = .5)) 


# next explore association between auc/outcome and cmax/outcome: 
pk_matches = stud_adjust_responders %>% filter(good_id %in% pkData$good_id)
pk_and_study_join = inner_join(pkData,pk_matches)
par(mfrow=c(1,1))
plot(pk_and_study_join$AUC, pk_and_study_join$mtb_l_on_days_sc) #auc_vs_response
plot(pk_and_study_join$CMAX, pk_and_study_join$mtb_l_on_days_sc) #cmax_vs_response
#ggplot() + geom_violin()
auc_vs_CATresponse = ggplot() + geom_boxplot(data = pk_and_study_join, mapping = aes(x =CAT_response, y = AUC, colour = CAT_response))
cmax_vs_CATresponse = ggplot() + geom_boxplot(data = pk_and_study_join, mapping = aes(x =CAT_response, y = CMAX, colour = CAT_response))
#higher auc and cmax are associated with faster response

# stat learning:#######################################
# possible vars to drop: Obs,age, study_id,id,birthregion, cxr_bilat, cd4cells, cd4cells_, all the mtb_COLUMNS that aren't our outcome?
# unknown vars: ast, rptdose?, weight_ipt,weight_pwk

stud_adjust_responders %>% filter(weight_ipt == ".") %>% summarise(n = n()) #81 people
stud_adjust_responders %>% filter(weight_pwk == ".") %>% summarise(n = n()) #81 people

colnames(stud_adjust_responders)[grep("mtb_",colnames(stud_adjust_responders))] 
pk_and_study_join_filtered = pk_and_study_join %>% select(-c(Obs,age,id,ID,OID,birthregion, cxr_bilat, cd4cells, cd4cells_,
                                                                       mtb_s_on_days_te,mtb_s_all_days_te,mtb_s_on_days_sc,
                                                                       mtb_l_on_days_sc,rptdose,karnofsky, platelet_,weight_ipt,weight_pwk,
                                                             Gene1,Gene2,Gene3,Biolmarker.1,Biomarker2,Biomarker3,cav_enroll,cav_bi))#

cont = pk_and_study_join_filtered %>% select(age_,ast,ast_uln,hgb,platelet,wbc,height,weight_enroll,AUC, CMAX)


quals = pk_and_study_join_filtered %>% select(-c(age_,ast,ast_uln,hgb,platelet,wbc,height,weight_enroll,good_id,AUC, CMAX))
quals <- as.data.frame(sapply(quals, as.factor))
stud_adjust_clean = cbind(quals,cont)
stud_adjust_clean$good_id = pk_and_study_join_filtered$good_id
stud_adjust_clean$CAT_response = pk_and_study_join_filtered$CAT_response
# there are some weird missing values in the first row
stud_adjust_clean$homeless[1] = 0 # table(stud_adjust_clean$homeless); way more 0's
stud_adjust_clean$unemployed[1] = 0# table(stud_adjust_clean$unemployed)
stud_adjust_clean$drug_inj[1] = 0
stud_adjust_clean$drug_noninj[1] = 0
stud_adjust_clean$drug_use[1] = 0
stud_adjust_clean$alcohol_use[1] = 0
########### enter traain#########
library(caret)
set.seed(42)
#split the training set into a train/test set (which I'm calling validate), so that testSet is a true hold out set
inTrain <- createDataPartition(stud_adjust_clean$CAT_response, p = 0.8, list = F)
# weight_ipt, and weight_pwk have lots of NAs, removing them for the moment, we could impute later
train_set <- stud_adjust_clean[inTrain,]  
training = train_set %>% select(-c(good_id,CAT_response)); prop.table(table(train_set$CAT_response)) # fast = .75% of samples
test_set <- stud_adjust_clean[-inTrain,]; prop.table(table(test_set$CAT_response)) # fast = .75% of samples
rownames(test_set) = NULL
train_response = train_set$CAT_response 

myControl <- trainControl(method = "repeatedcv", number = 10, repeats = 5 )
binary_rf_1 = train(training,train_response,method = "rf", metric = "Kappa",tuneLength =5,ntree = 500, 
                    proximity = T, importance = T, type = 1, trControl = myControl)
# in sample performance sucks: Kappa = 0.08
save(binary_rf_1, file = "binary_rf_1.RData")

########### do some rfe ##############
predictors_df <- training # Make sure you're only doing this on your training data!!!!
labels <- train_response
set.seed(1)
rfe(
  predictors_df[complete.cases(predictors_df),],
  labels[complete.cases(predictors_df)],
  sizes = seq(2, 10, by = 2),
  rfeControl = rfeControl(functions = rfFuncs, method = 'cv', number = 5) #rfFuncs means use a randomForest
)

#The top 5 variables (out of 53): site, CMAX, cxrclass, AUC, smear_who
#############
rfe_train_vars = training %>% select(c(site, CMAX, cxrclass, AUC, smear_who))
rfe_test = test_set %>% select(c(site, CMAX, cxrclass, AUC, smear_who,CAT_response))
rfe_rf_5vars = train(rfe_train_vars,train_response,method = "rf", metric = "Kappa",tuneLength =4,ntree = 500, 
                    proximity = T, importance = T, type = 1, trControl = myControl)
#This is marginally better than all vars:
#mtry 5     Acc: 0.7125123  Kapp: 0.1735163  

# does maximizing ROC help?
myControl_roc <- trainControl(method = "repeatedcv", number = 10, repeats = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
rfe_rf_5vars_roc = train(rfe_train_vars,train_response,method = "rf", metric = "ROC",tuneLength =4,ntree = 500, 
                     proximity = T, importance = T, type = 1, trControl = myControl_roc)
###############
#test out of sample
testy = rfe_test %>% select(-c(CAT_response))
any(is.na(testy),arr.ind = T); sum(is.na(testy),arr.ind = T); which(is.na(testy),arr.ind = T);
#row col
#441  91  15
#make predictions
rfe_rf_5vars_predictions = predict(rfe_rf_5vars,testy)
rfe_rf_5vars_confus = confusionMatrix(rfe_rf_5vars_predictions,test_set$CAT_response) 
#Acc = .68, Kapp = .07 # Got most fast people correct, got most slow people wrong

#try again this time with lots of variables, and removeing cav_enroll,cav_bi from the first join
set.seed(66)
rfe_2 = rfe(
  predictors_df[complete.cases(predictors_df),],
  labels[complete.cases(predictors_df)],
  sizes = seq(2, length(names(training)), by = 4),
  rfeControl = rfeControl(functions = rfFuncs, method = 'cv', number = 5) #rfFuncs means use a randomForest
)
#The top 5 variables (out of 40): site, AUC

########## feature selection / Next steps ################
# permute random seeds and take the top5 vars to see if that changes anything
library(doMC)
registerDoMC(cores=3)

model_list = list()
for (i in 1:5){
  set.seed(i)
  nam = paste('rf',i,sep="_")
  assign(nam,train(training,train_response, method = "rf", metric = "Kappa", tuneLength = 3, ntree = 1500, importance = T))
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
my_avgTop_vars = all_vars %>% select(Vars = myVars_1_1,contains("fast")) %>% mutate(avg_imp = rowSums(.[, -1])/length(model_list)) %>% arrange(desc(avg_imp)) %>% select(Vars,avg_imp) %>% top_n(n = n,wt=avg_imp)
# "site"         "cxrclass"     "CMAX"         "AUC"          "birth_region"

top_5 = train_set %>% select(c(site,cxrclass,CMAX,AUC,birth_region))
myControl_roc <- trainControl(method = "repeatedcv", number = 10, repeats = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
top_5_rf = train(top_5,train_response,method = "rf", metric = "ROC",tuneLength =4,ntree = 1500, 
                         proximity = T, importance = T, type = 1, trControl = myControl_roc)
#mtry 5     ROC: 0.6231633  Sens: 0.8466667  Spec: 0.2514286

png("imp_vars_plot.png")
plot(varImp(top_5_rf)); ggsave(imp_vars_plot,file = "imp_vars_plot.png")
dev.off()
# becareful about how you interpret this. Ask Beau if you have questions.

auc_vs_CATresponse = ggplot() + geom_boxplot(data = pk_and_study_join, mapping = aes(x =CAT_response, y = AUC, colour = CAT_response)); ggsave(auc_vs_CATresponse,file = "auc_vs_CATresponse.png")
cmax_vs_CATresponse = ggplot() + geom_boxplot(data = pk_and_study_join, mapping = aes(x =CAT_response, y = CMAX, colour = CAT_response)); ggsave(cmax_vs_CATresponse,file = "cmax_vs_CATresponse.png")
site_vs_CATresponse = ggplot() + geom_boxplot(data = pk_and_study_join, mapping = aes(x =CAT_response, y = site, colour = CAT_response)); ggsave(site_vs_CATresponse,file = "site_vs_CATresponse.png")
cxrclass_vs_CATresponse = ggplot() + geom_boxplot(data = pk_and_study_join, mapping = aes(x =CAT_response, y = cxrclass, colour = CAT_response)); ggsave(cxrclass_vs_CATresponse,file = "cxrclass_vs_CATresponse.png")
birth_region_vs_CATresponse = ggplot() + geom_boxplot(data = pk_and_study_join, mapping = aes(x =CAT_response, y = birth_region, colour = CAT_response)); ggsave(birth_region_vs_CATresponse,file = "birth_region_vs_CATresponse.png")

site_vs_CATresponse_hist = ggplot() + geom_histogram(data = pk_and_study_join, mapping = aes(x =site, fill = CAT_response)); ggsave(site_vs_CATresponse_hist,file = "site_vs_CATresponse_hist.png")