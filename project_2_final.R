setwd("/Users/beaunorgeot/systemsPharma")

library(caret)
library(dplyr)
library(tidyr)
library(doMC)
registerDoMC(cores = 3)

final_df = read.csv("df_beau_final2_lessfactors.csv",na.strings=c("",".","NA"))
final_df = final_df %>% mutate(education = ifelse(is.na(education),"1",education))                                                                                                        
final_df = final_df %>% mutate(study = ifelse(is.na(study),"29x",study))
final_df = final_df %>% filter(!is.na(mtb_l_on_days_sc))

final_df$mtb_l_on_days_sc = as.numeric(as.character(final_df$mtb_l_on_days_sc))
final_df = final_df %>% dplyr::mutate(CAT_response = ifelse(mtb_l_on_days_sc <= 70,"fast","slow")) 
final_df$CAT_response = as.factor(final_df$CAT_response)

final_df = final_df %>% mutate(fever_cl = ifelse(is.na(fever),1,fever), bmi_c = ifelse(is.na(bmi_),2,bmi_), education_c = ifelse(is.na(education),1,education))
final_df = final_df %>% select(-c(education,bmi_,fever,FOODCAT,mtb_l_on_days_sc,rifrpt,smear_who,ethn_hisp,weight_ipt,weight_pwk,homeless,unemployed,drug_inj,drug_noninj,drug_use,alcohol_use,cigarette_use))


inTrain <- createDataPartition(final_df$CAT_response, p = 0.8, list = F)
trainSet <- final_df[inTrain,]
train_response = trainSet$CAT_response
training = trainSet %>% select(-CAT_response)
testSet <- final_df[-inTrain,]

myControl_f = trainControl(method = "repeatedcv", number = 10, repeats = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
set.seed(17)
final_all_rf = train(training, train_response, method= "rf", metric = "ROC",tuneLength =4,ntree = 1500, 
                     importance = T, trControl = myControl_f)
# mtry = 2 ROC: 0.6528994  Sens:0.9570575  Spec:0.1085455
### out of sample predictions
preds = predict(final_all_rf, testSet)
confused = confusionMatrix(preds,testSet$CAT_response)
#Sens = .97, Spec = .12, Kapp = .12

varImp(final_all_rf)
# top 5 vars
#real_predicted_auc, AUC_using_penetration,cxrclass,antitb,education_c

#plot them up
auc_vs_CATresponse = ggplot() + geom_boxplot(data = final_df, mapping = aes(x =CAT_response, y = real_predicted_auc, colour = CAT_response)); ggsave(auc_vs_CATresponse,file = "auc_vs_CATresponse2.png")
aucPENETRATION_vs_CATresponse = ggplot() + geom_boxplot(data = final_df, mapping = aes(x =CAT_response, y = AUC_using_penetration, colour = CAT_response)); ggsave(aucPENETRATION_vs_CATresponse,file = "aucPENETRATION_vs_CATresponse.png")
cxrclass_vs_CATresponse = ggplot() + geom_boxplot(data = final_df, mapping = aes(x =CAT_response, y = cxrclass, colour = CAT_response)); ggsave(cxrclass_vs_CATresponse,file = "cxrclass_vs_CATresponse.png")
education_vs_CATresponse_hist = ggplot() + geom_histogram(data = final_df, mapping = aes(x =education_c, fill = CAT_response)); ggsave(education_vs_CATresponse_hist,file = "education_vs_CATresponse_hist.png")
antitb_vs_CATresponse = ggplot() + geom_boxplot(data = final_df, mapping = aes(x =CAT_response, y = antitb, colour = CAT_response)); ggsave(antitb_vs_CATresponse,file = "antitb_vs_CATresponse.png")
# I don't know how to interpret this antib plot
table(final_df$antitb) # there were 357 0's and 136 1's 

someshit = final_df %>% 
  group_by(CAT_response,antitb) %>%
  summarise(n = n())

# all Basically, nearly all of the slow responders had antitb as 0