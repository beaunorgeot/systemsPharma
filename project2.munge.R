
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
stud_adjust_responders_filtered = stud_adjust_responders %>% select(-c(Obs,age,id,birthregion, cxr_bilat, cd4cells, cd4cells_,
                                                                       mtb_s_on_days_te,mtb_s_all_days_te,mtb_s_on_days_sc,
                                                                       mtb_l_on_days_sc,rptdose,karnofsky, platelet_))

cont = stud_adjust_responders_filtered %>% select(age_,ast,ast_uln,hgb,platelet,wbc,height,weight_enroll,weight_ipt,weight_pwk,
                                                  Biolmarker.1,Biomarker2,Biomarker3)

cont$weight_ipt = as.numeric(as.character(cont$weight_ipt))
cont$weight_pwk = as.numeric(as.character(cont$weight_pwk))

quals = stud_adjust_responders_filtered %>% select(-c(age_,ast,ast_uln,hgb,platelet,wbc,height,weight_enroll,weight_ipt,weight_pwk,
                                                     Biolmarker.1,Biomarker2,Biomarker3,good_id))
quals <- as.data.frame(sapply(quals, as.factor))
stud_adjust_clean = cbind(quals,cont)
stud_adjust_clean$good_id = stud_adjust_responders_filtered$good_id
stud_adjust_clean$CAT_response = stud_adjust_responders_filtered$CAT_response
library(caret)
set.seed(42)
#split the training set into a train/test set (which I'm calling validate), so that testSet is a true hold out set
inTrain <- createDataPartition(stud_adjust_clean$CAT_response, p = 0.8, list = F)
train_set <- stud_adjust_clean[inTrain,]
training = train_set %>% select(-c(good_id,CAT_response))
test_set <- stud_adjust_clean[-inTrain,]
train_response = train_set$CAT_response

myControl <- trainControl(method = "repeatedcv", number = 10, repeats = 5 )
binary_rf_1 = train(training,train_response,method = "rf", metric = "Kappa",tuneLength =5,ntree = 500, 
                    proximity = T, importance = T, type = 1, trControl = myControl)

testy = test_set %>% select(-c(good_id,CAT_response))
#make predictions
binary_rf_1_predictions = predict(binary_rf_1,testy)
binary_rf_1_confus = confusionMatrix(binary_rf_1_predictions,test_set$CAT_response) 