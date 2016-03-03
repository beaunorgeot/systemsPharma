
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
#set a threshold for high/low response
fast_folks = stud_adjust_responders %>% filter(mtb_l_on_days_sc <= 70) %>% summarise(n = n()) # 554
stud_adjust_responders= stud_adjust_responders %>% mutate(CAT_response = ifelse(mtb_l_on_days_sc <= 70,"fast","slow")) 
stud_adjust_responders$CAT_response = as.factor(stud_adjust_responders$CAT_response)
pk_matches = stud_adjust_responders %>% filter(good_id %in% pkData$good_id)
pk_and_study_join = inner_join(pkData,pk_matches)

library(tidyr)
stud_adjust_responders = separate(stud_adjust_responders, into = c("a1","a2","a3") ,col = study_id, sep= "-", remove = T) %>% mutate(good_id = paste(a1,a2,a3, sep = ""))
pkData = pkData %>% mutate(good_id = as.character(pkData$OID))
# next calc r2 between auc/outcome and cmax/outcome: cor()

#mtb_l_on_days_sc :: days from first study dose to stable conversion in liquid medium
summary(dataStud$mtb_l_on_days_sc)
par(mfrow=c(2,1))
hist(as.numeric(dataStud$mtb_l_on_days_sc, breaks = 20))
hist(as.numeric(dataStud$mtb_l_on_days_sc, breaks = 10))
#hist(as.numeric(dataStud$mtb_s_on_days_sc))
barchart(dataStud$mtb_s_on_days_sc)
library(ggplot2)
ggplot() + geom_histogram(data=dataStud, mapping=aes(x = mtb_l_on_days_sc, binwidth = .5)) 

dataStud$mtb_l_on_days_sc = as.numeric(dataStud$mtb_l_on_days_sc)
#where are the hard responders
fast_folks = dataStud %>% filter(mtb_l_on_days_sc <= 70) %>% summarise(n = n())
dataStud$mtb_l_on_days_sc = as.numeric(as.character(dataStud$mtb_l_on_days_sc))
