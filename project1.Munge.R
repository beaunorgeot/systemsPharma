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
allFeatures = cbind(cnv,geneX)
save(allFeatures, file = "allFeatures.RData")

sum(is.na(allFeatures)) #68530/1610549 ~ 4% of the data is missing

#prep for heat map
only_drugs = drugs %>% select(-cLine)
corrs = cor(only_drugs, use = "pairwise.complete.obs")
sigs = which(corrs > .5 & corrs != 1,arr.ind = T)

# Beau: start a RF preprocesing w/KNN. How to choose value of K?
