######################################################################################
## Project: Validation of Maternal Vaginal Microbiota Signatures in Pregnancy: The ECHO Vaginal Microbiome Consortium - Dr Kimberly McKee
## Script name: Random Forests.R
## Script purpose: Creating random forest plots using complete-case dataset
## Date: March 2023
## Author: Beatrice Palazzolo
## Organization: Department of Family Medicine, University of Michigan
######################################################################################

# Libraries ---------------------------------------------------------------

pkgs <- c(
  "readxl", "tidyverse", "ggplot2", "knitr", "readr", "haven", "gtsummary", 
  "kable", "kableExtra", "writexl", "webshot2", "randomForest", "caret", "ROCR", "ggpubr", "stringr"
)

pkg_check <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "https://cran.rstudio.com")
    library(x, character.only = TRUE)
  }
}
lapply(pkgs, pkg_check)

# Read main datasets ------------------------------------------------------

# Set directory

all_data_cst <- readRDS("all_data_cst.rds")

all_data <- all_data_cst %>%
  dplyr::select(3:5229,5237,5241,5249,5255,5447)

all_data <- all_data %>%
  mutate(
    white=as.factor(white),
    AbxEver=as.factor(AbxEver),
    educ=as.factor(educ),
    Parity_cat=ifelse(Parity_cat=="No prior", 0,
                      ifelse(Parity_cat=="1", "1+",
                             ifelse(Parity_cat>1, "1+", NA))),
    Parity_cat=as.factor(Parity_cat)
  ) 

# Random forest -----------------------------------------------------------

# All cohorts ####

# White, Maternal age, Education, Antibiotics Ever in Pregnancy, Parity Category

# white ####

# Set directory

white_data <- all_data %>%
  dplyr::select(1:5227,5232)
white_data <- white_data[,colSums(is.na(white_data))<nrow(white_data)]
white_data <- white_data[,rowSums(is.na(white_data))<ncol(white_data)]
white_data<-white_data[complete.cases(white_data), ]
# n=675

#Set random seed before splitting for repeatability 
set.seed(12345) 

#Set a portion of data aside for testing (80/20)
samp <- sample(nrow(white_data), 0.8 * nrow(white_data))
train <- white_data[samp, ]
test <- white_data[-samp, ]
#Check dimensions
dim(train)
dim(test)
# n=540 in train
# n=135 in test

#Run randomForest & tune
library(randomForest)
set.seed(12345) # need to reset seed for reproducibility 
# Start with model with default settings (ntrees=500,mtry=72)
model1 <- randomForest(white ~ ., data = train, proximity=TRUE)
model1
# OOB error 5% so accuracy 95%.
set.seed(12345) # need to reset seed for reproducibility 
model2 <- randomForest(white ~ ., data = train, ntree=800,proximity=TRUE)
model2
# OOB error 4.63% so accuracy around 95.4%.
set.seed(12345) # need to reset seed for reproducibility 
model3 <- randomForest(white ~ ., data = train, ntree=1500,proximity=TRUE)
model3
# OOB error 4.44% so accuracy around 95.6%.
set.seed(12345) # need to reset seed for reproducibility 
model4 <- randomForest(white ~ ., data = train, ntree=1000,proximity=TRUE)
model4
# OOB error 4.44% so accuracy around 95.6%. 1000 best.
set.seed(12345) # need to reset seed for reproducibility 
model5 <- randomForest(white ~ ., data = train, ntree=750,proximity=TRUE)
model5
# OOB error 4.63% so accuracy around 95.4%.
set.seed(12345) # need to reset seed for reproducibility 
features <- setdiff(names(train),"white")
model6 <- tuneRF(
  x = train[features],
  y = train$white,
  ntreeTry=1000,
  mtryStart=5,
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE
)
# 73 mtry for lowest OOB error--close to default value of 72
set.seed(12345) # need to reset seed for reproducibility 
model4 <- randomForest(white ~ ., data = train, ntree = 1000, mtry=73, proximity=TRUE)
model4
model=model4
model

#Confusion matrix
model$confusion 

#Validate model using test data
prediction <- predict(model, newdata = test)
table(prediction, test$white)
prediction

#Display the predicted vs. the actual values
results<-cbind(prediction,test$white)
results
colnames(results)<-c('pred','real')
results<-as.data.frame(results)
View(results)

#Calculate model accuracy
sum(prediction==test$white) / nrow(test) # 95.6% accuracy

#Plot random forest
png(file = "Rf plot.png", res=90)
plot(model)
dev.off()

#Number of nodes for trees
png(file = "Number of nodes.png", res=90)
hist(treesize(model),
     main = "No. of Nodes for the Trees ('White')")
#col = "red")
dev.off()

#Predictor importance
imp <- importance(model)
imp <- data.frame(predictors = rownames(imp), imp)
# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
# Select the top predictors
imp.20 <- imp.sort[1:20, ]

png(file = "varImp.png", res=90, width = 800, height = 480)
varImp <- varImpPlot(model,
                     sort = T,
                     n.var = 20,
                     main = "Top 20 - Variable Importance ('white')")
dev.off()

# Higher the value of mean decrease accuracy or mean decrease Gini score, higher 
# the importance of the variable in the model

# ggplot

imp.20$predictors=factor(imp.20$predictors,levels=c("pt__00013",
                                                    "pt__00406",
                                                    "pt__00008",
                                                    "pt__00652",
                                                    "pt__00127",
                                                    "pt__00006",
                                                    "pt__00014",
                                                    "pt__00018",
                                                    "pt__00004",
                                                    "pt__00009",
                                                    "pt__00112",
                                                    "pt__00002",
                                                    "pt__00199",
                                                    "pt__00001",
                                                    "pt__00271",
                                                    "pt__00041",
                                                    "pt__00539",
                                                    "pt__00005",
                                                    "pt__00042",
                                                    "pt__00122"),
                         labels=c("Streptococcus oralis",
                                  "Chlamydia trachomatis / Elusimicrobium minutum",
                                  "Lactobacillus gallinarum",
                                  "Peptoniphilus rhinitidis / Peptoniphilus harei",
                                  "Anaerococcus vaginalis",
                                  "Sneathia amnii_1",
                                  "Dialister micraerophilus_1",
                                  "Haemophilus parainfluenzae_1",
                                  "Atopobium vaginae_1",
                                  "Veillonella dispar",
                                  "Peptoniphilus harei / Peptoniphilus grossensis / Peptoniphilus vaginalis",
                                  "Lactobacillus iners_1",
                                  "Anaerococcus octavius_1",
                                  "Gardnerella vaginalis",
                                  "Peptoniphilus pacaensis",
                                  "Finegoldia magna",
                                  "Anaerococcus prevotii / Anaerococcus tetradius",
                                  "Megasphaera genomosp. type_1_1",
                                  "Streptococcus salivarius_1",
                                  "Varibaculum cambriense"))

predimp1 <- ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most Important Taxa for 'White'") +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face="bold"),    
        axis.title.y = element_text(face="bold")) +
  xlab("Taxa") + ylab("Mean Decrease Gini")
ggsave("Preds importance.pdf")

# ROC in ggplot
library(pROC)
predicted <- predict(model, test, type="prob")
rocobj <- roc(test$white, predicted[,2])
# auc <- round(auc(test$white, predicted),4)
auc <- round(auc(rocobj),4)
roc_gg1<- ggroc(rocobj, colour = 'steelblue', size = 1) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) +
  theme(plot.title = element_text(face="bold")) +
  xlab("Specificity") + ylab("Sensitivity") +
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="red", linetype="dashed")
ggsave("ROC.pdf")

#Multi-dimensional scaling plot of proximity matrix
png(file = "MDSplot.png", res=90)
MDSplot(model, train$white)
dev.off()

# SHAP analysis ####
install.packages("xgboost")
install.packages("SHAPforxgboost")
install.packages("caret")
library(xgboost)
library(SHAPforxgboost)
library(caret)

mtrain <- train %>%
  mutate(
    white=as.numeric(white)
  ) 
colnames(mtrain)[colnames(mtrain) %in% c("pt__00013", "pt__00652", "pt__00020","pt__00406","pt__00153")] <- c("Streptococcus oralis", "Peptoniphilus rhinitidis / Peptoniphilus harei","Megasphaera genomosp. type_1_2","[Bacteroides] coagulans_1","Parvibacter caecicola")
labels <- mtrain$white
labels <- as.numeric(labels)-1
mtrain=as.matrix(mtrain[,-5228])
bstDMatrix <- xgboost(data = mtrain, 
                      max.depth = 2, 
                      eta = 1, 
                      nthread = 2, 
                      nrounds = 2, 
                      label=labels,
                      objective = "binary:logistic")
shap_values <- shap.values(xgb_model = bstDMatrix, X_train = mtrain)
shap_values$mean_shap_score
shap_values_white <- shap_values$shap_score
shap_long_white <- shap.prep(xgb_model = bstDMatrix, X_train = mtrain)
# **SHAP summary plot**
# Shows top 20 most important predictors
shap.plot.summary.wrap1(bstDMatrix, X = as.matrix(mtrain[,-5228]), top_n = 20)
# The y-axis indicates the variable name, in order of importance from top to bottom. 
# The value next to them is the mean SHAP value.
# On the x-axis is the SHAP value. Indicates how much is the change in log-odds. 
# Low Streptococcus oralis associated with higher log-odds of being white on average (though some data points with lower).
# Low Parvibacter caecicola more associated with lower log-odds of being white (though mean SHAP is close to zero and sparse data points associated with higher log-odds of being white). 
# Low Peptoniphilus rhinitidis / Peptoniphilus harei associated with lower log-odds of being white (though some data points have odds just to the right of zero). 
# Some evidence of higher [Bacteroides] coagulans 1 associated higher log-odds of being white. Mostly though, lower [Bacteroides] coagulans 1 associated with lower log-odds of being white.    
# Low Megasphaera genomosp. type 1 2 mostly associated with higher (but small) log-odds of being white. 
png(file = "SHAP plot.png", res=100, width=800)
shap.plot.summary.wrap1(bstDMatrix, X = as.matrix(mtrain[,-5228]), top_n = 5)
dev.off()

# educ ####
# Set directory

ed_data <- all_data %>%
  dplyr::select(1:5228)
ed_data <- ed_data[,colSums(is.na(ed_data))<nrow(ed_data)]
ed_data <- ed_data[,rowSums(is.na(ed_data))<ncol(ed_data)]
ed_data<-ed_data[complete.cases(ed_data), ]
# n=678

#Set random seed before splitting for repeatability 
set.seed(12345) 

#Set a portion of data aside for testing (80/20)
samp <- sample(nrow(ed_data), 0.8 * nrow(ed_data))
train <- ed_data[samp, ]
test <- ed_data[-samp, ]
#Check dimensions
dim(train)
dim(test)
# n=542 in train
# n=136 in test

#Run randomForest & tune
library(randomForest)
set.seed(12345) # need to reset seed for reproducibility 
# Start with model with default settings (ntrees=500,mtry=72)
model1 <- randomForest(educ ~ ., data = train, proximity=TRUE)
model1
# OOB error 48.34% so accuracy 51.7%.
set.seed(12345) # need to reset seed for reproducibility 
model2 <- randomForest(educ ~ ., data = train, ntree=800,proximity=TRUE)
model2
# OOB error 47.6% so accuracy around 52.4%.
set.seed(12345) # need to reset seed for reproducibility 
model3 <- randomForest(educ ~ ., data = train, ntree=1500,proximity=TRUE)
model3
# OOB error 46.86% so accuracy around 53.14%. 1500 best.
set.seed(12345) # need to reset seed for reproducibility 
model4 <- randomForest(educ ~ ., data = train, ntree=1000,proximity=TRUE)
model4
# OOB error 47.97% so accuracy around 52%. 
set.seed(12345) # need to reset seed for reproducibility 
model5 <- randomForest(educ ~ ., data = train, ntree=2000,proximity=TRUE)
model5
# OOB error 47.23% so accuracy around 52.8%. 2000 is best. 
set.seed(12345) # need to reset seed for reproducibility 
model6 <- randomForest(educ ~ ., data = train, ntree=3000,proximity=TRUE)
model6
# OOB error 47.23% so accuracy around 52.8%.
set.seed(12345) # need to reset seed for reproducibility 
features <- setdiff(names(train),"educ")
model7 <- tuneRF(
  x = train[features],
  y = train$educ,
  ntreeTry=2000,
  mtryStart=5,
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE
)
# 15 mtry for lowest OOB error--not close to default value of 72, but 72 performs better so keep it
set.seed(12345) # need to reset seed for reproducibility 
model8 <- randomForest(educ ~ ., data = train, ntree = 2000, mtry=72, proximity=TRUE)
model8
model=model8
model

#Confusion matrix
model$confusion 

#Validate model using test data
prediction <- predict(model, newdata = test)
table(prediction, test$educ)
prediction

#Display the predicted vs. the actual values
results<-cbind(prediction,test$educ)
results
colnames(results)<-c('pred','real')
results<-as.data.frame(results)
View(results)

#Calculate model accuracy
sum(prediction==test$educ) / nrow(test) # 45.6% accuracy

#Plot random forest
png(file = "Rf plot.png", res=90)
plot(model)
dev.off()

#Number of nodes for trees
png(file = "Number of nodes.png", res=90)
hist(treesize(model),
     main = "No. of Nodes for the Trees ('Education Level')")
#col = "red")
dev.off()

#Predictor importance
imp <- importance(model)
imp <- data.frame(predictors = rownames(imp), imp)
# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
# Select the top predictors
imp.20 <- imp.sort[1:20, ]

png(file = "varImp.png", res=90, width = 800, height = 480)
varImp <- varImpPlot(model,
                     sort = T,
                     n.var = 20,
                     main = "Top 20 - Variable Importance ('Education Level')")
dev.off()

# Higher the value of mean decrease accuracy or mean decrease Gini score, higher 
# the importance of the variable in the model

# ggplot
imp.20$predictors=factor(imp.20$predictors,levels=c("pt__00015",
                                                    "pt__00002",
                                                    "pt__00001",
                                                    "pt__00013",
                                                    "pt__00014",
                                                    "pt__00004",
                                                    "pt__00011",
                                                    "pt__00006",
                                                    "pt__00005",
                                                    "pt__00034",
                                                    "pt__00018",
                                                    "pt__00019",
                                                    "pt__00025",
                                                    "pt__00020",
                                                    "pt__00008",
                                                    "pt__00016",
                                                    "pt__00009",
                                                    "pt__00041",
                                                    "pt__00035",
                                                    "pt__00010"),
                         labels=c("Aerococcus christensenii_1",
                                  "Lactobacillus iners_1",
                                  "Gardnerella vaginalis_1",
                                  "Streptococcus oralis",
                                  "Dialister micraerophilus_1",
                                  "Atopobium vaginae_1",
                                  "Lactobacillus crispatus",
                                  "Sneathia amnii_1",
                                  "Megasphaera genomosp. type_1_1",
                                  "Parvibacter caecicola / Adlercreutzia equolifaciens",
                                  "Haemophilus parainfluenzae_1",
                                  "Ureaplasma parvum",
                                  "Prevotella disiens_1",
                                  "Megasphaera genomosp. type_1_2",
                                  "Lactobacillus gallinarum",
                                  "Lactobacillus gasseri",
                                  "Veillonella dispar",
                                  "Finegoldia magna",
                                  "Staphylococcus epidermidis_1",
                                  "Lactobacillus jensenii_1"))
predimp2 <- ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most Important Taxa for 'Education Level'") +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face="bold"),    
        axis.title.y = element_text(face="bold")) +
  xlab("Taxa") + ylab("Mean Decrease Gini")
ggsave("Preds importance.pdf")

# ROC in ggplot
library(pROC)
predicted <- predict(model, test, type="prob")
rocobj <- roc(test$educ, predicted[,2])
# auc <- round(auc(test$white, predicted),4)
auc <- round(auc(rocobj),4)
roc_gg2<- ggroc(rocobj, colour = 'steelblue', size = 1) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) +
  theme(plot.title = element_text(face="bold")) +
  xlab("Specificity") + ylab("Sensitivity") +
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="red", linetype="dashed")
ggsave("ROC.pdf")

#Multi-dimensional scaling plot of proximity matrix
png(file = "MDSplot.png", res=90)
MDSplot(model, train$educ)
dev.off()

# AbxEver ####
# Set directory

abx_data <- all_data %>%
  dplyr::select(1:5227,5229)
abx_data <- abx_data[,colSums(is.na(abx_data))<nrow(abx_data)]
abx_data <- abx_data[,rowSums(is.na(abx_data))<ncol(abx_data)]
abx_data<-abx_data[complete.cases(abx_data), ]
# n=515

all_data_cst %>%
  dplyr::select(cohort,AbxEver) %>%
  filter(cohort=="Emory" & is.na(AbxEver)) %>%
  count()
all_data_cst %>%
  dplyr::select(cohort,AbxEver) %>%
  filter(cohort=="MARCH" & is.na(AbxEver)) %>%
  count()
# Only 1 missing data from Emory (all missing from Wisconsin)

#Set random seed before splitting for repeatability 
set.seed(12345) 

#Set a portion of data aside for testing (80/20)
samp <- sample(nrow(abx_data), 0.8 * nrow(abx_data))
train <- abx_data[samp, ]
test <- abx_data[-samp, ]
#Check dimensions
dim(train)
dim(test)
# n=412 in train
# n=103 in test

#Run randomForest & tune
library(randomForest)
set.seed(12345) # need to reset seed for reproducibility 
# Start with model with default settings (ntrees=500,mtry=72)
model1 <- randomForest(AbxEver ~ ., data = train, proximity=TRUE)
model1
# OOB error 38.83% so accuracy 61.2%.
set.seed(12345) # need to reset seed for reproducibility 
model2 <- randomForest(AbxEver ~ ., data = train, ntree=800,proximity=TRUE)
model2
# OOB error 38.83% so accuracy around 61.2%.
set.seed(12345) # need to reset seed for reproducibility 
model3 <- randomForest(AbxEver ~ ., data = train, ntree=1500,proximity=TRUE)
model3
# OOB error 37.62% so accuracy around 62.3%. 1500 best.
set.seed(12345) # need to reset seed for reproducibility 
model4 <- randomForest(AbxEver ~ ., data = train, ntree=2000,proximity=TRUE)
model4
# OOB error 37.86% so accuracy around 62.14%. 
set.seed(12345) # need to reset seed for reproducibility 
features <- setdiff(names(train),"AbxEver")
model5 <- tuneRF(
  x = train[features],
  y = train$AbxEver,
  ntreeTry=1500,
  mtryStart=5,
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE
)
# 15 mtry for lowest OOB error--not close to default value of 72, but 72 performs better so keep it
set.seed(12345) # need to reset seed for reproducibility 
model6 <- randomForest(AbxEver ~ ., data = train, ntree = 1500, mtry=4, proximity=TRUE)
model6
model=model6
model
set.seed(12345) # need to reset seed for reproducibility
model7 <- randomForest(AbxEver ~ ., data = train, ntree = 1500, mtry=72, proximity=TRUE)
model7
model=model7
model

#Confusion matrix
model$confusion 

#Validate model using test data
prediction <- predict(model, newdata = test)
table(prediction, test$AbxEver)
prediction

#Display the predicted vs. the actual values
results<-cbind(prediction,test$AbxEver)
results
colnames(results)<-c('pred','real')
results<-as.data.frame(results)
View(results)

#Calculate model accuracy
sum(prediction==test$AbxEver) / nrow(test) # 68% accuracy

#Plot random forest
png(file = "Rf plot.png", res=90)
plot(model)
dev.off()

#Number of nodes for trees
png(file = "Number of nodes.png", res=90)
hist(treesize(model),
     main = "No. of Nodes for the Trees ('Antibiotics Ever in Pregnancy')")
#col = "red")
dev.off()

#Predictor importance
imp <- importance(model)
imp <- data.frame(predictors = rownames(imp), imp)
# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
# Select the top predictors
imp.20 <- imp.sort[1:20, ]

png(file = "varImp.png", res=90, width = 800, height = 480)
varImp <- varImpPlot(model,
                     sort = T,
                     n.var = 20,
                     main = "Top 20 - Variable Importance ('Antibiotics Ever in Pregnancy')")
dev.off()

# Higher the value of mean decrease accuracy or mean decrease gini score , higher 
# the importance of the variable in the model

# ggplot
imp.20$predictors=factor(imp.20$predictors,levels=c("pt__00002",
                                                    "pt__00001",
                                                    "pt__00004",
                                                    "pt__00006",
                                                    "pt__00011",
                                                    "pt__00015",
                                                    "pt__00044",
                                                    "pt__00034",
                                                    "pt__00020",
                                                    "pt__00014",
                                                    "pt__00013",
                                                    "pt__00008",
                                                    "pt__00005",
                                                    "pt__00035",
                                                    "pt__00058",
                                                    "pt__00003",
                                                    "pt__00018",
                                                    "pt__00030",
                                                    "pt__00010",
                                                    "pt__00017"),
                         labels=c("Lactobacillus iners_1",
                                  "Gardnerella vaginalis_1",
                                  "Atopobium vaginae_1",
                                  "Sneathia amnii_1",
                                  "Lactobacillus crispatus",
                                  "Aerococcus christensenii_1",
                                  "Mycoplasma hominis_1",
                                  "Parvibacter caecicola / Adlercreutzia equolifaciens",
                                  "Megasphaera genomosp. type_1_2",
                                  "Dialister micraerophilus_1",
                                  "Streptococcus oralis",
                                  "Lactobacillus gallinarum",
                                  "Megasphaera genomosp. type_1_1",
                                  "Staphylococcus epidermidis_1",
                                  "Corynebacterium tuberculostearicum",
                                  "Prevotella buccalis_1",
                                  "Haemophilus parainfluenzae_1",
                                  "Dialister micraerophilus_2",
                                  "Lactobacillus jensenii_1",
                                  "Prevotella bivia_1"))
predimp3 <- ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most Important Taxa for 'Antibiotics Ever in Pregnancy'") +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face="bold"),    
        axis.title.y = element_text(face="bold")) +
  xlab("Taxa") + ylab("Mean Decrease Gini")
ggsave("Preds importance.pdf",width=9)

# ROC in ggplot
library(pROC)
predicted <- predict(model, test, type="prob")
rocobj <- roc(test$AbxEver, predicted[,2])
# auc <- round(auc(test$white, predicted),4)
auc <- round(auc(rocobj),4)
roc_gg3<- ggroc(rocobj, colour = 'steelblue', size = 1) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) +
  theme(plot.title = element_text(face="bold")) +
  xlab("Specificity") + ylab("Sensitivity") +
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="red", linetype="dashed")
ggsave("ROC.pdf")

#Multi-dimensional scaling plot of proximity matrix
png(file = "MDSplot.png", res=90)
MDSplot(model, train$AbxEver)
dev.off()

# Parity_cat ####
# Set directory

p_data <- all_data %>%
  dplyr::select(1:5227,5231)
p_data <- p_data[,colSums(is.na(p_data))<nrow(p_data)]
p_data <- p_data[,rowSums(is.na(p_data))<ncol(p_data)]
p_data<-p_data[complete.cases(p_data), ]
# n=679

#Set random seed before splitting for repeatability 
set.seed(12345) 

#Set a portion of data aside for testing (80/20)
samp <- sample(nrow(p_data), 0.8 * nrow(p_data))
train <- p_data[samp, ]
test <- p_data[-samp, ]
#Check dimensions
dim(train)
dim(test)
# n=543 in train
# n=136 in test

#Run randomForest & tune
library(randomForest)
set.seed(12345) # need to reset seed for reproducibility 
# Start with model with default settings (ntrees=500,mtry=72)
model1 <- randomForest(Parity_cat ~ ., data = train, proximity=TRUE)
model1
# OOB error 42.54% so accuracy 57.46%.
set.seed(12345) # need to reset seed for reproducibility 
model2 <- randomForest(Parity_cat ~ ., data = train, ntree=800,proximity=TRUE)
model2
# OOB error 41.25% so accuracy around 58.75%.
set.seed(12345) # need to reset seed for reproducibility 
model3 <- randomForest(Parity_cat ~ ., data = train, ntree=1500,proximity=TRUE)
model3
# OOB error 41.07% so accuracy around 58.93%. 
set.seed(12345) # need to reset seed for reproducibility 
model4 <- randomForest(Parity_cat ~ ., data = train, ntree=1000,proximity=TRUE)
model4
# OOB error 40.88% so accuracy around 59.12%. 1000 best.
set.seed(12345) # need to reset seed for reproducibility 
features <- setdiff(names(train),"Parity_cat")
model5 <- tuneRF(
  x = train[features],
  y = train$Parity_cat,
  ntreeTry=1000,
  mtryStart=5,
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE
)
# 7 mtry for lowest OOB error
set.seed(12345) # need to reset seed for reproducibility 
model6 <- randomForest(Parity_cat ~ ., data = train, ntree = 1000, mtry=7, proximity=TRUE)
model6
model=model6
model
#Confusion matrix
model$confusion 

#Validate model using test data
prediction <- predict(model, newdata = test)
table(prediction, test$Parity_cat)
prediction

#Display the predicted vs. the actual values
results<-cbind(prediction,test$Parity_cat)
results
colnames(results)<-c('pred','real')
results<-as.data.frame(results)
View(results)

#Calculate model accuracy
sum(prediction==test$Parity_cat) / nrow(test) # 60.3% accuracy

#Plot random forest
png(file = "Rf plot.png", res=90)
plot(model)
dev.off()

#Number of nodes for trees
png(file = "Number of nodes.png", res=90)
hist(treesize(model),
     main = "No. of Nodes for the Trees ('Parity Category')")
#col = "red")
dev.off()

#Predictor importance
imp <- importance(model)
imp <- data.frame(predictors = rownames(imp), imp)
# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
# Select the top predictors
imp.20 <- imp.sort[1:20, ]

png(file = "varImp.png", res=90, width = 800, height = 480)
varImp <- varImpPlot(model,
                     sort = T,
                     n.var = 20,
                     main = "Top 20 - Variable Importance ('Parity Category')")
dev.off()

# Higher the value of mean decrease accuracy or mean decrease Gini score, higher 
# the importance of the variable in the model

# ggplot
imp.20$predictors=factor(imp.20$predictors,levels=c("pt__00011",
                                                    "pt__00001",
                                                    "pt__00009",
                                                    "pt__00070",
                                                    "pt__00010",
                                                    "pt__00008",
                                                    "pt__00007",
                                                    "pt__00031",
                                                    "pt__00075",
                                                    "pt__00016",
                                                    "pt__00150",
                                                    "pt__00036",
                                                    "pt__00100",
                                                    "pt__00092",
                                                    "pt__00018",
                                                    "pt__00003",
                                                    "pt__00015",
                                                    "pt__00042",
                                                    "pt__00326",
                                                    "pt__00012"),
                         labels=c("Lactobacillus crispatus",
                                  "Gardnerella vaginalis",
                                  "Veillonella dispar",
                                  "Lactobacillus coleohominis",
                                  "Lactobacillus jensenii_1",
                                  "Lactobacillus gallinarum",
                                  "Prevotella timonensis_1",
                                  "Porphyromonas uenonis",
                                  "Bifidobacterium longum",
                                  "Lactobacillus gasseri",
                                  "Corynebacterium amycolatum_1",
                                  "Fenollaria timonensis / Fenollaria massiliensis",
                                  "Lactobacillus reuteri_1",
                                  "Mageeibacillus indolicus_1",
                                  "Haemophilus parainfluenzae_1",
                                  "Prevotella buccalis_1",
                                  "Aerococcus christensenii_1",
                                  "Streptococcus salivarius_1",
                                  "Megasphaera micronuciformis_2",
                                  "Prevotella timonensis_2"))
predimp4 <- ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most Important Taxa for 'Parity Category'") +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face="bold"),    
        axis.title.y = element_text(face="bold")) +
  xlab("Taxa") + ylab("Mean Decrease Gini")
ggsave("Preds importance.pdf",width=7)

# ROC in ggplot
library(pROC)
predicted <- predict(model, test, type="prob")
rocobj <- roc(test$Parity_cat, predicted[,2])
# auc <- round(auc(test$white, predicted),4)
auc <- round(auc(rocobj),4)
roc_gg4<- ggroc(rocobj, colour = 'steelblue', size = 1) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) +
  theme(plot.title = element_text(face="bold")) +
  xlab("Specificity") + ylab("Sensitivity") +
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="red", linetype="dashed")
ggsave("ROC.pdf")

#Multi-dimensional scaling plot of proximity matrix
png(file = "MDSplot.png", res=90)
MDSplot(model, train$Parity_cat)
dev.off()

# Momage ####
# Set directory

mom_data <- all_data_cst %>%
  dplyr::select(3:5229,5236)
mom_data <- mom_data[,colSums(is.na(mom_data))<nrow(mom_data)]
mom_data <- mom_data[,rowSums(is.na(mom_data))<ncol(mom_data)]
mom_data<-mom_data[complete.cases(mom_data), ]
# n=678

#Set random seed before splitting for repeatability 
set.seed(12345) 

#Set a portion of data aside for testing (80/20)
samp <- sample(nrow(mom_data), 0.8 * nrow(mom_data))
train <- mom_data[samp, ]
test <- mom_data[-samp, ]
#Check dimensions
dim(train)
dim(test)
# n=542 in train
# n=136 in test

#Run randomForest & tune
library(randomForest)
set.seed(12345) # need to reset seed for reproducibility 
# Start with model with default settings (ntrees=500,mtry=72)
model1 <- randomForest(Momage ~ ., data = train, proximity=TRUE)
model1
mom_pred <- predict(model1, test)
print(mean((mom_pred-test$Momage)^2))
# MSE 18.6.
set.seed(12345) # need to reset seed for reproducibility 
model2 <- randomForest(Momage ~ ., data = train, ntree=800,proximity=TRUE)
model2
mom_pred <- predict(model2, test)
print(mean((mom_pred-test$Momage)^2))
# MSE 18.5.
set.seed(12345) # need to reset seed for reproducibility 
model3 <- randomForest(Momage ~ ., data = train, ntree=1500,proximity=TRUE)
model3
mom_pred <- predict(model3, test)
print(mean((mom_pred-test$Momage)^2))
# MSE 18.6.
set.seed(12345) # need to reset seed for reproducibility 
model4 <- randomForest(Momage ~ ., data = train, ntree=1000,proximity=TRUE)
model4
mom_pred <- predict(model4, test)
print(mean((mom_pred-test$Momage)^2))
# MSE 18.5. best is 1,000.
set.seed(12345) # need to reset seed for reproducibility 
features <- setdiff(names(train),"Momage")
model5 <- tuneRF(
  x = train[features],
  y = train$Momage,
  ntreeTry=1000,
  mtryStart=5,
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE
)
# 33 mtry for lowest OOB error
set.seed(12345) # need to reset seed for reproducibility 
model6 <- randomForest(Momage ~ ., data = train, ntree = 1000, mtry=72, proximity=TRUE)
model6
mom_pred <- predict(model6, test)
print(mean((mom_pred-test$Momage)^2))
model=model4 # 4 was best
model

#Confusion matrix
model$confusion 

#Validate model using test data
prediction <- predict(model, newdata = test)
table(prediction, test$Momage)
prediction

#Display the predicted vs. the actual values
results<-cbind(prediction,test$Momage)
results
colnames(results)<-c('pred','real')
results<-as.data.frame(results)
View(results)

#Calculate model accuracy
sum(prediction==test$Momage) / nrow(test) 

#Plot random forest
png(file = "Rf plot.png", res=90)
plot(model)
dev.off()

#Number of nodes for trees
png(file = "Number of nodes.png", res=90)
hist(treesize(model),
     main = "No. of Nodes for the Trees ('Participant Age')")
#col = "red")
dev.off()

#Predictor importance
imp <- importance(model)
imp <- data.frame(predictors = rownames(imp), imp)
# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(IncNodePurity))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
# Select the top predictors
imp.20 <- imp.sort[1:20, ]

colnames(imp.20) <- c("predictors","Mean Decrease Gini")
png(file = "varImp.png", res=90, width = 800, height = 480)
varImp <- varImpPlot(model,
                     sort = T,
                     n.var = 20,
                     ylim=c(0,1400),
                     main = "Top 20 - Variable Importance ('Participant Age')")
dev.off()

# Higher the value of mean decrease accuracy or mean decrease Gini score, higher 
# the importance of the variable in the model

# ggplot
imp.20$predictors=factor(imp.20$predictors,levels=c("pt__00013",
                                                    "pt__00014",
                                                    "pt__00015",
                                                    "pt__00008",
                                                    "pt__00006",
                                                    "pt__00004",
                                                    "pt__00034",
                                                    "pt__00652",
                                                    "pt__00011",
                                                    "pt__00002",
                                                    "pt__00005",
                                                    "pt__00010",
                                                    "pt__00019",
                                                    "pt__00041",
                                                    "pt__00122",
                                                    "pt__00020",
                                                    "pt__00001",
                                                    "pt__00023",
                                                    "pt__00191",
                                                    "pt__00058"),
                         labels=c("Streptococcus oralis",
                                  "Dialister micraerophilus",
                                  "Aerococcus christensenii_1",
                                  "Lactobacillus gallinarum",
                                  "Sneathia amnii_1",
                                  "Atopobium vaginae_1",
                                  "Parvibacter caecicola / Adlercreutzia equolifaciens",
                                  "Peptoniphilus rhinitidis / Peptoniphilus harei",
                                  "Lactobacillus crispatus",
                                  "Lactobacillus iners_1",
                                  "Megasphaera genomosp. type_1_1",
                                  "Lactobacillus jensenii_1",
                                  "Ureaplasma parvum",
                                  "Finegoldia magna",
                                  "Varibaculum cambriense",
                                  "Megasphaera genomosp. type_1",
                                  "Gardnerella vaginalis",
                                  "Escherichia coli",
                                  "Actinotignum schaalii_1",
                                  "Corynebacterium tuberculostearicum"))
imp.20$IncNodePurity <- as.numeric(imp.20$IncNodePurity)
predimp5 <- ggplot(imp.20, aes(x = predictors, y = IncNodePurity)) +
  geom_bar(stat = "identity", fill = "indianred") +
  # geom_point() +
  scale_y_continuous(limits = c(0, 1400)) +
  # coord_flip() +
  ggtitle("Most Important Taxa for 'Participant Age'") +
  coord_flip() +
  xlab("Taxa") + ylab("Increase in Node Purity") +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face="bold"),    
        axis.title.y = element_text(face="bold")) 
ggsave("Preds importance.pdf",width=7)

#Multi-dimensional scaling plot of proximity matrix
png(file = "MDSplot.png", res=90)
MDSplot(model, train$Momage)
dev.off()

# Combine plots ####

# Set directory

predimp1 <- predimp1 +
  labs(title = "Self-reported White Race") +
  theme(
    axis.title.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y = element_text(vjust=0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 23))
predimp1
# predimp1 <- predimp1 +
#   geom_label_repel(aes("Taxa"))
predimp2 <- predimp2 +
  labs(title = "Education Level") +
  theme(
    axis.title.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y = element_text(vjust=0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 23))
predimp2
predimp3 <- predimp3 +
  labs(title = "Antibiotics Ever in Pregnancy") +
  theme(
    axis.title.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y = element_text(vjust=0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 23))
predimp4 <- predimp4 +
  labs(title = "Parity Category") +
  theme(
    axis.title.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y = element_text(vjust=0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 23))
predimp5 <- predimp5 +
  labs(title = "Participant Age") +
  theme(
    axis.title.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y = element_text(vjust=0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 23))
ggarrange(predimp1, predimp5, predimp2, predimp3, predimp4, 
          # labels = c("White", "Participant Age","Education Level", "Antibiotics Ever in Pregnancy", "Parity Category"),
          # font.label = list(size = 10),
          vjust = 0.3,
          hjust = 0,
          align = c("hv"),
          ncol = 3, nrow = 2)
ggsave("Panel Preds Importance.pdf", height=17, width=20)

roc_gg1 <- roc_gg1+
  xlab("Specificity") + ylab("Sensitivity") +
  theme(plot.title = element_text(face = "bold"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)) +
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="red", linetype="dashed") +
  labs(title = "Self-reported White Race (AUC = 0.973)") 
roc_gg2 <- roc_gg2 +
  xlab("Specificity") + ylab("Sensitivity") +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="red", linetype="dashed") +
  labs(title = "Education Level (AUC = 0.9025)")
roc_gg3 <- roc_gg3 +
  xlab("Specificity") + ylab("Sensitivity") +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="red", linetype="dashed") +
  labs(title = "Antibiotics Ever in Pregnancy (AUC = 0.7078)")
roc_gg4 <- roc_gg4 +
  xlab("Specificity") + ylab("Sensitivity") +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="red", linetype="dashed") +
  labs(title = "Parity Category (AUC = 0.6114)")
ggarrange(roc_gg1, roc_gg2, roc_gg3, roc_gg4, 
          # labels = c("White", "Participant Age","Education Level", "Antibiotics Ever in Pregnancy", "Parity Category"),
          # font.label = list(size = 10),
          align = c("hv"),
          ncol = 2, nrow = 2)
ggsave("Panel ROCs.pdf", width=9)
