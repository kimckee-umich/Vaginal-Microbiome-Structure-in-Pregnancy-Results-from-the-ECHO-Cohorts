######################################################################################
## Project: Validation of Maternal Vaginal Microbiota Signatures in Pregnancy: The ECHO Vaginal Microbiome Consortium - Dr Kimberly McKee
## Script name: Random Forests_Imputed Data.R
## Script purpose: Creating random forest plots with imputed data
## Date: May 2023
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
  dplyr::select(3:5229,5237,5241,5249,5255,5447,Momage,cohort,specimen)

all_data <- all_data %>%
  mutate(
    white=as.factor(white),
    AbxEver=as.factor(AbxEver),
    educ=as.factor(educ),
    Parity_cat=ifelse(Parity_cat=="No prior", 0,
                      ifelse(Parity_cat=="1", "1+",
                             ifelse(Parity_cat>1, "1+", NA))),
    Parity_cat=as.factor(Parity_cat),
    Momage=as.numeric(Momage)
  ) 

march <- all_data %>%
  filter(cohort=="MARCH")
emory <- all_data %>% 
  filter(cohort=="Emory")
wisc <- all_data %>% 
  filter(cohort=="Wisconsin")


# Data imputation ---------------------------------------------------------

march$Momage[is.na(march$Momage)] <- median(march$Momage,na.rm = TRUE)
emory$Momage[is.na(emory$Momage)] <- median(emory$Momage,na.rm = TRUE)
# Create mode() function to calculate mode
mode <- function(x, na.rm = FALSE) {
  if(na.rm){ #if na.rm is TRUE, remove NA values from input x
    x = x[!is.na(x)]
  }
  val <- unique(x)
  return(val[which.max(tabulate(match(x, val)))])
}
mode(march$white)
mode(wisc$white)
march$white[is.na(march$white)] <- 1
wisc$white[is.na(wisc$white)] <- 1
march$white <- as.factor(march$white)
wisc$white <- as.factor(wisc$white)
mode(march$educ) # BA or Higher
mode(wisc$educ) # BA or Higher
march$educ[is.na(march$educ)] <- "BA or Higher"
wisc$educ[is.na(wisc$educ)] <- "BA or Higher"
mode(march$Parity_cat)
march$Parity_cat[is.na(march$Parity_cat)] <- "1+"
march$Parity_cat <- as.factor(march$Parity_cat)
mode(emory$AbxEver) # No 
emory$AbxEver[is.na(emory$AbxEver)] <- "No"
emory$AbxEver <- as.factor(emory$AbxEver)
march_emory <- full_join(march,emory)
all_data_imp <- full_join(march_emory,wisc)
which(is.na(all_data_imp)) 
all_data = all_data_imp

all_data <- all_data %>%
  dplyr::select(-c(cohort,specimen))

# Random forest -----------------------------------------------------------

# All cohorts ####

# White, Maternal age, Education, Antibiotics Ever in Pregnancy, Parity Category

# white ####
# Set directory

white_data <- all_data %>%
  dplyr::select(1:5227,5232)

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
# n=544 in train
# n=136 in test

#Run randomForest & tune
library(randomForest)
set.seed(12345) # need to reset seed for reproducibility 
# Start with model with default settings (ntrees=500,mtry=72)
model1 <- randomForest(white ~ ., data = train, proximity=TRUE)
model1
# OOB error 4.78% so accuracy 95.22%.
set.seed(12345) # need to reset seed for reproducibility 
model2 <- randomForest(white ~ ., data = train, ntree=800,proximity=TRUE)
model2
# OOB error 4.96% so accuracy around 95.04%.
set.seed(12345) # need to reset seed for reproducibility 
model3 <- randomForest(white ~ ., data = train, ntree=700,proximity=TRUE)
model3
# OOB error 4.96% so accuracy around 95.04%.
# OOB error 4.63% so accuracy around 95.4%.
set.seed(12345) # need to reset seed for reproducibility 
features <- setdiff(names(train),"white")
model4 <- tuneRF(
  x = train[features],
  y = train$white,
  ntreeTry=500,
  mtryStart=5,
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE
)
# 33 mtry for lowest OOB error, but not performing so well
set.seed(12345) # need to reset seed for reproducibility 
model5 <- randomForest(white ~ ., data = train, ntree = 500, mtry=33, proximity=TRUE)
model5
model=model1
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
sum(prediction==test$white) / nrow(test) # 96.3% accuracy

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
                                                    "pt__00008",
                                                    "pt__00652",
                                                    "pt__00127",
                                                    "pt__00018",
                                                    "pt__00406",
                                                    "pt__00009",
                                                    "pt__00112",
                                                    "pt__00271",
                                                    "pt__00014",
                                                    "pt__00089",
                                                    "pt__00006",
                                                    "pt__00004",
                                                    "pt__00199",
                                                    "pt__00005",
                                                    "pt__00001",
                                                    "pt__00011",
                                                    "pt__00539",
                                                    "pt__00020",
                                                    "pt__00002"),
                         labels=c("Streptococcus oralis",
                                  "Lactobacillus gallinarum",
                                  "Peptoniphilus rhinitidis / Peptoniphilus harei",
                                  "Anaerococcus vaginalis",
                                  "Haemophilus parainfluenzae_1",
                                  "Chlamydia trachomatis / Elusimicrobium minutum",
                                  "Veillonella dispar",
                                  "Peptoniphilus harei / Peptoniphilus grossensis / Peptoniphilus vaginalis",
                                  "Peptoniphilus pacaensis",
                                  "Dialister micraerophilus_1",
                                  "Dialister propionicifaciens_2",
                                  "Sneathia amnii_1",
                                  "Atopobium vaginae_1",
                                  "Anaerococcus octavius_1",
                                  "Megasphaera genomosp. type_1_1",
                                  "Gardnerella vaginalis_1",
                                  "Lactobacillus crispatus",
                                  "Anaerococcus prevotii / Anaerococcus tetradius",
                                  "Megasphaera genomosp. type_1_2",
                                  "Lactobacillus iners_1"))

predimp1 <- ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most Important Taxa for 'White'") +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face="bold"),    
        axis.title.y = element_text(face="bold")) +
  xlab("Taxa") + ylab("Mean Decrease Gini")
ggsave("Preds importance.pdf", width=10)

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
library(xgboost)
library(SHAPforxgboost)
library(caret)
mtrain <- train %>%
  mutate(
    white=as.numeric(white)
  ) 
colnames(mtrain)[colnames(mtrain) %in% 
                   c("pt__00013", "pt__00652", "pt__00014","pt__00260","pt__04455")] <- c("Streptococcus oralis", "Peptoniphilus rhinitidis / Peptoniphilus harei","Dialister micraerophilus_1","Lachnoclostridium edouardi_2","Muribaculaceae")

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
# Low Peptoniphilus rhinitidis / Peptoniphilus harei associated with higher odds of being white. 
# Low Dialister micraerophilus_1 associated with lower odds of being white. 
# Some evidence of higher Lachnoclostridium edouardi_2 associated higher odds of being white. Mostly though, lower Lachnoclostridium edouardi_2 associated with lower odds of being white.    
# Low Muribaculaceae mostly associated with higher (but small) odds of being white. 
png(file = "SHAP plot.png", res=100, width=800)
shap.plot.summary.wrap1(bstDMatrix, X = as.matrix(mtrain[,-5228]), top_n = 5)
dev.off()

# educ ####
# Set directory

ed_data <- all_data %>%
  dplyr::select(1:5228)
# n=680

#Set random seed before splitting for repeatability 
set.seed(12345) 

#Set a portion of data aside for testing (80/20)
samp <- sample(nrow(ed_data), 0.8 * nrow(ed_data))
train <- ed_data[samp, ]
test <- ed_data[-samp, ]
#Check dimensions
dim(train)
dim(test)
# n=544 in train
# n=136 in test

#Run randomForest & tune
library(randomForest)
set.seed(12345) # need to reset seed for reproducibility 
# Start with model with default settings (ntrees=500,mtry=72)
model1 <- randomForest(educ ~ ., data = train, proximity=TRUE)
model1
# OOB error 47.98% so accuracy 52.02%.
set.seed(12345) # need to reset seed for reproducibility 
model2 <- randomForest(educ ~ ., data = train, ntree=800,proximity=TRUE)
model2
# OOB error 47.61% so accuracy around 52.39%. 800 best.
set.seed(12345) # need to reset seed for reproducibility 
model3 <- randomForest(educ ~ ., data = train, ntree=1500,proximity=TRUE)
model3
# OOB error 48.35% so accuracy around 51.65%.
set.seed(12345) # need to reset seed for reproducibility 
model4 <- randomForest(educ ~ ., data = train, ntree=1000,proximity=TRUE)
model4
# OOB error 47.79% so accuracy around 52.21%. 
set.seed(12345) # need to reset seed for reproducibility 
features <- setdiff(names(train),"educ")
model5 <- tuneRF(
  x = train[features],
  y = train$educ,
  ntreeTry=800,
  mtryStart=5,
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE
)
# 10 mtry for lowest OOB error--not close to default value of 72, but 72 performs better so keep it
set.seed(12345) # need to reset seed for reproducibility 
model6 <- randomForest(educ ~ ., data = train, ntree = 800, mtry=10, proximity=TRUE)
model6
model=model2
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
sum(prediction==test$educ) / nrow(test) # 50.7% accuracy

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
imp.20$predictors=factor(imp.20$predictors,levels=c("pt__00002",
                                                    "pt__00015",
                                                    "pt__00013",
                                                    "pt__00001",
                                                    "pt__00014",
                                                    "pt__00004",
                                                    "pt__00006",
                                                    "pt__00011",
                                                    "pt__00041",
                                                    "pt__00008",
                                                    "pt__00034",
                                                    "pt__00005",
                                                    "pt__00018",
                                                    "pt__00020",
                                                    "pt__00012",
                                                    "pt__00025",
                                                    "pt__00019",
                                                    "pt__00022",
                                                    "pt__00016",
                                                    "pt__00029"),
                         labels=c("Lactobacillus iners_1",
                                  "Aerococcus christensenii_1",
                                  "Streptococcus oralis",
                                  "Gardnerella vaginalis_1",
                                  "Dialister micraerophilus_1",
                                  'Atopobium vaginae_1',
                                  "Sneathia amnii_1",
                                  "Lactobacillus crispatus",
                                  "Finegoldia magna",
                                  "Lactobacillus gallinarum",
                                  "Parvibacter caecicola / Adlercreutzia equolifaciens",
                                  "Megasphaera genomosp. type_1_1",
                                  "Haemophilus parainfluenzae_1",
                                  "Megasphaera genomosp. type_1_2",
                                  "Prevotella timonensis_2",
                                  "Prevotella disiens_1",
                                  "Ureaplasma parvum",
                                  "Prevotella corporis_1",
                                  "Lactobacillus gasseri",
                                  "Haemophilus haemolyticus_1"))
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
# n=516

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
# n=104 in test

#Run randomForest & tune
library(randomForest)
set.seed(12345) # need to reset seed for reproducibility 
# Start with model with default settings (ntrees=500,mtry=72)
model1 <- randomForest(AbxEver ~ ., data = train, proximity=TRUE)
model1
# OOB error 35.19% so accuracy 64.81%. 500 best.
set.seed(12345) # need to reset seed for reproducibility 
model2 <- randomForest(AbxEver ~ ., data = train, ntree=800,proximity=TRUE)
model2
# OOB error 36.41% so accuracy around 63.59%.
set.seed(12345) # need to reset seed for reproducibility 
model3 <- randomForest(AbxEver ~ ., data = train, ntree=700,proximity=TRUE)
model3
# OOB error 35.92% so accuracy around 64.08%. 
set.seed(12345) # need to reset seed for reproducibility 
features <- setdiff(names(train),"AbxEver")
model4 <- tuneRF(
  x = train[features],
  y = train$AbxEver,
  ntreeTry=500,
  mtryStart=5,
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE
)
# 4 mtry for lowest OOB error--not close to default value of 72, but 72 performs better so keep it
set.seed(12345) # need to reset seed for reproducibility 
model5 <- randomForest(AbxEver ~ ., data = train, ntree = 500, mtry=4, proximity=TRUE)
model5
model=model1
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

# Higher the value of mean decrease accuracy or mean decrease Gini score, higher 
# the importance of the variable in the model

# ggplot
imp.20$predictors=factor(imp.20$predictors,levels=c("pt__00004",
                                                    "pt__00015",
                                                    "pt__00002",
                                                    "pt__00001",
                                                    "pt__00011",
                                                    "pt__00013",
                                                    "pt__00014",
                                                    "pt__00044",
                                                    "pt__00006",
                                                    "pt__00020",
                                                    "pt__00008",
                                                    "pt__00005",
                                                    "pt__00034",
                                                    "pt__00010",
                                                    "pt__00035",
                                                    "pt__00003",
                                                    "pt__00028",
                                                    "pt__00018",
                                                    "pt__00007",
                                                    "pt__00022"),
                         labels=c("Atopobium vaginae_1",
                                  "Aerococcus christensenii_1",
                                  "Lactobacillus iners_1",
                                  "Gardnerella vaginalis_1",
                                  "Lactobacillus crispatus",
                                  "Streptococcus oralis",
                                  "Dialister micraerophilus_1",
                                  "Mycoplasma hominis_1",
                                  "Sneathia amnii_1",
                                  "Megasphaera genomosp. type_1_2",
                                  "Lactobacillus gallinarum",
                                  "Megasphaera genomosp. type_1_1",
                                  "Parvibacter caecicola / Adlercreutzia equolifaciens",
                                  "Lactobacillus jensenii_1",
                                  "Staphylococcus epidermidis_1",
                                  "Prevotella buccalis_1",
                                  "Prevotella melaninogenica_1",
                                  "Haemophilus parainfluenzae_1",
                                  "Prevotella timonensis_1",
                                  "Prevotella corporis_1"))
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
# n=544 in train
# n=136 in test

#Run randomForest & tune
library(randomForest)
set.seed(12345) # need to reset seed for reproducibility 
# Start with model with default settings (ntrees=500,mtry=72)
model1 <- randomForest(Parity_cat ~ ., data = train, proximity=TRUE)
model1
# OOB error 44.67% so accuracy 55.33%.
set.seed(12345) # need to reset seed for reproducibility 
model2 <- randomForest(Parity_cat ~ ., data = train, ntree=800,proximity=TRUE)
model2
# OOB error 44.49% so accuracy around 55.51%.800 best.
set.seed(12345) # need to reset seed for reproducibility 
model3 <- randomForest(Parity_cat ~ ., data = train, ntree=1500,proximity=TRUE)
model3
# OOB error 45.59% so accuracy around 54.41%. 
set.seed(12345) # need to reset seed for reproducibility 
model4 <- randomForest(Parity_cat ~ ., data = train, ntree=1000,proximity=TRUE)
model4
# OOB error 44.85% so accuracy around 55.15%.
set.seed(12345) # need to reset seed for reproducibility 
features <- setdiff(names(train),"Parity_cat")
model5 <- tuneRF(
  x = train[features],
  y = train$Parity_cat,
  ntreeTry=800,
  mtryStart=5,
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE
)
# 7 mtry for lowest OOB error
set.seed(12345) # need to reset seed for reproducibility 
model6 <- randomForest(Parity_cat ~ ., data = train, ntree = 800, mtry=7, proximity=TRUE)
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
sum(prediction==test$Parity_cat) / nrow(test) # 65.4% accuracy

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
imp.20$predictors=factor(imp.20$predictors,levels=c("pt__00035",
                                                    "pt__00037",
                                                    "pt__00002",
                                                    "pt__00009",
                                                    "pt__00075",
                                                    "pt__00016",
                                                    "pt__00036",
                                                    "pt__00042",
                                                    "pt__00041",
                                                    "pt__00007",
                                                    "pt__00008",
                                                    "pt__00011",
                                                    "pt__00072",
                                                    "pt__00001",
                                                    "pt__00013",
                                                    "pt__00010",
                                                    "pt__00018",
                                                    "pt__00169",
                                                    "pt__00209",
                                                    "pt__02511"),
                         labels=c("Staphylococcus epidermidis_1",
                                  "Parvimonas micra_1",
                                  "Lactobacillus iners_1",
                                  "Veillonella dispar",
                                  "Bifidobacterium longum_1",
                                  "Lactobacillus gasseri",
                                  "Fenollaria timonensis / Fenollaria massiliensis",
                                  "Streptococcus salivarius_1",
                                  "Finegoldia magna",
                                  "Prevotella timonensis_1",
                                  "Lactobacillus gallinarum",
                                  "Lactobacillus crispatus",
                                  "Prevotella corporis_2",
                                  "Gardnerella vaginalis_1",
                                  "Streptococcus oralis",
                                  "Lactobacillus jensenii_1",
                                  "Haemophilus parainfluenzae_1",
                                  "Corynebacterium mycetoides",
                                  "Anaerococcus hydrogenalis",
                                  "Shuttleworthia satelles / Lachnobacterium bovis"))
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

mom_data <- all_data %>%
  dplyr::select(1:5227,5233)
# n=680

#Set random seed before splitting for repeatability 
set.seed(12345) 

#Set a portion of data aside for testing (80/20)
samp <- sample(nrow(mom_data), 0.8 * nrow(mom_data))
train <- mom_data[samp, ]
test <- mom_data[-samp, ]
#Check dimensions
dim(train)
dim(test)
# n=544 in train
# n=136 in test

#Run randomForest & tune
library(randomForest)
set.seed(12345) # need to reset seed for reproducibility 
# Start with model with default settings (ntrees=500,mtry=72)
model1 <- randomForest(Momage ~ ., data = train, proximity=TRUE)
model1
mom_pred <- predict(model1, test)
print(mean((mom_pred-test$Momage)^2))
# MSE 20.9.
set.seed(12345) # need to reset seed for reproducibility 
model2 <- randomForest(Momage ~ ., data = train, ntree=800,proximity=TRUE)
model2
mom_pred <- predict(model2, test)
print(mean((mom_pred-test$Momage)^2))
# MSE 20.8. Best is 800.
set.seed(12345) # need to reset seed for reproducibility 
model3 <- randomForest(Momage ~ ., data = train, ntree=1500,proximity=TRUE)
model3
mom_pred <- predict(model3, test)
print(mean((mom_pred-test$Momage)^2))
# MSE 21. Worst.
set.seed(12345) # need to reset seed for reproducibility 
model4 <- randomForest(Momage ~ ., data = train, ntree=1000,proximity=TRUE)
model4
mom_pred <- predict(model4, test)
print(mean((mom_pred-test$Momage)^2))
# MSE 20.9. 
set.seed(12345) # need to reset seed for reproducibility 
features <- setdiff(names(train),"Momage")
model5 <- tuneRF(
  x = train[features],
  y = train$Momage,
  ntreeTry=800,
  mtryStart=5,
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE
)
# 49 mtry for lowest OOB error
set.seed(12345) # need to reset seed for reproducibility 
model6 <- randomForest(Momage ~ ., data = train, ntree = 800, mtry=49, proximity=TRUE)
model6
mom_pred <- predict(model6, test)
print(mean((mom_pred-test$Momage)^2))
model=model6 
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
imp.20$predictors=factor(imp.20$predictors,levels=c("pt__00014",
                                                    "pt__00013",
                                                    "pt__00015",
                                                    "pt__00008",
                                                    "pt__00004",
                                                    "pt__00034",
                                                    "pt__00002",
                                                    "pt__00001",
                                                    "pt__00011",
                                                    "pt__00020",
                                                    "pt__00042",
                                                    "pt__00005",
                                                    "pt__00006",
                                                    "pt__00170",
                                                    "pt__00041",
                                                    "pt__00122",
                                                    "pt__00652",
                                                    "pt__00033",
                                                    "pt__00040",
                                                    "pt__00018"),
                         labels=c("Dialister micraerophilus",
                                  "Streptococcus oralis",
                                  "Aerococcus christensenii_1",
                                  "Lactobacillus gallinarum",
                                  "Atopobium vaginae_1",
                                  "Parvibacter caecicola / Adlercreutzia equolifaciens",
                                  "Lactobacillus iners_1",
                                  "Gardnerella vaginalis_1",
                                  "Lactobacillus crispatus",
                                  "Megasphaera genomosp. type_1",
                                  "Streptococcus salivarius_1",
                                  "Megasphaera genomosp. type_1_1",
                                  "Sneathia amnii_1",
                                  "Actinomyces neuii",
                                  "Finegoldia magna",
                                  "Varibaculum cambriense",
                                  "Peptoniphilus rhinitidis / Peptoniphilus harei",
                                  "Sneathia sanguinegens_1",
                                  "Prevotella buccalis_3",
                                  "Haemophilus parainfluenzae_1"))
imp.20$IncNodePurity <- as.numeric(imp.20$IncNodePurity)
predimp5 <- ggplot(imp.20, aes(x = predictors, y = IncNodePurity)) +
  geom_bar(stat = "identity", fill = "indianred") +
  # geom_point() +
  scale_y_continuous(limits = c(0, 1400)) +
  # coord_flip() +
  ggtitle("Most Important Taxa for 'Participant Age'") +
  coord_flip() +
  xlab("Taxa") + ylab("Increase in Node Purity ") +
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

library(ggpubr)
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