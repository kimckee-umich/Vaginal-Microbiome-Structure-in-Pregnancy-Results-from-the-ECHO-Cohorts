######################################################################################
## Project: Validation of Maternal Vaginal Microbiota Signatures in Pregnancy: The ECHO Vaginal Microbiome Consortium - Dr Kimberly McKee
## Script name: Data preparation_and_PERMANOVA.R
## Script purpose: Merging the harmonized data with the taxonomy data and conducting PERMANOVA
## Date: March 2023
## Author: Beatrice Palazzolo
## Organization: Department of Family Medicine, University of Michigan
######################################################################################

# Libraries ---------------------------------------------------------------

pkgs <- c(
  "readxl", "tidyverse", "nnet", "foreign", "tidyverse",
  "stargazer", "ggplot2", "poLCA", "vegan", "moderndive", "dslabs",
  "infer", "janitor", "remotes", "knitr", "ggplot2",
  "usethis", "readr", "skimr", "printr", "tab", "summarytools",
  "table1","lubridate", "vegan", "writexl", "kableExtra", "data.table", 
  "writexl", "ggpubr", "scales", "RColorBrewer", "ggpubr"
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

taxonomy <- read.csv("ECHO.phylotype.0.1.rel.csv")
harm_data <- read.csv("HarmonizedData.csv")

taxonomy <- taxonomy %>%
  dplyr::rename("specimen"="X")

#### Link MARCH ####

# Create new dataframe that includes only MARCH from harmonized data
harm_data %>% tabyl(cohort)
MARCH=harm_data %>% filter(cohort=="MARCH")
MARCH=MARCH %>% mutate(Studyid=as.numeric(substring(SAMPLEID,first=2)))

#Generate link file for UM
linkfile=read_excel("MARCH_VaginalFecalMetadata_id_timepoint.xlsx")
linkfile %>% tabyl(VagSample_Y)
linkfile=linkfile %>% filter(VagSample_Y==1)
linkfile=linkfile %>% filter(Time==3)
MARCH2=left_join(MARCH,linkfile,by="Studyid")

# Add MSU
MSUlink=read_excel("MSU_linkinfo_fromdietdata.xlsx")
MSUlink=MSUlink %>% dplyr::rename(SAMPLEID=CHARMID)
MARCH3=left_join(MARCH2,MSUlink,by="SAMPLEID")

# Generate link for MSU
MARCH3=MARCH3 %>% mutate(linkvar=ifelse(KIT=="77","MS007",
                                        ifelse(KIT=="80","MS008",
                                               ifelse(KIT=="78","MS014",
                                                      ifelse(KIT=="116","MS043",KIT)))))
MARCH3=MARCH3 %>% mutate(linkvar=ifelse(is.na(linkvar),seqID,linkvar))
#How many obs do not have a link ID?
sum(is.na(MARCH3$linkvar)) # 17 do not have a link ID

# Filter complete taxonomy to only include MARCH and rename specimens that are inconsistently named between MARCH taxonomy and MARCH harmonized data
taxo_MARCH <- taxonomy %>%
  filter(str_detect(specimen, "^M|^GUT")) %>%
  mutate(type=substring(specimen,1,3),
         linkvar=ifelse(type=="GUT",substring(specimen,first=4),
                        ifelse(specimen=="M_8139a","M_81390513192a",
                               ifelse(specimen=="M_8397a", "M_83970516193a",
                                      ifelse(specimen=="M_8411a", "M_84110514192a",
                                             ifelse(specimen=="M_8722a","M_87220514193a",specimen))))))

# Final MARCH taxonomy and MARCH harmonized data join
all_MARCH=taxo_MARCH %>% left_join(MARCH3,by="linkvar")
sum(is.na(all_MARCH$Lactobacillus.iners))
all_MARCH %>% tabyl(cohort)
all_MARCH = all_MARCH %>% filter(!is.na(cohort))

#### Link Emory ####

harm_Emory=harm_data %>% filter(cohort=="Emory")

#Generate link file for Emory
Emorylink=read.csv("Atlanta ECHO Metadata_Vaginal V1_McKEE OIF_ns.csv", na.strings=c("","NA"))
Emorylink %>% tabyl(VaginalVisit1SampleAvailable)
Emorylink=Emorylink %>% filter(VaginalVisit1SampleAvailable==1)

# Change dot in specimen to underscore
Emorylink$specimen <- str_replace_all(Emorylink$specimen,"[.]","_")

# Join linkfile and taxonomy
taxo_Emory <- taxonomy %>%
  filter(str_detect(specimen, "VagM1|E"))

taxo_Emory_2=left_join(Emorylink, taxo_Emory)

# Join with harmonized data
all_Emory=left_join(harm_Emory, taxo_Emory_2)
sum(is.na(all_Emory$Lactobacillus.iners)) # 3 NAs
all_Emory %>% tabyl(cohort)
all_Emory = all_Emory %>% filter(!is.na(cohort))

# Join all MARCH and Emory data
MARCH_Emory <- full_join(all_MARCH, all_Emory)
# Remove columns that have all NA values
MARCH_Emory <- MARCH_Emory[,colSums(is.na(MARCH_Emory))<nrow(MARCH_Emory)]
# Remove any rows that only have NA values
MARCH_Emory %>% filter(rowSums(is.na(.)) != ncol(.))

#### Link Wisconsin ####

harm_Wisconsin <- harm_data %>%
  filter(cohort=="Wisconsin")

# CREW-WISC ####
CREW_mapping=read_excel("WISC_CREW_key_mapping.xlsx")
CREW_mapping <- CREW_mapping %>%
  filter(SID_num!="NA") %>%
  dplyr::select(specimen=Specimen, SID_NUM=SID_num)

# Join CREW-WISC metadata to taxonomy
taxo_CREW <- taxonomy %>%
  filter(str_detect(specimen, "^V")) 
CREW_taxo <- full_join(taxo_CREW,CREW_mapping, by="specimen")
CREW_taxo <- na.omit(CREW_taxo)  

# Join CREW-WISC metadata to harmonized data
CREW_taxo <- CREW_taxo %>%
  mutate(SID_NUM=as.numeric(SID_NUM))
CREWlink <- read.csv("McKee_eLab_data_pull_WISC_040121.csv", na.strings=c("","NA"))
CREWlink <- CREWlink %>%
  dplyr::select(SID_NUM,SUBJID)
linkfile_2=left_join(harm_Wisconsin,CREWlink,by="SUBJID")
# Remove columns that have all NA values
linkfile_2 <- linkfile_2[,colSums(is.na(linkfile_2))<nrow(linkfile_2)]
# Join
WISC1=left_join(linkfile_2,CREW_taxo,by="SID_NUM")
# Check that cohort data are complete and remove if no ID or specimen
WISC1 <- WISC1 %>%
  filter(!is.na(cohort))
WISC1 <- WISC1[complete.cases(WISC1$SID_NUM),]
WISC1 <- WISC1[complete.cases(WISC1$specimen),]

# MAAP ####

# Join MAAP metadata to harmonized data 
Wisclink_2 <- read.csv("McKee_eLab_data_pull_MAAP_040121.csv", na.strings=c("","NA"))
Wisclink_2 <- Wisclink_2 %>% dplyr::rename(famID=FAMILYID)
WISC_2=left_join(harm_Wisconsin,Wisclink_2)
WISC_2 <- WISC_2 %>%
  filter(!is.na(cohort))

# Dataset where MAAP specimen data are
MAAP_mapping=read.csv("McKee_mapping.csv", na.strings=c("","NA"))
# Filter out specimens sequenced for a separate analysis but not needed for the PERMANOVA analysis //
# These specimens are recorded in columns 3:5 as E+digits codes, so we keep the NAs in those columns
MAAP_nna <- MAAP_mapping %>% filter_at(vars(lynchLabID_mom_stool_pd,lynchLabID_mom_stool_7d,
                                            lynchLabID_mom_stool_6m),all_vars(is.na(.)))
# Generate MAAP specimens by collating IDs from 4 variables 
id_1 <- MAAP_nna %>%
  dplyr::select(specimen=lynchLabID_mom_vagrect_pd, famID) %>%
  filter(!is.na(specimen))
id_2 <- MAAP_nna %>%
  dplyr::select(specimen=lynchLabID_mom_vagrect_del, famID) %>%
  filter(!is.na(specimen))
id_3 <- MAAP_nna %>%
  dplyr::select(specimen=lynchLabID_mom_rect_pd, famID) %>%
  filter(!is.na(specimen))
id_4 <- MAAP_nna %>%
  dplyr::select(specimen=lynchLabID_mom_rect_del, famID) %>%
  filter(!is.na(specimen))
Wisc_specimen <- rbind(id_1,id_2,id_3,id_4)
Wisc_complete <- full_join(MAAP_nna, Wisc_specimen, by="famID") 
Wisc_complete <- Wisc_complete %>%
  dplyr::select(1,27)

# Join mapping & specimen data with MAAP metadata & harmonized data 
Wisc_complete_final = left_join(Wisc_complete,WISC_2) 
# Remove columns that have all NA values
Wisc_complete_final <- Wisc_complete_final[,colSums(is.na(Wisc_complete_final))<nrow(Wisc_complete_final)]
# And with taxonomy
taxo_MAAP <- taxonomy %>%
  filter(str_detect(specimen, "^F")) 
all_MAAP=left_join(Wisc_complete_final,taxo_MAAP)
# Remove columns that have all NA values
all_MAAP <- all_MAAP[,colSums(is.na(all_MAAP))<nrow(all_MAAP)]
all_MAAP <- all_MAAP %>%
  filter(!is.na(cohort))

# Join all Wisconsin data

all_Wisconsin=full_join(WISC1,all_MAAP)

# Join with MARCH & EMORY

all_data_complete=full_join(MARCH_Emory,all_Wisconsin)
all_data_complete <- all_data_complete[complete.cases(all_data_complete[,c(2:5233)]),]

# TAXO contains only the taxonomy data
TAXO = all_data_complete[,2:5233]
dim(TAXO)

#Verify no taxonomy with ALL zero reads
allzero = 0
# Find which columns have all zeros
for(i in 1:(ncol(TAXO))){ 
  if(isTRUE(sum(TAXO[,i])==0)){allzero==c(allzero,1)}
  else{allzero=c(allzero, 0)}
}
# Remove the 0 at the start of allzero vector
allzero = allzero[2:length(allzero)]
allzero

# Save data --------------------------------------------------------------

all_MARCH <- all_MARCH %>%
  dplyr::select(-c(pt__02970,
                   pt__02985,
                   pt__03080,
                   pt__03110,
                   pt__03203))
saveRDS(all_MARCH, file = "all_MARCH.rds")
write_xlsx(all_MARCH, "all_MARCH.xlsx")

all_Emory <- all_Emory %>%
  dplyr::select(-c(pt__02970,
                   pt__02985,
                   pt__03080,
                   pt__03110,
                   pt__03203))
saveRDS(all_Emory, file = "all_Emory.rds")
write_xlsx(all_Emory, "all_Emory.xlsx")

MARCH_Emory <- MARCH_Emory %>%
  dplyr::select(-c(pt__02970,
                   pt__02985,
                   pt__03080,
                   pt__03110,
                   pt__03203))
saveRDS(MARCH_Emory, file = "MARCH_Emory.rds")
write_xlsx(MARCH_Emory, "MARCH_Emory.xlsx")

WISC1 <- WISC1 %>%
  dplyr::select(-c(pt__02970,
                   pt__02985,
                   pt__03080,
                   pt__03110,
                   pt__03203))
saveRDS(WISC1, file = "WISC1.rds")
write_xlsx(WISC1, "WISC1.xlsx")

all_MAAP <- all_MAAP %>%
  dplyr::select(-c(pt__02970,
                   pt__02985,
                   pt__03080,
                   pt__03110,
                   pt__03203))
saveRDS(all_MAAP, file = "all_MAAP.rds")
write_xlsx(all_MAAP, "all_MAAP.xlsx")

all_Wisconsin <- all_Wisconsin %>%
  dplyr::select(-c(pt__02970,
                   pt__02985,
                   pt__03080,
                   pt__03110,
                   pt__03203))
saveRDS(all_Wisconsin, file = "all_Wisconsin.rds")
write_xlsx(all_Wisconsin, "all_Wisconsin.xlsx")

all_data_complete <- all_data_complete %>%
  dplyr::select(-c(pt__02970,
                   pt__02985,
                   pt__03080,
                   pt__03110,
                   pt__03203))
saveRDS(all_data_complete, file = "all_data_complete.rds")
write_xlsx(all_data_complete, "all_data_complete.xlsx")

TAXO <- TAXO %>%
  dplyr::select(-c(pt__02970,
                   pt__02985,
                   pt__03080,
                   pt__03110,
                   pt__03203))
saveRDS(TAXO, file = "TAXO.rds")
write_xlsx(TAXO, "TAXO.xlsx")

# Dispersion diagnostics --------------------------------------------------

# Set directory

# Cohort ####
# Dispersion diagnostics performed using the three cohorts as the groups 
march <- all_data_complete %>%
  dplyr::select(2:5228,cohort) %>%
  filter(cohort=="MARCH")
emory <- all_data_complete %>%
  dplyr::select(2:5228,cohort) %>%
  filter(cohort=="Emory")
wisconsin <- all_data_complete %>%
  dplyr::select(2:5228,cohort) %>%
  filter(cohort=="Wisconsin")
data <- rbind(march,emory,wisconsin)
table(data$cohort)
groups <- factor(c(rep("MARCH",123),rep("Emory",393),rep("Wisconsin",164)))
data <- data %>%
  dplyr::select(-cohort)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.8733, so do not reject the null hypothesis.
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-values are all non-significant.
jpeg(file="PCoA plot.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# white ####
# Dispersion diagnostics performed using white as the group
data_complete=all_data_complete %>% mutate(white=ifelse(race_recode=="White",1,
                                                        ifelse(is.na(race_recode),NA,0)),
                                           white=as.factor(white))
data <- data_complete %>%
  dplyr::select(2:5228,white) 
white <- data_complete %>%
  dplyr::select(2:5228,white) %>%
  filter(white=="1")
nonwhite <- data_complete %>%
  dplyr::select(2:5228,white) %>%
  filter(white=="0")
data <- rbind(white,nonwhite)
data=data %>% mutate(white=ifelse(white=="1","White",
                                  ifelse(is.na(white),NA,"Non-white")),
                     white=as.factor(white))
table(data$white)
groups <- factor(c(rep("White",250),rep("Non-white",425)))
data <- data %>%
  dplyr::select(-white)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.9273, so do not reject the null hypothesis.
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-values are all non-significant.
jpeg(file="PCoA plot_white.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# hispanic ####
# Dispersion diagnostics performed using hispanic as the group
data <- all_data_complete %>%
  dplyr::select(2:5228,hispanic) %>%
  mutate(hispanic=as.factor(hispanic))
hisp <- data %>%
  dplyr::select(2:5228,hispanic) %>%
  filter(hispanic=="Hispanic")
nonhisp <- data %>%
  dplyr::select(2:5228,hispanic) %>%
  filter(hispanic=="Non-hispanic")
data <- rbind(hisp,nonhisp)
table(data$hispanic)
groups <- factor(c(rep("Hispanic",13),rep("Non-hispanic",658)))
data <- data %>%
  dplyr::select(-hispanic)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.6439, so do not reject the null hypothesis.
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-values are all non-significant.
jpeg(file="PCoA plot_hispanic.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# educ ####
# Dispersion diagnostics performed using education levels as the groups 
data <- all_data_complete %>%
  dplyr::select(2:5228,educ) %>%
  mutate(educ=as.factor(educ))
ba <- data %>%
  dplyr::select(2:5228,educ) %>%
  filter(educ=="BA or Higher")
hsged <- data %>%
  dplyr::select(2:5228,educ) %>%
  filter(educ=="HS/GED")
hs <- data %>%
  dplyr::select(2:5228,educ) %>%
  filter(educ=="Less than HS")
coll <- data %>%
  dplyr::select(2:5228,educ) %>%
  filter(educ=="Some College/Assoc.")
data <- rbind(ba,hsged,hs,coll)
table(data$educ)
groups <- factor(c(rep("BA or Higher",252),rep("HS/GED",171),rep("Less than HS",63),rep("Some College/Assoc.",192)))
data <- data %>%
  dplyr::select(-educ)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.04812, so reject the null hypothesis.
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-values are all non-significant.
jpeg(file="PCoA plot_educ.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# PublicIns ####
# Dispersion diagnostics performed using the PublicIns as the group
data <- all_data_complete %>%
  dplyr::select(2:5228,PublicIns) %>%
  mutate(PublicIns=ifelse(PublicIns=="Yes","Public Insurance",
                          ifelse(PublicIns=="No", "No Public Insurance",NA)),
         PublicIns=as.factor(PublicIns))
pubins <- data %>%
  dplyr::select(2:5228,PublicIns) %>%
  filter(PublicIns=="Public Insurance")
nopubins <- data %>%
  dplyr::select(2:5228,PublicIns) %>%
  filter(PublicIns=="No Public Insurance")
data <- rbind(pubins,nopubins)
table(data$PublicIns)
groups <- factor(c(rep("Public Insurance",332),rep("No Public Insurance",228)))
data <- data %>%
  dplyr::select(-PublicIns)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.07621, so do not reject the null hypothesis (at 5%).
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-value is non-significant (at 5%).
jpeg(file="PCoA plot_PublicIns.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# AbxEver ####
# Dispersion diagnostics performed using the AbxEver as the group
data <- all_data_complete %>%
  dplyr::select(2:5228,AbxEver) %>%
  mutate(AbxEver=ifelse(AbxEver=="Yes","Antibiotics in Pregnancy",
                        ifelse(AbxEver=="No", "No Antibiotics in Pregnancy",NA)),
         AbxEver=as.factor(AbxEver))
abx <- data %>%
  dplyr::select(2:5228,AbxEver) %>%
  filter(AbxEver=="Antibiotics in Pregnancy")
noabx <- data %>%
  dplyr::select(2:5228,AbxEver) %>%
  filter(AbxEver=="No Antibiotics in Pregnancy")
data <- rbind(abx,noabx)
table(data$AbxEver)
groups <- factor(c(rep("Antibiotics in Pregnancy",195),rep("No Antibiotics in Pregnancy",320)))
data <- data %>%
  dplyr::select(-AbxEver)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.1517, so do not reject the null hypothesis.
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-value is non-significant.
jpeg(file="PCoA plot_abxever.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# SexB ####
# Dispersion diagnostics performed using the SexB as the group
data <- all_data_complete %>%
  dplyr::select(2:5228,SexB) %>%
  mutate(SexB=as.factor(SexB))
female <- data %>%
  dplyr::select(2:5228,SexB) %>%
  filter(SexB=="Female")
male <- data %>%
  dplyr::select(2:5228,SexB) %>%
  filter(SexB=="Male")
data <- rbind(female,male)
table(data$SexB)
groups <- factor(c(rep("Female",346),rep("Male",333)))
data <- data %>%
  dplyr::select(-SexB)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.7023, so do not reject the null hypothesis.
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-value is non-significant.
jpeg(file="PCoA plot_SexB.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# smoking ####
# Dispersion diagnostics performed using smoking as the group 
data <- all_data_complete %>%
  dplyr::select(2:5228,smoking) %>%
  mutate(smoking=ifelse(smoking=="No","Non-smoking",
                        ifelse(smoking=="Yes","Smoking",NA)),
         smoking=as.factor(smoking))
smoking <- data %>%
  dplyr::select(2:5228,smoking) %>%
  filter(smoking=="Smoking")
nonsmoking <- data %>%
  dplyr::select(2:5228,smoking) %>%
  filter(smoking=="Non-smoking")
data <- rbind(smoking,nonsmoking)
table(data$smoking)
groups <- factor(c(rep("Smoking",82),rep("Non-smoking",587)))
data <- data %>%
  dplyr::select(-smoking)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.03933, so reject the null hypothesis.
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-value is significant.
jpeg(file="PCoA plot_smoking.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# gest_diab ####
# Dispersion diagnostics performed using gest_diab as the group
data <- all_data_complete %>%
  dplyr::select(2:5228,gest_diab) %>%
  mutate(gest_diab=ifelse(gest_diab=="Yes","Gestational Diabetes",
                          ifelse(gest_diab=="No","No Gestational Diabetes",NA)),
         gest_diab=as.factor(gest_diab))
gestdiab <- data %>%
  dplyr::select(2:5228,gest_diab) %>%
  filter(gest_diab=="Gestational Diabetes")
nogestdiab <- data %>%
  dplyr::select(2:5228,gest_diab) %>%
  filter(gest_diab=="No Gestational Diabetes")
data <- rbind(gestdiab,nogestdiab)
table(data$gest_diab)
groups <- factor(c(rep("Gestational Diabetes",23),rep("No Gestational Diabetes",602)))
data <- data %>%
  dplyr::select(-gest_diab)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.9877, so do not reject the null hypothesis.
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-value is not significant.
jpeg(file="PCoA plot_gestdiab.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# htn ####
# Dispersion diagnostics performed using htn as the group
data <- all_data_complete %>%
  dplyr::select(2:5228,htn) %>%
  mutate(htn=ifelse(htn=="Yes","Hypertension",
                    ifelse(htn=="No","No Hypertension",NA)),
         htn=as.factor(htn))
htn <- data %>%
  dplyr::select(2:5228,htn) %>%
  filter(htn=="Hypertension")
nohtn <- data %>%
  dplyr::select(2:5228,htn) %>%
  filter(htn=="No Hypertension")
data <- rbind(htn,nohtn)
table(data$htn)
groups <- factor(c(rep("Hypertension",77),rep("No Hypertension",548)))
data <- data %>%
  dplyr::select(-htn)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.04475, so reject the null hypothesis.
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-value is significant.
jpeg(file="PCoA plot_htn.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# BMI_cat ####
# Dispersion diagnostics performed using BMI_cat as the groups 
data <- all_data_complete %>%
  dplyr::select(2:5228,BMI_cat) %>%
  mutate(BMI_cat=as.factor(BMI_cat))
normal <- data %>%
  dplyr::select(2:5228,BMI_cat) %>%
  filter(BMI_cat=="Normal(>=18.5 & <25)")
under <- data %>%
  dplyr::select(2:5228,BMI_cat) %>%
  filter(BMI_cat=="Underweight(<18.5)")
over <- data %>%
  dplyr::select(2:5228,BMI_cat) %>%
  filter(BMI_cat=="Overweight(>=25 & <30)")
obese <- data %>%
  dplyr::select(2:5228,BMI_cat) %>%
  filter(BMI_cat=="Obese(>=30)")
data <- rbind(normal,under,over,obese)
table(data$BMI_cat)
groups <- factor(c(rep("Normal(>=18.5 & <25",271),rep("Underweight(<18.5)",22),rep("Overweight(>=25 & <30)",142),rep("Obese(>=30)",191)))
data <- data %>%
  dplyr::select(-BMI_cat)
dis <- vegdist(data,method="bray")
mod <- betadisper(dis, groups)
anova(mod)
# Null hypothesis: no difference in dispersion between groups.
# p-value is 0.11, so do not reject the null hypothesis.
TukeyHSD(mod)
# Tukey's test verifies if and which groups differ in relation to their variances.
# Adjusted p-values are all non-significant.
jpeg(file="PCoA plot_BMI.jpeg")
plot(mod, main="Groups and Distances to Centroids on the First Two PCoA Axes")
dev.off()
# Plotting the groups and distances to centroids on the first two PCoA axes
# Interpretation: objects ordinated closer to one another are more similar //
# than those ordinated further away. 

# Parity_cat ####
# Analysis ---------------------------------------------------------------

# All cohorts merged ####

set.seed(2922)

all_data_complete=all_data_complete %>% mutate(white=ifelse(race_recode=="White",1,
                                                            ifelse(is.na(race_recode),NA,0)))

saveRDS(all_data_complete, file = "all_data_complete.rds")
write_xlsx(all_data_complete, "all_data_complete.xlsx")

#Save the columns that are covariates
covs=c("white", "hispanic", "Momage", "educ", "PublicIns", "AbxEver",
       "SexB", "smoking", "gest_diab","htn", "BMI_cat", "Parity_cat")
pvals=1:length(covs)
rsq=1:length(covs)

all_data_complete <- all_data_complete %>%
  dplyr::select(1:5228,5234:5237,5240,5245,5247,5249,5250,5253,5254,5446)

for(i in 1:length(covs)){
  
  all_data_complete %>% tabyl(covs[i])
  
  TAXO_temp=TAXO[complete.cases(all_data_complete[,covs[i]]),]
  alldata_temp=all_data_complete[complete.cases(all_data_complete[,covs[i]]),]
  modelfit=adonis2(TAXO_temp~alldata_temp[,covs[i]], data=alldata_temp, method="bray", perm = 10000)
  pvals[i] = modelfit$`Pr(>F)`[1]
  rsq[i] = modelfit$R2[1]
}

# Unadjusted model ####

TAXO_sub=TAXO[complete.cases(all_data_complete[,covs]),]

all_sub=all_data_complete[complete.cases(all_data_complete[,covs]),]

permanova=as.data.frame(cbind(covs, pvals, rsq))
permanova$pvaladj = p.adjust(permanova$pvals,"fdr")

permanova$pval_cat[which(permanova$pvaladj<0.001)] = 1 
permanova$pval_cat[which(permanova$pvaladj>=0.001 & permanova$pvaladj<0.01)] = 2 
permanova$pval_cat[which(permanova$pvaladj>=0.01 & permanova$pvaladj<0.05)] = 3 
permanova$pval_cat[which(permanova$pvaladj>=0.05 & permanova$pvaladj<0.1)] = 4 
permanova$pval_cat[which(permanova$pvaladj>=0.1 & permanova$pvaladj<0.2)] = 5 
permanova$pval_cat[which(permanova$pvaladj>=0.2)] = 6
permanova$pval_cat = factor(permanova$pval_cat,levels=c(1,2,3,4,5,6),
                            labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

permanova$covs = factor(permanova$covs,
                        levels=c( "white", "hispanic", "Momage", "educ", "PublicIns", "AbxEver",
                                  "SexB", "smoking", "gest_diab","htn", "BMI_cat", "Parity_cat"),
                        labels=c("Self-reported White Race", "Self-reported Hispanic", "Participant Age", "Education Level", "Public Insurance",
                                 "Antibiotics Ever in Pregnancy",
                                 "Infant Sex", "Smoking in Pregnancy",
                                 "Gestational Diabetes", "Hypertension",
                                 "BMI Category", "Parity Category"))
permanova$temp=rep(1, dim(permanova)[1])
permanova$rsq=as.numeric(permanova$rsq)
permanova$pvals=as.numeric(permanova$pvals)

knitr::kable(permanova[,-c(5:6)],digits=4, caption="Multiple Factor Model Results", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/.../Multiple Factor Model Results_all cohorts.png") 

permanova$temp=factor(permanova$temp, levels=1, labels="Multiple Factor Model")

p_unadj <- ggplot(permanova, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red", "lightgoldenrod"))+  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_unadj

# Adjusted model ####

# Keep variables with FDR-adjusted p-values <0.2.
# Drop any variables that are redundant (like private insurance and race).

TAXO_sub_2=TAXO[complete.cases(all_data_complete[,c("white","Momage", "educ","PublicIns","AbxEver","htn","BMI_cat", "Parity_cat")]),]
all_sub_2=all_data_complete[complete.cases(all_data_complete[,c("white","Momage", "educ","PublicIns","AbxEver","htn","BMI_cat", "Parity_cat")]),]

fullmodel=adonis2(TAXO_sub_2~white+Momage+educ+PublicIns+AbxEver+htn+BMI_cat+Parity_cat, data=all_sub_2, method="bray", perm = 10000)

newcovs=c("white","Momage", "educ","PublicIns","AbxEver","htn","BMI_cat", "Parity_cat")

fullmodel1 <- fullmodel %>% dplyr::select(5,3) 
fullmodel1 <- fullmodel1 [-c(9,10), ]
adj_permanova=as.data.frame(cbind(newcovs, fullmodel1))
names(adj_permanova)=c("newcovs", "pvals","rsq")
adj_permanova$pvals=as.numeric(adj_permanova$pvals)
adj_permanova$pvaladj = p.adjust(adj_permanova$pvals,"fdr")
adj_permanova$rsq=as.numeric(adj_permanova$rsq)

adj_permanova$newcovs=factor(adj_permanova$newcovs,levels=c("white","Momage", "educ","PublicIns","AbxEver","htn","BMI_cat", "Parity_cat"),
                             labels=c("Self-reported White Race", "Participant Age", "Education Level", "Public Insurance", "Antibiotics Ever in Pregnancy", "Hypertension", "BMI Category", "Parity Category"))

knitr::kable(adj_permanova[,c(1:4)], digits=4, caption="Adjusted Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable(".../Adjusted Model Results_all cohorts.png") 

adj_permanova$pval_cat = rep(0,length(adj_permanova$pvaladj))
adj_permanova$pval_cat[which(adj_permanova$pvaladj<0.001)] = 1 
adj_permanova$pval_cat[which(adj_permanova$pvaladj>=0.001 & adj_permanova$pvaladj<0.01)] = 2 
adj_permanova$pval_cat[which(adj_permanova$pvaladj>=0.01 & adj_permanova$pvaladj<0.05)] = 3 
adj_permanova$pval_cat[which(adj_permanova$pvaladj>=0.05 & adj_permanova$pvaladj<0.1)] = 4 
adj_permanova$pval_cat[which(adj_permanova$pvaladj>=0.1 & adj_permanova$pvaladj<0.2)] = 5 
adj_permanova$pval_cat[which(adj_permanova$pvaladj>=0.2)] = 6
adj_permanova$pval_cat = factor(adj_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

adj_permanova$temp=rep(3,length(newcovs))
adj_permanova <- adj_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Adjusted Model", NA),
    temp=as.factor(temp)
  )

p_adj <- ggplot(adj_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red3", "red", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  guides(colour = guide_legend(order = 1),size = guide_legend(order = 2))
p_adj 

full_data=bind_rows(permanova,adj_permanova)

full_plot <- ggplot(full_data, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)")
full_plot

ggsave(".../Full Model Plot.png") 

# White ####

TAXO_white=TAXO[complete.cases(all_data_complete[,c("white")]),]
all_white=all_data_complete[complete.cases(all_data_complete[,c("white")]),]

white=adonis2(TAXO_white~white, data=all_white, method="bray", perm = 10000)

newcovs=c("white")

white1 <- white %>% dplyr::select(5,3) 
white1 <- white1 [-c(2,3), ]
white_permanova=as.data.frame(cbind(newcovs, white1))
names(white_permanova)=c("newcovs", "pvals","rsq")
white_permanova$pvals=as.numeric(white_permanova$pvals)
white_permanova$pvaladj = p.adjust(white_permanova$pvals,"fdr")
white_permanova$rsq=as.numeric(white_permanova$rsq)

white_permanova$newcovs=factor(white_permanova$newcovs,levels=c("white"),
                             labels=c("Self-reported White Race"))

knitr::kable(white_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

white_permanova$pval_cat = rep(0,length(white_permanova$pvaladj))
white_permanova$pval_cat[which(white_permanova$pvaladj<0.001)] = 1 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.001 & white_permanova$pvaladj<0.01)] = 2 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.01 & white_permanova$pvaladj<0.05)] = 3 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.05 & white_permanova$pvaladj<0.1)] = 4 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.1 & white_permanova$pvaladj<0.2)] = 5 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.2)] = 6
white_permanova$pval_cat = factor(white_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

white_permanova$temp=rep(3,length(newcovs))
white_permanova <- white_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_white <- ggplot(white_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red3", "red", "coral1", "sandybrown", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_white

# Hispanic ####

TAXO_hisp=TAXO[complete.cases(all_data_complete[,c("hispanic")]),]
all_hisp=all_data_complete[complete.cases(all_data_complete[,c("hispanic")]),]

hisp=adonis2(TAXO_hisp~hispanic, data=all_hisp, method="bray", perm = 10000)

newcovs=c("hispanic")

hisp1 <- hisp %>% dplyr::select(5,3) 
hisp1 <- hisp1 [-c(2,3), ]
hisp_permanova=as.data.frame(cbind(newcovs, hisp1))
names(hisp_permanova)=c("newcovs", "pvals","rsq")
hisp_permanova$pvals=as.numeric(hisp_permanova$pvals)
hisp_permanova$pvaladj = p.adjust(hisp_permanova$pvals,"fdr")
hisp_permanova$rsq=as.numeric(hisp_permanova$rsq)

hisp_permanova$newcovs=factor(hisp_permanova$newcovs,levels=c("hispanic"),
                               labels=c("Self-reported Hispanic Race"))

knitr::kable(hisp_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

hisp_permanova$pval_cat = rep(0,length(hisp_permanova$pvaladj))
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj<0.001)] = 1 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.001 & hisp_permanova$pvaladj<0.01)] = 2 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.01 & hisp_permanova$pvaladj<0.05)] = 3 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.05 & hisp_permanova$pvaladj<0.1)] = 4 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.1 & hisp_permanova$pvaladj<0.2)] = 5 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.2)] = 6
hisp_permanova$pval_cat = factor(hisp_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                  labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

hisp_permanova$temp=rep(3,length(newcovs))
hisp_permanova <- hisp_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_hisp <- ggplot(hisp_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_hisp

# Momage ####

TAXO_age=TAXO[complete.cases(all_data_complete[,c("Momage")]),]
all_age=all_data_complete[complete.cases(all_data_complete[,c("Momage")]),]

age=adonis2(TAXO_age~Momage, data=all_age, method="bray", perm = 10000)

newcovs=c("Momage")

age1 <- age %>% dplyr::select(5,3) 
age1 <- age1 [-c(2,3), ]
age_permanova=as.data.frame(cbind(newcovs, age1))
names(age_permanova)=c("newcovs", "pvals","rsq")
age_permanova$pvals=as.numeric(age_permanova$pvals)
age_permanova$pvaladj = p.adjust(age_permanova$pvals,"fdr")
age_permanova$rsq=as.numeric(age_permanova$rsq)

age_permanova$newcovs=factor(age_permanova$newcovs,levels=c("Momage"),
                              labels=c("Participant Age"))

knitr::kable(age_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

age_permanova$pval_cat = rep(0,length(age_permanova$pvaladj))
age_permanova$pval_cat[which(age_permanova$pvaladj<0.001)] = 1 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.001 & age_permanova$pvaladj<0.01)] = 2 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.01 & age_permanova$pvaladj<0.05)] = 3 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.05 & age_permanova$pvaladj<0.1)] = 4 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.1 & age_permanova$pvaladj<0.2)] = 5 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.2)] = 6
age_permanova$pval_cat = factor(age_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

age_permanova$temp=rep(3,length(newcovs))
age_permanova <- age_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_age <- ggplot(age_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_age

# educ ####

TAXO_educ=TAXO[complete.cases(all_data_complete[,c("educ")]),]
all_educ=all_data_complete[complete.cases(all_data_complete[,c("educ")]),]

educ=adonis2(TAXO_educ~educ, data=all_educ, method="bray", perm = 10000)

newcovs=c("educ")

educ1 <- educ %>% dplyr::select(5,3) 
educ1 <- educ1 [-c(2,3), ]
educ_permanova=as.data.frame(cbind(newcovs, educ1))
names(educ_permanova)=c("newcovs", "pvals","rsq")
educ_permanova$pvals=as.numeric(educ_permanova$pvals)
educ_permanova$pvaladj = p.adjust(educ_permanova$pvals,"fdr")
educ_permanova$rsq=as.numeric(educ_permanova$rsq)

educ_permanova$newcovs=factor(educ_permanova$newcovs,levels=c("educ"),
                             labels=c("Education Level"))

knitr::kable(educ_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

educ_permanova$pval_cat = rep(0,length(educ_permanova$pvaladj))
educ_permanova$pval_cat[which(educ_permanova$pvaladj<0.001)] = 1 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.001 & educ_permanova$pvaladj<0.01)] = 2 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.01 & educ_permanova$pvaladj<0.05)] = 3 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.05 & educ_permanova$pvaladj<0.1)] = 4 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.1 & educ_permanova$pvaladj<0.2)] = 5 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.2)] = 6
educ_permanova$pval_cat = factor(educ_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

educ_permanova$temp=rep(3,length(newcovs))
educ_permanova <- educ_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_educ <- ggplot(educ_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_educ

# PublicIns ####

TAXO_ins=TAXO[complete.cases(all_data_complete[,c("PublicIns")]),]
all_ins=all_data_complete[complete.cases(all_data_complete[,c("PublicIns")]),]

ins=adonis2(TAXO_ins~PublicIns, data=all_ins, method="bray", perm = 10000)

newcovs=c("PublicIns")

ins1 <- ins %>% dplyr::select(5,3) 
ins1 <- ins1 [-c(2,3), ]
ins_permanova=as.data.frame(cbind(newcovs, ins1))
names(ins_permanova)=c("newcovs", "pvals","rsq")
ins_permanova$pvals=as.numeric(ins_permanova$pvals)
ins_permanova$pvaladj = p.adjust(ins_permanova$pvals,"fdr")
ins_permanova$rsq=as.numeric(ins_permanova$rsq)

ins_permanova$newcovs=factor(ins_permanova$newcovs,levels=c("PublicIns"),
                              labels=c("Public Insurance"))

knitr::kable(ins_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

ins_permanova$pval_cat = rep(0,length(ins_permanova$pvaladj))
ins_permanova$pval_cat[which(ins_permanova$pvaladj<0.001)] = 1 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.001 & ins_permanova$pvaladj<0.01)] = 2 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.01 & ins_permanova$pvaladj<0.05)] = 3 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.05 & ins_permanova$pvaladj<0.1)] = 4 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.1 & ins_permanova$pvaladj<0.2)] = 5 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.2)] = 6
ins_permanova$pval_cat = factor(ins_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

ins_permanova$temp=rep(3,length(newcovs))
ins_permanova <- ins_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_ins <- ggplot(ins_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_ins

# AbxEver ####

TAXO_abx=TAXO[complete.cases(all_data_complete[,c("AbxEver")]),]
all_abx=all_data_complete[complete.cases(all_data_complete[,c("AbxEver")]),]

abx=adonis2(TAXO_abx~AbxEver, data=all_abx, method="bray", perm = 10000)

newcovs=c("AbxEver")

abx1 <- abx %>% dplyr::select(5,3) 
abx1 <- abx1 [-c(2,3), ]
abx_permanova=as.data.frame(cbind(newcovs, abx1))
names(abx_permanova)=c("newcovs", "pvals","rsq")
abx_permanova$pvals=as.numeric(abx_permanova$pvals)
abx_permanova$pvaladj = p.adjust(abx_permanova$pvals,"fdr")
abx_permanova$rsq=as.numeric(abx_permanova$rsq)

abx_permanova$newcovs=factor(abx_permanova$newcovs,levels=c("AbxEver"),
                             labels=c("Antibiotics Ever in Pregnancy"))

knitr::kable(abx_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

abx_permanova$pval_cat = rep(0,length(abx_permanova$pvaladj))
abx_permanova$pval_cat[which(abx_permanova$pvaladj<0.001)] = 1 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.001 & abx_permanova$pvaladj<0.01)] = 2 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.01 & abx_permanova$pvaladj<0.05)] = 3 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.05 & abx_permanova$pvaladj<0.1)] = 4 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.1 & abx_permanova$pvaladj<0.2)] = 5 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.2)] = 6
abx_permanova$pval_cat = factor(abx_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

abx_permanova$temp=rep(3,length(newcovs))
abx_permanova <- abx_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_abx <- ggplot(abx_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_abx

# SexB ####

TAXO_sex=TAXO[complete.cases(all_data_complete[,c("SexB")]),]
all_sex=all_data_complete[complete.cases(all_data_complete[,c("SexB")]),]

sex=adonis2(TAXO_sex~SexB, data=all_sex, method="bray", perm = 10000)

newcovs=c("SexB")

sex1 <- sex %>% dplyr::select(5,3) 
sex1 <- sex1 [-c(2,3), ]
sex_permanova=as.data.frame(cbind(newcovs, sex1))
names(sex_permanova)=c("newcovs", "pvals","rsq")
sex_permanova$pvals=as.numeric(sex_permanova$pvals)
sex_permanova$pvaladj = p.adjust(sex_permanova$pvals,"fdr")
sex_permanova$rsq=as.numeric(sex_permanova$rsq)

sex_permanova$newcovs=factor(sex_permanova$newcovs,levels=c("SexB"),
                             labels=c("Infant Sex"))

knitr::kable(sex_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

sex_permanova$pval_cat = rep(0,length(sex_permanova$pvaladj))
sex_permanova$pval_cat[which(sex_permanova$pvaladj<0.001)] = 1 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.001 & sex_permanova$pvaladj<0.01)] = 2 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.01 & sex_permanova$pvaladj<0.05)] = 3 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.05 & sex_permanova$pvaladj<0.1)] = 4 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.1 & sex_permanova$pvaladj<0.2)] = 5 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.2)] = 6
sex_permanova$pval_cat = factor(sex_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

sex_permanova$temp=rep(3,length(newcovs))
sex_permanova <- sex_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_sex <- ggplot(sex_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_sex

# Smoking ####

TAXO_smoking=TAXO[complete.cases(all_data_complete[,c("smoking")]),]
all_smoking=all_data_complete[complete.cases(all_data_complete[,c("smoking")]),]

smoking=adonis2(TAXO_smoking~smoking, data=all_smoking, method="bray", perm = 10000)

newcovs=c("smoking")

smoking1 <- smoking %>% dplyr::select(5,3) 
smoking1 <- smoking1 [-c(2,3), ]
smoking_permanova=as.data.frame(cbind(newcovs, sex1))
names(smoking_permanova)=c("newcovs", "pvals","rsq")
smoking_permanova$pvals=as.numeric(smoking_permanova$pvals)
smoking_permanova$pvaladj = p.adjust(smoking_permanova$pvals,"fdr")
smoking_permanova$rsq=as.numeric(smoking_permanova$rsq)

smoking_permanova$newcovs=factor(smoking_permanova$newcovs,levels=c("smoking"),
                             labels=c("Smoking in Pregnancy"))

knitr::kable(smoking_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

smoking_permanova$pval_cat = rep(0,length(smoking_permanova$pvaladj))
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj<0.001)] = 1 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.001 & smoking_permanova$pvaladj<0.01)] = 2 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.01 & smoking_permanova$pvaladj<0.05)] = 3 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.05 & smoking_permanova$pvaladj<0.1)] = 4 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.1 & smoking_permanova$pvaladj<0.2)] = 5 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.2)] = 6
smoking_permanova$pval_cat = factor(smoking_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

smoking_permanova$temp=rep(3,length(newcovs))
smoking_permanova <- smoking_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_smoking <- ggplot(smoking_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_smoking

# gest_diab ####

TAXO_diab=TAXO[complete.cases(all_data_complete[,c("gest_diab")]),]
all_diab=all_data_complete[complete.cases(all_data_complete[,c("gest_diab")]),]

gest_diab=adonis2(TAXO_diab~gest_diab, data=all_diab, method="bray", perm = 10000)

newcovs=c("gest_diab")

gest_diab1 <- gest_diab %>% dplyr::select(5,3) 
gest_diab1 <- gest_diab1 [-c(2,3), ]
diab_permanova=as.data.frame(cbind(newcovs, gest_diab1))
names(diab_permanova)=c("newcovs", "pvals","rsq")
diab_permanova$pvals=as.numeric(diab_permanova$pvals)
diab_permanova$pvaladj = p.adjust(diab_permanova$pvals,"fdr")
diab_permanova$rsq=as.numeric(diab_permanova$rsq)

diab_permanova$newcovs=factor(diab_permanova$newcovs,levels=c("gest_diab"),
                                 labels=c("Gestational Diabetes"))

knitr::kable(diab_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

diab_permanova$pval_cat = rep(0,length(diab_permanova$pvaladj))
diab_permanova$pval_cat[which(diab_permanova$pvaladj<0.001)] = 1 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.001 & diab_permanova$pvaladj<0.01)] = 2 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.01 & diab_permanova$pvaladj<0.05)] = 3 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.05 & diab_permanova$pvaladj<0.1)] = 4 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.1 & diab_permanova$pvaladj<0.2)] = 5 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.2)] = 6
diab_permanova$pval_cat = factor(diab_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                    labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

diab_permanova$temp=rep(3,length(newcovs))
diab_permanova <- diab_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_diab <- ggplot(diab_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_diab

# htn ####

TAXO_htn=TAXO[complete.cases(all_data_complete[,c("htn")]),]
all_htn=all_data_complete[complete.cases(all_data_complete[,c("htn")]),]

htn=adonis2(TAXO_htn~htn, data=all_htn, method="bray", perm = 10000)

newcovs=c("htn")

htn1 <- htn %>% dplyr::select(5,3) 
htn1 <- htn1 [-c(2,3), ]
htn_permanova=as.data.frame(cbind(newcovs, htn1))
names(htn_permanova)=c("newcovs", "pvals","rsq")
htn_permanova$pvals=as.numeric(htn_permanova$pvals)
htn_permanova$pvaladj = p.adjust(htn_permanova$pvals,"fdr")
htn_permanova$rsq=as.numeric(htn_permanova$rsq)

htn_permanova$newcovs=factor(htn_permanova$newcovs,levels=c("htn"),
                              labels=c("Hypertension"))

knitr::kable(htn_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

htn_permanova$pval_cat = rep(0,length(htn_permanova$pvaladj))
htn_permanova$pval_cat[which(htn_permanova$pvaladj<0.001)] = 1 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.001 & htn_permanova$pvaladj<0.01)] = 2 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.01 & htn_permanova$pvaladj<0.05)] = 3 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.05 & htn_permanova$pvaladj<0.1)] = 4 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.1 & htn_permanova$pvaladj<0.2)] = 5 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.2)] = 6
htn_permanova$pval_cat = factor(htn_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

htn_permanova$temp=rep(3,length(newcovs))
htn_permanova <- htn_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_htn <- ggplot(htn_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_htn

# BMI_cat ####

TAXO_bmi=TAXO[complete.cases(all_data_complete[,c("BMI_cat")]),]
all_bmi=all_data_complete[complete.cases(all_data_complete[,c("BMI_cat")]),]

bmi=adonis2(TAXO_bmi~BMI_cat, data=all_bmi, method="bray", perm = 10000)

newcovs=c("BMI_cat")

bmi1 <- bmi %>% dplyr::select(5,3) 
bmi1 <- bmi1 [-c(2,3), ]
bmi_permanova=as.data.frame(cbind(newcovs, bmi1))
names(bmi_permanova)=c("newcovs", "pvals","rsq")
bmi_permanova$pvals=as.numeric(bmi_permanova$pvals)
bmi_permanova$pvaladj = p.adjust(bmi_permanova$pvals,"fdr")
bmi_permanova$rsq=as.numeric(bmi_permanova$rsq)

bmi_permanova$newcovs=factor(bmi_permanova$newcovs,levels=c("BMI_cat"),
                             labels=c("BMI Category"))

knitr::kable(bmi_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

bmi_permanova$pval_cat = rep(0,length(bmi_permanova$pvaladj))
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj<0.001)] = 1 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.001 & bmi_permanova$pvaladj<0.01)] = 2 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.01 & bmi_permanova$pvaladj<0.05)] = 3 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.05 & bmi_permanova$pvaladj<0.1)] = 4 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.1 & bmi_permanova$pvaladj<0.2)] = 5 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.2)] = 6
bmi_permanova$pval_cat = factor(bmi_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

bmi_permanova$temp=rep(3,length(newcovs))
bmi_permanova <- bmi_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_bmi <- ggplot(bmi_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_bmi

# Parity_cat ####

TAXO_parity=TAXO[complete.cases(all_data_complete[,c("Parity_cat")]),]
all_parity=all_data_complete[complete.cases(all_data_complete[,c("Parity_cat")]),]

parity=adonis2(TAXO_parity~Parity_cat, data=all_parity, method="bray", perm = 10000)

newcovs=c("Parity_cat")

parity1 <- parity %>% dplyr::select(5,3) 
parity1 <- parity1 [-c(2,3), ]
parity_permanova=as.data.frame(cbind(newcovs, parity1))
names(parity_permanova)=c("newcovs", "pvals","rsq")
parity_permanova$pvals=as.numeric(parity_permanova$pvals)
parity_permanova$pvaladj = p.adjust(parity_permanova$pvals,"fdr")
parity_permanova$rsq=as.numeric(parity_permanova$rsq)

parity_permanova$newcovs=factor(parity_permanova$newcovs,levels=c("Parity_cat"),
                             labels=c("Parity Category"))

knitr::kable(parity_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

parity_permanova$pval_cat = rep(0,length(parity_permanova$pvaladj))
parity_permanova$pval_cat[which(parity_permanova$pvaladj<0.001)] = 1 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.001 & parity_permanova$pvaladj<0.01)] = 2 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.01 & parity_permanova$pvaladj<0.05)] = 3 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.05 & parity_permanova$pvaladj<0.1)] = 4 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.1 & parity_permanova$pvaladj<0.2)] = 5 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.2)] = 6
parity_permanova$pval_cat = factor(parity_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

parity_permanova$temp=rep(3,length(newcovs))
parity_permanova <- parity_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_parity <- ggplot(parity_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_parity

# Single/Multiple Factor Model Plot ####

sfm_data=bind_rows(white_permanova,hisp_permanova,age_permanova,
                    educ_permanova, ins_permanova, abx_permanova,
                    sex_permanova, smoking_permanova, diab_permanova,
                    htn_permanova, bmi_permanova, parity_permanova)

sfm_plot <- ggplot(sfm_data, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)")
sfm_plot

ggsave("/...Single Factor Model Plot.png") 

full_data=bind_rows(sfm_data,adj_permanova)

full_plot <- ggplot(full_data, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red3", "red", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)")
full_plot

ggsave("/...Full Model Plot_option1.png") 

full_plot2 <- ggarrange(sfm_plot,p_adj,
                       ncol = 2, nrow = 1, common.legend=FALSE,
                       align = c("hv"),
                       font.label = list(size = 13)
                       # legend="bottom"
)
full_plot2

ggsave("/...Full Model Plot_option2.png",
       width=9, height=6) 

# MARCH ####
all_MARCH=all_MARCH %>% mutate(white=ifelse(race_recode=="White",1,
                                            ifelse(is.na(race_recode),NA,0)))

# Remove columns/rows that have all NA values
all_MARCH <- all_MARCH[,colSums(is.na(all_MARCH))<nrow(all_MARCH)]
all_MARCH <- all_MARCH[,rowSums(is.na(all_MARCH))<ncol(all_MARCH)]
dim(all_MARCH)

# TAXO contains only the taxonomy data
TAXO_MARCH = all_MARCH[,2:5228]
dim(TAXO_MARCH)
# Remove columns/rows that have all NA values
TAXO_MARCH <- TAXO_MARCH[,colSums(is.na(TAXO_MARCH))<nrow(TAXO_MARCH)]
TAXO_MARCH <- TAXO_MARCH[,rowSums(is.na(TAXO_MARCH))<ncol(TAXO_MARCH)]
dim(TAXO_MARCH)
# # Remove obs 73, 172, and 217 because they have NAs
# TAXO_MARCH <- TAXO_MARCH[-c(73, 172, 217), ]
# Check no NAs
is.na(TAXO_MARCH)
apply(is.na(TAXO_MARCH), 2, which)

#Verify no taxonomy with ALL zero reads
allzero_MARCH = 0
# Find which columns have all zeros
for(i in 1:(ncol(TAXO_MARCH))){ 
  if(isTRUE(sum(TAXO_MARCH[,i])==0)){allzero_MARCH==c(allzero_MARCH,1)}
  else{allzero_MARCH=c(allzero_MARCH, 0)}
}
# Remove the 0 at the start of allzero vector
allzero_MARCH = allzero_MARCH[2:length(allzero_MARCH)]
allzero_MARCH

# Clean up dataset with only vars needed
all_MARCH <- all_MARCH %>%
  dplyr::select(1:5228,5234:5237,5240,5245,5247,5249,5250,5252,5253,5263)

for(i in 1:length(covs)){
  
  all_MARCH %>% tabyl(covs[i])
  
  TAXO_temp_MARCH=TAXO_MARCH[complete.cases(all_MARCH[,covs[i]]),]
  alldata_temp_MARCH=all_MARCH[complete.cases(all_MARCH[,covs[i]]),]
  modelfit=adonis2(TAXO_temp_MARCH~alldata_temp_MARCH[,covs[i]], data=alldata_temp_MARCH, method="bray", perm = 10000)
  pvals[i] = modelfit$`Pr(>F)`[1]
  rsq[i] = modelfit$R2[1]
}

# Unadjusted model ####

permanova_MARCH=as.data.frame(cbind(covs, pvals, rsq))
permanova_MARCH$pvaladj = p.adjust(permanova_MARCH$pvals,"fdr")
permanova_MARCH$pval_cat = rep(0,length(permanova_MARCH$pvaladj))
permanova_MARCH$pval_cat[which(permanova_MARCH$pvaladj<0.001)] = 1 
permanova_MARCH$pval_cat[which(permanova_MARCH$pvaladj>=0.001 & permanova_MARCH$pvaladj<0.01)] = 2 
permanova_MARCH$pval_cat[which(permanova_MARCH$pvaladj>=0.01 & permanova_MARCH$pvaladj<0.05)] = 3 
permanova_MARCH$pval_cat[which(permanova_MARCH$pvaladj>=0.05 & permanova_MARCH$pvaladj<0.1)] = 4 
permanova_MARCH$pval_cat[which(permanova_MARCH$pvaladj>=0.1 & permanova_MARCH$pvaladj<0.2)] = 5 
permanova_MARCH$pval_cat[which(permanova_MARCH$pvaladj>=0.2)] = 6
permanova_MARCH$pval_cat = factor(permanova_MARCH$pval_cat,levels=c(1,2,3,4,5,6),
                                  labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

permanova_MARCH$covs = factor(permanova_MARCH$covs,
                              levels=c( "white", "hispanic", "Momage", "educ", "PublicIns",
                                        "AbxEver", "SexB", "smoking", "gest_diab",
                                        "htn", "BMI_cat", "Parity_cat"),
                              labels=c("Self-reported White Race", "Self-reported Hispanic Race", "Participant Age", "Education Level", "Public Insurance",
                                       "Antibiotics Ever in Pregnancy",
                                       "Infant Sex", "Smoking in Pregnancy",
                                       "Gestational Diabetes", "Hypertension",
                                       "BMI Category", "Parity Category"))
permanova_MARCH$temp=rep(1, dim(permanova_MARCH)[1])
permanova_MARCH$rsq=as.numeric(permanova_MARCH$rsq)
permanova_MARCH$pvals=as.numeric(permanova_MARCH$pvals)

knitr::kable(permanova_MARCH[,c(1:4)],digits=4, caption="Single Factor Model Results - MARCH", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results_MARCH.png") 

permanova_MARCH$temp <- "Single Factor Model"
permanova_MARCH$temp <- as.factor(permanova_MARCH$temp)

p_unadj_MARCH <- ggplot(permanova_MARCH, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "coral1", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)")+
  guides(size = guide_legend(order = 2),colour = guide_legend(order = 1))
p_unadj_MARCH

# Adjusted model ####

# Keep variables with FDR-adjusted p-values <0.2.

TAXO_sub_2_MARCH=TAXO_MARCH[complete.cases(all_MARCH[,c("Momage", "educ","PublicIns", "AbxEver","Parity_cat")]),]
all_sub_2_MARCH=all_MARCH[complete.cases(all_MARCH[,c("Momage", "educ","PublicIns", "AbxEver","Parity_cat")]),]

fullmodel_MARCH=adonis2(TAXO_sub_2_MARCH~Momage+educ+PublicIns+AbxEver+Parity_cat, data=all_sub_2_MARCH, method="bray", perm = 10000)

newcovs_MARCH=c("Momage", "educ","PublicIns", "AbxEver","Parity_cat")

fullmodel1_MARCH <- fullmodel_MARCH %>% dplyr::select(5,3) 
fullmodel1_MARCH <- fullmodel1_MARCH [-c(6,7), ]
adj_permanova_MARCH=as.data.frame(cbind(newcovs_MARCH, fullmodel1_MARCH))
names(adj_permanova_MARCH)=c("newcovs", "pvals","rsq")
adj_permanova_MARCH$pvals=as.numeric(adj_permanova_MARCH$pvals)
adj_permanova_MARCH$pvaladj = p.adjust(adj_permanova_MARCH$pvals,"fdr")
adj_permanova_MARCH$rsq=as.numeric(adj_permanova_MARCH$rsq)

adj_permanova_MARCH$newcovs=factor(adj_permanova_MARCH$newcovs,levels=c("Momage", "educ","PublicIns", "AbxEver","Parity_cat"),
                                   labels=c("Participant Age", "Education Level", "Public Insurance",
                                            "Antibiotics Ever in Pregnancy",
                                            "Parity Category"))

knitr::kable(adj_permanova_MARCH[,c(1:4)], digits=4, caption="Adjusted Model - MARCH", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Adjusted Model_MARCH.png") 

adj_permanova_MARCH$pval_cat = rep(0,length(adj_permanova_MARCH$pvaladj))
adj_permanova_MARCH$pval_cat[which(adj_permanova_MARCH$pvaladj<0.001)] = 1 
adj_permanova_MARCH$pval_cat[which(adj_permanova_MARCH$pvaladj>=0.001 & adj_permanova_MARCH$pvaladj<0.01)] = 2 
adj_permanova_MARCH$pval_cat[which(adj_permanova_MARCH$pvaladj>=0.01 & adj_permanova_MARCH$pvaladj<0.05)] = 3 
adj_permanova_MARCH$pval_cat[which(adj_permanova_MARCH$pvaladj>=0.05 & adj_permanova_MARCH$pvaladj<0.1)] = 4 
adj_permanova_MARCH$pval_cat[which(adj_permanova_MARCH$pvaladj>=0.1 & adj_permanova_MARCH$pvaladj<0.2)] = 5 
adj_permanova_MARCH$pval_cat[which(adj_permanova_MARCH$pvaladj>=0.2)] = 6
adj_permanova_MARCH$pval_cat = factor(adj_permanova_MARCH$pval_cat,levels=c(1,2,3,4,5,6),
                                  labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))
adj_permanova_MARCH$temp=rep(2,length(adj_permanova_MARCH$newcovs))

adj_permanova_MARCH <- adj_permanova_MARCH %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Adjusted Model", NA),
    temp=as.factor(temp)
  )
adj_permanova_MARCH$temp <- "Adjusted Model"
adj_permanova_MARCH$temp <- as.factor(adj_permanova_MARCH$temp)

p_adj_MARCH <- ggplot(adj_permanova_MARCH, aes(x=temp,y =covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red3", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value",size="Variance Explained (%)") +
  guides(size = guide_legend(order = 2),colour = guide_legend(order = 1))
p_adj_MARCH

permanova_MARCH$temp <- as.factor(permanova_MARCH$temp)
adj_permanova_MARCH$temp <- as.factor(adj_permanova_MARCH$temp)
full_data_MARCH=bind_rows(permanova_MARCH,adj_permanova_MARCH)

full_plot_MARCH <- ggplot(full_data_MARCH, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "coral1", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)")
full_plot_MARCH

full_plot_MARCH <- ggarrange(p_unadj_MARCH, p_adj_MARCH,
                          ncol = 2, nrow = 1, common.legend=FALSE
                          # legend="bottom"
                          )
full_plot_MARCH

ggsave("/...full_plot_MARCH.png", plot = full_plot_MARCH)

# White ####

TAXO_white=TAXO_MARCH[complete.cases(all_MARCH[,c("white")]),]
all_white=all_MARCH[complete.cases(all_MARCH[,c("white")]),]

white=adonis2(TAXO_white~white, data=all_white, method="bray", perm = 10000)

newcovs=c("white")

white1 <- white %>% dplyr::select(5,3) 
white1 <- white1 [-c(2,3), ]
white_permanova=as.data.frame(cbind(newcovs, white1))
names(white_permanova)=c("newcovs", "pvals","rsq")
white_permanova$pvals=as.numeric(white_permanova$pvals)
white_permanova$pvaladj = p.adjust(white_permanova$pvals,"fdr")
white_permanova$rsq=as.numeric(white_permanova$rsq)

white_permanova$newcovs=factor(white_permanova$newcovs,levels=c("white"),
                               labels=c("Self-reported White Race"))

knitr::kable(white_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

white_permanova$pval_cat = rep(0,length(white_permanova$pvaladj))
white_permanova$pval_cat[which(white_permanova$pvaladj<0.001)] = 1 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.001 & white_permanova$pvaladj<0.01)] = 2 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.01 & white_permanova$pvaladj<0.05)] = 3 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.05 & white_permanova$pvaladj<0.1)] = 4 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.1 & white_permanova$pvaladj<0.2)] = 5 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.2)] = 6
white_permanova$pval_cat = factor(white_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                  labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

white_permanova$temp=rep(3,length(newcovs))
white_permanova <- white_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_white <- ggplot(white_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_white

# Hispanic ####

TAXO_hisp=TAXO_MARCH[complete.cases(all_MARCH[,c("hispanic")]),]
all_hisp=all_MARCH[complete.cases(all_MARCH[,c("hispanic")]),]

hisp=adonis2(TAXO_hisp~hispanic, data=all_hisp, method="bray", perm = 10000)

newcovs=c("hispanic")

hisp1 <- hisp %>% dplyr::select(5,3) 
hisp1 <- hisp1 [-c(2,3), ]
hisp_permanova=as.data.frame(cbind(newcovs, hisp1))
names(hisp_permanova)=c("newcovs", "pvals","rsq")
hisp_permanova$pvals=as.numeric(hisp_permanova$pvals)
hisp_permanova$pvaladj = p.adjust(hisp_permanova$pvals,"fdr")
hisp_permanova$rsq=as.numeric(hisp_permanova$rsq)

hisp_permanova$newcovs=factor(hisp_permanova$newcovs,levels=c("hispanic"),
                              labels=c("Self-reported Hispanic Race"))

knitr::kable(hisp_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

hisp_permanova$pval_cat = rep(0,length(hisp_permanova$pvaladj))
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj<0.001)] = 1 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.001 & hisp_permanova$pvaladj<0.01)] = 2 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.01 & hisp_permanova$pvaladj<0.05)] = 3 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.05 & hisp_permanova$pvaladj<0.1)] = 4 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.1 & hisp_permanova$pvaladj<0.2)] = 5 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.2)] = 6
hisp_permanova$pval_cat = factor(hisp_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

hisp_permanova$temp=rep(3,length(newcovs))
hisp_permanova <- hisp_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_hisp <- ggplot(hisp_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_hisp

# Momage ####

TAXO_age=TAXO_MARCH[complete.cases(all_MARCH[,c("Momage")]),]
all_age=all_MARCH[complete.cases(all_MARCH[,c("Momage")]),]

age=adonis2(TAXO_age~Momage, data=all_age, method="bray", perm = 10000)

newcovs=c("Momage")

age1 <- age %>% dplyr::select(5,3) 
age1 <- age1 [-c(2,3), ]
age_permanova=as.data.frame(cbind(newcovs, age1))
names(age_permanova)=c("newcovs", "pvals","rsq")
age_permanova$pvals=as.numeric(age_permanova$pvals)
age_permanova$pvaladj = p.adjust(age_permanova$pvals,"fdr")
age_permanova$rsq=as.numeric(age_permanova$rsq)

age_permanova$newcovs=factor(age_permanova$newcovs,levels=c("Momage"),
                             labels=c("Participant Age"))

knitr::kable(age_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

age_permanova$pval_cat = rep(0,length(age_permanova$pvaladj))
age_permanova$pval_cat[which(age_permanova$pvaladj<0.001)] = 1 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.001 & age_permanova$pvaladj<0.01)] = 2 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.01 & age_permanova$pvaladj<0.05)] = 3 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.05 & age_permanova$pvaladj<0.1)] = 4 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.1 & age_permanova$pvaladj<0.2)] = 5 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.2)] = 6
age_permanova$pval_cat = factor(age_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

age_permanova$temp=rep(3,length(newcovs))
age_permanova <- age_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_age <- ggplot(age_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_age

# educ ####

TAXO_educ=TAXO_MARCH[complete.cases(all_MARCH[,c("educ")]),]
all_educ=all_MARCH[complete.cases(all_MARCH[,c("educ")]),]

educ=adonis2(TAXO_educ~educ, data=all_educ, method="bray", perm = 10000)

newcovs=c("educ")

educ1 <- educ %>% dplyr::select(5,3) 
educ1 <- educ1 [-c(2,3), ]
educ_permanova=as.data.frame(cbind(newcovs, educ1))
names(educ_permanova)=c("newcovs", "pvals","rsq")
educ_permanova$pvals=as.numeric(educ_permanova$pvals)
educ_permanova$pvaladj = p.adjust(educ_permanova$pvals,"fdr")
educ_permanova$rsq=as.numeric(educ_permanova$rsq)

educ_permanova$newcovs=factor(educ_permanova$newcovs,levels=c("educ"),
                              labels=c("Education Level"))

knitr::kable(educ_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

educ_permanova$pval_cat = rep(0,length(educ_permanova$pvaladj))
educ_permanova$pval_cat[which(educ_permanova$pvaladj<0.001)] = 1 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.001 & educ_permanova$pvaladj<0.01)] = 2 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.01 & educ_permanova$pvaladj<0.05)] = 3 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.05 & educ_permanova$pvaladj<0.1)] = 4 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.1 & educ_permanova$pvaladj<0.2)] = 5 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.2)] = 6
educ_permanova$pval_cat = factor(educ_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

educ_permanova$temp=rep(3,length(newcovs))
educ_permanova <- educ_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_educ <- ggplot(educ_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_educ

# PublicIns ####

TAXO_ins=TAXO_MARCH[complete.cases(all_MARCH[,c("PublicIns")]),]
all_ins=all_MARCH[complete.cases(all_MARCH[,c("PublicIns")]),]

ins=adonis2(TAXO_ins~PublicIns, data=all_ins, method="bray", perm = 10000)

newcovs=c("PublicIns")

ins1 <- ins %>% dplyr::select(5,3) 
ins1 <- ins1 [-c(2,3), ]
ins_permanova=as.data.frame(cbind(newcovs, ins1))
names(ins_permanova)=c("newcovs", "pvals","rsq")
ins_permanova$pvals=as.numeric(ins_permanova$pvals)
ins_permanova$pvaladj = p.adjust(ins_permanova$pvals,"fdr")
ins_permanova$rsq=as.numeric(ins_permanova$rsq)

ins_permanova$newcovs=factor(ins_permanova$newcovs,levels=c("PublicIns"),
                             labels=c("Public Insurance"))

knitr::kable(ins_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

ins_permanova$pval_cat = rep(0,length(ins_permanova$pvaladj))
ins_permanova$pval_cat[which(ins_permanova$pvaladj<0.001)] = 1 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.001 & ins_permanova$pvaladj<0.01)] = 2 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.01 & ins_permanova$pvaladj<0.05)] = 3 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.05 & ins_permanova$pvaladj<0.1)] = 4 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.1 & ins_permanova$pvaladj<0.2)] = 5 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.2)] = 6
ins_permanova$pval_cat = factor(ins_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

ins_permanova$temp=rep(3,length(newcovs))
ins_permanova <- ins_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_ins <- ggplot(ins_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_ins

# AbxEver ####

TAXO_abx=TAXO_MARCH[complete.cases(all_MARCH[,c("AbxEver")]),]
all_abx=all_MARCH[complete.cases(all_MARCH[,c("AbxEver")]),]

abx=adonis2(TAXO_abx~AbxEver, data=all_abx, method="bray", perm = 10000)

newcovs=c("AbxEver")

abx1 <- abx %>% dplyr::select(5,3) 
abx1 <- abx1 [-c(2,3), ]
abx_permanova=as.data.frame(cbind(newcovs, abx1))
names(abx_permanova)=c("newcovs", "pvals","rsq")
abx_permanova$pvals=as.numeric(abx_permanova$pvals)
abx_permanova$pvaladj = p.adjust(abx_permanova$pvals,"fdr")
abx_permanova$rsq=as.numeric(abx_permanova$rsq)

abx_permanova$newcovs=factor(abx_permanova$newcovs,levels=c("AbxEver"),
                             labels=c("Antibiotics Ever in Pregnancy"))

knitr::kable(abx_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

abx_permanova$pval_cat = rep(0,length(abx_permanova$pvaladj))
abx_permanova$pval_cat[which(abx_permanova$pvaladj<0.001)] = 1 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.001 & abx_permanova$pvaladj<0.01)] = 2 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.01 & abx_permanova$pvaladj<0.05)] = 3 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.05 & abx_permanova$pvaladj<0.1)] = 4 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.1 & abx_permanova$pvaladj<0.2)] = 5 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.2)] = 6
abx_permanova$pval_cat = factor(abx_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

abx_permanova$temp=rep(3,length(newcovs))
abx_permanova <- abx_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_abx <- ggplot(abx_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_abx

# SexB ####

TAXO_sex=TAXO_MARCH[complete.cases(all_MARCH[,c("SexB")]),]
all_sex=all_MARCH[complete.cases(all_MARCH[,c("SexB")]),]

sex=adonis2(TAXO_sex~SexB, data=all_sex, method="bray", perm = 10000)

newcovs=c("SexB")

sex1 <- sex %>% dplyr::select(5,3) 
sex1 <- sex1 [-c(2,3), ]
sex_permanova=as.data.frame(cbind(newcovs, sex1))
names(sex_permanova)=c("newcovs", "pvals","rsq")
sex_permanova$pvals=as.numeric(sex_permanova$pvals)
sex_permanova$pvaladj = p.adjust(sex_permanova$pvals,"fdr")
sex_permanova$rsq=as.numeric(sex_permanova$rsq)

sex_permanova$newcovs=factor(sex_permanova$newcovs,levels=c("SexB"),
                             labels=c("Infant Sex"))

knitr::kable(sex_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

sex_permanova$pval_cat = rep(0,length(sex_permanova$pvaladj))
sex_permanova$pval_cat[which(sex_permanova$pvaladj<0.001)] = 1 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.001 & sex_permanova$pvaladj<0.01)] = 2 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.01 & sex_permanova$pvaladj<0.05)] = 3 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.05 & sex_permanova$pvaladj<0.1)] = 4 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.1 & sex_permanova$pvaladj<0.2)] = 5 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.2)] = 6
sex_permanova$pval_cat = factor(sex_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

sex_permanova$temp=rep(3,length(newcovs))
sex_permanova <- sex_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_sex <- ggplot(sex_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_sex

# Smoking ####

TAXO_smoking=TAXO_MARCH[complete.cases(all_MARCH[,c("smoking")]),]
all_smoking=all_MARCH[complete.cases(all_MARCH[,c("smoking")]),]

smoking=adonis2(TAXO_smoking~smoking, data=all_smoking, method="bray", perm = 10000)

newcovs=c("smoking")

smoking1 <- smoking %>% dplyr::select(5,3) 
smoking1 <- smoking1 [-c(2,3), ]
smoking_permanova=as.data.frame(cbind(newcovs, sex1))
names(smoking_permanova)=c("newcovs", "pvals","rsq")
smoking_permanova$pvals=as.numeric(smoking_permanova$pvals)
smoking_permanova$pvaladj = p.adjust(smoking_permanova$pvals,"fdr")
smoking_permanova$rsq=as.numeric(smoking_permanova$rsq)

smoking_permanova$newcovs=factor(smoking_permanova$newcovs,levels=c("smoking"),
                                 labels=c("Smoking in Pregnancy"))

knitr::kable(smoking_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

smoking_permanova$pval_cat = rep(0,length(smoking_permanova$pvaladj))
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj<0.001)] = 1 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.001 & smoking_permanova$pvaladj<0.01)] = 2 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.01 & smoking_permanova$pvaladj<0.05)] = 3 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.05 & smoking_permanova$pvaladj<0.1)] = 4 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.1 & smoking_permanova$pvaladj<0.2)] = 5 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.2)] = 6
smoking_permanova$pval_cat = factor(smoking_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                    labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

smoking_permanova$temp=rep(3,length(newcovs))
smoking_permanova <- smoking_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_smoking <- ggplot(smoking_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_smoking

# gest_diab ####

TAXO_diab=TAXO_MARCH[complete.cases(all_MARCH[,c("gest_diab")]),]
all_diab=all_MARCH[complete.cases(all_MARCH[,c("gest_diab")]),]

gest_diab=adonis2(TAXO_diab~gest_diab, data=all_diab, method="bray", perm = 10000)

newcovs=c("gest_diab")

gest_diab1 <- gest_diab %>% dplyr::select(5,3) 
gest_diab1 <- gest_diab1 [-c(2,3), ]
diab_permanova=as.data.frame(cbind(newcovs, gest_diab1))
names(diab_permanova)=c("newcovs", "pvals","rsq")
diab_permanova$pvals=as.numeric(diab_permanova$pvals)
diab_permanova$pvaladj = p.adjust(diab_permanova$pvals,"fdr")
diab_permanova$rsq=as.numeric(diab_permanova$rsq)

diab_permanova$newcovs=factor(diab_permanova$newcovs,levels=c("gest_diab"),
                              labels=c("Gestational Diabetes"))

knitr::kable(diab_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

diab_permanova$pval_cat = rep(0,length(diab_permanova$pvaladj))
diab_permanova$pval_cat[which(diab_permanova$pvaladj<0.001)] = 1 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.001 & diab_permanova$pvaladj<0.01)] = 2 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.01 & diab_permanova$pvaladj<0.05)] = 3 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.05 & diab_permanova$pvaladj<0.1)] = 4 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.1 & diab_permanova$pvaladj<0.2)] = 5 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.2)] = 6
diab_permanova$pval_cat = factor(diab_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

diab_permanova$temp=rep(3,length(newcovs))
diab_permanova <- diab_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_diab <- ggplot(diab_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_diab

# htn ####

TAXO_htn=TAXO_MARCH[complete.cases(all_MARCH[,c("htn")]),]
all_htn=all_MARCH[complete.cases(all_MARCH[,c("htn")]),]

htn=adonis2(TAXO_htn~htn, data=all_htn, method="bray", perm = 10000)

newcovs=c("htn")

htn1 <- htn %>% dplyr::select(5,3) 
htn1 <- htn1 [-c(2,3), ]
htn_permanova=as.data.frame(cbind(newcovs, htn1))
names(htn_permanova)=c("newcovs", "pvals","rsq")
htn_permanova$pvals=as.numeric(htn_permanova$pvals)
htn_permanova$pvaladj = p.adjust(htn_permanova$pvals,"fdr")
htn_permanova$rsq=as.numeric(htn_permanova$rsq)

htn_permanova$newcovs=factor(htn_permanova$newcovs,levels=c("htn"),
                             labels=c("Hypertension"))

knitr::kable(htn_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

htn_permanova$pval_cat = rep(0,length(htn_permanova$pvaladj))
htn_permanova$pval_cat[which(htn_permanova$pvaladj<0.001)] = 1 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.001 & htn_permanova$pvaladj<0.01)] = 2 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.01 & htn_permanova$pvaladj<0.05)] = 3 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.05 & htn_permanova$pvaladj<0.1)] = 4 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.1 & htn_permanova$pvaladj<0.2)] = 5 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.2)] = 6
htn_permanova$pval_cat = factor(htn_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

htn_permanova$temp=rep(3,length(newcovs))
htn_permanova <- htn_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_htn <- ggplot(htn_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_htn

# BMI_cat ####

TAXO_bmi=TAXO_MARCH[complete.cases(all_MARCH[,c("BMI_cat")]),]
all_bmi=all_MARCH[complete.cases(all_MARCH[,c("BMI_cat")]),]

bmi=adonis2(TAXO_bmi~BMI_cat, data=all_bmi, method="bray", perm = 10000)

newcovs=c("BMI_cat")

bmi1 <- bmi %>% dplyr::select(5,3) 
bmi1 <- bmi1 [-c(2,3), ]
bmi_permanova=as.data.frame(cbind(newcovs, bmi1))
names(bmi_permanova)=c("newcovs", "pvals","rsq")
bmi_permanova$pvals=as.numeric(bmi_permanova$pvals)
bmi_permanova$pvaladj = p.adjust(bmi_permanova$pvals,"fdr")
bmi_permanova$rsq=as.numeric(bmi_permanova$rsq)

bmi_permanova$newcovs=factor(bmi_permanova$newcovs,levels=c("BMI_cat"),
                             labels=c("BMI Category"))

knitr::kable(bmi_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

bmi_permanova$pval_cat = rep(0,length(bmi_permanova$pvaladj))
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj<0.001)] = 1 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.001 & bmi_permanova$pvaladj<0.01)] = 2 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.01 & bmi_permanova$pvaladj<0.05)] = 3 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.05 & bmi_permanova$pvaladj<0.1)] = 4 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.1 & bmi_permanova$pvaladj<0.2)] = 5 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.2)] = 6
bmi_permanova$pval_cat = factor(bmi_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

bmi_permanova$temp=rep(3,length(newcovs))
bmi_permanova <- bmi_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_bmi <- ggplot(bmi_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_bmi

# Parity_cat ####

TAXO_parity=TAXO_MARCH[complete.cases(all_MARCH[,c("Parity_cat")]),]
all_parity=all_MARCH[complete.cases(all_MARCH[,c("Parity_cat")]),]

parity=adonis2(TAXO_parity~Parity_cat, data=all_parity, method="bray", perm = 10000)

newcovs=c("Parity_cat")

parity1 <- parity %>% dplyr::select(5,3) 
parity1 <- parity1 [-c(2,3), ]
parity_permanova=as.data.frame(cbind(newcovs, parity1))
names(parity_permanova)=c("newcovs", "pvals","rsq")
parity_permanova$pvals=as.numeric(parity_permanova$pvals)
parity_permanova$pvaladj = p.adjust(parity_permanova$pvals,"fdr")
parity_permanova$rsq=as.numeric(parity_permanova$rsq)

parity_permanova$newcovs=factor(parity_permanova$newcovs,levels=c("Parity_cat"),
                                labels=c("Parity Category"))

knitr::kable(parity_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

parity_permanova$pval_cat = rep(0,length(parity_permanova$pvaladj))
parity_permanova$pval_cat[which(parity_permanova$pvaladj<0.001)] = 1 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.001 & parity_permanova$pvaladj<0.01)] = 2 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.01 & parity_permanova$pvaladj<0.05)] = 3 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.05 & parity_permanova$pvaladj<0.1)] = 4 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.1 & parity_permanova$pvaladj<0.2)] = 5 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.2)] = 6
parity_permanova$pval_cat = factor(parity_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                   labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

parity_permanova$temp=rep(3,length(newcovs))
parity_permanova <- parity_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_parity <- ggplot(parity_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_parity

# Single/Multiple Factor Model Plot ####

sfm_data=bind_rows(white_permanova,hisp_permanova,age_permanova,
                   educ_permanova, ins_permanova, abx_permanova,
                   sex_permanova, smoking_permanova, diab_permanova,
                   htn_permanova, bmi_permanova, parity_permanova)

sfm_plot <- ggplot(sfm_data, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)")
sfm_plot

ggsave("/...Single Factor Model Plot.png") 

full_data=bind_rows(sfm_data,adj_permanova_MARCH)

full_plot_MARCH <- ggplot(full_data, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red3", "red", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)")
full_plot_MARCH

ggsave("/...Full Model Plot_option1.png") 

full_plotMARCH2 <- ggarrange(sfm_plot,p_adj_MARCH,
                        ncol = 2, nrow = 1, common.legend=FALSE,
                        align = c("hv"),
                        font.label = list(size = 13)
                        # legend="bottom"
)
full_plotMARCH2

ggsave("/...Full Model Plot_option2.png",
       width=9, height=6) 

# Emory ####

all_Emory=all_Emory %>% mutate(white=ifelse(race_recode=="White",1,
                                            ifelse(is.na(race_recode),NA,0)))

# Remove columns/rows that have all NA values
all_Emory <- all_Emory[,colSums(is.na(all_Emory))<nrow(all_Emory)]
all_Emory <- all_Emory[,rowSums(is.na(all_Emory))<ncol(all_Emory)]
dim(all_Emory)

# TAXO contains only the taxonomy data
TAXO_Emory = all_Emory[,112:5338]
dim(TAXO_Emory)
# Remove columns/rows that have all NA values
TAXO_Emory <- TAXO_Emory[,colSums(is.na(TAXO_Emory))<nrow(TAXO_Emory)]
TAXO_Emory <- TAXO_Emory[,rowSums(is.na(TAXO_Emory))<ncol(TAXO_Emory)]
dim(TAXO_Emory)
# Remove obs 73, 172, and 217 because they have NAs
TAXO_Emory <- TAXO_Emory[-c(73, 172, 217), ]
# Check no NAs
is.na(TAXO_Emory)
apply(is.na(TAXO_Emory), 2, which)

#Verify no taxonomy with ALL zero reads
allzero_Emory = 0
# Find which columns have all zeros
for(i in 1:(ncol(TAXO_Emory))){ 
  if(isTRUE(sum(TAXO_Emory[,i])==0)){allzero_Emory==c(allzero_Emory,1)}
  else{allzero_Emory=c(allzero_Emory, 0)}
}
# Remove the 0 at the start of allzero vector
allzero_Emory = allzero_Emory[2:length(allzero_Emory)]
allzero_Emory

# Remove obs 73, 172, and 217 because they have NAs
all_Emory <- all_Emory[-c(73, 172, 217), ]

# Save covs
covs=c("Momage", "educ", "PublicIns", "AbxEver", "SexB", "smoking", "gest_diab","htn", "BMI_cat", "Parity_cat")
pvals=1:length(covs)
rsq=1:length(covs)

# Clean up dataset
all_Emory <- all_Emory %>%
  dplyr::select(112:5338,4:6,8,13,15,17,18,21,22)

for(i in 1:length(covs)){
  
  all_Emory %>% tabyl(covs[i])
  
  TAXO_temp_Emory=TAXO_Emory[complete.cases(all_Emory[,covs[i]]),]
  alldata_temp_Emory=all_Emory[complete.cases(all_Emory[,covs[i]]),]
  modelfit=adonis2(TAXO_temp_Emory~alldata_temp_Emory[,covs[i]], data=alldata_temp_Emory, method="bray", perm = 10000)
  pvals[i] = modelfit$`Pr(>F)`[1]
  rsq[i] = modelfit$R2[1]
}

# Unadjusted model ####

permanova_Emory=as.data.frame(cbind(covs, pvals, rsq))
permanova_Emory$pvaladj = p.adjust(permanova_Emory$pvals,"fdr")

permanova_Emory$pval_cat[which(permanova_Emory$pvaladj<0.001)] = 1 
permanova_Emory$pval_cat[which(permanova_Emory$pvaladj>=0.001 & permanova_Emory$pvaladj<0.01)] = 2 
permanova_Emory$pval_cat[which(permanova_Emory$pvaladj>=0.01 & permanova_Emory$pvaladj<0.05)] = 3 
permanova_Emory$pval_cat[which(permanova_Emory$pvaladj>=0.05 & permanova_Emory$pvaladj<0.1)] = 4 
permanova_Emory$pval_cat[which(permanova_Emory$pvaladj>=0.1 & permanova_Emory$pvaladj<0.2)] = 5 
permanova_Emory$pval_cat[which(permanova_Emory$pvaladj>=0.2)] = 6
permanova_Emory$pval_cat = factor(permanova_Emory$pval_cat,levels=c(1,2,3,4,5,6),
                                      labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

permanova_Emory$covs = factor(permanova_Emory$covs, 
                              levels=c( "Momage", "educ", "PublicIns", "AbxEver", 
                                        "SexB", "smoking", 
                                        "gest_diab","htn", "BMI_cat", "Parity_cat"),
                              labels=c("Participant Age", "Education Level", "Public Insurance",
                                       "Antibiotics Ever in Pregnancy",
                                       "Infant Sex", "Smoking in Pregnancy",
                                       "Gestational Diabetes", "Hypertension",
                                       "BMI Category", "Parity Category"))
permanova_Emory$temp=rep(1, dim(permanova_Emory)[1])
permanova_Emory$rsq=as.numeric(permanova_Emory$rsq)
permanova_Emory$pvals=as.numeric(permanova_Emory$pvals)

knitr::kable(permanova_Emory[,c(1:4)],digits=4, caption="Single Factor Model Results - Emory", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results_Emory.png") 

permanova_Emory$temp=factor(permanova_Emory$temp, levels=1, labels="Single Factor Model")

p_unadj_Emory <- ggplot(permanova_Emory, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red3", "red", "coral1", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)")
p_unadj_Emory

# Adjusted model ####

# Keep variables with FDR-adjusted p-values <0.2.

TAXO_sub_2_Emory=TAXO_Emory[complete.cases(all_Emory[,c("Momage","educ","PublicIns","AbxEver","BMI_cat","Parity_cat")]),]
all_sub_2_Emory=all_Emory[complete.cases(all_Emory[,c("Momage","educ","PublicIns","AbxEver","BMI_cat","Parity_cat")]),]

fullmodel_Emory=adonis2(TAXO_sub_2_Emory~Momage+educ+PublicIns+AbxEver+BMI_cat+Parity_cat, data=all_sub_2_Emory, method="bray", perm = 10000)

newcovs_Emory=c("Momage","educ","PublicIns","AbxEver","BMI_cat","Parity_cat")

fullmodel1_Emory <- fullmodel_Emory %>% dplyr::select(5,3) 
fullmodel1_Emory <- fullmodel1_Emory [-c(7,8), ]
adj_permanova_Emory=as.data.frame(cbind(newcovs_Emory, fullmodel1_Emory))
names(adj_permanova_Emory)=c("newcovs", "pvals","rsq")
adj_permanova_Emory$pvals=as.numeric(adj_permanova_Emory$pvals)
adj_permanova_Emory$pvaladj = p.adjust(adj_permanova_Emory$pvals,"fdr")
adj_permanova_Emory$rsq=as.numeric(adj_permanova_Emory$rsq)

adj_permanova_Emory$newcovs=factor(adj_permanova_Emory$newcovs,levels=c("Momage","educ","PublicIns","AbxEver","BMI_cat","Parity_cat"),
                                   labels=c("Participant Age", "Education Level", "Public Insurance",
                                            "Antibiotics Ever in Pregnancy", 
                                            "BMI Category","Parity Category"))

knitr::kable(adj_permanova_Emory[,c(1:4)], digits=4, caption="Adjusted Model - Emory", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Adjusted Model_Emory.png") 

adj_permanova_Emory$pval_cat = rep(0,length(adj_permanova_Emory$pvaladj))

adj_permanova_Emory$pval_cat[which(adj_permanova_Emory$pvaladj<0.001)] = 1 
adj_permanova_Emory$pval_cat[which(adj_permanova_Emory$pvaladj>=0.001 & adj_permanova_Emory$pvaladj<0.01)] = 2 
adj_permanova_Emory$pval_cat[which(adj_permanova_Emory$pvaladj>=0.01 & adj_permanova_Emory$pvaladj<0.05)] = 3 
adj_permanova_Emory$pval_cat[which(adj_permanova_Emory$pvaladj>=0.05 & adj_permanova_Emory$pvaladj<0.1)] = 4 
adj_permanova_Emory$pval_cat[which(adj_permanova_Emory$pvaladj>=0.1 & adj_permanova_Emory$pvaladj<0.2)] = 5 
adj_permanova_Emory$pval_cat[which(adj_permanova_Emory$pvaladj>=0.2)] = 6
adj_permanova_Emory$pval_cat = factor(adj_permanova_Emory$pval_cat,levels=c(1,2,3,4,5,6),
                                  labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

adj_permanova_Emory$temp=rep(3,length(adj_permanova_Emory$newcovs))

adj_permanova_Emory <- adj_permanova_Emory %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Adjusted Model", NA),
    temp=as.factor(temp)
  )

p_adj_Emory <- ggplot(adj_permanova_Emory, aes(x=temp,y =covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red", "coral1","sandybrown", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)")
p_adj_Emory

full_data_Emory=bind_rows(permanova_Emory,adj_permanova_Emory)

full_plot_Emory <- ggplot(full_data_Emory, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red3", "red", "coral1", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)")
full_plot_Emory

full_plot_Emory <- ggarrange(p_unadj_Emory, p_adj_Emory,
                             ncol = 2, nrow = 1, common.legend=FALSE
                             # legend="bottom"
)
full_plot_Emory

ggsave("/...full_plot_Emory.png", plot = full_plot_Emory)

# Momage ####

TAXO_age=TAXO_Emory[complete.cases(all_Emory[,c("Momage")]),]
all_age=all_Emory[complete.cases(all_Emory[,c("Momage")]),]

age=adonis2(TAXO_age~Momage, data=all_age, method="bray", perm = 10000)

newcovs=c("Momage")

age1 <- age %>% dplyr::select(5,3) 
age1 <- age1 [-c(2,3), ]
age_permanova=as.data.frame(cbind(newcovs, age1))
names(age_permanova)=c("newcovs", "pvals","rsq")
age_permanova$pvals=as.numeric(age_permanova$pvals)
age_permanova$pvaladj = p.adjust(age_permanova$pvals,"fdr")
age_permanova$rsq=as.numeric(age_permanova$rsq)

age_permanova$newcovs=factor(age_permanova$newcovs,levels=c("Momage"),
                             labels=c("Participant Age"))

knitr::kable(age_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

age_permanova$pval_cat = rep(0,length(age_permanova$pvaladj))
age_permanova$pval_cat[which(age_permanova$pvaladj<0.001)] = 1 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.001 & age_permanova$pvaladj<0.01)] = 2 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.01 & age_permanova$pvaladj<0.05)] = 3 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.05 & age_permanova$pvaladj<0.1)] = 4 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.1 & age_permanova$pvaladj<0.2)] = 5 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.2)] = 6
age_permanova$pval_cat = factor(age_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

age_permanova$temp=rep(3,length(newcovs))
age_permanova <- age_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_age <- ggplot(age_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_age

# educ ####

TAXO_educ=TAXO_Emory[complete.cases(all_Emory[,c("educ")]),]
all_educ=all_Emory[complete.cases(all_Emory[,c("educ")]),]

educ=adonis2(TAXO_educ~educ, data=all_educ, method="bray", perm = 10000)

newcovs=c("educ")

educ1 <- educ %>% dplyr::select(5,3) 
educ1 <- educ1 [-c(2,3), ]
educ_permanova=as.data.frame(cbind(newcovs, educ1))
names(educ_permanova)=c("newcovs", "pvals","rsq")
educ_permanova$pvals=as.numeric(educ_permanova$pvals)
educ_permanova$pvaladj = p.adjust(educ_permanova$pvals,"fdr")
educ_permanova$rsq=as.numeric(educ_permanova$rsq)

educ_permanova$newcovs=factor(educ_permanova$newcovs,levels=c("educ"),
                              labels=c("Education Level"))

knitr::kable(educ_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

educ_permanova$pval_cat = rep(0,length(educ_permanova$pvaladj))
educ_permanova$pval_cat[which(educ_permanova$pvaladj<0.001)] = 1 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.001 & educ_permanova$pvaladj<0.01)] = 2 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.01 & educ_permanova$pvaladj<0.05)] = 3 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.05 & educ_permanova$pvaladj<0.1)] = 4 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.1 & educ_permanova$pvaladj<0.2)] = 5 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.2)] = 6
educ_permanova$pval_cat = factor(educ_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

educ_permanova$temp=rep(3,length(newcovs))
educ_permanova <- educ_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_educ <- ggplot(educ_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_educ

# PublicIns ####

TAXO_ins=TAXO_Emory[complete.cases(all_Emory[,c("PublicIns")]),]
all_ins=all_Emory[complete.cases(all_Emory[,c("PublicIns")]),]

ins=adonis2(TAXO_ins~PublicIns, data=all_ins, method="bray", perm = 10000)

newcovs=c("PublicIns")

ins1 <- ins %>% dplyr::select(5,3) 
ins1 <- ins1 [-c(2,3), ]
ins_permanova=as.data.frame(cbind(newcovs, ins1))
names(ins_permanova)=c("newcovs", "pvals","rsq")
ins_permanova$pvals=as.numeric(ins_permanova$pvals)
ins_permanova$pvaladj = p.adjust(ins_permanova$pvals,"fdr")
ins_permanova$rsq=as.numeric(ins_permanova$rsq)

ins_permanova$newcovs=factor(ins_permanova$newcovs,levels=c("PublicIns"),
                             labels=c("Public Insurance"))

knitr::kable(ins_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

ins_permanova$pval_cat = rep(0,length(ins_permanova$pvaladj))
ins_permanova$pval_cat[which(ins_permanova$pvaladj<0.001)] = 1 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.001 & ins_permanova$pvaladj<0.01)] = 2 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.01 & ins_permanova$pvaladj<0.05)] = 3 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.05 & ins_permanova$pvaladj<0.1)] = 4 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.1 & ins_permanova$pvaladj<0.2)] = 5 
ins_permanova$pval_cat[which(ins_permanova$pvaladj>=0.2)] = 6
ins_permanova$pval_cat = factor(ins_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

ins_permanova$temp=rep(3,length(newcovs))
ins_permanova <- ins_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_ins <- ggplot(ins_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red3"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_ins

# AbxEver ####

TAXO_abx=TAXO_Emory[complete.cases(all_Emory[,c("AbxEver")]),]
all_abx=all_Emory[complete.cases(all_Emory[,c("AbxEver")]),]

abx=adonis2(TAXO_abx~AbxEver, data=all_abx, method="bray", perm = 10000)

newcovs=c("AbxEver")

abx1 <- abx %>% dplyr::select(5,3) 
abx1 <- abx1 [-c(2,3), ]
abx_permanova=as.data.frame(cbind(newcovs, abx1))
names(abx_permanova)=c("newcovs", "pvals","rsq")
abx_permanova$pvals=as.numeric(abx_permanova$pvals)
abx_permanova$pvaladj = p.adjust(abx_permanova$pvals,"fdr")
abx_permanova$rsq=as.numeric(abx_permanova$rsq)

abx_permanova$newcovs=factor(abx_permanova$newcovs,levels=c("AbxEver"),
                             labels=c("Antibiotics Ever in Pregnancy"))

knitr::kable(abx_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

abx_permanova$pval_cat = rep(0,length(abx_permanova$pvaladj))
abx_permanova$pval_cat[which(abx_permanova$pvaladj<0.001)] = 1 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.001 & abx_permanova$pvaladj<0.01)] = 2 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.01 & abx_permanova$pvaladj<0.05)] = 3 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.05 & abx_permanova$pvaladj<0.1)] = 4 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.1 & abx_permanova$pvaladj<0.2)] = 5 
abx_permanova$pval_cat[which(abx_permanova$pvaladj>=0.2)] = 6
abx_permanova$pval_cat = factor(abx_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

abx_permanova$temp=rep(3,length(newcovs))
abx_permanova <- abx_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_abx <- ggplot(abx_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red3"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_abx

# SexB ####

TAXO_sex=TAXO_Emory[complete.cases(all_Emory[,c("SexB")]),]
all_sex=all_Emory[complete.cases(all_Emory[,c("SexB")]),]

sex=adonis2(TAXO_sex~SexB, data=all_sex, method="bray", perm = 10000)

newcovs=c("SexB")

sex1 <- sex %>% dplyr::select(5,3) 
sex1 <- sex1 [-c(2,3), ]
sex_permanova=as.data.frame(cbind(newcovs, sex1))
names(sex_permanova)=c("newcovs", "pvals","rsq")
sex_permanova$pvals=as.numeric(sex_permanova$pvals)
sex_permanova$pvaladj = p.adjust(sex_permanova$pvals,"fdr")
sex_permanova$rsq=as.numeric(sex_permanova$rsq)

sex_permanova$newcovs=factor(sex_permanova$newcovs,levels=c("SexB"),
                             labels=c("Infant Sex"))

knitr::kable(sex_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

sex_permanova$pval_cat = rep(0,length(sex_permanova$pvaladj))
sex_permanova$pval_cat[which(sex_permanova$pvaladj<0.001)] = 1 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.001 & sex_permanova$pvaladj<0.01)] = 2 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.01 & sex_permanova$pvaladj<0.05)] = 3 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.05 & sex_permanova$pvaladj<0.1)] = 4 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.1 & sex_permanova$pvaladj<0.2)] = 5 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.2)] = 6
sex_permanova$pval_cat = factor(sex_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

sex_permanova$temp=rep(3,length(newcovs))
sex_permanova <- sex_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_sex <- ggplot(sex_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_sex

# Smoking ####

TAXO_smoking=TAXO_Emory[complete.cases(all_Emory[,c("smoking")]),]
all_smoking=all_Emory[complete.cases(all_Emory[,c("smoking")]),]

smoking=adonis2(TAXO_smoking~smoking, data=all_smoking, method="bray", perm = 10000)

newcovs=c("smoking")

smoking1 <- smoking %>% dplyr::select(5,3) 
smoking1 <- smoking1 [-c(2,3), ]
smoking_permanova=as.data.frame(cbind(newcovs, sex1))
names(smoking_permanova)=c("newcovs", "pvals","rsq")
smoking_permanova$pvals=as.numeric(smoking_permanova$pvals)
smoking_permanova$pvaladj = p.adjust(smoking_permanova$pvals,"fdr")
smoking_permanova$rsq=as.numeric(smoking_permanova$rsq)

smoking_permanova$newcovs=factor(smoking_permanova$newcovs,levels=c("smoking"),
                                 labels=c("Smoking in Pregnancy"))

knitr::kable(smoking_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

smoking_permanova$pval_cat = rep(0,length(smoking_permanova$pvaladj))
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj<0.001)] = 1 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.001 & smoking_permanova$pvaladj<0.01)] = 2 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.01 & smoking_permanova$pvaladj<0.05)] = 3 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.05 & smoking_permanova$pvaladj<0.1)] = 4 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.1 & smoking_permanova$pvaladj<0.2)] = 5 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.2)] = 6
smoking_permanova$pval_cat = factor(smoking_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                    labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

smoking_permanova$temp=rep(3,length(newcovs))
smoking_permanova <- smoking_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_smoking <- ggplot(smoking_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_smoking

# gest_diab ####

TAXO_diab=TAXO_Emory[complete.cases(all_Emory[,c("gest_diab")]),]
all_diab=all_Emory[complete.cases(all_Emory[,c("gest_diab")]),]

gest_diab=adonis2(TAXO_diab~gest_diab, data=all_diab, method="bray", perm = 10000)

newcovs=c("gest_diab")

gest_diab1 <- gest_diab %>% dplyr::select(5,3) 
gest_diab1 <- gest_diab1 [-c(2,3), ]
diab_permanova=as.data.frame(cbind(newcovs, gest_diab1))
names(diab_permanova)=c("newcovs", "pvals","rsq")
diab_permanova$pvals=as.numeric(diab_permanova$pvals)
diab_permanova$pvaladj = p.adjust(diab_permanova$pvals,"fdr")
diab_permanova$rsq=as.numeric(diab_permanova$rsq)

diab_permanova$newcovs=factor(diab_permanova$newcovs,levels=c("gest_diab"),
                              labels=c("Gestational Diabetes"))

knitr::kable(diab_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

diab_permanova$pval_cat = rep(0,length(diab_permanova$pvaladj))
diab_permanova$pval_cat[which(diab_permanova$pvaladj<0.001)] = 1 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.001 & diab_permanova$pvaladj<0.01)] = 2 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.01 & diab_permanova$pvaladj<0.05)] = 3 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.05 & diab_permanova$pvaladj<0.1)] = 4 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.1 & diab_permanova$pvaladj<0.2)] = 5 
diab_permanova$pval_cat[which(diab_permanova$pvaladj>=0.2)] = 6
diab_permanova$pval_cat = factor(diab_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

diab_permanova$temp=rep(3,length(newcovs))
diab_permanova <- diab_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_diab <- ggplot(diab_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_diab

# htn ####

TAXO_htn=TAXO_Emory[complete.cases(all_Emory[,c("htn")]),]
all_htn=all_Emory[complete.cases(all_Emory[,c("htn")]),]

htn=adonis2(TAXO_htn~htn, data=all_htn, method="bray", perm = 10000)

newcovs=c("htn")

htn1 <- htn %>% dplyr::select(5,3) 
htn1 <- htn1 [-c(2,3), ]
htn_permanova=as.data.frame(cbind(newcovs, htn1))
names(htn_permanova)=c("newcovs", "pvals","rsq")
htn_permanova$pvals=as.numeric(htn_permanova$pvals)
htn_permanova$pvaladj = p.adjust(htn_permanova$pvals,"fdr")
htn_permanova$rsq=as.numeric(htn_permanova$rsq)

htn_permanova$newcovs=factor(htn_permanova$newcovs,levels=c("htn"),
                             labels=c("Hypertension"))

knitr::kable(htn_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

htn_permanova$pval_cat = rep(0,length(htn_permanova$pvaladj))
htn_permanova$pval_cat[which(htn_permanova$pvaladj<0.001)] = 1 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.001 & htn_permanova$pvaladj<0.01)] = 2 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.01 & htn_permanova$pvaladj<0.05)] = 3 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.05 & htn_permanova$pvaladj<0.1)] = 4 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.1 & htn_permanova$pvaladj<0.2)] = 5 
htn_permanova$pval_cat[which(htn_permanova$pvaladj>=0.2)] = 6
htn_permanova$pval_cat = factor(htn_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

htn_permanova$temp=rep(3,length(newcovs))
htn_permanova <- htn_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_htn <- ggplot(htn_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_htn

# BMI_cat ####

TAXO_bmi=TAXO_Emory[complete.cases(all_Emory[,c("BMI_cat")]),]
all_bmi=all_Emory[complete.cases(all_Emory[,c("BMI_cat")]),]

bmi=adonis2(TAXO_bmi~BMI_cat, data=all_bmi, method="bray", perm = 10000)

newcovs=c("BMI_cat")

bmi1 <- bmi %>% dplyr::select(5,3) 
bmi1 <- bmi1 [-c(2,3), ]
bmi_permanova=as.data.frame(cbind(newcovs, bmi1))
names(bmi_permanova)=c("newcovs", "pvals","rsq")
bmi_permanova$pvals=as.numeric(bmi_permanova$pvals)
bmi_permanova$pvaladj = p.adjust(bmi_permanova$pvals,"fdr")
bmi_permanova$rsq=as.numeric(bmi_permanova$rsq)

bmi_permanova$newcovs=factor(bmi_permanova$newcovs,levels=c("BMI_cat"),
                             labels=c("BMI Category"))

knitr::kable(bmi_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

bmi_permanova$pval_cat = rep(0,length(bmi_permanova$pvaladj))
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj<0.001)] = 1 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.001 & bmi_permanova$pvaladj<0.01)] = 2 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.01 & bmi_permanova$pvaladj<0.05)] = 3 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.05 & bmi_permanova$pvaladj<0.1)] = 4 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.1 & bmi_permanova$pvaladj<0.2)] = 5 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.2)] = 6
bmi_permanova$pval_cat = factor(bmi_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

bmi_permanova$temp=rep(3,length(newcovs))
bmi_permanova <- bmi_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_bmi <- ggplot(bmi_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("sandybrown"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_bmi

# Parity_cat ####

TAXO_parity=TAXO_Emory[complete.cases(all_Emory[,c("Parity_cat")]),]
all_parity=all_Emory[complete.cases(all_Emory[,c("Parity_cat")]),]

parity=adonis2(TAXO_parity~Parity_cat, data=all_parity, method="bray", perm = 10000)

newcovs=c("Parity_cat")

parity1 <- parity %>% dplyr::select(5,3) 
parity1 <- parity1 [-c(2,3), ]
parity_permanova=as.data.frame(cbind(newcovs, parity1))
names(parity_permanova)=c("newcovs", "pvals","rsq")
parity_permanova$pvals=as.numeric(parity_permanova$pvals)
parity_permanova$pvaladj = p.adjust(parity_permanova$pvals,"fdr")
parity_permanova$rsq=as.numeric(parity_permanova$rsq)

parity_permanova$newcovs=factor(parity_permanova$newcovs,levels=c("Parity_cat"),
                                labels=c("Parity Category"))

knitr::kable(parity_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

parity_permanova$pval_cat = rep(0,length(parity_permanova$pvaladj))
parity_permanova$pval_cat[which(parity_permanova$pvaladj<0.001)] = 1 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.001 & parity_permanova$pvaladj<0.01)] = 2 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.01 & parity_permanova$pvaladj<0.05)] = 3 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.05 & parity_permanova$pvaladj<0.1)] = 4 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.1 & parity_permanova$pvaladj<0.2)] = 5 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.2)] = 6
parity_permanova$pval_cat = factor(parity_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                   labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

parity_permanova$temp=rep(3,length(newcovs))
parity_permanova <- parity_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_parity <- ggplot(parity_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_parity

# Single/Multiple Factor Model Plot ####

sfm_data=bind_rows(age_permanova,
                   educ_permanova, ins_permanova, abx_permanova,
                   sex_permanova, smoking_permanova, diab_permanova,
                   htn_permanova, bmi_permanova, parity_permanova)

sfm_plot <- ggplot(sfm_data, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red3", "red", "sandybrown", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)") +
  guides(size = guide_legend(order = 2),colour = guide_legend(order = 1))
sfm_plot

ggsave("/...Single Factor Model Plot.png") 

full_data=bind_rows(sfm_data,adj_permanova_Emory)

full_plot_Emory <- ggplot(full_data, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("firebrick4", "red3", "red", "coral1", "sandybrown","lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)")
full_plot_Emory

ggsave("/...Full Model Plot_option1.png") 

full_plotEmory2 <- ggarrange(sfm_plot,p_adj_Emory,
                        ncol = 2, nrow = 1, common.legend=FALSE,
                        align = c("hv"),
                        font.label = list(size = 13)
                        # legend="bottom"
)
full_plotEmory2

ggsave("/...Full Model Plot_option2.png",
       width=9, height=6) 

# All Wisconsin ####

all_Wisconsin=all_Wisconsin %>% mutate(white=ifelse(race_recode=="White",1,
                                                    ifelse(is.na(race_recode),NA,0)))

# Remove columns/rows that have all NA values
all_Wisconsin <- all_Wisconsin[,colSums(is.na(all_Wisconsin))<nrow(all_Wisconsin)]
all_Wisconsin <- all_Wisconsin[,rowSums(is.na(all_Wisconsin))<ncol(all_Wisconsin)]
dim(all_Wisconsin)

# TAXO contains only the taxonomy data
TAXO_Wisc = all_Wisconsin[,23:5249]
dim(TAXO_Wisc)
# Remove columns/rows that have all NA values
TAXO_Wisc <- TAXO_Wisc[,colSums(is.na(TAXO_Wisc))<nrow(TAXO_Wisc)]
TAXO_Wisc <- TAXO_Wisc[,rowSums(is.na(TAXO_Wisc))<ncol(TAXO_Wisc)]
dim(TAXO_Wisc)
# Check no NAs
is.na(TAXO_Wisc)
apply(is.na(TAXO_Wisc), 2, which)

#Verify no taxonomy with ALL zero reads
allzero_Wisc = 0
# Find which columns have all zeros
for(i in 1:(ncol(TAXO_Wisc))){ 
  if(isTRUE(sum(TAXO_Wisc[,i])==0)){allzero_Wisc==c(allzero_Wisc,1)}
  else{allzero_Wisc=c(allzero_Wisc, 0)}
}
# Remove the 0 at the start of allzero vector
allzero_Wisc = allzero_Wisc[2:length(allzero_Wisc)]
allzero_Wisc

#Save the columns that are covariates
covs_Wisc=c("white", "hispanic", "Momage", "educ",
            "SexB", "smoking", "BMI_cat", "Parity_cat")
pvals_Wisc=1:length(covs_Wisc)
rsq_Wisc=1:length(covs_Wisc)

#Clean up dataset
all_Wisconsin <- all_Wisconsin %>%
  dplyr::select(23:5249,3:5,5341,10,12,19:20)

for(i in 1:length(covs_Wisc)){
  
  all_Wisconsin %>% tabyl(covs_Wisc[i])
  
  TAXO_temp_Wisc=TAXO_Wisc[complete.cases(all_Wisconsin[,covs_Wisc[i]]),]
  alldata_temp_Wisc=all_Wisconsin[complete.cases(all_Wisconsin[,covs_Wisc[i]]),]
  modelfit=adonis2(TAXO_temp_Wisc~alldata_temp_Wisc[,covs_Wisc[i]], data=alldata_temp_Wisc, method="bray", perm = 10000)
  pvals_Wisc[i] = modelfit$`Pr(>F)`[1]
  rsq_Wisc[i] = modelfit$R2[1]
}

# Unadjusted model ####

permanova_Wisc=as.data.frame(cbind(covs_Wisc, pvals_Wisc, rsq_Wisc))
permanova_Wisc$pvaladj = p.adjust(permanova_Wisc$pvals_Wisc,"fdr")

permanova_Wisc$pval_cat = rep(0,length(permanova_Wisc$pvaladj))

permanova_Wisc$pval_cat[which(permanova_Wisc$pvaladj<0.001)] = 1 
permanova_Wisc$pval_cat[which(permanova_Wisc$pvaladj>=0.001 & permanova_Wisc$pvaladj<0.01)] = 2 
permanova_Wisc$pval_cat[which(permanova_Wisc$pvaladj>=0.01 & permanova_Wisc$pvaladj<0.05)] = 3 
permanova_Wisc$pval_cat[which(permanova_Wisc$pvaladj>=0.05 & permanova_Wisc$pvaladj<0.1)] = 4 
permanova_Wisc$pval_cat[which(permanova_Wisc$pvaladj>=0.1 & permanova_Wisc$pvaladj<0.2)] = 5 
permanova_Wisc$pval_cat[which(permanova_Wisc$pvaladj>=0.2)] = 6
permanova_Wisc$pval_cat = factor(permanova_Wisc$pval_cat,levels=c(1,2,3,4,5,6),
                                      labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

permanova_Wisc$covs_Wisc = factor(permanova_Wisc$covs_Wisc,
                                  levels=c("white", "hispanic", "Momage", "educ", "SexB", "smoking", 
                                           "BMI_cat", "Parity_cat"),
                                  labels=c("Self-reported White Race", "Self-reported Hispanic Race", "Participant Age", "Education Level", "Baby Sex", "Smoking in Pregnancy",
                                           "BMI Category", "Parity Category"))

permanova_Wisc$temp=rep(1, dim(permanova_Wisc)[1])
permanova_Wisc$rsq_Wisc=as.numeric(permanova_Wisc$rsq_Wisc)
permanova_Wisc$pvals_Wisc=as.numeric(permanova_Wisc$pvals_Wisc)

knitr::kable(permanova_Wisc[,c(1:4)],digits=4, caption="Single Factor Model Results - all Wisconsin", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results_all Wisconsin.png") 

permanova_Wisc$temp=factor(permanova_Wisc$temp, levels=1, labels="Single Factor Model")

p_unadj_Wisc <- ggplot(permanova_Wisc, aes(x=temp,y = covs_Wisc))+
  geom_point(aes(size =rsq_Wisc*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red", "gold"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)")
p_unadj_Wisc

# Adjusted model ####

# Keep variables with FDR-adjusted p-values <0.2.

TAXO_sub_2_Wisc=TAXO_Wisc[complete.cases(all_Wisconsin[,c("white", "Parity_cat")]),]
all_sub_2_Wisc=all_Wisconsin[complete.cases(all_Wisconsin[,c("white", "Parity_cat")]),]

fullmodel_Wisc=adonis2(TAXO_sub_2_Wisc~white+Parity_cat, data=all_sub_2_Wisc, method="bray", perm = 10000)

newcovs_Wisc=c("white", "Parity_cat")

fullmodel1_Wisc <- fullmodel_Wisc %>% dplyr::select(5,3) 
fullmodel1_Wisc <- fullmodel1_Wisc [-c(3,4), ]
adj_permanova_Wisc=as.data.frame(cbind(newcovs_Wisc, fullmodel1_Wisc))
names(adj_permanova_Wisc)=c("newcovs_Wisc", "pvals","rsq")
adj_permanova_Wisc$pvals=as.numeric(adj_permanova_Wisc$pvals)
adj_permanova_Wisc$pvaladj = p.adjust(adj_permanova_Wisc$pvals,"fdr")
adj_permanova_Wisc$rsq=as.numeric(adj_permanova_Wisc$rsq)

adj_permanova_Wisc$newcovs_Wisc=factor(adj_permanova_Wisc$newcovs_Wisc,levels=c("white", "Parity_cat"),
                               labels=c("Self-reported White Race", 
                                 "Parity Category"))
adj_permanova_Wisc$covs=adj_permanova_Wisc$newcovs_Wisc
adj_permanova_Wisc <- adj_permanova_Wisc %>%
  dplyr::select(-newcovs_Wisc)
adj_permanova_Wisc <- adj_permanova_Wisc %>%
  dplyr::select(covs,pvals,rsq,pvaladj)

knitr::kable(adj_permanova_Wisc[,c(1:4)], digits=4, caption="Adjusted Model - all Wisconsin", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Adjusted Model_all Wisconsin.png") 

adj_permanova_Wisc$pval_cat = rep(0,length(adj_permanova_Wisc$pvaladj))

adj_permanova_Wisc$pval_cat[which(adj_permanova_Wisc$pvaladj<0.001)] = 1 
adj_permanova_Wisc$pval_cat[which(adj_permanova_Wisc$pvaladj>=0.001 & adj_permanova_Wisc$pvaladj<0.01)] = 2 
adj_permanova_Wisc$pval_cat[which(adj_permanova_Wisc$pvaladj>=0.01 & adj_permanova_Wisc$pvaladj<0.05)] = 3 
adj_permanova_Wisc$pval_cat[which(adj_permanova_Wisc$pvaladj>=0.05 & adj_permanova_Wisc$pvaladj<0.1)] = 4 
adj_permanova_Wisc$pval_cat[which(adj_permanova_Wisc$pvaladj>=0.1 & adj_permanova_Wisc$pvaladj<0.2)] = 5 
adj_permanova_Wisc$pval_cat[which(adj_permanova_Wisc$pvaladj>=0.2)] = 6
adj_permanova_Wisc$pval_cat = factor(adj_permanova_Wisc$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

adj_permanova_Wisc$temp=rep(3,length(adj_permanova_Wisc$covs))

adj_permanova_Wisc <- adj_permanova_Wisc %>%
  mutate(
    temp=ifelse(temp==3, "Adjusted Model", NA),
    temp=as.factor(temp)
  )

permanova_Wisc <- permanova_Wisc %>%
  dplyr::rename(covs=covs_Wisc,
         pvals=pvals_Wisc,
         rsq=rsq_Wisc)

p_adj_Wisc <- ggplot(adj_permanova_Wisc, aes(x=temp,y =covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red3","sandybrown"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)")+
  guides(size = guide_legend(order = 2),colour = guide_legend(order = 1))
p_adj_Wisc

fulldata_Wisc=bind_rows(permanova_Wisc,adj_permanova_Wisc)

full_plot_Wisc <- ggplot(fulldata_Wisc, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red3", "red", "gold"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)")
full_plot_Wisc

full_plot_Wisc <- ggarrange(p_unadj_Wisc, p_adj_Wisc,
                             ncol = 2, nrow = 1, common.legend=FALSE,
                            align = c("hv")
                             # legend="bottom"
)
full_plot_Wisc

ggsave("/...full_plot_all Wisconsin.png", plot = full_plot_Wisc)

# White ####

TAXO_white=TAXO_Wisc[complete.cases(all_Wisconsin[,c("white")]),]
all_white=all_Wisconsin[complete.cases(all_Wisconsin[,c("white")]),]

white=adonis2(TAXO_white~white, data=all_white, method="bray", perm = 10000)

newcovs=c("white")

white1 <- white %>% dplyr::select(5,3) 
white1 <- white1 [-c(2,3), ]
white_permanova=as.data.frame(cbind(newcovs, white1))
names(white_permanova)=c("newcovs", "pvals","rsq")
white_permanova$pvals=as.numeric(white_permanova$pvals)
white_permanova$pvaladj = p.adjust(white_permanova$pvals,"fdr")
white_permanova$rsq=as.numeric(white_permanova$rsq)

white_permanova$newcovs=factor(white_permanova$newcovs,levels=c("white"),
                               labels=c("Self-reported White Race"))

knitr::kable(white_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

white_permanova$pval_cat = rep(0,length(white_permanova$pvaladj))
white_permanova$pval_cat[which(white_permanova$pvaladj<0.001)] = 1 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.001 & white_permanova$pvaladj<0.01)] = 2 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.01 & white_permanova$pvaladj<0.05)] = 3 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.05 & white_permanova$pvaladj<0.1)] = 4 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.1 & white_permanova$pvaladj<0.2)] = 5 
white_permanova$pval_cat[which(white_permanova$pvaladj>=0.2)] = 6
white_permanova$pval_cat = factor(white_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                  labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

white_permanova$temp=rep(3,length(newcovs))
white_permanova <- white_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_white <- ggplot(white_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("sandybrown"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_white

# Hispanic ####

TAXO_hisp=TAXO_Wisc[complete.cases(all_Wisconsin[,c("hispanic")]),]
all_hisp=all_Wisconsin[complete.cases(all_Wisconsin[,c("hispanic")]),]

hisp=adonis2(TAXO_hisp~hispanic, data=all_hisp, method="bray", perm = 10000)

newcovs=c("hispanic")

hisp1 <- hisp %>% dplyr::select(5,3) 
hisp1 <- hisp1 [-c(2,3), ]
hisp_permanova=as.data.frame(cbind(newcovs, hisp1))
names(hisp_permanova)=c("newcovs", "pvals","rsq")
hisp_permanova$pvals=as.numeric(hisp_permanova$pvals)
hisp_permanova$pvaladj = p.adjust(hisp_permanova$pvals,"fdr")
hisp_permanova$rsq=as.numeric(hisp_permanova$rsq)

hisp_permanova$newcovs=factor(hisp_permanova$newcovs,levels=c("hispanic"),
                              labels=c("Self-reported Hispanic Race"))

knitr::kable(hisp_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

hisp_permanova$pval_cat = rep(0,length(hisp_permanova$pvaladj))
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj<0.001)] = 1 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.001 & hisp_permanova$pvaladj<0.01)] = 2 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.01 & hisp_permanova$pvaladj<0.05)] = 3 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.05 & hisp_permanova$pvaladj<0.1)] = 4 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.1 & hisp_permanova$pvaladj<0.2)] = 5 
hisp_permanova$pval_cat[which(hisp_permanova$pvaladj>=0.2)] = 6
hisp_permanova$pval_cat = factor(hisp_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

hisp_permanova$temp=rep(3,length(newcovs))
hisp_permanova <- hisp_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_hisp <- ggplot(hisp_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_hisp

# Momage ####

TAXO_age=TAXO_Wisc[complete.cases(all_Wisconsin[,c("Momage")]),]
all_age=all_Wisconsin[complete.cases(all_Wisconsin[,c("Momage")]),]

age=adonis2(TAXO_age~Momage, data=all_age, method="bray", perm = 10000)

newcovs=c("Momage")

age1 <- age %>% dplyr::select(5,3) 
age1 <- age1 [-c(2,3), ]
age_permanova=as.data.frame(cbind(newcovs, age1))
names(age_permanova)=c("newcovs", "pvals","rsq")
age_permanova$pvals=as.numeric(age_permanova$pvals)
age_permanova$pvaladj = p.adjust(age_permanova$pvals,"fdr")
age_permanova$rsq=as.numeric(age_permanova$rsq)

age_permanova$newcovs=factor(age_permanova$newcovs,levels=c("Momage"),
                             labels=c("Participant Age"))

knitr::kable(age_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

age_permanova$pval_cat = rep(0,length(age_permanova$pvaladj))
age_permanova$pval_cat[which(age_permanova$pvaladj<0.001)] = 1 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.001 & age_permanova$pvaladj<0.01)] = 2 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.01 & age_permanova$pvaladj<0.05)] = 3 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.05 & age_permanova$pvaladj<0.1)] = 4 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.1 & age_permanova$pvaladj<0.2)] = 5 
age_permanova$pval_cat[which(age_permanova$pvaladj>=0.2)] = 6
age_permanova$pval_cat = factor(age_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

age_permanova$temp=rep(3,length(newcovs))
age_permanova <- age_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_age <- ggplot(age_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_age

# educ ####

TAXO_educ=TAXO_Wisc[complete.cases(all_Wisconsin[,c("educ")]),]
all_educ=all_Wisconsin[complete.cases(all_Wisconsin[,c("educ")]),]

educ=adonis2(TAXO_educ~educ, data=all_educ, method="bray", perm = 10000)

newcovs=c("educ")

educ1 <- educ %>% dplyr::select(5,3) 
educ1 <- educ1 [-c(2,3), ]
educ_permanova=as.data.frame(cbind(newcovs, educ1))
names(educ_permanova)=c("newcovs", "pvals","rsq")
educ_permanova$pvals=as.numeric(educ_permanova$pvals)
educ_permanova$pvaladj = p.adjust(educ_permanova$pvals,"fdr")
educ_permanova$rsq=as.numeric(educ_permanova$rsq)

educ_permanova$newcovs=factor(educ_permanova$newcovs,levels=c("educ"),
                              labels=c("Education Level"))

knitr::kable(educ_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

educ_permanova$pval_cat = rep(0,length(educ_permanova$pvaladj))
educ_permanova$pval_cat[which(educ_permanova$pvaladj<0.001)] = 1 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.001 & educ_permanova$pvaladj<0.01)] = 2 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.01 & educ_permanova$pvaladj<0.05)] = 3 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.05 & educ_permanova$pvaladj<0.1)] = 4 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.1 & educ_permanova$pvaladj<0.2)] = 5 
educ_permanova$pval_cat[which(educ_permanova$pvaladj>=0.2)] = 6
educ_permanova$pval_cat = factor(educ_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                 labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

educ_permanova$temp=rep(3,length(newcovs))
educ_permanova <- educ_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_educ <- ggplot(educ_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_educ

# SexB ####

TAXO_sex=TAXO_Wisc[complete.cases(all_Wisconsin[,c("SexB")]),]
all_sex=all_Wisconsin[complete.cases(all_Wisconsin[,c("SexB")]),]

sex=adonis2(TAXO_sex~SexB, data=all_sex, method="bray", perm = 10000)

newcovs=c("SexB")

sex1 <- sex %>% dplyr::select(5,3) 
sex1 <- sex1 [-c(2,3), ]
sex_permanova=as.data.frame(cbind(newcovs, sex1))
names(sex_permanova)=c("newcovs", "pvals","rsq")
sex_permanova$pvals=as.numeric(sex_permanova$pvals)
sex_permanova$pvaladj = p.adjust(sex_permanova$pvals,"fdr")
sex_permanova$rsq=as.numeric(sex_permanova$rsq)

sex_permanova$newcovs=factor(sex_permanova$newcovs,levels=c("SexB"),
                             labels=c("Infant Sex"))

knitr::kable(sex_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

sex_permanova$pval_cat = rep(0,length(sex_permanova$pvaladj))
sex_permanova$pval_cat[which(sex_permanova$pvaladj<0.001)] = 1 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.001 & sex_permanova$pvaladj<0.01)] = 2 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.01 & sex_permanova$pvaladj<0.05)] = 3 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.05 & sex_permanova$pvaladj<0.1)] = 4 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.1 & sex_permanova$pvaladj<0.2)] = 5 
sex_permanova$pval_cat[which(sex_permanova$pvaladj>=0.2)] = 6
sex_permanova$pval_cat = factor(sex_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

sex_permanova$temp=rep(3,length(newcovs))
sex_permanova <- sex_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_sex <- ggplot(sex_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_sex

# Smoking ####

TAXO_smoking=TAXO_Wisc[complete.cases(all_Wisconsin[,c("smoking")]),]
all_smoking=all_Wisconsin[complete.cases(all_Wisconsin[,c("smoking")]),]

smoking=adonis2(TAXO_smoking~smoking, data=all_smoking, method="bray", perm = 10000)

newcovs=c("smoking")

smoking1 <- smoking %>% dplyr::select(5,3) 
smoking1 <- smoking1 [-c(2,3), ]
smoking_permanova=as.data.frame(cbind(newcovs, sex1))
names(smoking_permanova)=c("newcovs", "pvals","rsq")
smoking_permanova$pvals=as.numeric(smoking_permanova$pvals)
smoking_permanova$pvaladj = p.adjust(smoking_permanova$pvals,"fdr")
smoking_permanova$rsq=as.numeric(smoking_permanova$rsq)

smoking_permanova$newcovs=factor(smoking_permanova$newcovs,levels=c("smoking"),
                                 labels=c("Smoking in Pregnancy"))

knitr::kable(smoking_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

smoking_permanova$pval_cat = rep(0,length(smoking_permanova$pvaladj))
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj<0.001)] = 1 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.001 & smoking_permanova$pvaladj<0.01)] = 2 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.01 & smoking_permanova$pvaladj<0.05)] = 3 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.05 & smoking_permanova$pvaladj<0.1)] = 4 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.1 & smoking_permanova$pvaladj<0.2)] = 5 
smoking_permanova$pval_cat[which(smoking_permanova$pvaladj>=0.2)] = 6
smoking_permanova$pval_cat = factor(smoking_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                    labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

smoking_permanova$temp=rep(3,length(newcovs))
smoking_permanova <- smoking_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_smoking <- ggplot(smoking_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_smoking

# BMI_cat ####

TAXO_bmi=TAXO_Wisc[complete.cases(all_Wisconsin[,c("BMI_cat")]),]
all_bmi=all_Wisconsin[complete.cases(all_Wisconsin[,c("BMI_cat")]),]

bmi=adonis2(TAXO_bmi~BMI_cat, data=all_bmi, method="bray", perm = 10000)

newcovs=c("BMI_cat")

bmi1 <- bmi %>% dplyr::select(5,3) 
bmi1 <- bmi1 [-c(2,3), ]
bmi_permanova=as.data.frame(cbind(newcovs, bmi1))
names(bmi_permanova)=c("newcovs", "pvals","rsq")
bmi_permanova$pvals=as.numeric(bmi_permanova$pvals)
bmi_permanova$pvaladj = p.adjust(bmi_permanova$pvals,"fdr")
bmi_permanova$rsq=as.numeric(bmi_permanova$rsq)

bmi_permanova$newcovs=factor(bmi_permanova$newcovs,levels=c("BMI_cat"),
                             labels=c("BMI Category"))

knitr::kable(bmi_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

bmi_permanova$pval_cat = rep(0,length(bmi_permanova$pvaladj))
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj<0.001)] = 1 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.001 & bmi_permanova$pvaladj<0.01)] = 2 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.01 & bmi_permanova$pvaladj<0.05)] = 3 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.05 & bmi_permanova$pvaladj<0.1)] = 4 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.1 & bmi_permanova$pvaladj<0.2)] = 5 
bmi_permanova$pval_cat[which(bmi_permanova$pvaladj>=0.2)] = 6
bmi_permanova$pval_cat = factor(bmi_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

bmi_permanova$temp=rep(3,length(newcovs))
bmi_permanova <- bmi_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_bmi <- ggplot(bmi_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_bmi

# Parity_cat ####

TAXO_parity=TAXO_Wisc[complete.cases(all_Wisconsin[,c("Parity_cat")]),]
all_parity=all_Wisconsin[complete.cases(all_Wisconsin[,c("Parity_cat")]),]

parity=adonis2(TAXO_parity~Parity_cat, data=all_parity, method="bray", perm = 10000)

newcovs=c("Parity_cat")

parity1 <- parity %>% dplyr::select(5,3) 
parity1 <- parity1 [-c(2,3), ]
parity_permanova=as.data.frame(cbind(newcovs, parity1))
names(parity_permanova)=c("newcovs", "pvals","rsq")
parity_permanova$pvals=as.numeric(parity_permanova$pvals)
parity_permanova$pvaladj = p.adjust(parity_permanova$pvals,"fdr")
parity_permanova$rsq=as.numeric(parity_permanova$rsq)

parity_permanova$newcovs=factor(parity_permanova$newcovs,levels=c("Parity_cat"),
                                labels=c("Parity Category"))

knitr::kable(parity_permanova[,c(1:4)], digits=4, caption="Single Factor Model", col.names = c("Covariates", "P-value", "R-squared", "FDR Adjusted P-value")) %>% kable_styling() %>%
  save_kable("/...Single Factor Model Results.png") 

parity_permanova$pval_cat = rep(0,length(parity_permanova$pvaladj))
parity_permanova$pval_cat[which(parity_permanova$pvaladj<0.001)] = 1 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.001 & parity_permanova$pvaladj<0.01)] = 2 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.01 & parity_permanova$pvaladj<0.05)] = 3 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.05 & parity_permanova$pvaladj<0.1)] = 4 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.1 & parity_permanova$pvaladj<0.2)] = 5 
parity_permanova$pval_cat[which(parity_permanova$pvaladj>=0.2)] = 6
parity_permanova$pval_cat = factor(parity_permanova$pval_cat,levels=c(1,2,3,4,5,6),
                                   labels=c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.2",">=0.2"))

parity_permanova$temp=rep(3,length(newcovs))
parity_permanova <- parity_permanova %>%
  #dplyr::select(-1) %>%
  dplyr::rename(covs=newcovs) %>%
  mutate(
    temp=ifelse(temp==3, "Single Factor Model", NA),
    temp=as.factor(temp)
  )

p_parity <- ggplot(parity_permanova, aes(x=temp,y=covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red3"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted P-value", size="Variance Explained (%)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
p_parity

# Single/Multiple Factor Model Plot ####

sfm_data=bind_rows(white_permanova, hisp_permanova, age_permanova,
                   educ_permanova, 
                   sex_permanova, smoking_permanova, 
                   bmi_permanova, parity_permanova)

sfm_plot <- ggplot(sfm_data, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red3", "sandybrown", "lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)") +
  guides(size = guide_legend(order = 2),colour = guide_legend(order = 1))
sfm_plot

ggsave("/...Single Factor Model Plot.png") 

full_data=bind_rows(sfm_data,adj_permanova_Wisc)

full_plot_Wisc <- ggplot(full_data, aes(x=temp,y = covs))+
  geom_point(aes(size =rsq*100,colour=pval_cat),alpha=0.90)+
  scale_colour_manual(values = c("red3", "sandybrown","lightgoldenrod"))+
  xlab("")+ylab("")+
  labs(colour="FDR Adjusted\np-value", size="Variance Explained (%)")
full_plot_Wisc

ggsave("/...Full Model Plot_option1.png") 

full_plotWisc2 <- ggarrange(sfm_plot,p_adj_Wisc,
                        ncol = 2, nrow = 1, common.legend=FALSE,
                        align = c("hv"),
                        font.label = list(size = 13)
                        # legend="bottom"
)
full_plotWisc2

ggsave("/...Full Model Plot_option2.png",
       width=9, height=6) 

# Full cohort panel ####
full_plot_MARCH <- full_plot_MARCH +
  guides(size = guide_legend(order = 2),colour = guide_legend(order = 1))

full_plot <- ggarrange(full_plot_MARCH, full_plot_Emory, full_plot_Wisc,
                            ncol = 3, nrow = 1, common.legend=FALSE,
                       labels = c("MARCH", "Atlanta", "WISC and MAAP"),
                       align = c("hv"),
                       font.label = list(size = 11)
                            # legend="bottom"
)
full_plot

ggsave("/...full_plot_option1.png", 
       plot = full_plot, width=20, height=10, dpi=300)

full_plot <- ggarrange(full_plotMARCH2, full_plotEmory2, full_plotWisc2,
                       ncol = 1, nrow = 3, common.legend=FALSE,
                       labels = c("MARCH", "Atlanta", "WISC and MAAP"),
                       # vjust=0.08,
                       hjust=-0.1,
                       align = c("hv"),
                       font.label = list(size = 11)
                       # legend="bottom"
)
full_plot

ggsave("/...full_plot_option2.png", 
       plot = full_plot, width=15, height=25, dpi=300)
