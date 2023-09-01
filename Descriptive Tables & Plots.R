######################################################################################
## Project: Validation of Maternal Vaginal Microbiota Signatures in Pregnancy: The ECHO Vaginal Microbiome Consortium - Dr Kimberly McKee
## Script name: Descriptive Tables & Plots.R
## Script purpose: Creating CST plots / by cohort
## Date: May 2023
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
  "gt","webshot2", "gtsummary"
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

cst <- read.csv("valenciaCST_echo.20230302.csv")
WISC <- readRDS("WISC1.rds")
MAAP <- readRDS("all_MAAP.rds")
march_emory <- readRDS("MARCH_Emory.rds")

# Data preparation -------------------------------------------------------

WISC$cohort <- "WISC"
MAAP$cohort <- "MAAP"
all_Wisconsin=full_join(WISC,MAAP)
all_data_complete=full_join(march_emory,all_Wisconsin)
all_data_complete=all_data_complete %>% mutate(white=ifelse(race_recode=="White",1,
                                                            ifelse(is.na(race_recode),NA,0)))

cst2 <- cst %>%
  dplyr::select(1,CST)

# Check if any missing data for partition
which(is.na(cst2$CST))

cst2 <- cst2 %>%
  rename(specimen=sampleID)

all_data_cst <- cst2 %>%
  right_join(all_data_complete)

all_data <- all_data_cst %>%
  dplyr::select(CST,cohort,white,hispanic,Momage,educ,PublicIns,PrivateIns,gest_diab,
                AbxEver,AbxTri1,AbxTri2,AbxTri3,SexB,smoking,htn,BMI_cat,Parity_cat) %>%
  rename(Cohort=cohort,
         White=white,
         Hispanic=hispanic,
         `Maternal Age`=Momage,
         Education=educ,
         `Public Insurance`=PublicIns,
         `Private Insurance`=PrivateIns,
         `Antibiotics Ever in Pregnancy`=AbxEver,
         `Antibiotics in First Trimester`=AbxTri1,
         `Antibiotics in Second Trimester`=AbxTri2,
         `Antibiotics in Third Trimester`=AbxTri3,
         `Birth Sex` = SexB,
         `Smoking in Pregnancy`=smoking,
         `Gestational Diabetes`=gest_diab,
         Hypertension=htn,
         `BMI Category`=BMI_cat,
         `Parity Category`=Parity_cat) %>%
  mutate(
    `Maternal Age`=as.numeric(`Maternal Age`),
    White=as.factor(White),
    Hispanic=as.factor(Hispanic),
    Hypertension=ifelse(Hypertension=="Yes",1,
                        ifelse(Hypertension=="No",0,NA)),
    Hypertension=as.factor(Hypertension),
    `Public Insurance`=ifelse(`Public Insurance`=="Yes",1,
                              ifelse(`Public Insurance`=="No",0,NA)),
    `Public Insurance`=as.factor(`Public Insurance`),
    `Private Insurance`=ifelse(`Private Insurance`=="Yes",1,
                               ifelse(`Private Insurance`=="No",0,NA)),
    `Private Insurance`=as.factor(`Private Insurance`),
    `Antibiotics Ever in Pregnancy`=ifelse(`Antibiotics Ever in Pregnancy`=="Yes",1,
                                           ifelse(`Antibiotics Ever in Pregnancy`=="No",0,NA)),
    `Antibiotics Ever in Pregnancy`=as.factor(`Antibiotics Ever in Pregnancy`),
    `Antibiotics in First Trimester`=ifelse(`Antibiotics in First Trimester`=="Yes",1,
                                            ifelse(`Antibiotics in First Trimester`=="No",0,NA)),
    `Antibiotics in First Trimester`=as.factor(`Antibiotics in First Trimester`),
    `Antibiotics in Second Trimester`=ifelse(`Antibiotics in Second Trimester`=="Yes",1,
                                             ifelse(`Antibiotics in Second Trimester`=="No",0,NA)),
    `Antibiotics in Second Trimester`=as.factor(`Antibiotics in Second Trimester`),
    `Antibiotics in Third Trimester`=ifelse(`Antibiotics in Third Trimester`=="Yes",1,
                                            ifelse(`Antibiotics in Third Trimester`=="No",0,NA)),
    `Antibiotics in Third Trimester`=as.factor(`Antibiotics in Third Trimester`),
    `Smoking in Pregnancy`=ifelse(`Smoking in Pregnancy`=="Yes",1,
                                  ifelse(`Smoking in Pregnancy`=="No",0,NA)),
    `Smoking in Pregnancy`=as.factor(`Smoking in Pregnancy`),
    `Gestational Diabetes`=ifelse(`Gestational Diabetes`=="Yes",1,
                                  ifelse(`Gestational Diabetes`=="No",0,NA)),
    `Gestational Diabetes`=as.factor(`Gestational Diabetes`)
  ) 

all_data <- mutate(all_data,
                   CST = factor(CST, levels = c("I", "II", "III", "IV-B", "IV-C", "V")))

# Contingency Table -------------------------------------------------------

all_data_clps <- all_data %>%
  mutate(
    CST=ifelse(CST=="I" | CST=="II" | CST=="V", "Non-iners Lactobacillus (I,II,V)",
               ifelse(CST=="III", "Lactobacillus iners (III)",
                      ifelse(CST=="IV" | CST=="IV-B" | CST=="IV-C", "Diverse (IV-B, IV-C)", NA)))
  )

all_data_clps$CST <- factor(all_data_clps$CST,
                            levels = c("Non-iners Lactobacillus (I,II,V)", "Lactobacillus iners (III)", "Diverse (IV-B, IV-C)"))

# Table 1
table1_clps <- all_data_clps %>%
  dplyr::select(1:18) %>%
  tbl_summary(
    # type = list(all_categorical() ~ "categorical"),
    type = list(where(is.numeric) ~ "continuous2"),
    statistic = all_continuous() ~ c("{mean} ({sd})"),
    missing_text = "(Missing)") %>%
  add_n() %>%
  modify_header(label = "**Variable**") %>% 
  bold_labels() %>%
  modify_caption("**Table 1 - Collapsed CST**") 
#as_hux_table() 
table1_clps

cohort_only_clps <- all_data_clps %>%
  dplyr::select(1:18) %>%
  tbl_summary(by = Cohort,
              # type = list(all_categorical() ~ "categorical"),
              type = list(where(is.numeric) ~ "continuous2"),
              statistic = all_continuous() ~ c("{mean} ({sd})"),
              missing_text = "(Missing)") %>%
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value=TRUE)) %>% 
  add_overall()%>%
  modify_header(label = "**Variable**") %>% 
  bold_labels() %>%
  modify_caption("**Contingency table by cohort**") 
#as_hux_table() 
cohort_only_clps

# By cohort, stratified by partition 
cohort_part_clps <-
  all_data_clps %>% 
  dplyr::select(1:18) %>%
  tbl_strata(
    strata = Cohort,
    .tbl_fun =
      ~ .x %>%
      tbl_summary(by = CST,
                  type = list(where(is.numeric) ~ "continuous2"),
                  statistic = all_continuous() ~ c("{mean} ({sd})"),
                  missing_text = "(Missing)") %>%
      #add_overall(), #not possible to get overall
      add_p(test.args = all_tests("fisher.test") ~ list(workspace=2e9)),
    .header = "**{strata}**, N = {n}") %>%
  bold_labels %>%
  modify_caption("**Contingency table by cohort, stratified by partition**")
cohort_part_clps

# By partition
part_only_clps <- all_data_clps %>%
  dplyr::select(1:18) %>%
  tbl_summary(by = CST,
              type = list(where(is.numeric) ~ "continuous2"),
              statistic = all_continuous() ~ c("{mean} ({sd})"),
              missing_text = "(Missing)") %>%
  add_p() %>%
  add_overall()%>%
  modify_header(label = "**Variable**") %>% 
  bold_labels() %>%
  modify_caption("**Contingency table by partition**")
part_only_clps

# By partition, stratified by cohort
part_cohort_clps <-
  all_data_clps %>% 
  dplyr::select(1:18) %>%
  tbl_strata(
    strata = CST,
    .tbl_fun =
      ~ .x %>%
      tbl_summary(by = Cohort,
                  type = list(where(is.numeric) ~ "continuous2"),
                  statistic = all_continuous() ~ c("{mean} ({sd})"),
                  missing_text = "(Missing)") %>%
      #add_overall(), #not possible to get overall
      add_p(test.args = all_tests("fisher.test") ~ list(workspace=2e9)),
    .header = "**{strata}**, N = {n}") %>%
  bold_labels %>%
  modify_caption("**Contingency table by partition, stratified by cohort**")
part_cohort_clps

# Set directory

as_hux_xlsx(table1_clps, file="Table1.xlsx", include = everything(), bold_header_rows = TRUE)
as_hux_xlsx(cohort_only_clps, file="Contingency table by cohort.xlsx", include = everything(), bold_header_rows = TRUE)
as_hux_xlsx(cohort_part_clps, file="Contingency table by cohort_stratified by partition.xlsx", include = everything(), bold_header_rows = TRUE)
as_hux_xlsx(part_only_clps, file="Contingency table by partition.xlsx", include = everything(), bold_header_rows = TRUE)
as_hux_xlsx(part_cohort_clps, file="Contingency table by partition_stratified by cohort.xlsx", include = everything(), bold_header_rows = TRUE)
