######################################################################################
## Project: Validation of Maternal Vaginal Microbiota Signatures in Pregnancy: The ECHO Vaginal Microbiome Consortium - Dr Kimberly McKee
## Script name: Univariate Analysis.R
## Script purpose: Conducting statistical analysis using logistic models by cohort - 5 partitions, only adjusted model PERMANOVA vars
## Date: March 2023
## Author: Beatrice Palazzolo
## Organization: Department of Family Medicine, University of Michigan
######################################################################################

# Libraries ---------------------------------------------------------------

pkgs <- c(
  "readxl", "tidyverse", "nnet", "tidyverse",
  "stargazer", "ggplot2", "janitor", "knitr", "ggplot2", 
  "readr", "printr", "tab", "lubridate", "effects", "lme4", "MASS", 
  "miscTools", "plm", "naniar", "haven", "GGally", "stringr", "gee", "gtsummary", "rlang",
  "geepack", "multgee", "kable", "kableExtra", "writexl", "elrm", "magrittr", "performance", "gt", 
  "webshot2", "nnet", "finalfit", "DescTools", "car", "broom", "generalhoslem", "gofcat"
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

all_data <- readRDS("all_data_cst.rds")

all_data <- all_data %>%
  dplyr::select(
    CST,
    specimen,
    cohort,
    white,
    Momage,
    educ,
    AbxEver,
    BMI_cat,
    Parity_cat
  ) %>%
  mutate(
    CST=ifelse(CST=="I" | CST=="II" | CST=="V", "Non-iners Lactobacillus (I,II,V)",
               ifelse(CST=="III", "Lactobacillus iners (III)",
                      ifelse(CST=="IV" | CST=="IV-B" | CST=="IV-C", "Diverse (IV)", NA))),
    CST=as.factor(CST),
    cohort=as.factor(cohort),
    white=as.factor(white),
    educ=as.factor(educ),
    Momage=as.numeric(Momage),
    AbxEver=as.factor(AbxEver),
    Parity_cat=ifelse(Parity_cat=="No prior", 0,
                      ifelse(Parity_cat=="1", "1+",
                             ifelse(Parity_cat>1, "1+", NA))),
    Parity_cat=as.factor(Parity_cat),
    BMI_cat=as.factor(BMI_cat)) %>%
  filter(!is.na(CST))

# MARCH
march <- all_data %>%
  filter(cohort=="MARCH")

# Emory
emory <- all_data %>%
  filter(cohort=="Emory")

# Wisconsin
wisc <- all_data %>%
  filter(cohort=="Wisconsin")

tab1 <- table(all_data$cohort)

# Functions ---------------------------------------------------------------

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# Statistical analysis ----------------------------------------------------

# Multinomial fixed effects model ####

# All cohorts ####
# Set directory

# White, Maternal age, Education, Antibiotics Ever in Pregnancy, Parity Category

## White 
all_data$CST <- relevel(all_data$CST, ref = "Non-iners Lactobacillus (I,II,V)")
all_data$cohort <- relevel(all_data$cohort, ref = "MARCH")
all_data$Parity_cat <- relevel(all_data$Parity_cat, ref = "1+")

race <- multinom(CST~white,
                 data=all_data,
                 na.action=na.exclude) 
summary(race)

stargazer(race, type="html", out="race.htm",
          add.lines = list(c("n", nrow(race$residuals), nrow(race$residuals))))

# Results table - odds ratios
# Set up wider table
multinom_pivot_wider <- function(x) {
  # check inputs match expectatations
  if (!inherits(x, "tbl_regression") || !inherits(x$model_obj, "multinom")) {
    stop("`x=` must be class 'tbl_regression' summary of a `nnet::multinom()` model.")
  }
  
  # create tibble of results
  df <- tibble::tibble(outcome_level = unique(x$table_body$groupname_col))
  df$tbl <- 
    purrr::map(
      df$outcome_level,
      function(lvl) {
        gtsummary::modify_table_body(
          x, 
          ~dplyr::filter(.x, .data$groupname_col %in% lvl) %>%
            dplyr::ungroup() %>%
            dplyr::select(-.data$groupname_col)
        )
      }
    )
  
  tbl_merge(df$tbl, tab_spanner = paste0("**", df$outcome_level, "**"))
}

race <- multinom(CST~white, all_data, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>%
  modify_header(label = "**Variable**") %>% 
  multinom_pivot_wider() 
race

gt::gtsave(as_gt(race), file = "race.png")

as_hux_xlsx(race, file="race.xlsx", include = everything(), bold_header_rows = TRUE)

## White & Momage
race_ma <- multinom(CST~white+Momage,
                    data=all_data,
                    na.action=na.exclude) 
summary(race_ma)

stargazer(race_ma, type="html", out="race_ma.htm",
          add.lines = list(c("n", nrow(race_ma$residuals), nrow(race_ma$residuals))))

race_ma <-
  multinom(CST~white+Momage, all_data, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>%
  modify_header(label = "**Variable**") %>% 
  multinom_pivot_wider() 
race_ma
gt::gtsave(as_gt(race_ma), file = file.path("race_ma.png"))

as_hux_xlsx(race_ma, file="race_ma.xlsx", include = everything(), bold_header_rows = TRUE)

## White & Momage & educ
race_ma_educ <- multinom(CST~white+Momage+educ,
                         data=all_data,
                         na.action=na.exclude) 
summary(race_ma_educ)

stargazer(race_ma_educ, type="html", out="race_ma_educ.htm",
          add.lines = list(c("n", nrow(race_ma_educ$residuals), nrow(race_ma_educ$residuals))))

race_ma_educ <-
  multinom(CST~white+Momage+educ, all_data, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>%
  modify_header(label = "**Variable**") %>% 
  multinom_pivot_wider() 
race_ma_educ
gt::gtsave(as_gt(race_ma_educ), file = file.path("race_ma_educ.png"))

as_hux_xlsx(race_ma_educ, file="race_ma_educ.xlsx", include = everything(), bold_header_rows = TRUE)

## White & Momage & educ & Parity_cat
race_ma_educ_p <- multinom(CST~white+Momage+educ+Parity_cat,
                         data=all_data,
                         na.action=na.exclude) 
summary(race_ma_educ_p)

stargazer(race_ma_educ_p, type="html", out="race_ma_educ_p.htm",
          add.lines = list(c("n", nrow(race_ma_educ_p$residuals), nrow(race_ma_educ_p$residuals))))

race_ma_educ_p <-
  multinom(CST~white+Momage+educ+Parity_cat, all_data, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  modify_header(label = "**Variable**") %>% 
  bold_labels() %>%
  multinom_pivot_wider() 
race_ma_educ_p
gt::gtsave(as_gt(race_ma_educ_p), file = file.path("race_ma_educ_p.png"))

as_hux_xlsx(race_ma_educ_p, file="race_ma_educ_p.xlsx", include = everything(), bold_header_rows = TRUE)

# All vars
all <- multinom(CST~white+
                  educ+
                  Momage+
                  Parity_cat+
                  AbxEver,
                data=all_data,
                na.action=na.exclude) 
summary(all)

stargazer(all, type="html", out="all.htm",
          add.lines = list(c("n", nrow(all$residuals), nrow(all$residuals))))

all <-
  multinom(CST ~ white+
             educ+
             Momage+
             Parity_cat+
             AbxEver, all_data, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>%
  modify_header(label = "**Variable**") %>% 
  multinom_pivot_wider() 
all
gt::gtsave(as_gt(all), file = file.path("all.png"))

as_hux_xlsx(all, file="all.xlsx", include = everything(), bold_header_rows = TRUE)

merge1 <-
  tbl_merge(
    tbls = list(race, race_ma, race_ma_educ, race_ma_educ_p, all),
    # tab_spanner = str_glue("**{tbls$bdy_style}**")
    tab_spanner = c("**Model 1**", "**Model 2**", "**Model 3**", "**Model 4**", "**Final Model**")
  ) 
merge1
gt::gtsave(as_gt(merge1), file = "all models.png")

as_hux_xlsx(merge1, file="all models.xlsx", include = everything(), bold_header_rows = TRUE)

# Diagnostics -------------------------------------------------------------
# Goodness-of-fit ##
all_data$CST <- relevel(all_data$CST, ref = "Non-iners Lactobacillus (I,II,V)")
all_data$cohort <- relevel(all_data$cohort, ref = "MARCH")
all_data$Parity_cat <- relevel(all_data$Parity_cat, ref = "1+")
all <- multinom(CST~white+
                  educ+
                  Momage+
                  AbxEver+
                  Parity_cat,
                data=all_data,
                na.action=na.exclude) 
hosmerlem(all, group = 10)
# As per Hosmer-Lemeshow test, model could be a good fit (though HL is not too reliable)

# MARCH ####
# Set directory

# Maternal Age, Education, Parity Category

## White 
march$CST <- relevel(march$CST, ref = "Non-iners Lactobacillus (I,II,V)")
march$cohort <- relevel(march$cohort, ref = "MARCH")
march$Parity_cat <- relevel(march$Parity_cat, ref = "1+")

## Momage
ma <- multinom(CST~Momage,
                    data=march,
                    na.action=na.exclude) 
summary(ma)

stargazer(ma, type="html", out="ma.htm",
          add.lines = list(c("n", nrow(ma$residuals), nrow(ma$residuals))))

ma <-
  multinom(CST~Momage, march, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider() 
ma
gt::gtsave(as_gt(ma), file = "ma.png")

as_hux_xlsx(ma, file="ma.xlsx", include = everything(), bold_header_rows = TRUE)

## Momage & educ
ma_educ <- multinom(CST~Momage+educ,
                         data=march,
                         na.action=na.exclude) 
summary(ma_educ)
colSums(model.matrix(ma_educ))

stargazer(ma_educ, type="html", out="ma_educ.htm",
          add.lines = list(c("n", nrow(ma_educ$residuals), nrow(ma_educ$residuals))))

ma_educ <-
  multinom(CST~Momage+educ, march, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider() 
ma_educ
gt::gtsave(as_gt(ma_educ), file = file.path("ma_educ.png"))

as_hux_xlsx(ma_educ, file="ma_educ.xlsx", include = everything(), bold_header_rows = TRUE)

## Momage & educ & Parity_cat
ma_educ_p <- multinom(CST~Momage+educ+Parity_cat,
                           data=march,
                           na.action=na.exclude) 
summary(ma_educ_p)

stargazer(ma_educ_p, type="html", out="ma_educ_p.htm",
          add.lines = list(c("n", nrow(ma_educ_p$residuals), nrow(ma_educ_p$residuals))))

ma_educ_p <-
  multinom(CST~Momage+educ+Parity_cat, march, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider() 
ma_educ_p
gt::gtsave(as_gt(ma_educ_p), file = file.path("ma_educ_p.png"))

as_hux_xlsx(ma_educ_p, file="ma_educ_p.xlsx", include = everything(), bold_header_rows = TRUE)

merge1 <-
  tbl_merge(
    tbls = list(ma, ma_educ, ma_educ_p),
    tab_spanner = c("**Model 1**", "**Model 2**", "**Final Model**")
  ) 
gt::gtsave(as_gt(merge1), file = file.path("all models.png"))

as_hux_xlsx(merge1, file="all models.xlsx", include = everything(), bold_header_rows = TRUE)

# Emory ####
# Set directory

emory$CST <- relevel(emory$CST, ref = "Non-iners Lactobacillus (I,II,V)")
emory$cohort <- relevel(emory$cohort, ref = "MARCH")
emory$Parity_cat <- relevel(emory$Parity_cat, ref = "1+")

# Maternal Age, Education, Antibiotics in Pregnancy, BMI Category, Parity Category

## Momage
ma <- multinom(CST~Momage,
               data=emory,
               na.action=na.exclude) 
summary(ma)

stargazer(ma, type="html", out="ma.htm",
          add.lines = list(c("n", nrow(ma$residuals), nrow(ma$residuals))))

ma <-
  multinom(CST~Momage, emory, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider() 
ma
gt::gtsave(as_gt(ma), file = file.path("ma.png"))

as_hux_xlsx(ma, file="ma.xlsx", include = everything(), bold_header_rows = TRUE)

## Momage & educ
ma_educ <- multinom(CST~Momage+educ,
                    data=emory,
                    na.action=na.exclude) 
summary(ma_educ)

stargazer(ma_educ, type="html", out="ma_educ.htm",
          add.lines = list(c("n", nrow(ma_educ$residuals), nrow(ma_educ$residuals))))

ma_educ <-
  multinom(CST~Momage+educ, emory, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider() 
ma_educ
gt::gtsave(as_gt(ma_educ), file = file.path("ma_educ.png"))

as_hux_xlsx(ma_educ, file="ma_educ.xlsx", include = everything(), bold_header_rows = TRUE)

## Momage & educ & AbxEver
ma_educ_a <- multinom(CST~Momage+educ+AbxEver,
                      data=emory,
                      na.action=na.exclude) 
summary(ma_educ_a)

stargazer(ma_educ_a, type="html", out="ma_educ_a.htm",
          add.lines = list(c("n", nrow(ma_educ_a$residuals), nrow(ma_educ_a$residuals))))

ma_educ_a <-
  multinom(CST~Momage+educ+AbxEver, emory, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider() 
ma_educ_a
gt::gtsave(as_gt(ma_educ_a), file = file.path("ma_educ_a.png"))

as_hux_xlsx(ma_educ_a, file="ma_educ_a.xlsx", include = everything(), bold_header_rows = TRUE)

## Momage & educ & AbxEver & BMI_cat
ma_educ_a_bmi <- multinom(CST~Momage+educ+AbxEver+BMI_cat,
                      data=emory,
                      na.action=na.exclude) 
summary(ma_educ_a_bmi)

stargazer(ma_educ_a_bmi, type="html", out="ma_educ_a_bmi.htm",
          add.lines = list(c("n", nrow(ma_educ_a_bmi$residuals), nrow(ma_educ_a_bmi$residuals))))

ma_educ_a_bmi <-
  multinom(CST~Momage+educ+AbxEver+BMI_cat, emory, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider() 
ma_educ_a_bmi
gt::gtsave(as_gt(ma_educ_a_bmi), file = file.path("ma_educ_a_bmi.png"))

as_hux_xlsx(ma_educ_a_bmi, file="ma_educ_a_bmi.xlsx", include = everything(), bold_header_rows = TRUE)

## Momage & educ & AbxEver & BMI_cat & Parity_cat
ma_educ_a_bmi_p <- multinom(CST~Momage+educ+AbxEver+BMI_cat+Parity_cat,
                      data=emory,
                      na.action=na.exclude) 
summary(ma_educ_a_bmi_p)

stargazer(ma_educ_a_bmi_p, type="html", out="ma_educ_a_bmi_p.htm",
          add.lines = list(c("n", nrow(ma_educ_a_bmi_p$residuals), nrow(ma_educ_a_bmi_p$residuals))))

ma_educ_a_bmi_p <-
  multinom(CST~Momage+educ+AbxEver+BMI_cat+Parity_cat, emory, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider() 
ma_educ_a_bmi_p
gt::gtsave(as_gt(ma_educ_a_bmi_p), file = file.path("ma_educ_a_bmi_p.png"))

as_hux_xlsx(ma_educ_a_bmi_p, file="ma_educ_a_bmi_p.xlsx", include = everything(), bold_header_rows = TRUE)

merge1 <-
  tbl_merge(
    tbls = list(ma, ma_educ, ma_educ_a, ma_educ_a_bmi, ma_educ_a_bmi_p),
    tab_spanner = c("**Model 1**", "**Model 2**", "**Model 3**", "**Model 4**", "**Final Model**")
  ) 
gt::gtsave(as_gt(merge1), file = file.path("all models.png"))

as_hux_xlsx(merge1, file="all models.xlsx", include = everything(), bold_header_rows = TRUE)

# Wisconsin ####

# Set directory

## White 
wisc$CST <- relevel(wisc$CST, ref = "Non-iners Lactobacillus (I,II,V)")
wisc$cohort <- relevel(wisc$cohort, ref = "MARCH")
wisc$Parity_cat <- relevel(wisc$Parity_cat, ref = "1+")

race <- multinom(CST~white,
                 data=wisc,
                 na.action=na.exclude) 
summary(race)

stargazer(race, type="html", out="race.htm",
          add.lines = list(c("n", nrow(race$residuals), nrow(race$residuals))))

race <-
  multinom(CST~white, wisc, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider() 
race
gt::gtsave(as_gt(race), file = file.path("race.png"))

as_hux_xlsx(race, file="race.xlsx", include = everything(), bold_header_rows = TRUE)

all <- multinom(CST~white+Parity_cat,
                data=wisc,
                na.action=na.exclude) 
summary(all)

stargazer(all, type="html", out="all.htm",
          add.lines = list(c("n", nrow(all$residuals), nrow(all$residuals))))

all <-
  multinom(CST~white+Parity_cat, wisc, na.action=na.exclude) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider() 
all
gt::gtsave(as_gt(all), file = file.path("all.png"))

as_hux_xlsx(all, file="all.xlsx", include = everything(), bold_header_rows = TRUE)

merge1 <-
  tbl_merge(
    tbls = list(race,all),
    tab_spanner = c("**Model 1**", "**Final Model**")
  ) 
gt::gtsave(as_gt(merge1), file = file.path("all models.png"))

as_hux_xlsx(merge1, file="all models.xlsx", include = everything(), bold_header_rows = TRUE)
