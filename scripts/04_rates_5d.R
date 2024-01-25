## ----load_lib, include = FALSE, warning=FALSE, message=FALSE, echo=FALSE--------------------------------------------------------------------------------------
# remotes::install_github("tlverse/sl3@devel")
# remotes::install_github("tlverse/tmle3")

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(here)
library(future)
library(gridExtra)
library(corrplot)
library(tmle3)
library(sl3)
library(gtsummary)
library(knitr)



## ----setup, include = FALSE-----------------------------------------------------------------------------------------------------------------------------------
plotFolder <- here("results","images", "04_5d")
if(!file.exists(plotFolder)) dir.create(plotFolder,recursive=TRUE)

rdataFolder <- here("data", "rdata")

knitr::opts_chunk$set(
  cache=FALSE, autodep=FALSE, warning=FALSE, message=FALSE, echo=FALSE,
  results = 'markup', dev='png', dpi=150, fig.align = "center", fig.path=paste0(plotFolder, "/"),
  cache.path=".cache/",
  duplicate.label="allow"
)

set.seed(123)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
load(here(rdataFolder, "02_2_varFusion.RData"))

PreprocessDataFolder <- here("data", "DataPreprocessing")

df <- read.csv(here(PreprocessDataFolder,"02_3_data_all_blocks.csv"), row.names = 1, header = T)  
df_imp <- read.csv(here(PreprocessDataFolder,"02_3_data_imputed_all_blocks.csv"), row.names = 1, header = T)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# discharged within the first 24 hours
df_tmp <- df[df$half_day <= 2, ]
discharge_in24_id <- df_tmp$ICUSTAY_ID[df_tmp$Y_t == 2]

df <- df[! df$ICUSTAY_ID %in% discharge_in24_id, ]
df_imp <- df_imp[! df_imp$ICUSTAY_ID %in% discharge_in24_id, ]


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
data_ids <- read.csv(here(PreprocessDataFolder,"data_ids.csv"), row.names = 1, header = T)  

data_ids <- data_ids[data_ids$ICUSTAY_ID %in% df$ICUSTAY_ID, ]
length(unique(data_ids$SUBJECT_ID))
length(unique(data_ids$ICUSTAY_ID))


## ----results='asis'-------------------------------------------------------------------------------------------------------------------------------------------
df_tmp <- df_imp %>% 
  group_by(ICUSTAY_ID) %>% 
  mutate(outcome = max(Y_t)) %>%
  dplyr::select(ICUSTAY_ID, half_day, outcome, Y_t, Y_t_plus_1, everything()) %>%
  ungroup() %>%
  mutate(outcome = as.character(outcome),
         Y_t = as.character(Y_t))


df_5d_1 <- df_tmp %>%
  dplyr::select(ICUSTAY_ID, outcome) %>%
  unique()

df_5d_2 <- df_tmp %>% 
  dplyr::select(unname(c("ICUSTAY_ID", 
                         var_name_list$num_vitals, 
                         var_name_list$hours_missing_vitals, 
                         var_name_list$num_labs))) %>%
  group_by(ICUSTAY_ID) %>% 
  summarise(across(everything(), sum))

df_5d_3 <- df_tmp %>% 
  dplyr::select(unname(c("ICUSTAY_ID", 
                         var_name_list$num_baselines, 
                         var_name_list$baselines, 
                         var_name_list$vitals, 
                         var_name_list$labs))) %>%
  group_by(ICUSTAY_ID) %>% 
  summarise(across(everything(), mean))

df_5d <- merge(df_5d_1, df_5d_2, by="ICUSTAY_ID") %>%
  merge(df_5d_3, by="ICUSTAY_ID")
rm(df_5d_1, df_5d_2, df_5d_3, df_tmp)

df_5d$outcome[df_5d$outcome == "0"] = "Longer ICU Stay"
df_5d$outcome[df_5d$outcome == "1"] = "ICU Death"
df_5d$outcome[df_5d$outcome == "2"] = "Live Discharge"


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# standardize the measurement rates for patients that died within the first 24 hours
library(lubridate)
ADMISSIONS <- read.csv(here("data","raw_data","mimic-iii-clinical-database-1.4", "ADMISSIONS.csv")) %>% 
  select(SUBJECT_ID, HADM_ID, DEATHTIME)
PATIENTS <- read.csv(here("data","raw_data","mimic-iii-clinical-database-1.4", "PATIENTS.csv")) %>% 
  select(SUBJECT_ID, DOD_SSN)
ICUSTAYS <- read.csv(here("data","raw_data","mimic-iii-clinical-database-1.4", "ICUSTAYS.csv")) %>% 
  select(SUBJECT_ID, HADM_ID, ICUSTAY_ID, INTIME, OUTTIME)

df_Y <- ADMISSIONS %>% 
  merge(PATIENTS, by = "SUBJECT_ID", all = T) %>%
  merge(ICUSTAYS, by = c("SUBJECT_ID", "HADM_ID"), all = T) %>%
  unique()

df_Y$DOD_SSN[is.na(df_Y$DOD_SSN)] = df_Y$DEATHTIME[is.na(df_Y$DOD_SSN)]
df_Y$DEATHTIME[df_Y$DEATHTIME == ""] =  df_Y$DOD_SSN[df_Y$DEATHTIME == ""]

df_Y <- df_Y %>% 
  mutate(INTIME =  ymd_hms(INTIME),
         OUTTIME = ymd_hms(OUTTIME),
         DEATHTIME = ymd_hms(DEATHTIME),
         DISCHARGE_AFTER_ICUS_IN_hours = (as.numeric(difftime(OUTTIME, INTIME))/(60*60)),
         DEATH_AFTER_ICUS_IN_hours =  (as.numeric(difftime(DEATHTIME, INTIME))/60)
         ) 

df_Y <- df_Y %>%  rowwise() %>%
  mutate(ICU_length_hours = min(c(DISCHARGE_AFTER_ICUS_IN_hours, DEATH_AFTER_ICUS_IN_hours), na.rm = T)) %>%
  select(ICUSTAY_ID, ICU_length_hours) 


df_5d <- merge(df_5d, df_Y, by = "ICUSTAY_ID")

df_5d_1 <- df_5d[df_5d$outcome != "Longer ICU Stay", ] 
df_5d_1$less_5d <- 1

df_5d_2 <- df_5d %>% filter(outcome=="Longer ICU Stay")
df_5d_2$less_5d <- 0

for (var in c(var_name_list$hours_missing_vitals, var_name_list$num_vitals, var_name_list$num_labs)) {
  df_5d_1[var] = 120 * df_5d_1[var] /  df_5d_1$ICU_length_hours
}

df_5d <- rbind(df_5d_1, df_5d_2)


## ----get_severity_score---------------------------------------------------------------------------------------------------------------------------------------
sofa <- read.csv(here("data", "SeverityScores", "sofa.csv")) 
sapsii <- read.csv(here("data", "SeverityScores", "sapsii.csv")) 

df_5d <- df_5d %>% 
  merge(sofa %>% 
          dplyr::select(icustay_id, sofa) %>% 
          rename(ICUSTAY_ID = icustay_id), 
        by="ICUSTAY_ID",
        all.x = T) %>%
  merge(sapsii %>% 
          dplyr::select(icustay_id, sapsii) %>% 
          rename(ICUSTAY_ID = icustay_id), 
        by="ICUSTAY_ID",
        all.x = T)  




## ----prepare_function-----------------------------------------------------------------------------------------------------------------------------------------
run_tmle3_outcome_rates <- function(Y, A, confounders){
  
  df <- data.frame(Y, A, confounders)
  node_list <- list(
    W = names(confounders),
    A = "A",
    Y = "Y"
  )
  
  
  processed <- process_missing(df, node_list)
  df <- processed$data
  node_list <- processed$node_list
  
  
  ## Define ATE's outcome levels
  ate_spec <- tmle_ATE(
    treatment_level = "1",
    control_level = "0"
  )
  
  
  ## Define learners
  lrn_mean <- Lrnr_mean$new()
  
  lrn_glmnet <- Lrnr_glmnet$new()
  lrn_earth <- Lrnr_earth$new()
  lrn_tree <- Lrnr_glmtree$new()
  lrn_tree2 <- Lrnr_rpart$new()
  lrn_gbm <- Lrnr_lightgbm$new()
  lrn_rf <- Lrnr_ranger$new()
  
  
  # define metalearners appropriate to data types
  ls_metalearner <- make_learner( # leaner for binary outcome
    Lrnr_nnls
  ) 
  mn_metalearner <- make_learner(  # learner for multinomial outcome
    Lrnr_solnp,
    metalearner_linear_multinomial,
    loss_loglik_multinomial
  )
  c_metalearner <- make_learner(  # learner for continuout outcome
    Lrnr_solnp,
    metalearner_linear,
    loss_squared_error
  )
  
  
  sl_Y <- Lrnr_sl$new(
    learners = list(lrn_mean, lrn_glmnet, lrn_tree2, lrn_rf),
    metalearner = c_metalearner
  )
  
  if(length(unique(A)) == 2){
    sl_A <- Lrnr_sl$new(
      learners = list(lrn_glmnet, lrn_tree2, lrn_rf),
      metalearner = ls_metalearner
    )
  } else {
    sl_A <- Lrnr_sl$new(
      learners = list(lrn_glmnet, lrn_tree2, lrn_rf),
      metalearner = mn_metalearner
    )
  }
  
  
  learner_list <- list(A = sl_A, Y = sl_Y)
  
  
  # Do it for each level in one step
  tmle_task <- ate_spec$make_tmle_task(df, node_list)
  
  initial_likelihood <- ate_spec$make_initial_likelihood(
    tmle_task,
    learner_list
  )
  
  tsm_spec <- tmle_TSM_all()
  
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
  all_tsm_params <- tsm_spec$make_params(tmle_task, targeted_likelihood)
  
  ate_params <- lapply(all_tsm_params[-1],
                       function(tsm_param){define_param(
                         Param_delta,
                         targeted_likelihood,
                         delta_param_RR,
                         list(all_tsm_params[[1]], tsm_param)
                       )})
  
  all_params <- c(all_tsm_params, ate_params)
  
  tmle_fit_multiparam_base_cfd <- fit_tmle3(
    tmle_task, targeted_likelihood, all_params,
    targeted_likelihood$updater
  )
  
  return(tmle_fit_multiparam_base_cfd)
}

make_box_plots <- function(df_box, y_vars, x){
  meltData <- melt(df_box %>% dplyr::select(y_vars, x), id.vars = x) 
  names(meltData)[1] = "x"
  if(x == 'age_cat'){
    meltData$x <- factor(meltData$x, levels = c('18-30', '31-45', '46-65', '>65')) 
  } else {
    meltData$x <- as.factor(meltData$x)
  }
  meltData$item = sub("^(n_|h_m_)", "", meltData$variable)
  pi <- ggplot(meltData, aes(x = x, y = value, color = x)) +
    geom_boxplot(position = position_dodge(width = 1)) +  
    facet_wrap(~ item, scales = 'free', nrow = 1) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(color="x", x = "", y = "count")
  return(pi)
}


make_box_plots_no_outliers <- function(df_box, y_vars, x){
  meltData <- melt(df_box %>% dplyr::select(y_vars, x), id.vars = x) 
  names(meltData)[1] = "x"
  if(x == 'age_cat'){
    meltData$x <- factor(meltData$x, levels = c('18-30', '31-45', '46-65', '>65')) 
  } else {
    meltData$x <- as.factor(meltData$x)
  }
  meltData$item = sub("^(n_|h_m_)", "", meltData$variable)
  
  # Custom calculation for the whiskers
  whisker_data <- meltData %>%
    group_by(item, x) %>%
    summarize(lower = min(value),
              upper = max(value)) %>%
    ungroup()
  
  pi <- ggplot(meltData, aes(x = x, y = value, color = x)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1)) +
    geom_linerange(data = whisker_data, aes(x = x, ymin = lower, ymax = upper), color = "black") +
    facet_wrap(~ item, scales = 'free', nrow = 1) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(color="x", x = "", y = "count")
  return(pi)
}



make_box_plots_tmle_rlts <- function(df_box, y_type){
  pi <- ggplot(df_box %>% filter(Y_type == y_type), 
               aes(x = group)) +
    geom_errorbar(aes(ymin=lower, ymax=upper, color=group), width=0.7) +
    geom_point(aes(y=tmle_est, color = group))+
    facet_wrap(~ Y_name, scales = 'free', nrow = 1) +
    labs(y="TMLE estimated rates") + # subtitle = 'TMLE estimation with other baseline variables and SOFA score as confounders'
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))  # plot.subtitle = element_text(hjust = 0.5)
  
  return(pi)
}

make_box_plots_tmle_rlts_text <- function(df_box, y_type){
  
  df_tmp <- df_box %>% 
    filter(Y_type == y_type) %>%
    mutate(lower_2dgt = round(lower, 2),
           upper_2dgt = round(upper, 2),
           tmle_est_2dgt = round(tmle_est, 2))
  
  pi <- ggplot(df_tmp, aes(x = group)) +
    geom_errorbar(aes(ymin=lower, ymax=upper, color=group), width=0.7) +
    geom_point(aes(y=tmle_est, color = group)) +
    geom_text(aes(y=upper_2dgt, label=upper_2dgt), vjust=-0.5, color="black") +
    geom_text(aes(y=lower_2dgt, label=lower_2dgt), vjust=1.5, color="black") +
    geom_text(aes(y=tmle_est_2dgt, label=tmle_est_2dgt), vjust=1.5, color="black") +
    facet_wrap(~ Y_name, scales = 'free', nrow = 1) +
    labs(y="TMLE estimated rates") +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  
  return(pi)
}




## ----age_hist-------------------------------------------------------------------------------------------------------------------------------------------------
cat("Summary of age distribution: \n")
hist(df_5d$AGE)
summary(df_5d$AGE)
age_c = cut(df_5d$AGE,
            breaks = c(18,30,45,65,101),
            labels = c("18-30","31-45","46-65",">65"))
df_5d <- df_5d %>% mutate(age_cat = age_c)
rm(age_c)

df_5d_all <- df_5d




## ----fig.height=3, fig.width=16-------------------------------------------------------------------------------------------------------------------------------
x = "age_cat"


## ----echo=FALSE, results='asis'-------------------------------------------------------------------------------------------------------------------------------
baselines = var_name_list$baselines[var_name_list$baselines != "AGE"]
y_vars <- c(var_name_list$num_vitals, var_name_list$hours_missing_vitals, var_name_list$num_labs)
idx <- 1:length(y_vars)

glmfit_0_rlt_list <- list()
glmfit_sofa_rlt_list <- list()
glmfit_saps_rlt_list <- list()
tmle_fit_base_rlt_list <- list()
tmle_fit_base_sofa_rlt_list <- list()
tmle_fit_base_saps_rlt_list <- list()

for (i in 1:length(y_vars)) {
# for (i in 16:18) {
  y = y_vars[i]
  cat(paste0("\n## ", i, ". ", y, "\n\n"))

  source(here("scripts", "04_rates_5d_template.R"))

  glmfit_0_rlt_list[[i]] <- glmfit_0_rlt
  glmfit_sofa_rlt_list[[i]] <- glmfit_sofa_rlt
  glmfit_saps_rlt_list[[i]] <- glmfit_saps_rlt
  tmle_fit_base_rlt_list[[i]] <- tmle_fit_base_rlt
  tmle_fit_base_sofa_rlt_list[[i]] <- tmle_fit_base_sofa_rlt
  tmle_fit_base_saps_rlt_list[[i]] <- tmle_fit_base_saps_rlt

  print("haha\n\n")
}

save.image(here('04_5d_age.RData'))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
df_baseline <- read.csv(here("data", "DataPreprocessing", "df_baseline.csv"))  %>% dplyr::select(-X)
df_5d <- merge(df_5d, df_baseline %>% select(ICUSTAY_ID, ETHNICITY), by="ICUSTAY_ID", all.x = T )

df_5d$ETHNICITY[is.na(df_5d$ETHNICITY)] = "Missing"

df_5d$ethnicity = factor(df_5d$ETHNICITY , levels = c("WHITE", "ASIAN", "BLACK", "HISPANIC", "Other", "Missing"))
df_5d_all$ethnicity <- df_5d$ethnicity





## ----fig.height=3, fig.width=16-------------------------------------------------------------------------------------------------------------------------------
x = "ethnicity"


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
x = "ethnicity"
baselines = var_name_list$baselines[! var_name_list$baselines %in% c("ETHNICITY_W", "ETHNICITY_B", "ETHNICITY_A", "ETHNICITY_H", "ETHNICITY_O")]
y_vars <- c(var_name_list$num_vitals, var_name_list$hours_missing_vitals, var_name_list$num_labs)
idx <- 1:length(y_vars)

glmfit_0_rlt_list_eth <- list()
glmfit_sofa_rlt_list_eth  <- list()
glmfit_saps_rlt_list_eth  <- list()
tmle_fit_base_rlt_list_eth  <- list()
tmle_fit_base_sofa_rlt_list_eth  <- list()
tmle_fit_base_saps_rlt_list_eth  <- list()

for (i in 1:length(y_vars)) {
# for (i in 16:18) {
  y = y_vars[i]
  cat(paste0("\n## ", i, ". ", y, "\n\n"))

  source(here("scripts", "04_rates_5d_template.R"))

  glmfit_0_rlt_list_eth[[i]] <- glmfit_0_rlt
  glmfit_sofa_rlt_list_eth[[i]] <- glmfit_sofa_rlt
  glmfit_saps_rlt_list_eth[[i]] <- glmfit_saps_rlt
  tmle_fit_base_rlt_list_eth[[i]] <- tmle_fit_base_rlt
  tmle_fit_base_sofa_rlt_list_eth[[i]] <- tmle_fit_base_sofa_rlt
  tmle_fit_base_saps_rlt_list_eth[[i]] <- tmle_fit_base_saps_rlt

  print("haha\n\n")
}
save.image(here('04_5d_eth.RData'))



## ----gender_hist----------------------------------------------------------------------------------------------------------------------------------------------
df_baseline <- read.csv(here("data", "DataPreprocessing", "df_baseline.csv"), row.names = 1)
df_5d <- merge(df_5d, df_baseline %>% select(ICUSTAY_ID, GENDER), by="ICUSTAY_ID", all.x = T )
df_5d$gender = ifelse(df_5d$GENDER == "F", "Female", "Male")

df_5d$gender = factor(df_5d$gender)
df_5d_all$gender <- df_5d$gender



## ----gender_severity_Scores, fig.height=3, fig.width=8--------------------------------------------------------------------------------------------------------
p1 <- ggplot(df_5d) +
        geom_boxplot(aes(x = gender, y=sofa, col=gender)) +
        theme(legend.position = "none") +
        labs(x = "Gender", y = "score", title = "SOFA")
p2 <- ggplot(df_5d) +
        geom_boxplot(aes(x = gender, y=sapsii, col=gender)) +
        theme(legend.position = "none") +
        labs(x = "Gender", y = "score", title = "SAPS-II")
grid.arrange(p1, p2, nrow = 1)


## ----fig.height=3, fig.width=16-------------------------------------------------------------------------------------------------------------------------------
x = "gender"

## ----echo=FALSE, results='asis'-------------------------------------------------------------------------------------------------------------------------------
x = 'GENDER_F'
baselines = var_name_list$baselines[! var_name_list$baselines %in% c("GENDER_F")]
y_vars <- c(var_name_list$num_vitals, var_name_list$hours_missing_vitals, var_name_list$num_labs)

glmfit_0_rlt_list_gender <- list()
glmfit_sofa_rlt_list_gender  <- list()
glmfit_saps_rlt_list_gender  <- list()
tmle_fit_base_rlt_list_gender  <- list()
tmle_fit_base_sofa_rlt_list_gender  <- list()
tmle_fit_base_saps_rlt_list_gender  <- list()

for (i in 1:length(y_vars)) {
# for (i in 16:18) {
  y = y_vars[i]
  cat(paste0("\n## ", i, ". ", y, "\n\n"))

  source(here("scripts", "04_rates_5d_template.R"))

  glmfit_0_rlt_list_gender[[i]] <- glmfit_0_rlt
  glmfit_sofa_rlt_list_gender[[i]] <- glmfit_sofa_rlt
  glmfit_saps_rlt_list_gender[[i]] <- glmfit_saps_rlt
  tmle_fit_base_rlt_list_gender[[i]] <- tmle_fit_base_rlt
  tmle_fit_base_sofa_rlt_list_gender[[i]] <- tmle_fit_base_sofa_rlt
  tmle_fit_base_saps_rlt_list_gender[[i]] <- tmle_fit_base_saps_rlt

  print("\n haha\n\n")
}
save.image(here('04_5d_gender.RData'))





## ----insurance_hist-------------------------------------------------------------------------------------------------------------------------------------------
df_baseline <- read.csv(here("data", "DataPreprocessing", "df_baseline.csv"), row.names = 1)
df_5d <- merge(df_5d, df_baseline %>% select(ICUSTAY_ID, INSURANCE), by="ICUSTAY_ID", all.x = T )

df_5d$insurance = factor(df_5d$INSURANCE, levels = c('not_insuranced', 'insuranced'))
df_5d_all$insurance <- df_5d$insurance

cat("Summary of insurance: \n")

## ----fig.height=3, fig.width=16-------------------------------------------------------------------------------------------------------------------------------
x = "insurance"

## ----echo=FALSE, results='asis'-------------------------------------------------------------------------------------------------------------------------------
x = 'INSURANCE_Y'
baselines = var_name_list$baselines[! var_name_list$baselines %in% c("INSURANCE_Y")]
y_vars <- c(var_name_list$num_vitals, var_name_list$hours_missing_vitals, var_name_list$num_labs)

glmfit_0_rlt_list_insurance <- list()
glmfit_sofa_rlt_list_insurance  <- list()
glmfit_saps_rlt_list_insurance  <- list()
tmle_fit_base_rlt_list_insurance  <- list()
tmle_fit_base_sofa_rlt_list_insurance  <- list()
tmle_fit_base_saps_rlt_list_insurance  <- list()

for (i in 1:length(y_vars)) {
# for (i in 16:18) {
  y = y_vars[i]
  cat(paste0("\n## ", i, ". ", y, "\n\n"))

  source(here("scripts", "04_rates_5d_template.R"))

  glmfit_0_rlt_list_insurance[[i]] <- glmfit_0_rlt
  glmfit_sofa_rlt_list_insurance[[i]] <- glmfit_sofa_rlt
  glmfit_saps_rlt_list_insurance[[i]] <- glmfit_saps_rlt
  tmle_fit_base_rlt_list_insurance[[i]] <- tmle_fit_base_rlt
  tmle_fit_base_sofa_rlt_list_insurance[[i]] <- tmle_fit_base_sofa_rlt
  tmle_fit_base_saps_rlt_list_insurance[[i]] <- tmle_fit_base_saps_rlt

  print("\n haha\n\n")
}
save.image(here('04_5d_insurance.RData'))




## ----LANGUAGE_hist--------------------------------------------------------------------------------------------------------------------------------------------
df_baseline <- read.csv(here("data", "DataPreprocessing", "df_baseline.csv"), row.names = 1)
df_5d <- merge(df_5d, df_baseline %>% select(ICUSTAY_ID, LANGUAGE), by="ICUSTAY_ID", all.x = T )

df_5d$language <- ifelse(is.na(df_5d$LANGUAGE), "Missing", df_5d$LANGUAGE)

df_5d$language = factor(df_5d$language, levels = c("ENGL", "Other", "Missing"))
df_5d_all$language <- df_5d$language



## ----fig.height=3, fig.width=16-------------------------------------------------------------------------------------------------------------------------------
x = "language"

## ----echo=FALSE, results='asis'-------------------------------------------------------------------------------------------------------------------------------
x = "language"
baselines = var_name_list$baselines[! var_name_list$baselines %in% c("LANGUAGE_ENGL")]
y_vars <- c(var_name_list$num_vitals, var_name_list$hours_missing_vitals, var_name_list$num_labs)

glmfit_0_rlt_list_language <- list()
glmfit_sofa_rlt_list_language  <- list()
glmfit_saps_rlt_list_language  <- list()
tmle_fit_base_rlt_list_language  <- list()
tmle_fit_base_sofa_rlt_list_language  <- list()
tmle_fit_base_saps_rlt_list_language  <- list()

for (i in 1:length(y_vars)) {
  y = y_vars[i]
  cat(paste0("\n## ", i, ". ", y, "\n\n"))

  source(here("scripts", "04_rates_5d_template.R"))

  glmfit_0_rlt_list_language[[i]] <- glmfit_0_rlt
  glmfit_sofa_rlt_list_language[[i]] <- glmfit_sofa_rlt
  glmfit_saps_rlt_list_language[[i]] <- glmfit_saps_rlt
  tmle_fit_base_rlt_list_language[[i]] <- tmle_fit_base_rlt
  tmle_fit_base_sofa_rlt_list_language[[i]] <- tmle_fit_base_sofa_rlt
  tmle_fit_base_saps_rlt_list_language[[i]] <- tmle_fit_base_saps_rlt

  print("\n haha\n\n")
}
save.image(here('04_5d_language.RData'))




## ----marital_hist---------------------------------------------------------------------------------------------------------------------------------------------
df_baseline <- read.csv(here("data", "DataPreprocessing", "df_baseline.csv"), row.names = 1)
df_5d <- merge(df_5d, df_baseline %>% select(ICUSTAY_ID, MARITAL_STATUS), by="ICUSTAY_ID", all.x = T )

df_5d$marital_status = tolower(df_5d$MARITAL_STATUS)
df_5d$marital_status[is.na(df_5d$marital_status)] = 'Missing'
df_5d$marital_status = factor(df_5d$marital_status, levels = c("partnered", "not_partnered", "Missing"))

df_5d_all$marital_status <- df_5d$marital_status




## ----fig.height=3, fig.width=16-------------------------------------------------------------------------------------------------------------------------------
x = "marital_status"

## ----echo=FALSE, results='asis'-------------------------------------------------------------------------------------------------------------------------------
x = 'marital_status'
baselines = var_name_list$baselines[! var_name_list$baselines %in% c("MARITAL_PARTN")]
y_vars <- c(var_name_list$num_vitals, var_name_list$hours_missing_vitals, var_name_list$num_labs)

glmfit_0_rlt_list_partner <- list()
glmfit_sofa_rlt_list_partner  <- list()
glmfit_saps_rlt_list_partner  <- list()
tmle_fit_base_rlt_list_partner  <- list()
tmle_fit_base_sofa_rlt_list_partner  <- list()
tmle_fit_base_saps_rlt_list_partner  <- list()

for (i in 1:length(y_vars)) {
# for (i in 16:18) {
  y = y_vars[i]
  cat(paste0("\n## ", i, ". ", y, "\n\n"))

  source(here("scripts", "04_rates_5d_template.R"))

  glmfit_0_rlt_list_partner[[i]] <- glmfit_0_rlt
  glmfit_sofa_rlt_list_partner[[i]] <- glmfit_sofa_rlt
  glmfit_saps_rlt_list_partner[[i]] <- glmfit_saps_rlt
  tmle_fit_base_rlt_list_partner[[i]] <- tmle_fit_base_rlt
  tmle_fit_base_sofa_rlt_list_partner[[i]] <- tmle_fit_base_sofa_rlt
  tmle_fit_base_saps_rlt_list_partner[[i]] <- tmle_fit_base_saps_rlt

  print("\n haha\n\n")
}
save.image(here('04_5d_partner.RData'))




## ----religion_hist--------------------------------------------------------------------------------------------------------------------------------------------
df_baseline <- read.csv(here("data", "DataPreprocessing", "df_baseline.csv"), row.names = 1)
df_5d <- merge(df_5d, df_baseline %>% select(ICUSTAY_ID, RELIGION), by="ICUSTAY_ID", all.x = T )

df_5d$religion <- tolower(df_5d$RELIGION)
df_5d$religion[is.na(df_5d$religion)] <- "Missing"

df_5d$religion[df_5d$religion == 'christian'] = 'Christian'
df_5d$religion[df_5d$religion == 'other'] = 'Other'

df_5d$religion = factor(df_5d$religion, levels = c("Christian", "Other", "Missing"))
df_5d_all$religion <- df_5d$religion




## ----fig.height=3, fig.width=16-------------------------------------------------------------------------------------------------------------------------------
x = "religion"

## ----echo=FALSE, results='asis'-------------------------------------------------------------------------------------------------------------------------------
x = "religion"
baselines = var_name_list$baselines[! var_name_list$baselines %in% c("RELIGION_Chrt")]
y_vars <- c(var_name_list$num_vitals, var_name_list$hours_missing_vitals, var_name_list$num_labs).

glmfit_0_rlt_list_religion <- list()
glmfit_sofa_rlt_list_religion  <- list()
glmfit_saps_rlt_list_religion  <- list()
tmle_fit_base_rlt_list_religion  <- list()
tmle_fit_base_sofa_rlt_list_religion  <- list()
tmle_fit_base_saps_rlt_list_religion  <- list()

for (i in 1:length(y_vars)) {
# for (i in 16:18) {
  y = y_vars[i]
  cat(paste0("\n## ", i, ". ", y, "\n\n"))

  source(here("scripts", "04_rates_5d_template.R"))

  glmfit_0_rlt_list_religion[[i]] <- glmfit_0_rlt
  glmfit_sofa_rlt_list_religion[[i]] <- glmfit_sofa_rlt
  glmfit_saps_rlt_list_religion[[i]] <- glmfit_saps_rlt
  tmle_fit_base_rlt_list_religion[[i]] <- tmle_fit_base_rlt
  tmle_fit_base_sofa_rlt_list_religion[[i]] <- tmle_fit_base_sofa_rlt
  tmle_fit_base_saps_rlt_list_religion[[i]] <- tmle_fit_base_saps_rlt

  print("\n haha\n\n")
}
save.image(here('04_5d_religion.RData'))




## ----careunit_hist--------------------------------------------------------------------------------------------------------------------------------------------
ICUSTAYS <- read.csv(here("data", "raw_data", "mimic-iii-clinical-database-1.4", "ICUSTAYS.csv"))
df_5d <- merge(df_5d, ICUSTAYS %>% select(ICUSTAY_ID, FIRST_CAREUNIT), by="ICUSTAY_ID", all.x = T )

df_5d$FIRST_CAREUNIT <- ifelse(is.na(df_5d$FIRST_CAREUNIT), "Missing", df_5d$FIRST_CAREUNIT)
df_5d_all$FIRST_CAREUNIT <- df_5d$FIRST_CAREUNIT


## ----fig.height=3, fig.width=16-------------------------------------------------------------------------------------------------------------------------------
x = "FIRST_CAREUNIT"


## ----echo=FALSE, results='asis'-------------------------------------------------------------------------------------------------------------------------------
x = "FIRST_CAREUNIT"
baselines = var_name_list$baselines
y_vars <- c(var_name_list$num_vitals, var_name_list$hours_missing_vitals, var_name_list$num_labs)

glmfit_0_rlt_list_careunit <- list()
glmfit_sofa_rlt_list_careunit  <- list()
glmfit_saps_rlt_list_careunit  <- list()
tmle_fit_base_rlt_list_careunit  <- list()
tmle_fit_base_sofa_rlt_list_careunit  <- list()
tmle_fit_base_saps_rlt_list_careunit  <- list()

for (i in 1:length(y_vars)) {
# for (i in 16:18) {
  y = y_vars[i]
  cat(paste0("\n## ", i, ". ", y, "\n\n"))

  source(here("scripts", "04_rates_5d_template.R"))

  glmfit_0_rlt_list_careunit[[i]] <- glmfit_0_rlt
  glmfit_sofa_rlt_list_careunit[[i]] <- glmfit_sofa_rlt
  glmfit_saps_rlt_list_careunit[[i]] <- glmfit_saps_rlt
  tmle_fit_base_rlt_list_careunit[[i]] <- tmle_fit_base_rlt
  tmle_fit_base_sofa_rlt_list_careunit[[i]] <- tmle_fit_base_sofa_rlt
  tmle_fit_base_saps_rlt_list_careunit[[i]] <- tmle_fit_base_saps_rlt

  print("\n haha\n\n")
}
save.image(here('04_5d_careunit.RData'))

