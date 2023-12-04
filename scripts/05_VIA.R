## ----load_lib, include = FALSE, warning=FALSE, message=FALSE, echo=FALSE------------------------------------
# setwd("/home/seraphinashi/MissingDataInTheICU")
# install.packages("polspline")
library(dplyr)
library(data.table)
library(SuperLearner)
library(origami)
library(sl3)
library(ggplot2)
library(forcats)
library(here)
library(future)
library(ROCR) 


## ----setup, include = FALSE---------------------------------------------------------------------------------
plotFolder <- here("results","images","05_VIA")
if(!file.exists(plotFolder)) dir.create(plotFolder,recursive=TRUE)

rdataFolder <- here("data", "rdata")

knitr::opts_chunk$set(
  cache=FALSE, autodep=FALSE, warning=FALSE, message=FALSE, echo=FALSE,
  results = 'markup', dev='png', dpi=150, fig.align = "center", fig.path=plotFolder,
  cache.path=".cache/",
  duplicate.label="allow"
)

set.seed(123)



## -----------------------------------------------------------------------------------------------------------
load(here(rdataFolder, "02_2_varFusion.RData"))

PreprocessDataFolder <- here("data", "DataPreprocessing")

df <- read.csv(here(PreprocessDataFolder,"02_3_data_complt_blocks.csv"), row.names = 1, header = T)  
df_imp <- read.csv(here(PreprocessDataFolder,"02_3_data_imputed_complt_blocks.csv"), row.names = 1, header = T)
df_r_imp <- read.csv(here(PreprocessDataFolder,"02_3_data_random_imputed_complt_blocks.csv"), row.names = 1, header = T)


## -----------------------------------------------------------------------------------------------------------
cat("Missing value percentage for each column in the data: \n")
cat("  Baselines: \n")
print(round((colMeans(is.na(df %>% select(all_of(var_name_list$baselines)))))*100,4))

cat("  Vital signs: \n")
print(round((colMeans(is.na(df %>% select(var_name_list$vitals))))*100,4))
cat("  Lab tests: \n")
print(round((colMeans(is.na(df %>% select(var_name_list$labs))))*100,4))


## -----------------------------------------------------------------------------------------------------------
get_df_prev <- function(df){
    df_prev <- df %>% 
      select(c(var_name_list$ID, 
               var_name_list$vitals, var_name_list$hours_missing_vitals, var_name_list$num_vitals,
               var_name_list$labs, var_name_list$num_labs)) 
    names(df_prev)[-c(1:2)] <- paste0(names(df_prev)[-c(1:2)], "_prev")
    df_prev$half_day_prev <- df_prev$half_day 
    df_prev$half_day <- df_prev$half_day + 1
    
    df_prev  <- df_prev[df_prev$half_day <= 10, ]
    df <- merge(df, df_prev, by = c("ICUSTAY_ID", "half_day"), all.x = T)
    df <- df[df$half_day != 1, ]
    
    return(df)
}

set.seed(123)
smaller_id_set <- c(sample(unique(df$ICUSTAY_ID[df$Y_t_plus_1 == 0]), 3000), df$ICUSTAY_ID[df$Y_t_plus_1 == 1])
smaller_df_idx <- which(df$ICUSTAY_ID %in% smaller_id_set)

get_df_sl <- function(sets, smaller_df = T, prev_block = F, inputation = "median"){
  if(smaller_df){
    dat <- df[smaller_df_idx, ]
    dat_imp <- df_imp[smaller_df_idx,]
    dat_r_imp <- df_r_imp[smaller_df_idx, ]
  } else {
    dat <- df
    dat_imp <- df_imp
    dat_r_imp <- df_r_imp
  }

  
  if(prev_block){
    df_prev <- get_df_prev(dat)
    dat_imp_prev <- get_df_prev(dat_imp)
    dat_r_imp <- get_df_prev(df_r_imp)
  }
  
  vars <- c()
  
  if("VS" %in% sets)  vars <- c(vars, var_name_list$vitals)
  if("h_m_VS" %in% sets)  vars <- c(vars, var_name_list$hours_missing_vitals)
  if("n_VS" %in% sets)  vars <- c(vars, var_name_list$num_vitals)
  if("LT" %in% sets)  vars <- c(vars, var_name_list$labs)
  if("n_LT" %in% sets)  vars <- c(vars, var_name_list$num_labs)
  
  if(prev_block){
    if(length(vars) > 0){
      vars <- c(vars, paste0(vars, "_prev"))
    }
  }
  
  if("W" %in% sets)  vars <- c(vars, var_name_list$baselines)
  if("n_W" %in% sets)  vars <- c(vars, var_name_list$num_baselines)
  vars <- c(var_name_list$ID,
            var_name_list$outcome,
            vars)
  
  df_sl3 <- dat %>% dplyr::select(all_of(vars))
  
  if(inputation == "CCA"){
    df_sl3 <- df_sl3[complete.cases(df_sl3), ] 
  } else if(inputation == "random_unif"){
    df_sl3 <- dat_r_imp %>% dplyr::select(all_of(vars))
  } else {
    df_sl3 <- dat_imp %>% dplyr::select(all_of(vars))
  }
  return(df_sl3)
}

exp_var_lists <- list(c("W", "VS", "LT"),
                      c("n_W", "n_VS", "n_LT"),
                      c("n_W", "h_m_VS", "n_LT"),
                      c("n_W", "n_VS", "h_m_VS", "n_LT"),
                      c("W", "VS", "LT", "n_W", "n_VS", "n_LT"),
                      c("W", "VS", "LT", "n_W", "h_m_VS", "n_LT"),
                      c("W", "VS", "LT", "n_W", "n_VS", "h_m_VS", "n_LT")
)



## -----------------------------------------------------------------------------------------------------------
df_sl3 <- get_df_sl(sets = c("W", "VS", "LT", "n_W", "n_VS", "h_m_VS", "n_LT"), prev_block = F)



df_sl3 <- df_sl3 %>% select(-half_day)

# set.seed(123)
# df_sl3 <- df_sl3[c(sample(nrow(df_sl3), 100), which(df_sl3$Y_t_plus_1 == 1)), ]
## define sl3 task

# define the folds with clusters and stratification
set.seed(123)
folds <- make_folds(nrow(df_sl3), cluster_ids = df_sl3$ICUSTAY_ID)

# define the prediction task
sl3_task <- make_sl3_Task(
  data = df_sl3, 
  covariates = colnames(df_sl3)[-which(names(df_sl3) %in% c(var_name_list$outcome,
                                                              var_name_list$ID))],
  outcome = var_name_list$outcome,
  outcome_type = "binomial",
  id = "ICUSTAY_ID",
  folds = folds
)
  #====================
  
  ######################### build library for discrete SL ########################
  

# ,  , , , , , , , and the Ensemble SL itself

lrn_glm <- Lrnr_glm$new() # Generalized Linear Model (GLM)
lrn_earth <- Lrnr_earth$new() # Generalized Additive Model (GAM)
lrn_ranger <- Lrnr_ranger$new() # Random Forests, 
lrn_bayesglm <- Lrnr_pkg_SuperLearner$new("SL.bayesglm") # Bayesian GLM,
lrn_lasso <- make_learner(Lrnr_glmnet) # Lasso Regression
lrn_ridge <- Lrnr_glmnet$new(alpha = 0) # Ridge Regression
lrn_enet.5 <- make_learner(Lrnr_glmnet, alpha = 0.5) # Elastic-net Penalized Regression
lrn_polspline <- Lrnr_polspline$new() # Multivariate Adaptive Regression Splines
Lrn_lightgbm <- Lrnr_lightgbm$new() # Gradient Boosting Machine
Lrn_dbarts <- Lrnr_dbarts$new() # Bayesian Additive Regression Trees

base_lrnrs <- c(lrn_glm, lrn_earth, lrn_ranger, lrn_bayesglm, lrn_lasso, lrn_ridge, lrn_enet.5, Lrn_lightgbm, Lrn_dbarts) # Lrn_lightgbm, Lrn_dbarts
names(base_lrnrs) <- c("GLM", "Gam", "RF", "BayesianGLM", "Lassso", "Ridge", "ElasticNet", "XGBoost", "BayesianTrees") # "XGBoost", "BayesianTrees"

# base_lrnrs <- c(lrn_glm, lrn_earth, lrn_ranger, lrn_bayesglm, lrn_lasso, lrn_ridge, lrn_enet.5, lrn_polspline, lrn_polspline) # Lrn_lightgbm, Lrn_dbarts
# names(base_lrnrs) <- c("GLM", "Gam", "RF", "BayesianGLM", "Lassso", "Ridge", "ElasticNet", "Splines") # "XGBoost", "BayesianTrees"

base_lrnrs_stack <- Stack$new(base_lrnrs) 

sl_default <- Lrnr_sl$new(learners = base_lrnrs_stack)

candidate_lrnrs <- c(sl_default, base_lrnrs) 
names(candidate_lrnrs)[1:3] <- c("SL_defaultMeta", names(base_lrnrs))
candidate_lrnrs_stack <- Stack$new(candidate_lrnrs)
    


########################### Train discrete SL ##################################

discrete <- Lrnr_cv_selector$new(loss_loglik_binomial)

sl <- Lrnr_sl$new(learners = candidate_lrnrs_stack, metalearner = discrete)


############# how to run without parallelization: 
set.seed(45)
fit <- sl$train(sl3_task)
# save.image(file=here("data", "rdata", "05_VIA_sl3_all_0.RData"))



############################# summarize results ################################

# lines 71 and 72 should be the same
print(fit$coefficients)
print(fit$fit_object$cv_meta_fit$fit_object$coef)

print(fit$learner_fits$SL_defaultMeta)

cv_risk_table_NLL <- fit$cv_risk(eval_fun = loss_loglik_binomial)
cv_risk_table_NLL <- cv_risk_table_NLL %>% 
  rename(se_nll = se, fold_sd_nll = fold_sd) %>% 
  mutate(ci_l_nll = NLL - 1.96 * se_nll,
         ci_u_nll = NLL + 1.96 * se_nll,)

print(cv_risk_table_NLL)

auc_eval_function <- custom_ROCR_risk("auc")
cv_risk_table_AUC <- fit$cv_risk(eval_fun = auc_eval_function)
cv_risk_table_AUC <- cv_risk_table_AUC %>% rename(se_auc = se, fold_sd_auc = fold_sd) %>% 
  mutate(ci_l_auc = auc - 1.96 * se_auc,
         ci_u_auc = auc + 1.96 * se_auc,)

print(cv_risk_table_AUC)

# area under the precision recall curve is recommended for rare outcomes
aucpr_eval_function <- custom_ROCR_risk("aucpr")
cv_risk_table_AUCPR <- fit$cv_risk(eval_fun = aucpr_eval_function)
cv_risk_table_AUCPR <- cv_risk_table_AUCPR %>% rename(se_aucpr = se, fold_sd_aucpr = fold_sd) %>% 
  mutate(ci_l_aucpr = aucpr - 1.96 * se_aucpr,
         ci_u_aucpr = aucpr + 1.96 * se_aucpr,)

print(cv_risk_table_AUCPR)


save.image(file=here("data", "rdata", "05_VIA_sl3_all.RData"))

