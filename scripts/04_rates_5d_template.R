
#---------------------------
### `r i`.1. no confounders
set.seed(123)
glmfit_0 <- glm(formula = as.formula(paste0(y, "~ ", x, "+ less_5d")),
              data = df_5d)
glmfit_0_rlt <- summary(glmfit_0)$coefficients %>% as.data.frame() %>% mutate(across(where(is.numeric), round, 4))
print(summary(glmfit_0))

glmfit_0_rlt$Y = y

#---------------------------
### `r i`.2. SOFA score as a confounder
set.seed(123)
glmfit_sofa <- glm(formula =as.formula(paste0(y, "~ ", x, " + sofa + less_5d")),
              data = df_5d)
glmfit_sofa_rlt <- summary(glmfit_sofa)$coefficients %>% as.data.frame()
print(summary(glmfit_sofa))

glmfit_sofa_rlt$Y = y

#---------------------------
### `r i`.3. SAPS-II score as a confounder
set.seed(123)
glmfit_saps <- glm(formula =as.formula(paste0(y, "~ ", x, " + sapsii + less_5d")),
              data = df_5d)
glmfit_saps_rlt <- summary(glmfit_saps)$coefficients %>% as.data.frame()
print(summary(glmfit_saps))

glmfit_saps_rlt$Y = y

#---------------------------
### `r i`.4. Baseline confounders
Y = df_5d %>% pull(all_of(y))
A = df_5d %>% pull(all_of(x))

confounders <- df_5d %>% select(all_of(baselines), less_5d)

tmle_fit_base <- run_tmle3_outcome_rates(Y, A, confounders)
tmle_fit_base_rlt <- tmle_fit_base$summary %>% mutate(across(where(is.numeric), round, 4))
print(tmle_fit_base_rlt)
tmle_fit_base_rlt$Y = y

#---------------------------
### `r i`.5. Baseline & SOFA confounders
confounders <- df_5d %>% select(all_of(baselines), sofa, less_5d)

tmle_fit_base_sofa <- run_tmle3_outcome_rates(Y, A, confounders)
tmle_fit_base_sofa_rlt <- tmle_fit_base_sofa$summary %>% mutate(across(where(is.numeric), round, 4))
print(tmle_fit_base_sofa_rlt)
tmle_fit_base_sofa_rlt$Y = y

#---------------------------
### `r i`.6. Baseline & SAPS-II confounders
confounders <- df_5d %>% select(all_of(baselines), sapsii, less_5d)

tmle_fit_base_saps <- run_tmle3_outcome_rates(Y, A, confounders)
tmle_fit_base_saps_rlt <- tmle_fit_base_saps$summary %>% mutate(across(where(is.numeric), round, 4))
print(tmle_fit_base_saps_rlt)
tmle_fit_base_saps_rlt$Y = y


