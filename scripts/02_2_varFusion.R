## ----load_lib, include = FALSE, warning=FALSE, message=FALSE, echo=FALSE---------------------------------------------------------------------------------------------------
library(lubridate)
library(stringr)
library(tidyr)
library(R.utils)
library(data.table)
library(corrplot)
library(reshape)
library(ggplot2)
library(GGally)
library(dplyr)
library(cowplot)
library(ggfortify)
library(ggpubr)
library(mixOmics)
library(gplots)
library(here)


## ----setup, include = FALSE------------------------------------------------------------------------------------------------------------------------------------------------
plotFolder <- here("results","images", "DataPreprocessing")
if(!file.exists(plotFolder)) dir.create(plotFolder,recursive=TRUE)

dataFolder <- here("data", "DataPreprocessing")
if(!file.exists(dataFolder)) dir.create(dataFolder,recursive=TRUE)

rdataFolder <- here("data", "rdata")
if(!file.exists(rdataFolder)) dir.create(dataFolder,recursive=TRUE)

knitr::opts_chunk$set(
  cache=FALSE, autodep=FALSE, warning=FALSE, message=FALSE, echo=FALSE,
  results = 'markup', dev='png', dpi=150, fig.align = "center", fig.path=plotFolder,
  cache.path=".cache/",
  duplicate.label="allow"
)

set.seed(123)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load(here(rdataFolder, "01_DataPreprocessing.RData"))

df <- read.csv(file=here("data","DataPreprocessing","01_data_after_processing.csv"), row.names = 1, header= TRUE)
df_imp <- read.csv(here("data", "DataPreprocessing", "02_1_data_after_imputation.csv"), row.names = 1, header= TRUE)  
df_r_imp <- read.csv(here("data", "DataPreprocessing", "02_1_data_after_random_imputation.csv"), row.names = 1, header= TRUE)  


## ----correlation_matrix, fig.height = 20, fig.width = 20-------------------------------------------------------------------------------------------------------------------
corrplot(cor(df_imp %>% mutate(Y_t_plus_1 = as.numeric(Y_t_plus_1))), 
         method="color", 
         tl.col = unname(unlist(var_name_colors)))

png(file = sprintf("%s/correlation_matrix_original.png",plotFolder), width = 1800, height = 1800)
g <- corrplot(cor(df_imp %>% mutate(Y_t_plus_1 = as.numeric(Y_t_plus_1))), method="color", tl.col = unname(unlist(var_name_colors)))
dev.off()


## ----correlation_matrix_supplemental_vars, fig.height = 10, fig.width = 10-------------------------------------------------------------------------------------------------
cor_sup_var <- cor(df_imp %>% 
                     mutate(Y_t_plus_1 = as.numeric(Y_t_plus_1)) %>%
                      dplyr::select(c(var_name_list$num_baselines, 
                                      var_name_list$hours_missing_vitals, 
                                      var_name_list$num_vitals,
                                      var_name_list$num_outlier_vitals,
                                      var_name_list$num_labs)))
corrplot(cor_sup_var, 
         method="color", 
         tl.col = unname(c(var_name_colors$num_baselines, 
                                      var_name_colors$hours_missing_vitals, 
                                      var_name_colors$num_vitals,
                                      var_name_colors$num_outlier_vitals,
                                      var_name_colors$num_labs)))

png(file = sprintf("%s/correlation_supplemental_vars.png",plotFolder), width = 1800, height = 1800)
g <- corrplot(cor_sup_var, method="color", tl.col = unname(c(var_name_colors$num_baselines, 
                                      var_name_colors$hours_missing_vitals, 
                                      var_name_colors$num_vitals,
                                      var_name_colors$num_outlier_vitals,
                                      var_name_colors$num_labs)))
dev.off()


## ----correlation_matrix_supplemental_vars_vital, fig.height = 7, fig.width = 7---------------------------------------------------------------------------------------------
# cor_sup_var_vital <- cor(df_imp %>% 
#                      dplyr::select(c(var_name_list$hours_missing_vitals, 
#                                       var_name_list$num_vitals,
#                                       var_name_list$num_outlier_vitals)))
# corrplot(cor_sup_var_vital, 
#          method="color", 
#          tl.col = unname(c(var_name_colors$hours_missing_vitals, 
#                                       var_name_colors$num_vitals,
#                                       var_name_colors$num_outlier_vitals)))


## ----correlation_matrix_supplemental_vars_lab, fig.height = 8, fig.width = 8-----------------------------------------------------------------------------------------------
# cor_sup_var_lab <- cor(df_imp %>% 
#                       dplyr::select(c(var_name_list$num_labs)))
# corrplot(cor_sup_var_lab, 
#          method="color", 
#          tl.col = unname(c(var_name_colors$num_labs)))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# var_name_list$hours_missing_vitals <- c(var_name_list$hours_missing_vitals[! var_name_list$hours_missing_vitals %in% c("h_m_BPs", "h_m_BPd", "h_m_GCS_Eye", "h_m_GCS_Motor", "h_m_GCS_Verbal", "h_m_HR", "h_m_SpO2")], "h_m_BP", "h_m_GCS", "h_m_HR_SpO2")
var_name_list$hours_missing_vitals <- c("h_m_BP", "h_m_BPm", "h_m_GCS", "h_m_GCS_Total", "h_m_HR", "h_m_SpO2", "h_m_RR", "h_m_TempC")
var_name_list$num_vitals <-  c("n_BP", "n_BPm", "n_GCS", "n_GCS_Total", "n_HR", "n_SpO2", "n_RR", "n_TempC")
var_name_list$num_labs <- c(var_name_list$num_labs[! var_name_list$num_labs 
                                                   %in% c("n_ALT", "n_ALK", "n_AST", "n_TBil",
                                                          "n_pH", "n_PCO2", "n_PO2", "n_BE",
                                                          "n_Na", "n_K", "n_Cl", "n_HCO3", "n_AG", "n_BG", "n_BUN", "n_Cr",
                                                          "n_WBCs", "n_RBCs", "n_HGB", "n_MCV", "n_MCH", "n_MCHC", "n_RDW", "n_PLT",
                                                          "n_MO", "n_EO", "n_BA", "n_NE",
                                                          "n_Ca", "n_Mg", "n_Phos",
                                                          "n_PTT", "n_PT")], 
                            "n_lab_grp1", "n_lab_grp2", "n_lab_grp3", "n_lab_grp4", "n_lab_grp5", "n_lab_grp6", "n_lab_grp7")


var_name_colors <- list(ID = rep("gray", 2),
                        outcome = "black", 
                        baselines = rep("tan3", length(var_name_list$baselines)),
                        num_baselines = rep("tan4", length(var_name_list$num_baselines)),
                        vitals = rep("darkorchid1", length(var_name_list$vitals)),
                        hours_missing_vitals = rep("darkorchid4", length(var_name_list$hours_missing_vitals)),
                        num_vitals = rep("skyblue2", length(var_name_list$num_vitals)),
                        num_outlier_vitals = rep("steelblue4", length(var_name_list$num_outlier_vitals)),  
                        labs = rep("palegreen3", length(var_name_list$labs)),
                        num_labs = rep("palegreen4", length(var_name_list$num_labs)) )
                        #num_outlier_labs = rep("palegreen4", length(var_name_list$num_outlier_labs)))


var_fusion <- function(dat){
  dat <- dat %>%
    mutate(n_BP = rowMeans(dat %>% dplyr::select("n_BPs", "n_BPd")),
           h_m_BP = rowMeans(dat %>% dplyr::select("h_m_BPs", "h_m_BPd")),
           
           n_GCS = rowMeans(dat %>% dplyr::select("n_GCS_Eye", "n_GCS_Motor", "n_GCS_Verbal")),
           h_m_GCS = rowMeans(dat %>% dplyr::select("h_m_GCS_Eye", "h_m_GCS_Motor", "h_m_GCS_Verbal")),
           
           # n_HR_SpO2 = rowMeans(dat %>% dplyr::select("n_HR", "n_SpO2")),
           # h_m_HR_SpO2 = rowMeans(dat %>% dplyr::select("h_m_HR", "h_m_SpO2")),
           
           n_lab_grp1 = rowMeans(dat %>% dplyr::select("n_ALT", "n_ALK", "n_AST", "n_TBil")), #4
           n_lab_grp2 = rowMeans(dat %>% dplyr::select("n_pH", "n_PCO2", "n_PO2", "n_BE")), # 4
           n_lab_grp3 = rowMeans(dat %>% dplyr::select("n_Na", "n_K", "n_Cl", "n_HCO3", "n_AG", "n_BG", "n_BUN", "n_Cr")), # 8
           n_lab_grp4 = rowMeans(dat %>% dplyr::select("n_WBCs", "n_RBCs", "n_HGB", "n_MCV", "n_MCH", "n_MCHC", "n_RDW", "n_PLT")), # 8
           n_lab_grp5 = rowMeans(dat %>% dplyr::select("n_MO", "n_EO", "n_BA", "n_NE")), #4
           n_lab_grp6 = rowMeans(dat %>% dplyr::select("n_Ca", "n_Mg", "n_Phos")),
           n_lab_grp7 = rowMeans(dat %>% dplyr::select("n_PTT", "n_PT")),
           ) %>%
    dplyr::select(unname(unlist(var_name_list)))
  return(dat)
}



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df <- var_fusion(df)
df_imp <- var_fusion(df_imp)
df_r_imp <- var_fusion(df_r_imp)

## ----correlation_matrix_fusion, fig.height = 20, fig.width = 20------------------------------------------------------------------------------------------------------------
corrplot(cor(df_imp %>% mutate(Y_t_plus_1 = as.numeric(Y_t_plus_1))), 
         method="color", 
         tl.col = unname(unlist(var_name_colors)))

png(file = sprintf("%s/correlation_matrix_varfusion.png",plotFolder), width = 1800, height = 1800)
g <- corrplot(cor(df_imp %>% mutate(Y_t_plus_1 = as.numeric(Y_t_plus_1))), method="color", tl.col = unname(unlist(var_name_colors)))
dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(df, file=here(dataFolder,"02_2_data.csv"))
write.csv(df_imp, file=here(dataFolder,"02_2_data_imputed.csv"))
write.csv(df_r_imp, file=here(dataFolder,"02_2_data_random_imputed.csv"))


## ----blocks_outcome--------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(data = df) +
  geom_bar(aes(x = half_day, color = Y_t_plus_1, fill = Y_t_plus_1)) +
  scale_x_continuous(breaks= 1:10) +
  labs(x = "12h block")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list= ls()[!(ls() %in% c('var_name_list', 'var_name_colors', 'vitals_labs_items_election_info', 'vitals_labs_summary'))])
rdataFolder <- here("data", "rdata")

save.image(file=here(rdataFolder, "02_2_varFusion.RData"))

