## ----load_lib, include = FALSE, warning=FALSE, message=FALSE, echo=FALSE---------------------------------------------------------------------------------------------------
# setwd("/home/seraphinashi/MissingDataInTheICU")

library(dplyr)
library(here)

set.seed(123)
## ----setup, include = FALSE------------------------------------------------------------------------------------------------------------------------------------------------
# description:
# https://mimic.mit.edu/iv/datasets/core/
plotFolder <- here("results","images", "02_1_DataPreprocessing")
if(!file.exists(plotFolder)) dir.create(plotFolder,recursive=TRUE)

dataFolder <- here("data","DataPreprocessing")
if(!file.exists(dataFolder)) dir.create(dataFolder,recursive=TRUE)


rdataFolder <- here("data",  "rdata")

knitr::opts_chunk$set(
  cache=FALSE, autodep=FALSE, warning=FALSE, message=FALSE, echo=FALSE,
  results = 'markup', dev='png', dpi=150, fig.align = "center", fig.path=plotFolder,
  cache.path=".cache/",
  duplicate.label="allow"
)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load(here(rdataFolder, "01_DataPreprocessing.RData"))
df <- read.csv(file=here("data", "DataPreprocessing", "01_data_after_processing.csv"), row.names = 1, header= TRUE)

df_random_imp <- df
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# vital and lab values
for(col_i in var_name_list$vitals){
  median_i = vitals_labs_summary$vital_summary$median[vitals_labs_summary$vital_summary$ABB == col_i]
  min_i = vitals_labs_summary$vital_summary$min_clinical[vitals_labs_summary$vital_summary$ABB == col_i]
  max_i = vitals_labs_summary$vital_summary$max_clinical[vitals_labs_summary$vital_summary$ABB == col_i]
  
  df[[col_i]][is.na(df[[col_i]])] = median_i
  df_random_imp[[col_i]][is.na(df_random_imp[[col_i]])] = runif(1, min = min_i, max = max_i)
}

for(col_i in var_name_list$labs){
  median_i =  vitals_labs_summary$lab_summary$median[vitals_labs_summary$lab_summary$ABB == col_i]
  min_i = vitals_labs_summary$lab_summary$min[vitals_labs_summary$lab_summary$ABB == col_i]
  max_i = vitals_labs_summary$lab_summary$max[vitals_labs_summary$lab_summary$ABB == col_i]
  
  df[[col_i]][is.na(df[[col_i]])] = median_i
  df_random_imp[[col_i]][is.na(df_random_imp[[col_i]])] = runif(1, min = min_i, max = max_i)
}


## ----fig.height = 6, fig.width = 10----------------------------------------------------------------------------------------------------------------------------------------

Mode <- function(x, na.rm = T) {
  if(na.rm){
    x = x[!is.na(x)]
  }

  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

for(col_i in var_name_list$baselines[-1]){
  mode <- Mode(df[[col_i]])

  df[[col_i]][is.na(df[[col_i]])] = mode
  df_random_imp[[col_i]][is.na(df_random_imp[[col_i]])] = mode
}


## ----fig.height = 6, fig.width = 10----------------------------------------------------------------------------------------------------------------------------------------
# Imputing outliers to the closest upper or lower bound values. 
# for(col_i in var_name_list$vitals){
#   min <- vitals_labs_summary$vital_summary[vitals_labs_summary$vital_summary$ABB == col_i, "min_clinical"]
#   max <- vitals_labs_summary$vital_summary[vitals_labs_summary$vital_summary$ABB == col_i, "max_clinical"]
# 
#   df[[col_i]][df[[col_i]] > max] = max
#   df[[col_i]][df[[col_i]] < min] = min
# 
# }

for(col_i in var_name_list$vitals){
  median <- vitals_labs_summary$vital_summary$median[vitals_labs_summary$vital_summary$ABB == col_i]
  max_clinical <- vitals_labs_summary$vital_summary$max_clinical[vitals_labs_summary$vital_summary$ABB == col_i]
  min_clinical <- vitals_labs_summary$vital_summary$min_clinical[vitals_labs_summary$vital_summary$ABB == col_i]
  df[[col_i]][df[[col_i]] > max_clinical] = median
  df[[col_i]][df[[col_i]] < min_clinical] = median
  
  df_random_imp[[col_i]][df_random_imp[[col_i]] > max_clinical] = median
  df_random_imp[[col_i]][df_random_imp[[col_i]] < min_clinical] = median
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(df, file=here(dataFolder,"02_1_data_after_imputation.csv"))
write.csv(df_random_imp, file=here(dataFolder,"02_1_data_after_random_imputation.csv"))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list= ls()[!(ls() %in% c('df', 'var_name_list', 'var_name_colors', 'vitals_labs_items_election_info', 'vitals_labs_summary'))])

rdataFolder <- here("data", "rdata")
save.image(file=here(rdataFolder, "02_1_DataImputation.RData"))

summary(df)

# knitr::purl("scripts/02_DataImputation.Rmd")

