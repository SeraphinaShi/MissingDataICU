## ----load_lib, include = FALSE, warning=FALSE, message=FALSE, echo=FALSE-------------------------------------------------------------------------
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


## ----setup, include = FALSE----------------------------------------------------------------------------------------------------------------------
plotFolder <- here("results","images", "v2", "DataPreprocessing")
if(!file.exists(plotFolder)) dir.create(plotFolder,recursive=TRUE)

dataFolder <- here("data", "data_v2", "DataPreprocessing")
if(!file.exists(dataFolder)) dir.create(dataFolder,recursive=TRUE)

rdataFolder <- here("data", "data_v2", "rdata")
if(!file.exists(rdataFolder)) dir.create(dataFolder,recursive=TRUE)

knitr::opts_chunk$set(
  cache=FALSE, autodep=FALSE, warning=FALSE, message=FALSE, echo=FALSE,
  results = 'markup', dev='png', dpi=150, fig.align = "center", fig.path=plotFolder,
  cache.path=".cache/",
  duplicate.label="allow"
)

set.seed(123)


## ------------------------------------------------------------------------------------------------------------------------------------------------
load(here(rdataFolder, "02_2_varFusion.RData"))
df <- read.csv(here(dataFolder, "02_2_data.csv")) %>% dplyr::select(-X)
df_imp <- read.csv(here(dataFolder, "02_2_data_imputed.csv"))  %>% dplyr::select(-X)


## ------------------------------------------------------------------------------------------------------------------------------------------------
# get current clock info
df_Y <- read.csv(here("data", "data_v2", "DataPreprocessing", "df_Y.csv"), row.names = 1, header = T)

df <- merge(df_Y %>% dplyr::select(ICUSTAY_ID, half_day, Y_t, Y_t_plus_1), df, 
            by = c("ICUSTAY_ID", "half_day", "Y_t_plus_1"), all.y = T)
df_imp <- merge(df_Y %>% dplyr::select(ICUSTAY_ID, half_day, Y_t, Y_t_plus_1), df_imp, 
            by = c("ICUSTAY_ID", "half_day", "Y_t_plus_1"), all.y = T)

# with all blocks, including the incomplete ones
df_all <- df
df_imp_all <- df_imp

# with complete blocks
df <- df %>% 
  filter(! Y_t %in% c("1", "2", "Gone"))
df_imp <- df_imp %>% 
  filter(! Y_t %in% c("1", "2", "Gone"))


# change outcome to 0 or 1 for Y_t_plus_1
df$Y_t_plus_1 = ifelse(df$Y_t_plus_1 == "2","0", df$Y_t_plus_1)

# View(df[, c("ICUSTAY_ID", "half_day", "Y_t", "Y_t_plus_1")])


## ------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(df, file=here(dataFolder,"02_3_data_complt_blocks.csv"))
write.csv(df_imp, file=here(dataFolder,"02_3_data_imputed_complt_blocks.csv"))
write.csv(df_all, file=here(dataFolder,"02_3_data_all_blocks.csv"))
write.csv(df_imp_all, file=here(dataFolder,"02_3_data_imputed_all_blocks.csv"))


## ------------------------------------------------------------------------------------------------------------------------------------------------
rm(list= ls()[!(ls() %in% c('var_name_list', 'var_name_colors', 'vitals_labs_items_election_info', 'vitals_labs_summary'))])
rdataFolder <- here("data", "data_v2", "rdata")
save.image(file=here(rdataFolder, "02_3_complt_blocks.RData"))
# knitr::purl("scripts/02_2_varFusion.Rmd")

