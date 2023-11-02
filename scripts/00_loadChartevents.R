
library(lubridate)
library(dplyr)
library(stringr)
library(tidyr)
library(R.utils)
library(data.table)
library(here)

# nrow_CHARTEVENTS <- length(count.fields("mimic-iii-clinical-database-1.4/CHARTEVENTS.csv", sep = ",")) # 330712485
nrow_CHARTEVENTS <- 330712485
chunk_size <- (nrow_CHARTEVENTS - 1)/20 # 16535624

# create a empty list to store the chunks of data
CHARTEVENTS_lists <- list()

# read in the data 
CHARTEVENTS <- fread(here("data","raw_data","mimic-iii-clinical-database-1.4", "CHARTEVENTS.csv"), header = T, stringsAsFactors = FALSE)


CHARTEVENTS <- CHARTEVENTS %>% 
  filter(HADM_ID %in% all_adm_id$id) %>%
  select(ROW_ID, ICUSTAY_ID, ITEMID, CHARTTIME, VALUENUM, VALUEUOM, ERROR)  %>% 
  filter(ITEMID %in% D_ITEMS_vital_sub$ITEMID) %>%
  filter(!is.na(VALUENUM)) %>% 
  filter(VALUENUM != 0) %>%
  mutate(CHARTTIME  = ymd_hms(CHARTTIME))

CHARTEVENTS <- unique(CHARTEVENTS)

save.image(file=here(rdataFolder, "LoadCHARTEVENTS_lists_all.RData"))


