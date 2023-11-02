
library(data.table)
library(dplyr)
library(here)

D_ITEM_A_LINE <- c(1449, 51, 52, 53, 776, 777, 778, 779, 780, 2529, 8368, 8555, 6590, 6701, 6702, 6926, 6927, 
                   225210, 225556, 227292, 225575, 224289, 224291, 227543, 227546, 227547, 225698, 225722, 225737
                   , 225752, 226008, 228022, 228023, 228026, 225986, 226107, 220050, 220051, 220052, 220056, 220058, 
                   220224, 220227, 220235, 224828, 223830)

# nrow_CHARTEVENTS <- length(count.fields("mimic-iii-clinical-database-1.4/CHARTEVENTS.csv", sep = ",")) # 330712485
nrow_CHARTEVENTS <- 330712485
chunk_size <- (nrow_CHARTEVENTS - 1)/20 # 16535624


# read in the data 
CHARTEVENTS_i <- fread(here("data","raw_data","mimic-iii-clinical-database-1.4", "CHARTEVENTS.csv"),
                                header = T, nrows = chunk_size, stringsAsFactors = FALSE)

col_names_CHARTEVENTS <- colnames(CHARTEVENTS_i)
CHARTEVENTS_i <- CHARTEVENTS_i  %>% 
  select(SUBJECT_ID, HADM_ID, ICUSTAY_ID, ITEMID, VALUENUM) %>%
  filter(ITEMID %in% D_ITEM_A_LINE, 
         !is.na(VALUENUM), 
         VALUENUM != 0,
         !is.na(ICUSTAY_ID)) %>% 
  select(ICUSTAY_ID) %>% 
  unique()

ICUSTAY_ID_A_Line <- CHARTEVENTS_i$ICUSTAY_ID
  



for (i in 2:20){
  skip_rows = chunk_size * (i-1)
  CHARTEVENTS_i = fread(here("data","raw_data","mimic-iii-clinical-database-1.4", "CHARTEVENTS.csv"), 
                                 header = F, nrows = chunk_size, stringsAsFactors = FALSE, 
                                 skip = skip_rows+1)
  colnames(CHARTEVENTS_i) = col_names_CHARTEVENTS
  CHARTEVENTS_i <- CHARTEVENTS_i  %>% 
    select(SUBJECT_ID, HADM_ID, ICUSTAY_ID, ITEMID, VALUENUM) %>%
    filter(ITEMID %in% D_ITEM_A_LINE, 
           !is.na(VALUENUM), 
           VALUENUM != 0,
           !is.na(ICUSTAY_ID)) %>% 
    select(ICUSTAY_ID) %>% 
    unique()
  ICUSTAY_ID_A_Line <- c(ICUSTAY_ID_A_Line, CHARTEVENTS_i$ICUSTAY_ID)
}


ICUSTAY_ID_A_Line <- unique(ICUSTAY_ID_A_Line)

length(ICUSTAY_ID_A_Line)

rm(CHARTEVENTS_i, col_names_CHARTEVENTS, chunk_size, D_ITEM_A_LINE, i, nrow_CHARTEVENTS, skip_rows)

save.image(file = here("data", "rdata", "with_incp_b", "ICUSTAY_ID_A_Line.RData"))
