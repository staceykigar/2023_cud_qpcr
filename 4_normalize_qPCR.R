#' this script is used to calculate the normalized (delta delta Ct) values
#' for the PCR gene expression data. 

# Set-up ------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(naniar)

setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
all <- read_csv(file = "3_Addiction_clean.csv") #349 observations

# convert things to factors: 
all$cDNA_ID <- factor(all$cDNA_ID)
all$cDNA_batch <- factor(all$cDNA_batch)
all$plate <- factor(all$plate)

# only want rows with qPCR data:
df <- all %>% drop_na(cDNA_ID) #173 observations

# only want qPCR data where DRD2 amplified:
df %<>% drop_na(DRD2) #167 observations

# drop columns unrelated to pcr analysis:
df %<>% select(Key, FLAG, cDNA_ID, cDNA_batch, plate,
               `3xHSK`, DRD2, DRD3, DRD4, COMT)


# Calculate ddCt ----------------------------------------------------------
# dCt:
df %<>% mutate(dDRD2 = DRD2 - `3xHSK`,
               dDRD3 = DRD3 - `3xHSK`,
               dDRD4 = DRD4 - `3xHSK`,
               dCOMT = COMT - `3xHSK`) 

# create new object with non-numeric data & dCt data only:
ddCt <- df %>% select(Key, FLAG, cDNA_ID, cDNA_batch, plate,
                      dDRD2, dDRD3, dDRD4, dCOMT)

# Create 2^(-x) function: 
dCt.fun <- function(x){return(2^(-x))}

# Create the ddCt values
ddCt  %<>%  mutate(across(where(is.numeric),
                          function(x) x- min(x, na.rm = TRUE)))

# 2^(-ddCt): generates a calibrator
ddCt  %<>%  mutate(across(where(is.numeric), dCt.fun)) 

# round the numbers to 2 decimal points:
ddCt  %<>%  mutate(across(where(is.numeric), round, 2))

# rename 2^(-ddCt) columns to be more descriptive:
ddCt %<>% rename(final_ddDRD2 = dDRD2,
                 final_ddDRD3 = dDRD3,
                 final_ddDRD4 = dDRD4,
                 final_ddCOMT = dCOMT)

# can't work with the calibrator for ddCt values, so convert to NA
ddCt  %<>%  replace_with_na(replace = 
                                list( final_ddDRD2 = 1, 
                                      final_ddDRD3 = 1, 
                                      final_ddDRD4 = 1, 
                                      final_ddCOMT  =1))


# Save data (merge with primary database) ---------------------------------

# merge 
all <- full_join(all, ddCt, 
                 by = c("Key","FLAG","cDNA_ID", "cDNA_batch", "plate"))

# save data 
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
all %>% write_csv(file = "4_Addiction_clean.csv")
