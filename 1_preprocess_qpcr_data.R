#' run this script for each qPCR plate. The script may require inspection
#' of the raw data (2 QC checks), which is stored in Thermo's 
#' proprietary .eds format. This can be viewed through their cloud-based 
#' platform.

# Set-up ------------------------------------------------------------------

#' load libraries
library(tidyverse)
library(fs)
library(stringr)
library(magrittr)
library(data.table)

#' `enter current plate # here:`
curr_plate <- 12

#' set working directory
setwd("~/~r_projects/2023_kigar_cud_paper/raw_data/") 

#' there were 12 qPCR plates in the study. go through each individually
raw_data <- read_csv("210716_Addiction_plate12.csv") 

#' create new variable for file name at end based on above input
outfile <- paste0("plate", curr_plate, "_Addiction_qPCR_clean.csv")


# Wrangle -----------------------------------------------------------------
# drop extraneous columns
raw_data %<>% select(Well, Cq)

#' import plate template:
setwd("~/~r_projects/2023_kigar_cud_paper/data/")
template <- read_csv("Addiction_plate_key.csv")

#' join raw data & template
pcr <- inner_join(template, raw_data, by = "Well")

#' qPCR data outputs "Undetermined" when 
# no amplification detected in wells. Change to NA:
pcr$Cq <- replace(pcr$Cq, pcr$Cq=="Undetermined", NA)

#' coerce Cq column into numeric data. 
pcr$Cq <- as.numeric(pcr$Cq)

#' drop NAs:
pcr %<>% drop_na()

#' Organize rows by Sample:
pcr %<>% arrange(Sample)


# QC check #1 -------------------------------------------------------------

#' check no template controls did not amplify:
warning1 <- pcr %>% filter(Sample == "H2O" & Cq > 0)
print("Warning: amplification in negative control")
print(warning1)
#' if so, check that contamination is at least 2 orders of magnitude
#' less than sample amplification. if it is, proceed. if it isn't, exclude.


# QC check #2 -------------------------------------------------------------

#' Make new data table with aggregated means + sds for technical triplicates:
pcr.dt <- as.data.table(pcr)
pcr.dt <-  
  pcr.dt[, sapply(.SD, function(x) list(mean=mean(x), sd=sd(x))), 
         by=list(Sample, Target), .SDcols = 'Cq']

#' rename columns: 
pcr.dt  %<>% rename(Cq.avg = V1, Cq.sd = V2)

#' convert back to tibble
pcr.clean <- as_tibble(pcr.dt)

#' check for SD >= 0.5
warning2 <- pcr.clean %>% 
  filter(across(contains("Cq.sd")) >= 0.5 & Sample != "H2O") %>% 
  select(!contains("Cq.avg"))
print("Warning: standard deviation > 0.5")
print(warning2)
#' if so, visually inspect raw data to ensure no obvious technical errors
#' `if the sample is a plate or H2O control, ignore`
#' `if there is an obvious problem with one of the technical replicates:`
#' 1) duplicate the raw data file for the plate
#' 2) rename the duplicate file "210709_Addiction_plateX_edit.csv"
#' 3) rerun this script with the new data file
#' `if there is not an obvious problem with technical replicates, keep all`

# Create readable data table ----------------------------------------------

#' pivot data to convert targets to variables
wide <- pivot_wider(pcr.clean, 
                    names_from = Target, 
                    values_from = c(Cq.avg, Cq.sd))

#' create plate column to keep track of plate to plate variation
wide$plate <- curr_plate

#' move column to left
wide %<>% relocate(Sample, plate)

#' filter out plate controls 
wide %<>% filter(Sample != "H2O" & Sample !="plate_control")

#' Simplify column names:
wide %<>% rename_with(~ sub("Cq.avg_", "", .))
wide %<>% rename_with(~ sub("Cq.", "", .))

#' tidy 
wide$Sample <- as.numeric(wide$Sample)
wide %<>% arrange(Sample)

#' There were 193 samples total; samples were randomized and given a new #, 
#' 1-193, loaded in order across the 12 plates. This line will convert the 
#' Sample ID from being generic to the plate to matching the actual sample.
wide %<>% mutate_at(vars("Sample"), 
                    .funs = ~. + (curr_plate - 1) * 17)

# Create 3xHSK column -----------------------------------------------------
#' genes of interest will be normalized to an averaged housekeeping gene
#' in the next script.

#' average the three housekeeping genes (TBP, TMBIM4, GOLGA1) together:
avg <- wide %>% 
  mutate(`3xHSK` = select(.,TBP, TMBIM4, GOLGA1) %>% 
           rowMeans(na.rm = TRUE)) %>% relocate(Sample, plate,`3xHSK`)

# there were some extra controls on plate 12; this removes them.
avg %<>% filter(Sample < 194)


# Save data ---------------------------------------------------------------

#' write data to file:
setwd("~/~r_projects/2023_kigar_cud_paper/data/qPCR_data/")
avg %>% write_csv(file = outfile) 

# reset for next plate:
rm(list=ls())