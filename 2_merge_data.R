#' There are 3 major sections to this script; 
#' 1) stitch results from 1_preprocess_qpcr_data script together
#' 2) merge demographic data from 3 different studies (data from K. Ersche)
#' 3) merge molecular data with demographic data

# Set-up ------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(stringr) 
library(naniar)
library(data.table)

# 1) Merge qPCR files ----------------------------------------------------------

#' set working directory & merge files together:
merge <- list.files(path="~/~r_projects/2023_kigar_cud_paper/data/qPCR_data/", 
                    pattern = ".*plate.*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

#' tidy
merge %<>% arrange(Sample)


# Save merged qPCR data ------------------------------------------------------

#' write data to file:
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
merge %>% write_csv(file = "2_Addiction_qPCR_merge.csv")

#' reset for next step:
rm(list=ls())

#' manually merged this dataframe with SNP data collected by Lorinda Turner
#' (creates "SNP_qPCR_manual.csv")

# Import demographic data -------------------------------------------------

#' set working directory
setwd("~/~r_projects/2023_kigar_cud_paper/raw_data/patient_data/")

#' load data
endo <- read_csv(file = "Endo_231.csv")
aging <- read_csv(file = "Ageing_299.csv")
atx <- read_csv(file = "ATX_377.csv")
atx.ctq <- read_csv(file = "ATX_CTQabuse.csv")
atx.bis <- read_csv(file = "ATX_BIS11_subscores.csv")
age.onset <- read_csv(file = "Study_299_age_stim.csv")
age.depend <- read_csv(file = "Aging299_dependence.csv")


# Wrangle & tidy study data ---------------------------------------

#' patient data tables have missing data input as `#NULL!`; convert to NA
endo %<>% replace_with_na_all(condition = ~.x == "#NULL!")                         
aging %<>% replace_with_na_all(condition = ~.x == "#NULL!")
atx %<>% replace_with_na_all(condition = ~.x == "#NULL!")

#' atx dataset has commas in the thousands place for WBIC column. fix:
atx$WBIC <- gsub(",", "", atx$WBIC)

#' Temporarily remove character variables:
endo %<>% select(-Study_name)
aging %<>% select(-Study_name)
atx %<>% select(-Study_name)

#' NB: where there were `NULL!` entries, R thinks the column is character data
#' Coerce all remaining columns into numeric data:
endo %<>% mutate_at(1:43, as.numeric)
aging %<>% mutate_at(1:37, as.numeric)
atx %<>% mutate_at(1:36, as.numeric)

#' Reassign study names::
endo$Study_name <- "endo"
aging$Study_name <- "aging"
atx$Study_name <- "atx"

#' `endo` dataset has some duplicate variables. Remove:
endo %<>% select(-BIS11,-SSSVtot, -OCDUS, -STIMyrs)

#' R imported an empty row into the aging dataset. Remove:
aging %<>% dplyr::filter(!is.na(Patient_ID))

#' The atx data has the ID numbers copied at the end 2x. Remove:
atx %<>% select(-ID_No, -36)


# Merge 'atx' with psychometric data  -----------------------------------------

#' `CTQ` make column names consistent, then merge:
atx.ctq %<>% rename(Patient_ID = ID_ATX)
atx <- full_join(atx, atx.ctq, by = "Patient_ID")
rm(atx.ctq)

#' `BIS11` alter atx.bis column names to match endo df because will 
#' eventually merge:
atx.bis %<>% rename(Patient_ID = ID,
                    BIS11_attention = BIS11_Attention,
                    BIS11_motor = BIS11_Motor, 
                    BIS11_nonplan = BIS11_Nonplan)

#' slight differences in age (calculated as date of scan - DOB) for atx and
#' atx.bis dataframes. Karen says to go with atx dataset ages. 
#' Also drop study code from atx.bis
atx.bis %<>% select(-`Study Code`, -Age)

# replace all the atx bis information with atx.bis df info; the preexisting
# bis11 total numbers are similar in the "ATX_377" spreadsheet I received 
# from Karen, but they're off by a couple numbers. The new BIS11 numbers 
# can be tallied from the sub-scores, so are more believable. Karen confirms;
# drop the old BIS11 scores and go with the new ones from atx.bis
atx %<>% select(-BIS11_total)

# then merge the data frames.
atx <- full_join(atx, atx.bis, by = c("Patient_ID", "Group"))

# the atx.bis data set has a unique sample with no data (1999); remove:
atx %<>% dplyr::filter(Patient_ID != "1999")

#' I confirmed with Karen that it's okay to round the atx ages to an integer:
atx$Age <- round(atx$Age, digits = 0)

rm(atx.bis)


# Merge 'aging' with age.onset data ---------------------------------------
#' rename columns to match aging dataframe:
age.onset %<>% rename(Patient_ID = ID,
                      Age_stimulants = Age_onset_stim,
                      Stimulants_yrs = Yr_stim)

#' merge
aging <- full_join(aging, age.onset, by = "Patient_ID")

rm(age.onset)

# Merge 'aging' with age.depend data ---------------------------------------
#' rename columns to match aging dataframe:
age.depend %<>% rename(Patient_ID = ID,
                       Dep_Alcohol = Alcohol,
                       Dep_Opiates = Opiate,
                       Dep_Cocaine = Cocaine) 

#' remove columns from aging that are incomplete:
aging %<>% select(-Dependence_alcohol, 
                  -Dependence_opiates, 
                  - Dependence_cocaine)

#' remove 'group' from age.depend because different data format:
age.depend %<>% select(-Group)

#' merge
aging <- full_join(aging, age.depend, by = "Patient_ID")

rm(age.depend)

# 2) Merge all three studies together ----------------------------------------
#' `endo & aging`
col.endo <- endo %>% colnames()
col.aging <- aging %>% colnames()

#' find unique columns to each dataframe; if actually the same rename
#' for consistency
unique.endo <- col.endo[!col.endo %in% col.aging]
unique.aging <- col.aging[!col.aging %in% col.endo]

#' rename columns:
endo.fix <- endo
endo.fix %<>%  rename(No_cig = No_cigarettes,
                      STAI_S_total = STAI_State_total,
                      STAI_T_total = STAI_Trait_total,
                      SSSV_Es_nodrug = SSSV_Es_nondrug,
                      SSSV_Dis_nodrug = SSSV_Dis_nondrug,
                      No_cig = No_cigarettes)

aging.fix <- aging
aging.fix %<>%  rename(BIS11_total = BIS_11_total,
                       BIS11_attention = BIS_11_attention,
                       BIS11_motor = BIS_11_motor,
                       BIS11_nonplan = BIS_11_nonplan,
                       OCDUS_total = OCDUS)

#' now find true bespoke columns to each dataframe
col.endo.fix <- endo.fix %>% colnames()
col.aging.fix <- aging.fix %>% colnames()

unique.endo.fix <- col.endo.fix[!col.endo.fix %in% col.aging.fix]
unique.aging.fix <- col.aging.fix[!col.aging.fix %in% col.endo.fix]

#' IV, MMT, and Subutex all unique to aging + are incomplete so drop:
aging.fix %<>% select(-IV, -MMT, -Subutex)

#' Add missing OCI variable to endo data set:
endo.fix$OCI_total <- NA
endo.fix$THC <- NA

#' Add missing variables to aging data set:
aging.fix$CTQabuse <- NA
aging.fix$AUDIT_level <- NA
aging.fix$Smoking <- NA
aging.fix$SSSV_tot_nondrug <- NA

# recheck no unique columns
col.endo.fix <- endo.fix %>% colnames()
col.aging.fix <- aging.fix %>% colnames()

unique.endo.fix <- col.endo.fix[!col.endo.fix %in% col.aging.fix]
unique.aging.fix <- col.aging.fix[!col.aging.fix %in% col.endo.fix]

# merge these two dataframes:
merge <- rbind(endo.fix, aging.fix)

rm(col.aging, col.aging.fix, col.endo, col.endo.fix,
   unique.aging, unique.aging.fix, unique.endo, unique.endo.fix,
   aging.fix, endo.fix, endo, aging)

#' `endo/aging merge + atx`
# find unique columns to each dataframe; if actually the same rename
# for consistency
col.atx <- atx %>% colnames()
col.merge <- merge %>% colnames()
unique.merge <- col.merge[!col.merge %in% col.atx]
unique.atx <- col.atx[!col.atx %in% col.merge]

# rename columns:
atx.fix <- atx
atx.fix %<>%  rename(WBIC_No = WBIC,
                     Dep_Alcohol = Dependence_alcohol,
                     Dep_Opiates = Dependence_opiates,
                     Dep_Cocaine = Dependence_Cocaine,
                     THC = D3_THC,
                     Yrs_Ed = Edu_yrs,
                     Stimulants_yrs = Yrs_stim,
                     AUDIT = AUDIT_total,
                     Ethnic = Ethnicity,
                     SSSV_tot_nondrug = SSSV_total_nodrug)

#' recheck unique columns:
col.fix <- atx.fix %>% colnames()
col.merge <- merge %>% colnames()
unique.merge <- col.merge[!col.merge %in% col.fix]
unique.fix <- col.fix[!col.fix %in% col.merge]

#' rs36024 unique to atx, and I believe it's gone into
#' a separate publication, so exclude:
atx.fix %<>% select(-rs36024)

#' Add missing variables to atx data set:
atx.fix$AUDIT_level <- NA
atx.fix$DAST_20 <- NA
atx.fix$Age_stimulants <- NA
atx.fix$Smoking <- NA

#' recheck unique columns:
col.fix <- atx.fix %>% colnames()
col.merge <- merge %>% colnames()
unique.merge <- col.merge[!col.merge %in% col.fix]
unique.fix <- col.fix[!col.fix %in% col.merge]

#' merge dataframes:
merge <- rbind(merge, atx.fix)

rm(atx, atx.fix, col.atx, col.fix, col.merge, 
   unique.atx, unique.fix, unique.merge)

# Wrangle & tidy merged participant data ------------------------------------
#' `AUDIT_level` has lots of NA data, but can be determined from AUDIT score:
merge %<>% mutate(AUDIT_risk = if_else(AUDIT <= 7, 1, 
                                    if_else(AUDIT <= 15, 2, 
                                            if_else(AUDIT <= 19, 
                                                    3, 4))))

# AUDIT_risk and AUDIT_level are the same thing, so fix:
merge <- merge %>% select(-AUDIT_level) %>% 
  rename(AUDIT_level = AUDIT_risk)

#' `Age_stimulants` has lots of NA data, but can be determined 
#' from subtracting years of stimulant use from current age:
merge %<>% mutate(Age_stim_temp = Age - Stimulants_yrs)
#' not a perfect match for what's already in Age_stimulants column,
#' but should be more accurate and allows for more data inclusion.
#' replace pre-existing with new:
merge <- merge %>% select(-Age_stimulants) %>% 
  rename(Age_stimulants = Age_stim_temp)

#' `SSSV_tot_nodrug` has lots of NA data - can it be determined 
#' from adding SSSV_Tas, SSSV_Bs, SSSV_Es_nodrug, and SSSV_Dis_nodrug?:
test <- merge %>% 
  mutate(tempSSSV = SSSV_Tas + SSSV_Bs + SSSV_Es_nodrug + SSSV_Dis_nodrug)
merge$SSSV_tot_nondrug
# doesn't match perfectly. Will drop the _nodrug columns.
merge %<>% select(-SSSV_Es_nodrug, -SSSV_Dis_nodrug, -SSSV_tot_nondrug)

rm(test)

#' `Smoking` - missing half the data, not sure how this is different than 
#' Smoker. Drop:
merge %<>% select(-Smoking)

#' `Create Key`
merge$Key <- str_sub(merge$Patient_ID, 1, 1)

merge["Key"][merge["Key"] == "1"] <- "control"
merge["Key"][merge["Key"] == "2"] <- "addicted"
merge["Key"][merge["Key"] == "3"] <- "relative"
merge["Key"][merge["Key"] == "4"] <- "user"

#' for the analysis we want to combine first degree relatives and controls,
#' and we will also drop the recreational users. Store as a new key, called
#' `Key.combo`
merge$Key.combo <- merge$Key
merge[merge$Key.combo=="relative","Key.combo"] <- "control"
merge[merge$Key.combo=="user","Key.combo"] <- NA

#' there are two participants in the addicted group without cocaine dependence 
#' (Both are from study ID 377; the Patient IDs are 2020 and 2073). 
#' Karen says they should be excluded. Will flag:
merge$FLAG <- F

setDT(merge)[Study_ID=="377" & Patient_ID =="2020", FLAG := T]
merge[Study_ID=="377" & Patient_ID =="2073", FLAG := T]

# Save participant data: --------------------------------------------------
#' reorganize variables to be more logical (drop group):
merge %<>% select(Key, Key.combo, FLAG, Patient_ID, Study_ID, Study_name, 
                  WBIC_No, Age, Gender, Ethnic, Hand, Yrs_Ed, NART,
                  DAST_20, Smoker, No_cig, THC, AUDIT, AUDIT_level, 
                  Dep_Alcohol, Dep_Opiates, Dep_Cocaine, OCDUS_total,
                  Stimulant_choice, Age_stimulants, Stimulants_yrs,
                  OCI_total, BDI_total, CTQabuse, STAI_S_total, STAI_T_total,
                  BIS11_attention, BIS11_motor, BIS11_nonplan, BIS11_total,
                  SSSV_Tas, SSSV_Es, SSSV_Dis, SSSV_Bs, SSSV_total)

#' sort observations by key:
merge %<>% arrange(Key)

#' save file
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
merge %>% write_csv(file = "2_Addiction_demographics_complete.csv")


# 3) Import molecular data & merge with demographics: --------------------------

#' set wd
setwd("~/~r_projects/2023_kigar_cud_paper/data/")

#' I manually merged the genotyping data from Lori and qPCR data in excel:
qpcr.snp <- read_csv(file = "SNP_qPCR_manual.csv")

#' merge the demographic and molecular data 
df.merge <- full_join(merge, qpcr.snp, by = c("Patient_ID", "Study_ID")) 

# Tidy complete dataframe: ------------------------------------------------

# 377_2020 & 377_2073 should be excluded, according to Karen:
df.merge  %<>%  filter(is.na(DNA_ID) | DNA_ID != "266")
df.merge  %<>%  filter(is.na(DNA_ID) | DNA_ID != "25")

# some samples in the qpcr.snp df didn't match up with merge df. flag
df.merge[is.na(Key), FLAG := T]

#' these 12 samples apparently can't be matched up to a key and must be
#' dropped :(
df.merge %<>% drop_na(Key)

#' 'notes' column from qpcr.snp indicates some results are questionable.
#' Add to flags.
df.merge[!is.na(notes), FLAG := T]

#' no immediate reason to exclude so leave in. Delete notes column as it's 
#' redundant with the FLAG column: 
df.merge %<>% select(-notes)


# Save molecular+demographic data -----------------------------------------
#' reorganize variables to be more logical (drop RIN as hardly any data):
df.merge %<>% select(1:40, Date_drawn, DNA_ID, 
                     rs4680, rs6280, rs6277, rs1800497, rs12364283, rs686,
                     cDNA_ID, ng_ul, cDNA_batch, plate, 
                     `3xHSK`, TBP, TMBIM4, GOLGA1, DRD2, DRD3, DRD4, COMT)


#' save file
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
df.merge %>% write_csv(file = "2_Addiction_demograph_molecular_complete.csv")
