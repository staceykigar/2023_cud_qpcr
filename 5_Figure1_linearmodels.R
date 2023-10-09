#' this script is used to generate Figure 1, SF1, SF2, and accompanying
#' stats
#' 
#' NB: different genes have different calibrator samples used in the 
#' ddCt calculation. These have been converted to NA data, so each 
#' gene may have varying NA sample data.

# Set-up ------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(ggsignif)
library(ggpubr)
library(ggplot2)
library(dlookr)
library("PerformanceAnalytics")
library(flextable) # conflicts with 'font' in ggpubr
library(rcompanion)
library(car)
library(outliers)
library(report)

setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
all <- read_csv(file = "4_Addiction_clean.csv") 

# convert things to factors: 
all$Key <- factor(all$Key)
all$Patient_ID <- factor(all$Patient_ID)
all$cDNA_batch <- factor(all$cDNA_batch)
all$cDNA_ID <- factor(all$cDNA_ID)
all$plate <- factor(all$plate)
all$Study_ID <- factor(all$Study_ID)
all$Gender <- factor(all$Gender)
all$Ethnic <- factor(all$Ethnic)
all$Smoker <- factor(all$Smoker)
all$Dep_Alcohol <- factor(all$Dep_Alcohol)
all$Dep_Opiates <- factor(all$Dep_Opiates)

# create new df where variables will be decoded:
df <- all

## decode variables:
# gender, 1 = male, 2 = female
levels(df$Gender)[1:2] <- c("male", "female")

# ethnicity, 1 = white, 2 = black, 3 = asian, 4 = indian, 5 = hispanic, 6 = other
# BUT simplify into white/non-white:
levels(df$Ethnic)[1:6] <- c("white", 
                            "non-white", 
                            "non-white", 
                            "non-white", 
                            "non-white", 
                            "non-white")

# smoker, 0 = non-smoker, 1 = smoker, 2 = former smoker
# BUT at Karen's request, simplify former smokers and non-smokers:
levels(df$Smoker)[1:3] <- c("non-smoker", 
                            "smoker",
                            "non-smoker")

# alcohol dependence, 0 = no, 1 = yes
levels(df$Dep_Alcohol)[1:2] <- c("no", 
                            "yes")

# opiate dependence, 0 = no, 1 = yes
levels(df$Dep_Opiates)[1:2] <- c("no", 
                                 "yes")

# only want rows with qPCR data; use cDNA ID as proxy.
# only want variables related to qPCR data & regression.
target <- df %>% drop_na(cDNA_ID) %>% 
  dplyr::select(Key, FLAG, Patient_ID, Study_ID, 
         Age, Gender, Ethnic, Smoker, Dep_Alcohol, Dep_Opiates,
         cDNA_ID, cDNA_batch, plate, 
         final_ddDRD2, final_ddDRD3, final_ddDRD4, final_ddCOMT) 

# make function for looking at skewness and kurtosis
sk <- function(x)c(skew = skewness(x, type = 2),
                   kurt = kurtosis(x, type = 2))


# Plot raw data   -----------------------------------------------------------
# set wd for image output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/sup_fig1/")

# pivot to long format:
long <- target %>% dplyr::select(1:4,12:17) %>% 
  rename(DRD2 = final_ddDRD2, 
         DRD3 = final_ddDRD3,
         DRD4 = final_ddDRD4, 
         COMT = final_ddCOMT) %>% 
  pivot_longer(-c(Key, FLAG, Patient_ID, Study_ID, cDNA_batch, plate), 
               values_to = "mRNA", 
               names_to = "target") 


# make graph:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Key)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf1_addiction_violin_qpcr.pdf")

#' make graph - `highlight flags`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=FLAG)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sfX_addiction_violin_qpcr.pdf")
# nothing worrying

#' make graph - `visualize plate`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=plate)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) 

ggsave("sf1_addiction_violin_qpcr_plate.pdf")
# no clusters to my eye

#' make graph - `visualize cDNA batch`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=cDNA_batch)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) 

ggsave("sf1_addiction_violin_qpcr_cDNAbatch.pdf")
# there does appear to be some clustering here

#' make graph - `visualize study ID`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Study_ID)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) 

ggsave("sf1_addiction_violin_qpcr_studyID.pdf")
# maybe some clusters?

# Assess skew and kurtosis.  -----------------------------------
#' `check outliers`
hc <- target %>% filter(Key == "control")
cud <- target %>% filter(Key == "CUD")

grubbs.test(target$final_ddDRD2) #one outlier, 0.92*
grubbs.test(hc$final_ddDRD2) #one outlier, 0.78*
grubbs.test(cud$final_ddDRD2) #no outlier

grubbs.test(target$final_ddDRD3) #one outlier, 0.31**
grubbs.test(hc$final_ddDRD3) #one outlier, 0.3**
grubbs.test(cud$final_ddDRD3) #one outlier, 0.31**

grubbs.test(target$final_ddDRD4) #one outlier, 0.68***
grubbs.test(hc$final_ddDRD4) #no outlier
grubbs.test(cud$final_ddDRD4) #one outlier, 0.68****

grubbs.test(target$final_ddCOMT) #one outlier, 0.54****
grubbs.test(hc$final_ddCOMT) #one outlier, 0.54****
grubbs.test(cud$final_ddCOMT) #one outlier, 0.4*

#' `check normality:`
target %>% dplyr::select(Key,14:17) %>% group_by(Key) %>% 
  normality() %>% 
  mutate(across(where(is.numeric), round, 4)) %>% 
  regulartable() 
# everything is non-normal.

#' `DRD2`
# look at distribution of data:
hist(target$final_ddDRD2) # positive skew to data

# quantify skew/kurtosis
target %>% drop_na(final_ddDRD2) %>%
  aggregate(final_ddDRD2~Key, FUN = sk)
#         Key final_ddDRD2.skew final_ddDRD2.kurt
# 1 control        1.21692152        1.10407128
# 2     CUD        0.74056359        0.01861309
# control is highly skewed, CUD moderately skewed
# kurtosis looks okay

#' `DRD3`
# look at distribution of data:
hist(target$final_ddDRD3) # positive skew to data

# quantify skew/kurtosis
target %>% drop_na(final_ddDRD3) %>%
  aggregate(final_ddDRD3~Key, FUN = sk)
#       Key final_ddDRD3.skew final_ddDRD3.kurt
# 1 control          1.636290          3.391401
# 2     CUD          1.344041          1.831158
# skew is high
# kurtosis is high in control

#' `DRD4`
# look at distribution of data:
hist(target$final_ddDRD4) # positive skew to data

# quantify skew/kurtosis
target %>% drop_na(final_ddDRD4) %>%
  aggregate(final_ddDRD4~Key, FUN = sk)
#        Key final_ddDRD4.skew final_ddDRD4.kurt
# 1 control          1.388197          1.499898
# 2     CUD          1.497831          3.668313
# skew is high
# kurtosis is high in CUD

#' `COMT`
# look at distribution of data:
hist(target$final_ddCOMT) # positive skew to data

# quantify skew/kurtosis
target %>% drop_na(final_ddCOMT) %>%
  aggregate(final_ddCOMT~Key, FUN = sk)
#       Key final_ddCOMT.skew final_ddCOMT.kurt
# 1 control          2.694783          9.600482
# 2     CUD          1.052324          1.230659
# skew is very high
# kurtosis is very high


# Data transformation: sqrt -----------------------------------------------------
#' all but the CUD disorder group from DRD2 look highly skewed. Try 
#' transformations in order from least to most corrective until 
#' normal distribution (sqrt < ln < log10)

# first select the columns to transform:
cols <- target %>% dplyr::select(14:17) %>% colnames()
# store as new dataframe because original column data will be overwritten
target2 <- target %>% mutate_at(cols, sqrt)

# diagnostics
par(mfrow=c(1,2))

#' `check normality assumptions - sqrt:`
target2 %>% dplyr::select(Key,14:17) %>% group_by(Key) %>% 
  normality() %>% 
  mutate(across(where(is.numeric), round, 4)) %>% 
  regulartable() 
# DRD2 passes, otherwise only DRD4 CUD pass normality.

#' `DRD2 - sqrt`
# look at distribution of data before vs after:
hist(target$final_ddDRD2)
hist(target2$final_ddDRD2) # definitely improved

#' `DRD3 - sqrt`
# look at distribution of data before vs after:
hist(target$final_ddDRD3)
hist(target2$final_ddDRD3) # still a little right skewed

#' `DRD4 - sqrt`
# look at distribution of data before vs after:
hist(target$final_ddDRD4)
hist(target2$final_ddDRD4) # still a little right skewed

#' `COMT - sqrt`
# look at distribution of data before vs after:
hist(target$final_ddCOMT)
hist(target2$final_ddCOMT) # still pretty right skewed

# set wd for image output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/QC/")

# Reassess skew & kurtosis after sqrt transformation ---------------------------
#' `check outliers`
con <- target2 %>% filter(Key == "control")
cud <- target2 %>% filter(Key == "CUD")

grubbs.test(target2$final_ddDRD2) #no outlier
grubbs.test(hc$final_ddDRD2) #no outlier
grubbs.test(cud$final_ddDRD2) #no outlier

grubbs.test(target2$final_ddDRD3) #no outlier
grubbs.test(hc$final_ddDRD3) #one outlier, 0.547722557505166*
grubbs.test(cud$final_ddDRD3) #no outlier

grubbs.test(target2$final_ddDRD4) #no outlier
grubbs.test(hc$final_ddDRD4) #no outlier
grubbs.test(cud$final_ddDRD4) #one outlier, 0.824621125123532*

grubbs.test(target2$final_ddCOMT) #one outlier, 0.734846922834953***
grubbs.test(hc$final_ddCOMT) #one outlier, 0.734846922834953****
grubbs.test(cud$final_ddCOMT) #no outlier

# identify outliers:
target2 %>% filter(final_ddDRD3 == "0.547722557505166") 
#PatientID 1055, StudyID 377
target2 %>% filter(final_ddDRD4 == "0.824621125123532") 
#PatientID 2008, StudyID 377
target2 %>% filter(final_ddCOMT == "0.734846922834953") 
#PatientID 1055, StudyID 377

# diagnostics
par(mfrow=c(2,2))

#' `DRD2 - sqrt`
#' quantify skew/kurtosis 
target2 %>% dplyr::select(Key, final_ddDRD2) %>% 
  drop_na(final_ddDRD2) %>% 
  aggregate(final_ddDRD2~Key, FUN = sk)
#       Key final_ddDRD2.skew final_ddDRD2.kurt
# 1 control        0.47858579       -0.13550191
# 2     CUD        0.02129751       -0.45158314
# looks good

#' raw
cud <- target %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddDRD2) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # slightly skewed

con <- target %>% filter(Key == "control") %>% 
  dplyr::select(final_ddDRD2) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # quite skewed

#' transformed
cud <- target2 %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddDRD2) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # looks okay

con <- target2 %>% filter(Key == "control") %>% 
  dplyr::select(final_ddDRD2) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # looks okay
# manually saved as drd2_cw_rawCUD_rawCON_sqCON_sqCUD.pdf

#' `DRD3 - sqrt` outlier included
#' quantify skew/kurtosis 
target2 %>% dplyr::select(Key, final_ddDRD3) %>% 
  drop_na(final_ddDRD3) %>% 
  aggregate(final_ddDRD3~Key, FUN = sk)
#        Key final_ddDRD3.skew final_ddDRD3.kurt
# 1 control         0.7022237         0.9108573
# 2     CUD         0.7495966         0.1146639
# both are moderately skewed. kurtosis is okay.

#' raw
cud <- target %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddDRD3) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # definitely + skewed; lots of values that are =?

con <- target %>% filter(Key == "control") %>% 
  dplyr::select(final_ddDRD3) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # definitely + skewed; lots of values that are =?

#' transformed
cud <- target2 %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddDRD3) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # still a little wonky but better; still lots of values that are =?

con <- target2 %>% filter(Key == "control") %>% 
  dplyr::select(final_ddDRD3) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # still a little wonky but better; still lots of values that are =?
# manually saved as drd3_cw_rawCUD_rawCON_sqCON_sqCUD.pdf

#' `DRD4 - sqrt` outlier included
#' quantify skew/kurtosis 
target2 %>% dplyr::select(Key, final_ddDRD4) %>% 
  drop_na(final_ddDRD4) %>% 
  aggregate(final_ddDRD4~Key, FUN = sk)
# Key final_ddDRD4.skew final_ddDRD4.kurt
# 1 control         0.7172188         0.2808590
# 2     CUD         0.4488745         0.2635109
# moderately skewed, kurtosis fine

#' raw
cud <- target %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddDRD4) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # +skew

con <- target %>% filter(Key == "control") %>% 
  dplyr::select(final_ddDRD4) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # +skew

#' transformed
cud <- target2 %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddDRD4) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # much better except for extreme values on ends

con <- target2 %>% filter(Key == "control") %>% 
  dplyr::select(final_ddDRD4) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # wonky at the top but a bit better
# manually saved as drd4_cw_rawCUD_rawCON_sqCON_sqCUD.pdf

#' `COMT - sqrt` outlier included
#' quantify skew/kurtosis 
target2 %>% dplyr::select(Key, final_ddCOMT) %>% 
  drop_na(final_ddCOMT) %>% 
  aggregate(final_ddCOMT~Key, FUN = sk)
# Key final_ddCOMT.skew final_ddCOMT.kurt
# 1 control         1.8597733         5.4515573
# 2     CUD         0.5841473         0.2390007
# controls are highly skewed, CUD is good. kurtosis is okay

#' raw
cud <- target %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddCOMT) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # +skew, clumping

con <- target %>% filter(Key == "control") %>% 
  dplyr::select(final_ddCOMT) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # looks quite bad; nearly flat line.

#' transformed
cud <- target2 %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddCOMT) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # much better, still some clumping

con <- target2 %>% filter(Key == "control") %>% 
  dplyr::select(final_ddCOMT) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # still quite bad/nearly flat line
# manually saved as comt_cw_rawCUD_rawCON_sqCON_sqCUD.pdf

# Close double plotting window
dev.off()

# Data transformations: ln ------------------------------------------------
#' Try  harsher transformation for variables
#' (sqrt < *ln < log10)

# store as new dataframe because original column data will be overwritten
target3 <- target %>% mutate_at(cols, log) 

# diagnostics
par(mfrow=c(1,3))

#' `check normality assumptions - ln:`
target3 %>% dplyr::select(Key,14:17) %>% group_by(Key) %>% 
  normality() %>% 
  mutate(across(where(is.numeric), round, 4)) %>% 
  regulartable() # DRD2 fails now. DRD3 CUD fails. COMT ctrl fails.

#' `DRD2 - ln`
# look at distribution of data before vs after:
hist(target$final_ddDRD2)
hist(target2$final_ddDRD2) # definitely improved
hist(target3$final_ddDRD2) # flips to -skew
# manually saved as drd2_raw_sqrt_ln.pdf

#' `DRD3 - ln` 
# look at distribution of data before and after:
hist(target$final_ddDRD3)
hist(target2$final_ddDRD3) #better but still slightly +skewed
hist(target3$final_ddDRD3) # over-corrected/slightly -skewed
# manually saved as drd3_raw_sqrt_ln.pdf

#' `DRD4 - ln`
# look at distribution of data before and after:
hist(target$final_ddDRD4)
hist(target2$final_ddDRD4) # slightly under corrected/+skewed
hist(target3$final_ddDRD4) # looks better
# manually saved as drd4_raw_sqrt_ln.pdf

#' `COMT - ln`
# look at distribution of data before and after:
hist(target$final_ddCOMT)
hist(target2$final_ddCOMT) # slightly under corrected/+skewed
hist(target3$final_ddCOMT) # looks better
# manually saved as comt_raw_sqrt_ln.pdf

# Reassess skew & kurtosis after ln transformation ---------------------------
#' `check outliers`
con <- target3 %>% filter(Key == "control")
cud <- target3 %>% filter(Key == "CUD")

grubbs.test(target3$final_ddDRD2) #one outlier, -4.60517018598809**
grubbs.test(hc$final_ddDRD2) #one outlier, -4.60517018598809**
grubbs.test(cud$final_ddDRD2) #one outlier, -3.91202300542815*

grubbs.test(target3$final_ddDRD3) #one outlier, -4.60517018598809*
grubbs.test(hc$final_ddDRD3) #one outlier, -4.60517018598809*
grubbs.test(cud$final_ddDRD3) #no outlier

grubbs.test(target3$final_ddDRD4) #no outlier
grubbs.test(hc$final_ddDRD4) #no outlier
grubbs.test(cud$final_ddDRD4) #no outlier

grubbs.test(target3$final_ddCOMT) #one outlier, -0.616186139423817*
grubbs.test(hc$final_ddCOMT) #one outlier, -0.616186139423817**
grubbs.test(cud$final_ddCOMT) #no outlier

# identify outliers:
target3 %>% filter(final_ddDRD2 == "-4.60517018598809") # odd
#PatientID 1046, StudyID 377
target3 %>% filter(final_ddDRD2 == "-3.91202300542815") 
#PatientID 2051, StudyID 231
target3 %>% filter(final_ddDRD3 == "-4.60517018598809") # odd
#PatientID 1047, StudyID 231
target3 %>% filter(final_ddCOMT == "-0.616186139423817") 
#PatientID 1055, StudyID 377

# diagnostics
par(mfrow=c(1,2))

#' `DRD2 - ln` outliers included
#' quantify skew/kurtosis 
target3 %>% dplyr::select(Key, final_ddDRD2) %>% 
  drop_na(final_ddDRD2) %>% 
  aggregate(final_ddDRD2~Key, FUN = sk)
# Key final_ddDRD2.skew final_ddDRD2.kurt
# 1 control        -0.7052135         1.5850139
# 2     CUD        -1.0774408         1.5048878
# strong - skew

#' transformed
cud <- target3 %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddDRD2) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # okay other than the ends

con <- target3 %>% filter(Key == "control") %>% 
  dplyr::select(final_ddDRD2) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # looks compressed
# manually saved as drd2_lnCUD_lnCON.pdf

#' `DRD3 - ln` outlier included
#' quantify skew/kurtosis 
target3 %>% dplyr::select(Key, final_ddDRD3) %>% 
  drop_na(final_ddDRD3) %>% 
  aggregate(final_ddDRD3~Key, FUN = sk)
# Key final_ddDRD3.skew final_ddDRD3.kurt
# 1 control        -0.3944213         1.0347684
# 2     CUD         0.2267515        -0.7021493
# looks good

#' transformed
cud <- target3 %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddDRD3) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # still a little wonky but better; still lots of values that are =?

con <- target3 %>% filter(Key == "control") %>% 
  dplyr::select(final_ddDRD3) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # still a little wonky but better; still lots of values that are =?
# manually saved as drd3_lnCUD_lnCON.pdf

#' `DRD4 - ln` 
#' quantify skew/kurtosis 
target3 %>% dplyr::select(Key, final_ddDRD4) %>% 
  drop_na(final_ddDRD4) %>% 
  aggregate(final_ddDRD4~Key, FUN = sk)
# Key final_ddDRD4.skew final_ddDRD4.kurt
# 1 control       -0.11190563        0.07204777
# 2     CUD       -0.43425524       -0.22577416
# looks fine

#' transformed
cud <- target3 %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddDRD4) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # slight - skew, ends are wonky

con <- target3 %>% filter(Key == "control") %>% 
  dplyr::select(final_ddDRD4) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # ends are wonky and clumping of values
# manually saved as drd4_lnCUD_lnCON.pdf

#' `COMT - ln` outlier included
#' quantify skew/kurtosis 
target3 %>% dplyr::select(Key, final_ddCOMT) %>% 
  drop_na(final_ddCOMT) %>% 
  aggregate(final_ddCOMT~Key, FUN = sk)
# Key final_ddCOMT.skew final_ddCOMT.kurt
# 1 control         1.0589683         2.6561142
# 2     CUD         0.1354085        -0.2232968
# controls are still highly skewed, CUD is good. kurtosis is okay

#' transformed
cud <- target3 %>% filter(Key == "CUD") %>% 
  dplyr::select(final_ddCOMT) %>% drop_na() %>% pull()
qqnorm(cud)
qqline(cud) # much better, still some clumping

con <- target3 %>% filter(Key == "control") %>% 
  dplyr::select(final_ddCOMT) %>% drop_na() %>% pull()
qqnorm(con)
qqline(con) # still quite bad/clumping/wonky ends
# manually saved as comt_lnCUD_lnCON.pdf

# Close double plotting window
dev.off()

rm(con, cud)


# Assess independence of variables ----------------------------------------
#' want to regress out batch variables; check they aren't confounded.

table(target$plate, target$cDNA_batch) # THESE ARE CONFOUNDED
#     1  2  3
# 1  15  0  0
# 2  15  0  0
# 3  14  0  0
# 4  12  3  0
# 5   0 16  0
# 6   0 12  0
# 7   0 17  0
# 8   0 11  6
# 9   0  0 15
# 10  0  0 16
# 11  0  0 16
# 12  0  0  5

table(target$plate, target$Study_ID)
#    231 299 377
# 1    4   6   5
# 2    5   5   5
# 3    8   4   2
# 4    3   2  10
# 5    4   5   7
# 6    0   5   7
# 7    4   3  10
# 8    5   2  10
# 9    5   6   4
# 10   4   9   3
# 11   4   7   5
# 12   0   0   5

table(target$cDNA_batch, target$Study_ID)
#    231 299 377
# 1  19  16  21
# 2  11  15  33
# 3  16  23  19


#' `note: cDNA batch and plate are totally confounded/can't be co-adjusted`

# Descriptive stats: DRD2 -------------------------------------------
target %>% filter(Key == "control") %>% drop_na(final_ddDRD2) #86
target %>% filter(Key == "CUD") %>% drop_na(final_ddDRD2) #80

#' `raw`
#' check equal variance - data is non-normal
leveneTest(final_ddDRD2 ~ Key, data = target,
           na.action = na.exclude) #p-value = 0.01771 *

#' Welch's t test: raw data: normal-ish b/c robust to lots of samples, 
#' unequal variance
t.test(final_ddDRD2 ~ Key, data = target, 
            alternative = "two.sided", var.equal = FALSE,
       na.action = na.exclude)
# t = -3.2921, df = 150.22, p-value = 0.00124

# save p value for later
p.drd2.ttest.raw <- 0.00124

#' `sqrt transformed`
#' check equal variance: data is normal
bartlett.test(final_ddDRD2 ~ Key, data = target2,
              na.action = na.exclude) #p-value = 0.2064

#' Student's t test: transformed data: normal, equal variance
t.test(final_ddDRD2 ~ Key, data = target2, 
       alternative = "two.sided", var.equal = TRUE,
       na.action = na.exclude)
# t = -3.2432, df = 164, p-value = 0.001432

# save p value for later
p.drd2.ttest.sq <- 0.001432

#' `ln transformed`
#' check equal variance - data is non-normal
leveneTest(final_ddDRD2 ~ Key, data = target3,
           na.action = na.exclude) #p-value = 0.742

#' Student's t test: raw data: normal-ish b/c robust to lots of samples, 
#' equal variance
t.test(final_ddDRD2 ~ Key, data = target3, 
       alternative = "two.sided", var.equal = T,
       na.action = na.exclude)
# t = -2.8464, df = 164, p-value = 0.004987

# save p value for later
p.drd2.ttest.ln <- 0.004987


# Mixed linear regression: DRD2 ------------------------------------------
#' `sqrt`
batch.lm <- lm(final_ddDRD2 ~ Key + plate + Study_ID, 
                 data = target2, na.action = na.exclude)
summary(batch.lm) #explains 11.4% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.40462 -0.10382 -0.00512  0.09597  0.47686 
# 
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.375054   0.048264   7.771 1.11e-12 ***
# KeyCUD       0.070867   0.026834   2.641  0.00914 ** 
# plate2       0.055316   0.060218   0.919  0.35977    
# plate3       0.181895   0.063765   2.853  0.00494 ** 
# plate4       0.076251   0.060933   1.251  0.21273    
# plate5       0.131663   0.060603   2.173  0.03138 *  
# plate6      -0.002249   0.066024  -0.034  0.97287    
# plate7       0.139002   0.059027   2.355  0.01981 *  
# plate8       0.059061   0.060750   0.972  0.33251    
# plate9       0.082785   0.061223   1.352  0.17834    
# plate10      0.065853   0.059585   1.105  0.27084    
# plate11      0.074391   0.060796   1.224  0.22300    
# plate12     -0.009760   0.087186  -0.112  0.91102    
# Study_ID299 -0.027800   0.035453  -0.784  0.43418    
# Study_ID377  0.055175   0.034597   1.595  0.11284    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1647 on 151 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.1887,	Adjusted R-squared:  0.1135 
# F-statistic: 2.509 on 14 and 151 DF,  p-value: 0.003077
report(batch.lm)

# save p value for later
p.drd2.lm.sq <- 0.00914

# diagnostics on the residuals
par(mfrow=c(2,2))
plot(batch.lm) # looks pretty good; 4 points have more leverage tho.
# manually saved as drd2_sqrt_lm_batch.pdf

#' `ln`
batch.lm <- lm(final_ddDRD2 ~ Key + plate + Study_ID, 
               data = target3, na.action = na.exclude)
summary(batch.lm) #explains 11.4% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.97612 -0.35224  0.03616  0.47358  1.68446 
# 
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.00747    0.21744  -9.232 2.26e-16 ***
# KeyCUD       0.28196    0.12089   2.332   0.0210 *  
# plate2       0.21629    0.27130   0.797   0.4266    
# plate3       0.68999    0.28728   2.402   0.0175 *  
# plate4       0.31207    0.27452   1.137   0.2574    
# plate5       0.56817    0.27303   2.081   0.0391 *  
# plate6       0.01707    0.29745   0.057   0.9543    
# plate7       0.53148    0.26593   1.999   0.0475 *  
# plate8       0.19504    0.27370   0.713   0.4772    
# plate9       0.32797    0.27583   1.189   0.2363    
# plate10      0.23328    0.26845   0.869   0.3862    
# plate11      0.13667    0.27390   0.499   0.6185    
# plate12      0.01682    0.39280   0.043   0.9659    
# Study_ID299 -0.12049    0.15972  -0.754   0.4518    
# Study_ID377  0.24175    0.15587   1.551   0.1230    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.742 on 151 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.1609,	Adjusted R-squared:  0.08314 
# F-statistic: 2.069 on 14 and 151 DF,  p-value: 0.01645

# save p value for later
p.drd2.lm.ln <- 0.0210

# diagnostics on the residuals
par(mfrow=c(2,2))
plot(batch.lm) # residuals look okay but QQ is weird; 4 points have more leverage tho.
# manually saved as drd2_ln_lm_batch.pdf

# add residuals to new column:
target2$ddDRD2.res <- residuals(lm(final_ddDRD2 ~ 
       plate + Study_ID, data = target2, 
       na.action = na.exclude))

# add residuals to new column:
target3$ddDRD2.res <- residuals(lm(final_ddDRD2 ~ 
        plate + Study_ID, data = target3, 
        na.action = na.exclude))

# Close double plotting window
dev.off()

# Sensitivity analysis: DRD2 (sqrt) ----------------------------------------------
par(mfrow=c(2,2))

#' `regress out batch variables & other covariates (-Ethnic & -Smoker):`
batch.lm <- lm(final_ddDRD2 ~ 
                 Key + plate + Study_ID +
                 Gender + Age, data = target2, 
               na.action = na.exclude)
summary(batch.lm) # explains ~12% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.38618 -0.10596 -0.00651  0.10568  0.47655 
# 
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.4386981  0.0727866   6.027 1.25e-08 ***
# KeyCUD        0.0642195  0.0270325   2.376   0.0188 *  
# plate2        0.0287786  0.0620098   0.464   0.6433    
# plate3        0.1603301  0.0651842   2.460   0.0151 *  
# plate4        0.0649447  0.0610941   1.063   0.2895    
# plate5        0.1166302  0.0610394   1.911   0.0580 .  
# plate6       -0.0220684  0.0672688  -0.328   0.7433    
# plate7        0.1175682  0.0602735   1.951   0.0530 .  
# plate8        0.0394805  0.0616275   0.641   0.5227    
# plate9        0.0643921  0.0619574   1.039   0.3004    
# plate10       0.0586951  0.0600087   0.978   0.3296    
# plate11       0.0703054  0.0606508   1.159   0.2482    
# plate12      -0.0066593  0.0872086  -0.076   0.9392    
# Study_ID299  -0.0521395  0.0390563  -1.335   0.1839    
# Study_ID377   0.0373222  0.0379966   0.982   0.3276    
# Genderfemale -0.0870076  0.0511029  -1.703   0.0907 .  
# Age          -0.0006166  0.0014024  -0.440   0.6608    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1642 on 149 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.2045,	Adjusted R-squared:  0.1191 
# F-statistic: 2.395 on 16 and 149 DF,  p-value: 0.003233

# diagnostics on the residuals
plot(batch.lm) # looks okay, one point with the most leverage.
# manually saved as drd2_sqrt_lm_batch_age_sex.pdf

#' `regress out batch variables & other covariates (-Smoker):`
batch.lm <- lm(final_ddDRD2 ~ 
                 Key + plate + Study_ID +
                 Gender + Age + Ethnic, data = target2)
summary(batch.lm) # explains 12.8% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.36736 -0.11302 -0.00574  0.09771  0.44292 
# 
# Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.4353361  0.0722925   6.022 1.32e-08 ***
# KeyCUD           0.0489442  0.0278824   1.755   0.0813 .  
# plate2           0.0258611  0.0614472   0.421   0.6745    
# plate3           0.1682262  0.0646790   2.601   0.0102 *  
# plate4           0.0765807  0.0608399   1.259   0.2101    
# plate5           0.1213092  0.0605144   2.005   0.0468 *  
# plate6          -0.0193938  0.0666816  -0.291   0.7716    
# plate7           0.1108462  0.0598979   1.851   0.0662 .  
# plate8           0.0441668  0.0610877   0.723   0.4708    
# plate9           0.0874703  0.0626340   1.397   0.1647    
# plate10          0.0575592  0.0594527   0.968   0.3346    
# plate11          0.0713253  0.0600824   1.187   0.2371    
# plate12          0.0147543  0.0870913   0.169   0.8657    
# Study_ID299     -0.0418920  0.0389940  -1.074   0.2844    
# Study_ID377      0.0374854  0.0376906   0.995   0.3216    
# Genderfemale    -0.0870212  0.0506235  -1.719   0.0877 .  
# Age             -0.0008195  0.0013948  -0.588   0.5578    
# Ethnicnon-white  0.0585934  0.0339309   1.727   0.0863 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1626 on 147 degrees of freedom
# (8 observations deleted due to missingness)
# Multiple R-squared:  0.2179,	Adjusted R-squared:  0.1275 
# F-statistic: 2.409 on 17 and 147 DF,  p-value: 0.002522

# diagnostics on the residuals
par(mfrow=c(2,2))
plot(batch.lm) # residuals starting to creep up in scale-location
# manually saved as drd2_sqrt_lm_batch_age_sex_eth.pdf

#' `regress out batch variables & other covariates (-Ethnic):`
batch.lm <- lm(final_ddDRD2 ~ 
                 Key + plate + Study_ID +
                 Gender + Age + Smoker, data = target2)
summary(batch.lm) # explains 11.8% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.38081 -0.10641 -0.01166  0.10087  0.48501 
# 
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.4352209  0.0729500   5.966 1.72e-08 ***
# KeyCUD        0.0361323  0.0418894   0.863   0.3898    
# plate2        0.0234299  0.0623556   0.376   0.7076    
# plate3        0.1577759  0.0652991   2.416   0.0169 *  
# plate4        0.0593943  0.0614668   0.966   0.3355    
# plate5        0.1165176  0.0610864   1.907   0.0584 .  
# plate6       -0.0288606  0.0677633  -0.426   0.6708    
# plate7        0.1160133  0.0603457   1.922   0.0565 .  
# plate8        0.0315721  0.0623288   0.507   0.6132    
# plate9        0.0595103  0.0622537   0.956   0.3407    
# plate10       0.0552648  0.0601817   0.918   0.3600    
# plate11       0.0651924  0.0609760   1.069   0.2867    
# plate12      -0.0021323  0.0874278  -0.024   0.9806    
# Study_ID299  -0.0467799  0.0395598  -1.183   0.2389    
# Study_ID377   0.0422756  0.0384418   1.100   0.2732    
# Genderfemale -0.0874446  0.0511446  -1.710   0.0894 .  
# Age          -0.0006939  0.0014062  -0.493   0.6224    
# Smokersmoker  0.0359588  0.0409450   0.878   0.3812    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1643 on 148 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.2087,	Adjusted R-squared:  0.1178 
# F-statistic: 2.296 on 17 and 148 DF,  p-value: 0.004131

# diagnostics on the residuals
par(mfrow=c(2,2))
plot(batch.lm) # looks okay
# manually saved as drd2_sqrt_lm_batch_age_sex_smoke.pdf

# Close double plotting window
dev.off()

# Descriptive stats: DRD3 -------------------------------------------

#' `raw`
#' check equal variance - data is non-normal
leveneTest(final_ddDRD3 ~ Key, data = target,
           na.action = na.exclude) #p-value = 0.4976

#' Student's t test: raw data: normal-ish b/c robust to lots of samples, 
#' equal variance
t.test(final_ddDRD3 ~ Key, data = target, 
       alternative = "two.sided", var.equal = T,
       na.action = na.exclude)
# t = -1.4168, df = 164, p-value = 0.1584

# save p value for later
p.drd3.ttest.raw <- 0.1584

#' `sqrt transformed`
#' check equal variance: data is non-normal
leveneTest(final_ddDRD3 ~ Key, data = target2,
              na.action = na.exclude) #p-value = 0.7771

#' Student's t test: raw data: normal-ish b/c robust to lots of samples, 
#' equal variance
t.test(final_ddDRD3 ~ Key, data = target2, 
       alternative = "two.sided", var.equal = TRUE,
       na.action = na.exclude)
# t = -1.5748, df = 164, p-value = 0.1172

# save p value for later
p.drd3.ttest.sq <- 0.1172

#' `ln transformed`
#' check equal variance - data is non-normal
leveneTest(final_ddDRD3 ~ Key, data = target3,
           na.action = na.exclude) #p-value = 0.6742

#' Student's t test: raw data: normal-ish b/c robust to lots of samples, 
#' equal variance
t.test(final_ddDRD3 ~ Key, data = target3, 
       alternative = "two.sided", var.equal = T,
       na.action = na.exclude)
# t = -1.7338, df = 164, p-value = 0.08482

# save p value for later
p.drd3.ttest.ln <- 0.08482


# Mixed linear regression: DRD3 ------------------------------------------
#' `sqrt`
batch.lm <- lm(final_ddDRD3 ~ Key + plate + Study_ID, 
               data = target2, na.action = na.exclude)
summary(batch.lm) #explains 5.5% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.16956 -0.04844 -0.01110  0.04500  0.24091 
# 
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.269560   0.023617  11.414   <2e-16 ***
# KeyCUD       0.014116   0.013156   1.073   0.2850    
# plate2      -0.032149   0.029473  -1.091   0.2771    
# plate3       0.079181   0.031209   2.537   0.0122 *  
# plate4       0.035510   0.030290   1.172   0.2429    
# plate5       0.017499   0.029662   0.590   0.5561    
# plate6       0.045675   0.032314   1.414   0.1596    
# plate7       0.020658   0.028890   0.715   0.4757    
# plate8       0.030149   0.029326   1.028   0.3056    
# plate9       0.001323   0.029965   0.044   0.9648    
# plate10      0.009142   0.029163   0.313   0.7544    
# plate11      0.038849   0.029757   1.306   0.1937    
# plate12     -0.037918   0.042670  -0.889   0.3756    
# Study_ID299 -0.005842   0.017350  -0.337   0.7368    
# Study_ID377  0.007855   0.016911   0.464   0.6430    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0806 on 151 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.1355,	Adjusted R-squared:  0.05536 
# F-statistic: 1.691 on 14 and 151 DF,  p-value: 0.06286

# save p value for later
p.drd3.lm.sq <- 0.2850

# diagnostics on the residuals
par(mfrow=c(2,2))
plot(batch.lm) # QQ is okay, some curving at top, a little slanting in scale-loc
# manually saved as drd3_sqrt_lm_batch.pdf

#' `ln`
batch.lm <- lm(final_ddDRD3 ~ Key + plate + Study_ID, 
               data = target3, na.action = na.exclude)
summary(batch.lm) #explains 5.8% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.87459 -0.30891 -0.01122  0.36015  1.33552 
# 
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.73058    0.15878 -17.197   <2e-16 ***
# KeyCUD       0.10432    0.08845   1.179   0.2401    
# plate2      -0.21170    0.19815  -1.068   0.2870    
# plate3       0.53185    0.20982   2.535   0.0123 *  
# plate4       0.26486    0.20364   1.301   0.1954    
# plate5       0.12809    0.19942   0.642   0.5216    
# plate6       0.26680    0.21724   1.228   0.2213    
# plate7       0.13959    0.19423   0.719   0.4734    
# plate8       0.24757    0.19716   1.256   0.2112    
# plate9       0.05799    0.20145   0.288   0.7738    
# plate10      0.08573    0.19606   0.437   0.6625    
# plate11      0.30647    0.20006   1.532   0.1276    
# plate12     -0.24494    0.28687  -0.854   0.3945    
# Study_ID299 -0.01408    0.11664  -0.121   0.9041    
# Study_ID377  0.08413    0.11369   0.740   0.4604    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5419 on 151 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.1376,	Adjusted R-squared:  0.05765 
# F-statistic: 1.721 on 14 and 151 DF,  p-value: 0.05671

# save p value for later
p.drd3.lm.ln <- 0.2401

# diagnostics on the residuals
par(mfrow=c(2,2))
plot(batch.lm) # looks good; 4 points have more leverage tho.
# manually saved as drd3_ln_lm_batch.pdf

# add residuals to new column:
target2$ddDRD3.res <- residuals(lm(final_ddDRD3 ~ 
                                     plate + Study_ID, data = target2, 
                                   na.action = na.exclude))

# add residuals to new column:
target3$ddDRD3.res <- residuals(lm(final_ddDRD3 ~ 
                                     plate + Study_ID, data = target3, 
                                   na.action = na.exclude))

# Close double plotting window
dev.off()

# Descriptive stats: DRD4 -------------------------------------------

#' `raw`
#' check equal variance - data is non-normal
leveneTest(final_ddDRD4 ~ Key, data = target,
           na.action = na.exclude) #p-value = 0.5784

#' Student's t test: raw data: normal-ish b/c robust to lots of samples, 
#' equal variance
t.test(final_ddDRD4 ~ Key, data = target, 
       alternative = "two.sided", var.equal = T,
       na.action = na.exclude)
# t = 0.50576, df = 164, p-value = 0.6137

# save p value for later
p.drd4.ttest.raw <- 0.6137

#' `sqrt transformed`
#' check equal variance: one group is normal the other isn't
leveneTest(final_ddDRD4 ~ Key, data = target2,
           na.action = na.exclude) #p-value = 0.2582

#' Student's t test: raw data: normal-ish b/c robust to lots of samples, 
#' equal variance
t.test(final_ddDRD4 ~ Key, data = target2, 
       alternative = "two.sided", var.equal = TRUE,
       na.action = na.exclude)
# t = 0.74106, df = 164, p-value = 0.4597

# save p value for later
p.drd4.ttest.sq <- 0.4597

#' `ln transformed`
#' check equal variance - data is normal
bartlett.test(final_ddDRD4 ~ Key, data = target3,
           na.action = na.exclude) #p-value = 0.09181

#' Student's t test: raw data: normal-ish b/c robust to lots of samples, 
#' equal variance
t.test(final_ddDRD4 ~ Key, data = target3, 
       alternative = "two.sided", var.equal = T,
       na.action = na.exclude)
# t = 1.0468, df = 164, p-value = 0.2967

# save p value for later
p.drd4.ttest.ln <- 0.2967


# Mixed linear regression: DRD4 ------------------------------------------
#' `sqrt`
batch.lm <- lm(final_ddDRD4 ~ Key + plate + Study_ID, 
               data = target2, na.action = na.exclude)
summary(batch.lm) #explains 2.6% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.28144 -0.09396 -0.00465  0.06165  0.38471 
# 
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.389302   0.036442  10.683   <2e-16 ***
# KeyCUD      -0.020003   0.020250  -0.988   0.3248    
# plate2       0.044614   0.045482   0.981   0.3282    
# plate3       0.029504   0.048161   0.613   0.5411    
# plate4      -0.006893   0.046025  -0.150   0.8811    
# plate5       0.004626   0.046482   0.100   0.9209    
# plate6       0.077898   0.049866   1.562   0.1203    
# plate7       0.069377   0.044584   1.556   0.1218    
# plate8       0.081890   0.045258   1.809   0.0724 .  
# plate9      -0.043345   0.046241  -0.937   0.3501    
# plate10     -0.040213   0.045004  -0.894   0.3730    
# plate11      0.017266   0.045918   0.376   0.7074    
# plate12      0.020524   0.065855   0.312   0.7557    
# Study_ID299  0.015969   0.026778   0.596   0.5518    
# Study_ID377 -0.011280   0.026127  -0.432   0.6665    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1244 on 151 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.1085,	Adjusted R-squared:  0.02582 
# F-statistic: 1.312 on 14 and 151 DF,  p-value: 0.2061

# save p value for later
p.drd4.lm.sq <- 0.3248

# diagnostics on the residuals
par(mfrow=c(2,2))
plot(batch.lm) # QQ is wonky, some curving at top. some points w/leverage
# manually saved as drd4_sqrt_lm_batch.pdf

#' `ln`
batch.lm <- lm(final_ddDRD4 ~ Key + plate + Study_ID, 
               data = target3, na.action = na.exclude)
summary(batch.lm) #explains 3.1% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.87072 -0.38045  0.06628  0.39846  1.43589 
# 
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.95952    0.18757 -10.447   <2e-16 ***
# KeyCUD      -0.13436    0.10422  -1.289    0.199    
# plate2       0.23785    0.23409   1.016    0.311    
# plate3       0.13182    0.24788   0.532    0.596    
# plate4      -0.04734    0.23689  -0.200    0.842    
# plate5       0.03612    0.23924   0.151    0.880    
# plate6       0.37740    0.25666   1.470    0.144    
# plate7       0.36502    0.22947   1.591    0.114    
# plate8       0.36881    0.23294   1.583    0.115    
# plate9      -0.26254    0.23800  -1.103    0.272    
# plate10     -0.22722    0.23164  -0.981    0.328    
# plate11      0.09056    0.23634   0.383    0.702    
# plate12      0.16320    0.33896   0.481    0.631    
# Study_ID299  0.09302    0.13783   0.675    0.501    
# Study_ID377 -0.07596    0.13448  -0.565    0.573    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6402 on 151 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.1131,	Adjusted R-squared:  0.03088 
# F-statistic: 1.375 on 14 and 151 DF,  p-value: 0.1715

# save p value for later
p.drd4.lm.ln <- 0.199

# diagnostics on the residuals
par(mfrow=c(2,2))
plot(batch.lm) # looks ok; a few points have more leverage tho.
# manually saved as drd4_ln_lm_batch.pdf

# add residuals to new column:
target2$ddDRD4.res <- residuals(lm(final_ddDRD4 ~ 
                                     plate + Study_ID, data = target2, 
                                   na.action = na.exclude))

# add residuals to new column:
target3$ddDRD4.res <- residuals(lm(final_ddDRD4 ~ 
                                     plate + Study_ID, data = target3, 
                                   na.action = na.exclude))

# Close double plotting window
dev.off()

# Descriptive stats: COMT -------------------------------------------

#' `raw`
#' check equal variance - data is non-normal
leveneTest(final_ddCOMT ~ Key, data = target,
           na.action = na.exclude) #p-value = 0.2033

#' Student's t test: raw data: normal-ish b/c robust to lots of samples, 
#' equal variance
t.test(final_ddCOMT ~ Key, data = target, 
       alternative = "two.sided", var.equal = T,
       na.action = na.exclude)
# t = 0.14378, df = 164, p-value = 0.8859

# save p value for later
p.comt.ttest.raw <- 0.8859

#' `sqrt transformed`
#' check equal variance: non normal
leveneTest(final_ddCOMT ~ Key, data = target2,
           na.action = na.exclude) #p-value = 0.08462

#' Student's t test: raw data: normal-ish b/c robust to lots of samples, 
#' equal variance
t.test(final_ddCOMT ~ Key, data = target2, 
       alternative = "two.sided", var.equal = TRUE,
       na.action = na.exclude)
# t = 0.22505, df = 164, p-value = 0.8222

# save p value for later
p.comt.ttest.sq <- 0.8222

#' `ln transformed`
#' check equal variance - data is non-normal
leveneTest(final_ddCOMT ~ Key, data = target3,
              na.action = na.exclude) #p-value = 0.0329*

#' Welch's t test: raw data: normal-ish b/c robust to lots of samples, 
#' unequal variance
t.test(final_ddCOMT ~ Key, data = target3, 
       alternative = "two.sided", var.equal = F,
       na.action = na.exclude)
# t = 0.35845, df = 156.71, p-value = 0.7205

# save p value for later
p.comt.ttest.ln <- 0.7205


# Mixed linear regression: COMT ------------------------------------------
#' `sqrt`
batch.lm <- lm(final_ddCOMT ~ Key + plate + Study_ID, 
               data = target2, na.action = na.exclude)
summary(batch.lm) #explains 21% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.14670 -0.03796 -0.00412  0.03606  0.26786 
# 
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.4129647  0.0183011  22.565  < 2e-16 ***
# KeyCUD       0.0002221  0.0101949   0.022 0.982646    
# plate2       0.0813069  0.0228385   3.560 0.000496 ***
# plate3       0.0986292  0.0241840   4.078 7.32e-05 ***
# plate4       0.0499651  0.0234720   2.129 0.034901 *  
# plate5      -0.0094464  0.0229850  -0.411 0.681669    
# plate6       0.0012706  0.0250399   0.051 0.959596    
# plate7      -0.0121548  0.0223867  -0.543 0.587967    
# plate8      -0.0227272  0.0227249  -1.000 0.318861    
# plate9       0.0072441  0.0232197   0.312 0.755485    
# plate10      0.0233529  0.0225987   1.033 0.303082    
# plate11      0.0133254  0.0230586   0.578 0.564197    
# plate12     -0.0002120  0.0330649  -0.006 0.994893    
# Study_ID299 -0.0145949  0.0134446  -1.086 0.279403    
# Study_ID377  0.0040568  0.0131043   0.310 0.757307    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06246 on 151 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.281,	Adjusted R-squared:  0.2143 
# F-statistic: 4.215 on 14 and 151 DF,  p-value: 3.136e-06

# save p value for later
p.comt.lm.sq <- 0.982646

# diagnostics on the residuals
par(mfrow=c(2,2))
plot(batch.lm) # QQ a bit wonky, some curving at top. some points w/leverage
# manually saved as comt_sqrt_lm_batch.pdf

#' `ln`
batch.lm <- lm(final_ddCOMT ~ Key + plate + Study_ID, 
               data = target3, na.action = na.exclude)
summary(batch.lm) #explains 20% of variation
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.71270 -0.16597  0.00019  0.18186  0.95597 
# 
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.779547   0.081602 -21.808  < 2e-16 ***
# KeyCUD      -0.006708   0.045458  -0.148 0.882876    
# plate2       0.356463   0.101834   3.500 0.000611 ***
# plate3       0.420934   0.107833   3.904 0.000142 ***
# plate4       0.189663   0.104658   1.812 0.071941 .  
# plate5      -0.041924   0.102487  -0.409 0.683071    
# plate6       0.001600   0.111649   0.014 0.988588    
# plate7      -0.063715   0.099819  -0.638 0.524242    
# plate8      -0.113077   0.101327  -1.116 0.266213    
# plate9       0.034800   0.103533   0.336 0.737248    
# plate10      0.110746   0.100764   1.099 0.273493    
# plate11      0.067559   0.102815   0.657 0.512124    
# plate12     -0.007842   0.147432  -0.053 0.957648    
# Study_ID299 -0.064030   0.059947  -1.068 0.287178    
# Study_ID377  0.017724   0.058430   0.303 0.762056    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2785 on 151 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.2704,	Adjusted R-squared:  0.2027 
# F-statistic: 3.997 on 14 and 151 DF,  p-value: 7.611e-06

# save p value for later
p.comt.lm.ln <- 0.882876

# diagnostics on the residuals
par(mfrow=c(2,2))
plot(batch.lm) # looks ok; a few points have more leverage tho.
# manually saved as comt_ln_lm_batch.pdf

# add residuals to new column:
target2$ddCOMT.res <- residuals(lm(final_ddCOMT ~ 
                                     plate + Study_ID, data = target2, 
                                   na.action = na.exclude))

# add residuals to new column:
target3$ddCOMT.res <- residuals(lm(final_ddCOMT ~ 
                                     plate + Study_ID, data = target3, 
                                   na.action = na.exclude))

# Close double plotting window
dev.off()


# Multiple comparisons corrections ----------------------------------------
#' `t tests on raw data`
p.ttest.raw <- c(p.drd2.ttest.raw, 
                 p.drd3.ttest.raw, 
                 p.drd4.ttest.raw, 
                 p.comt.ttest.raw)
p.ttest.raw
# 0.00124 0.15840 0.61370 0.88590

p.adjust(p.ttest.raw, method="BH")
# 0.0049600 0.3168000 0.8182667 0.8859000

p.adjust(p.ttest.raw, method="bonferroni")
# 0.00496 0.63360 1.00000 1.00000


#' `t tests on sqrt transformed data`
p.ttest.sq <- c(p.drd2.ttest.sq, 
                  p.drd3.ttest.sq, 
                  p.drd4.ttest.sq, 
                  p.comt.ttest.sq)
p.ttest.sq
# 0.001432 0.117200 0.459700 0.822200

p.adjust(p.ttest.sq, method="BH")
# 0.0068520 0.1291000 0.2866667 0.9216000

p.adjust(p.ttest.sq, method="bonferroni")
# 0.006852 0.258200 0.860000 1.000000


#' `t tests on ln transformed data`
p.ttest.ln <- c(p.drd2.ttest.ln, 
                  p.drd3.ttest.ln, 
                  p.drd4.ttest.ln, 
                  p.comt.ttest.ln)
p.ttest.ln
# 0.004987 0.084820 0.296700 0.720500

p.adjust(p.ttest.ln, method="BH")
# 0.019948 0.169640 0.395600 0.720500

p.adjust(p.ttest.ln, method="bonferroni")
# 0.019948 0.339280 1.000000 1.000000


#' `linear mixed models on sqrt transformed data`
p.lm.sq <- c(p.drd2.lm.sq, 
               p.drd3.lm.sq, 
               p.drd4.lm.sq, 
               p.comt.lm.sq)
p.lm.sq
# 0.009140 0.285000 0.324800 0.982646

p.adjust(p.lm.sq, method="BH")
# 0.0365600 0.4330667 0.4330667 0.9826460

p.adjust(p.lm.sq, method="bonferroni")
# 0.03656 1.00000 1.00000 1.00000


#' `linear mixed models on ln transformed data`
p.lm.ln <- c(p.drd2.lm.ln, 
             p.drd3.lm.ln, 
             p.drd4.lm.ln, 
             p.comt.lm.ln)
p.lm.ln
# 0.021000 0.240100 0.199000 0.882876

p.adjust(p.lm.ln, method="BH")
# 0.0840000 0.3201333 0.3201333 0.8828760

p.adjust(p.lm.ln, method="bonferroni")
# 0.0840 0.9604 0.7960 1.0000


# Plot lm residuals ------------------------------------------------------

#' `sqrt`
# set wd for image output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/f1_sf2/sqrt/")

# pivot to long format:
long <- target2 %>% dplyr::select(1,3:10,12:13,18:21) %>% 
  rename(DRD2 = ddDRD2.res, 
         DRD3 = ddDRD3.res,
         DRD4 = ddDRD4.res, 
         COMT = ddCOMT.res,
         Ethnicity = Ethnic,
         `Smoking status` = Smoker,
         `Alcohol dependence` = Dep_Alcohol,
         `Opiate dependence` = Dep_Opiates) %>% 
  pivot_longer(-c(Key, Patient_ID, Study_ID, cDNA_batch, plate,
                  Age, Gender, Ethnicity, `Smoking status`, 
                  `Alcohol dependence`, `Opiate dependence`), 
               values_to = "mRNA", 
               names_to = "target")

# make graph:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Key)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("f1_addiction_violin_qpcr_lm_sq_res.pdf")

#' make graph `Age`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Age)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 14) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=25)

ggsave("sf2_addiction_violin_qpcr_lm_sq_res_age.pdf")

#' make graph `Gender`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Gender)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf2_addiction_violin_qpcr_lm_sq_res_sex.pdf")

#' make graph `Ethnicity`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Ethnicity)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf2_addiction_violin_qpcr_lm_sq_res_eth.pdf")

#' make graph `Smoking status`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=`Smoking status`)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf2_addiction_violin_qpcr_lm_sq_res_smoke.pdf")

#' make graph `Alcohol dependence`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=`Alcohol dependence`)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf2_addiction_violin_qpcr_lm_sq_res_alcohol.pdf")

#' make graph `Opiate dependence`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=`Opiate dependence`)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf2_addiction_violin_qpcr_lm_sq_res_opiate.pdf")


#' `ln`
# set wd for image output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/f1_sf2/ln/")

# pivot to long format:
long <- target3 %>% dplyr::select(1,3:10,12:13,18:21) %>% 
  rename(DRD2 = ddDRD2.res, 
         DRD3 = ddDRD3.res,
         DRD4 = ddDRD4.res, 
         COMT = ddCOMT.res,
         Ethnicity = Ethnic,
         `Smoking status` = Smoker,
         `Alcohol dependence` = Dep_Alcohol,
         `Opiate dependence` = Dep_Opiates) %>% 
  pivot_longer(-c(Key, Patient_ID, Study_ID, cDNA_batch, plate,
                  Age, Gender, Ethnicity, `Smoking status`, 
                  `Alcohol dependence`, `Opiate dependence`), 
               values_to = "mRNA", 
               names_to = "target")

# make graph:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Key)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("f1_addiction_violin_qpcr_lm_ln_res.pdf")

#' make graph `Age`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Age)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 14) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=25)

ggsave("sf2_addiction_violin_qpcr_lm_ln_res_age.pdf")

#' make graph `Gender`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Gender)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf2_addiction_violin_qpcr_lm_ln_res_sex.pdf")

#' make graph `Ethnicity`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Ethnicity)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf2_addiction_violin_qpcr_lm_ln_res_eth.pdf")

#' make graph `Smoking status`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=`Smoking status`)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf2_addiction_violin_qpcr_lm_ln_res_smoke.pdf")

#' make graph `Alcohol dependence`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=`Alcohol dependence`)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf2_addiction_violin_qpcr_lm_ln_res_alcohol.pdf")

#' make graph `Opiate dependence`:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = mRNA)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=`Opiate dependence`)) +
  theme_pubclean() +
  ggpubr::font("xy.text", size = 12) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("") + ylab("mRNA levels (batch corrected)\n") +
  facet_wrap(. ~ target, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("sf2_addiction_violin_qpcr_lm_ln_res_opiate.pdf")

# Merge residual columns together ------------------------------------------

#' want transformed values and residuals from target2.
#' first, rename columns, then drop extra columns in prep for merge:
target2 <- target2 %>% 
  rename(ddDRD2.sq = final_ddDRD2,
         ddDRD3.sq = final_ddDRD3,
         ddDRD4.sq = final_ddDRD4,
         ddCOMT.sq = final_ddCOMT,
         ddDRD2.res.sq = ddDRD2.res,
         ddDRD3.res.sq = ddDRD3.res,
         ddDRD4.res.sq = ddDRD4.res,
         ddCOMT.res.sq = ddCOMT.res) %>% 
  dplyr::select(Patient_ID, Study_ID, cDNA_ID,
                14:21)

# merge target2 with main df
merge <- full_join(df, target2, 
                   by = c("Patient_ID", "Study_ID", "cDNA_ID"))


#' want transformed values and residuals from target3.
#' first, rename columns, then drop extra columns in prep for merge:
target3 <- target3 %>% 
  rename(ddDRD2.ln = final_ddDRD2,
         ddDRD3.ln = final_ddDRD3,
         ddDRD4.ln = final_ddDRD4,
         ddCOMT.ln = final_ddCOMT,
         ddDRD2.res.ln = ddDRD2.res,
         ddDRD3.res.ln = ddDRD3.res,
         ddDRD4.res.ln = ddDRD4.res,
         ddCOMT.res.ln = ddCOMT.res) %>% 
  dplyr::select(Patient_ID, Study_ID, cDNA_ID,
                14:21)

# merge target2 with main df
merge <- full_join(merge, target3, 
                   by = c("Patient_ID", "Study_ID", "cDNA_ID"))


# Save data ---------------------------------------------------------------
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
merge %>% write_csv(file = "5_2023_Addiction_clean.csv") 





