#' this script is used to generate SF4 and accompanying
#' stats
#' 
# Set-up ------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(dlookr)
library("PerformanceAnalytics")
library(GGally)
library(ggplot2)
library(flextable)
library(ggcorrplot)
library(ggstatsplot)
library(ggpubr)
library(report)

# Import & tidy data 
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
df <- read_csv("7_2023_Addiction_clean.csv") #349 observations

# convert things to factors: 
df$Key <- factor(df$Key)

# convert things to characters: 
df$Patient_ID <- as.character(df$Patient_ID)
df$Study_ID <- as.character(df$Study_ID)

# set directory for output graphs: 
setwd("~/~r_projects/2023_kigar_cud_paper/results/sf4/")

#### Graph pcr-based PCs vs psychometric data ####
# select relevant variables, drop rows with missing values, rename
# variables to something more readable:
hyp.1 <- df %>% 
  select(Key, 3:4, BDI_total, STAI_S_total, 
         STAI_T_total, BIS11_total, SSSV_total, 
         pcr_PC1, pcr_PC2) %>% 
  drop_na() %>% 
  rename(BDI = BDI_total,
         `STAI-S` = STAI_S_total,
         `STAI-T` = STAI_T_total,
         BIS11 = BIS11_total,
         SSSV = SSSV_total,
         `PCR PC1` = pcr_PC1, 
         `PCR PC2` = pcr_PC2) #157 total

# check normality assumption:
hyp.1 %>% group_by(Key) %>% normality() %>% 
  mutate(across(is.numeric, ~round(., 4))) %>% 
  regulartable()

# many comparsions fail normality so use Spearman

# create the graph & save
ggpairs(hyp.1, columns = 4:10, aes(colour=Key, alpha = 0.5), 
        diag = list(continuous = wrap("densityDiag", alpha = 0.5)), 
        lower = list(continuous = "smooth"), 
        upper = list(continuous = wrap("cor", size = 5.5,
                                       method = "spearman"))) + 
  theme_classic() +
  scale_color_manual(values=c('#00BFC4', '#F8766D')) +
  scale_fill_manual(values=c('#00BFC4', '#F8766D')) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 12),
        axis.text.y = element_text(hjust=1, size = 12),
        strip.text = element_text(size = 18)) 

ggsave("sf4_scatterplot_correlations_pcrPC.pdf")

# create graph showing correlations after multiple comparisons
# correction & save:
grouped_ggcorrmat(hyp.1, type = "n", p.adjust.method = "fdr",
                  grouping.var = Key, 
                  colors = c("blue", "white", "red"))

ggsave("sf4_scatterplot_correlations_pcrPC_bygroup.pdf")

# Individual correlations for significant hits---------------------------------

con <- hyp.1 %>% filter(Key == "control")
cud <- hyp.1 %>% filter(Key == "CUD")

#' `PCR PC2 v BIS11`
# CUD samples:
cor.test(cud$`PCR PC2`, cud$BIS11, alternative = "two.sided",
         method = "spearman") %>% report()
# The Spearman's rank correlation rho between cud$`PCR PC2` and cud$BIS11 is negative,
# statistically significant, and medium (rho = -0.23, S = 86666.41, p = 0.044)


#### Graph pcr vs psychometric PCs ####
# select relevant variables, drop rows with missing values, rename
# variables to something more readable:
hyp.1 <- df %>% 
  select(Key, 3:4, 68:71, 
         psych_PC1, psych_PC2) %>% 
  drop_na() %>% 
  rename(`Psych PC1` = psych_PC1,
         `Psych PC2` = psych_PC2,
         DRD2 = ddDRD2.res.sq, 
         DRD3 = ddDRD3.res.sq,
         DRD4 = ddDRD4.res.sq, 
         COMT = ddCOMT.res.sq) #157 total

# check normality assumption:
hyp.1 %>% group_by(Key) %>% normality() %>% 
  mutate(across(is.numeric, ~round(., 4))) %>% 
  regulartable()

# many comparsions fail normality so use Spearman

# create the graph & save
ggpairs(hyp.1, columns = 4:9, aes(colour=Key, alpha = 0.5), 
        diag = list(continuous = wrap("densityDiag", alpha = 0.5)), 
        lower = list(continuous = "smooth"), 
        upper = list(continuous = wrap("cor", size = 5.5,
                                       method = "spearman"))) + 
  theme_classic() +
  scale_color_manual(values=c('#00BFC4', '#F8766D')) +
  scale_fill_manual(values=c('#00BFC4', '#F8766D')) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 12),
        axis.text.y = element_text(hjust=1, size = 12),
        strip.text = element_text(size = 18)) 

ggsave("sf4_scatterplot_correlations_psychPC.pdf")

# create graph showing correlations after multiple comparisons
# correction & save:
grouped_ggcorrmat(hyp.1, type = "n", p.adjust.method = "fdr",
                  grouping.var = Key, 
                  colors = c("blue", "white", "red"))

ggsave("sf4_scatterplot_correlations_psychPC_bygroup.pdf")

# Individual correlations for significant hits---------------------------------

con <- hyp.1 %>% filter(Key == "control")
cud <- hyp.1 %>% filter(Key == "CUD")

#' `Psych PC2 v DRD3`
# overall:
cor.test(hyp.1$`Psych PC2`, hyp.1$DRD3, alternative = "two.sided",
         method = "spearman") %>% report()
# The Spearman's rank correlation rho between hyp.1$`Psych PC2` and hyp.1$DRD3 is
# negative, statistically significant, and medium (rho = -0.20, S = 7.75e+05, p = 0.011)

# CUD samples:
cor.test(cud$`Psych PC2`, cud$DRD3, alternative = "two.sided",
         method = "spearman") %>% report()
# The Spearman's rank correlation rho between cud$`Psych PC2` and cud$DRD3 is negative,
# statistically significant, and medium (rho = -0.23, S = 86306.11, p = 0.049)

#' `Psych PC2 v COMT`
# CUD samples:
cor.test(cud$`Psych PC2`, cud$COMT, alternative = "two.sided",
         method = "spearman") %>% report()
# The Spearman's rank correlation rho between cud$`Psych PC2` and cud$COMT is positive,
# statistically significant, and medium (rho = 0.27, S = 51587.60, p = 0.021)

#### Graph pcr-based PCs vs psychometric PCs ####
# select relevant variables, drop rows with missing values, rename
# variables to something more readable:
hyp.1 <- df %>% 
  select(Key, 3:4, psych_PC1, psych_PC2,
         pcr_PC1, pcr_PC2) %>% 
  drop_na() %>% 
  rename(`Psych PC1` = psych_PC1,
         `Psych PC2` = psych_PC2,
         `PCR PC1` = pcr_PC1, 
         `PCR PC2` = pcr_PC2) #157 total

# check normality assumption:
hyp.1 %>% group_by(Key) %>% normality() %>% 
  mutate(across(is.numeric, ~round(., 4))) %>% 
  regulartable()

# some comparsions fail normality so use Spearman

# create the graph & save
ggpairs(hyp.1, columns = 4:7, aes(colour=Key, alpha = 0.5), 
        diag = list(continuous = wrap("densityDiag", alpha = 0.5)), 
        lower = list(continuous = "smooth"), 
        upper = list(continuous = wrap("cor", size = 5.5,
                                       method = "spearman"))) + 
  theme_classic() +
  scale_color_manual(values=c('#00BFC4', '#F8766D')) +
  scale_fill_manual(values=c('#00BFC4', '#F8766D')) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 12),
        axis.text.y = element_text(hjust=1, size = 12),
        strip.text = element_text(size = 18)) 

ggsave("sf4_scatterplot_correlations_PCvPC.pdf")


# create graph showing correlations after multiple comparisons
# correction & save:
grouped_ggcorrmat(hyp.1, type = "n", p.adjust.method = "fdr",
                  grouping.var = Key,
                  colors = c("blue", "white", "red"))

ggsave("sf4_scatterplot_correlations_PCvPC_bygroup.pdf")


# Individual correlations for significant hits---------------------------------

con <- hyp.1 %>% filter(Key == "control")
cud <- hyp.1 %>% filter(Key == "CUD")

#' `Psych PC2 v Psych PC1`
# overall:
cor.test(cud$`Psych PC1`, cud$`Psych PC2`, alternative = "two.sided",
         method = "spearman") %>% report()

# exploring other relationships -------------------------------------------

# subset data
hyp.1 <- df %>% 
  select(Key, 3:4, psych_PC1, psych_PC2,
         pcr_PC1, pcr_PC2, CTQabuse, Stimulants_yrs,
         OCDUS_total) %>% 
  drop_na(pcr_PC1) %>% 
  rename(`Psych PC1` = psych_PC1,
         `Psych PC2` = psych_PC2,
         `PCR PC1` = pcr_PC1, 
         `PCR PC2` = pcr_PC2,
         CTQ = CTQabuse,
         `Yrs stim use` = Stimulants_yrs,
         OCDUS = OCDUS_total)

# plot PCR PC2 vs CTQ
hyp.1 %>% 
  drop_na(CTQ) %>% 
  ggscatter(x = "PCR PC2", y = "CTQ",
          color = "Key",
          add = "reg.line", conf.int = T,
          cor.coef = F, cor.method = "pearson") +
  scale_color_manual(values=c('#00BFC4', '#F8766D')) +
  scale_fill_manual(values=c('#00BFC4', '#F8766D')) +
  stat_cor(aes(color = Key), label.x = 3) +
  theme(axis.text.x = element_text(hjust=1, size = 12),
        axis.text.y = element_text(hjust=1, size = 12)) # no relationship

# plot PCR PC2 vs OCDUS
hyp.1 %>% filter(Key == "CUD") %>% 
  drop_na(OCDUS) %>% 
  ggscatter(x = "PCR PC2", y = "OCDUS",
            add = "reg.line", conf.int = T,
            cor.coef = F, cor.method = "pearson") +
  stat_cor(label.x = 3) +
  theme(axis.text.x = element_text(hjust=1, size = 12),
        axis.text.y = element_text(hjust=1, size = 12)) # no relationship

# plot PCR PC2 v years of stimulant use
hyp.1 %>% filter(Key == "CUD") %>% 
  drop_na(`Yrs stim use`) %>% 
  ggscatter(x = "PCR PC2", y = "Yrs stim use",
            add = "reg.line", conf.int = T,
            cor.coef = F, cor.method = "pearson") +
  stat_cor(label.x = 3) +
  theme(axis.text.x = element_text(hjust=1, size = 12),
        axis.text.y = element_text(hjust=1, size = 12)) # no relationship


# Individual correlations for significant hits---------------------------------

con <- hyp.1 %>% filter(Key == "control")
cud <- hyp.1 %>% filter(Key == "CUD")

#' `DRD2 v BIS11`
# all samples:
cor.test(hyp.1$DRD2, hyp.1$BIS11, alternative = "two.sided",
         method = "spearman") %>% report()
# S = 670446, p-value = 0.6231, rho = -0.0395219 (Spearman)
# technically these variables both pass normality so could use pearson
# t = -0.37714, df = 155, p-value = 0.7066, cor = -0.03027865  (Pearson)


# CUD-specific relationships ----------------------------------------------

# subset
hyp.2 <- df %>% filter(Key == "CUD") %>% 
  select(Key, 3:4, psych_PC1, psych_PC2,
         pcr_PC1, pcr_PC2, Stimulants_yrs, Age_stimulants,
         OCDUS_total) %>% 
  drop_na(pcr_PC1) %>% 
  rename(`Psych PC1` = psych_PC1,
         `Psych PC2` = psych_PC2,
         `PCR PC1` = pcr_PC1, 
         `PCR PC2` = pcr_PC2,
         `Years using` = Stimulants_yrs,
         `Age onset` = Age_stimulants,
         OCDUS = OCDUS_total)


# plot
ggpairs(hyp.2, columns = 6:10, aes(alpha = 0.5), 
        diag = list(continuous = wrap("densityDiag", alpha = 0.5)), 
        lower = list(continuous = "smooth"), 
        upper = list(continuous = wrap("cor", size = 5.5,
                                       method = "spearman"))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 12),
        axis.text.y = element_text(hjust=1, size = 12),
        strip.text = element_text(size = 15)) 
