#' this script is used to generate Figure 2, SF3, and accompanying
#' stats
#' 
#' NB: different genes have different calibrator samples used in the 
#' ddCt calculation. These have been converted to NA data, so each 
#' gene may have varying NA sample data.

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
library(report)

# Import & tidy data 
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
df <- read_csv("5_2023_Addiction_clean.csv") #349 observations

# convert things to factors: 
df$Key <- factor(df$Key)

# convert things to characters: 
df$Patient_ID <- as.character(df$Patient_ID)
df$Study_ID <- as.character(df$Study_ID)

# Gene expression vs psychometric (sq) ------------------------------------

hyp.1 <- df %>% 
  select(Key, 3:4, BDI_total, STAI_S_total, 
         STAI_T_total, BIS11_total, SSSV_total, 64:67) %>% 
  drop_na() %>% 
  rename(BDI = BDI_total,
         `STAI-S` = STAI_S_total,
         `STAI-T` = STAI_T_total,
         BIS11 = BIS11_total,
         SSSV = SSSV_total,
         DRD2 = ddDRD2.sq, 
         DRD3 = ddDRD3.sq,
         DRD4 = ddDRD4.sq, 
         COMT = ddCOMT.sq) #157 total

# check normality assumption:
hyp.1 %>% group_by(Key) %>% normality() %>% 
  mutate(across(is.numeric, ~round(., 4))) %>% 
  regulartable()

# many fail normality; use spearman


# set working directory for graph output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/")

# create the graph & save
ggpairs(hyp.1, columns = 4:12, aes(colour=Key, alpha = 0.5), 
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

ggsave("f2_scatterplot_correlations_sq.pdf")

# create graph showing correlations after multiple comparisons
# correction & save:
grouped_ggcorrmat(hyp.1, type = "n", p.adjust.method = "fdr",
                  grouping.var = Key, 
                  colors = c("blue", "white", "red"))

ggsave("sf3_fdr_correlations_fdr_sq.pdf")


# Gene expression vs psychometric (sq res) ------------------------------------

hyp.1 <- df %>% 
  select(Key, 3:4, BDI_total, STAI_S_total, 
         STAI_T_total, BIS11_total, SSSV_total, 68:71) %>% 
  drop_na() %>% 
  rename(BDI = BDI_total,
         `STAI-S` = STAI_S_total,
         `STAI-T` = STAI_T_total,
         BIS11 = BIS11_total,
         SSSV = SSSV_total,
         DRD2 = ddDRD2.res.sq, 
         DRD3 = ddDRD3.res.sq,
         DRD4 = ddDRD4.res.sq, 
         COMT = ddCOMT.res.sq) #157 total

# check normality assumption:
hyp.1 %>% group_by(Key) %>% normality() %>% 
  mutate(across(is.numeric, ~round(., 4))) %>% 
  regulartable()

# many fail normality; use spearman


# set working directory for graph output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/")

# create the graph & save
ggpairs(hyp.1, columns = 4:12, aes(colour=Key, alpha = 0.5), 
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

ggsave("f2_scatterplot_correlations_sq_res.pdf")

# create graph showing correlations after multiple comparisons
# correction & save:
grouped_ggcorrmat(hyp.1, type = "n", p.adjust.method = "fdr",
                  grouping.var = Key, 
                  colors = c("blue", "white", "red"))

ggsave("sf3_fdr_correlations_fdr_sq_res.pdf")



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


# control samples:
cor.test(con$DRD2, con$BIS11, alternative = "two.sided",
         method = "spearman") %>% report()
# S = 113941, p-value = 0.02981*, rho = -0.2400895 (Spearman)
# technically these variables both pass normality so could use pearson
# t = -1.6516, df = 80, p-value = 0.1025, cor = -0.1815872 (Pearson)

# CUD samples:
cor.test(cud$DRD2, cud$BIS11, alternative = "two.sided",
         method = "spearman")
# S = 84868, p-value = 0.07444, rho = -0.2072227 (Spearman)
# technically these variables both pass normality so could use pearson
# t = -1.634, df = 73, p-value = 0.1066, cor = -0.1878407 (Pearson)


#' `DRD3 v BIS11`
# all samples:
cor.test(hyp.1$DRD3, hyp.1$BIS11, alternative = "two.sided",
         method = "spearman")
# S = 710962, p-value = 0.2021, rho = -0.1023425 (Spearman)

# control samples:
cor.test(con$DRD3, con$BIS11, alternative = "two.sided",
         method = "spearman")
# S = 99235, p-value = 0.4748, rho = -0.08003469 (Spearman)

# CUD samples:
cor.test(cud$DRD3, cud$BIS11, alternative = "two.sided",
         method = "spearman") %>% report()
# S = 90149, p-value = 0.01412*, rho = -0.2823405 (Spearman)


#' `DRD3 v SSSV`
# all samples:
cor.test(hyp.1$DRD3, hyp.1$SSSV, alternative = "two.sided",
         method = "spearman") %>% report()
# S = 755539, p-value = 0.03178*, rho = -0.1714575 (Spearman)

# control samples:
cor.test(con$DRD3, con$SSSV, alternative = "two.sided",
         method = "spearman")
# S = 111464, p-value = 0.05454, rho = -0.2131327  (Spearman)

# CUD samples:
cor.test(cud$DRD3, cud$SSSV, alternative = "two.sided",
         method = "spearman")
# S = 83795, p-value = 0.09896, rho = -0.1919646 (Spearman)


#' `COMT v SSSV`
# all samples:
cor.test(hyp.1$COMT, hyp.1$SSSV, alternative = "two.sided",
         method = "spearman")
# S = 586876, p-value = 0.262, rho = 0.09005283 (Spearman)

# control samples:
cor.test(con$COMT, con$SSSV, alternative = "two.sided",
         method = "spearman")
# S = 103591, p-value = 0.2538, rho = -0.127451  (Spearman)

# CUD samples:
cor.test(cud$COMT, cud$SSSV, alternative = "two.sided",
         method = "spearman") %>% report()
# S = 49905, p-value = 0.01158*, rho = 0.2901148  (Spearman)


# Gene expression vs psychometric (ln res) ------------------------------------

hyp.1 <- df %>% 
  select(Key, 3:4, BDI_total, STAI_S_total, 
         STAI_T_total, BIS11_total, SSSV_total, 76:79) %>% 
  drop_na() %>% 
  rename(BDI = BDI_total,
         `STAI-S` = STAI_S_total,
         `STAI-T` = STAI_T_total,
         BIS11 = BIS11_total,
         SSSV = SSSV_total,
         DRD2 = ddDRD2.res.ln, 
         DRD3 = ddDRD3.res.ln,
         DRD4 = ddDRD4.res.ln, 
         COMT = ddCOMT.res.ln) #157 total

# check normality assumption:
hyp.1 %>% group_by(Key) %>% normality() %>% 
  mutate(across(is.numeric, ~round(., 4))) %>% 
  regulartable()

# many fail normality; use spearman


# set working directory for graph output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/")

# create the graph & save
ggpairs(hyp.1, columns = 4:12, aes(colour=Key, alpha = 0.5), 
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

ggsave("f2_scatterplot_correlations_ln_res.pdf")

# create graph showing correlations after multiple comparisons
# correction & save:
grouped_ggcorrmat(hyp.1, type = "n", p.adjust.method = "fdr",
                  grouping.var = Key, 
                  colors = c("blue", "white", "red"))

ggsave("sf3_fdr_correlations_fdr_ln_res.pdf")



# CUD group only variables (sq) ------------------------------------------------
colnames(df)

hyp.2 <- df %>% 
  select(Key, 3:4, Stimulants_yrs, OCDUS_total, 
         64:67) %>% 
  drop_na() %>% 
  rename(`Years of use` = Stimulants_yrs,
         OCDUS = OCDUS_total,
         DRD2 = ddDRD2.sq, 
         DRD3 = ddDRD3.sq,
         DRD4 = ddDRD4.sq, 
         COMT = ddCOMT.sq) #75 total

# check normality assumption:
hyp.2 %>% group_by(Key) %>% normality() %>% 
  mutate(across(is.numeric, ~round(., 4))) %>% 
  regulartable()

# many fail normality; use Spearman

# set working directory for graph output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/")

# create the graph & save
ggpairs(hyp.2, columns = 4:9, aes(alpha = 0.5), 
        diag = list(continuous = wrap("densityDiag", alpha = 0.5)), 
        lower = list(continuous = "smooth"), 
        upper = list(continuous = wrap("cor", size = 5.5,
                                       method = "spearman"))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 12),
        axis.text.y = element_text(hjust=1, size = 12),
        strip.text = element_text(size = 18))

ggsave("sf3_cud_scatterplot_correlations_sq.pdf")


# CUD group only variables (sq res) ------------------------------------------------
colnames(df)

hyp.2 <- df %>% 
  select(Key, 3:4, Stimulants_yrs, Age_stimulants, OCDUS_total, 
         68:71) %>% 
  drop_na() %>% 
  rename(`Years using` = Stimulants_yrs,
         `Age onset` = Age_stimulants,
         OCDUS = OCDUS_total,
         DRD2 = ddDRD2.res.sq, 
         DRD3 = ddDRD3.res.sq,
         DRD4 = ddDRD4.res.sq, 
         COMT = ddCOMT.res.sq) #75 total

# check normality assumption:
hyp.2 %>% group_by(Key) %>% normality() %>% 
  mutate(across(is.numeric, ~round(., 4))) %>% 
  regulartable()

# many fail normality; spot check one
hyp.2 %>% filter(Key == "CUD") %>% 
  pull(DRD2) %>% shapiro.test()

# set working directory for graph output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/")

# create the graph & save
ggpairs(hyp.2, columns = 4:10, aes(alpha = 0.5), 
        diag = list(continuous = wrap("densityDiag", alpha = 0.5)), 
        lower = list(continuous = "smooth"), 
        upper = list(continuous = wrap("cor", size = 5.5,
                                       method = "spearman"))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 12),
        axis.text.y = element_text(hjust=1, size = 12),
        strip.text = element_text(size = 15))

ggsave("sf3_cud_scatterplot_correlations_sq_res.pdf")

# create graph showing correlations after multiple comparisons
# correction & save:
ggcorrmat(hyp.2, type = "n", p.adjust.method = "fdr", 
                  colors = c("blue", "white", "red"))
ggsave("sf3_cud_scatterplot_correlations_fdr.pdf")


# CUD group only variables (ln res)------------------------------------------------
colnames(df)

hyp.2 <- df %>% 
  select(Key, 3:4, Stimulants_yrs, OCDUS_total, 
         76:79) %>% 
  drop_na() %>% 
  rename(`Years of use` = Stimulants_yrs,
         OCDUS = OCDUS_total,
         DRD2 = ddDRD2.res.ln, 
         DRD3 = ddDRD3.res.ln,
         DRD4 = ddDRD4.res.ln, 
         COMT = ddCOMT.res.ln) #75 total

# check normality assumption:
hyp.2 %>% group_by(Key) %>% normality() %>% 
  mutate(across(is.numeric, ~round(., 4))) %>% 
  regulartable()

# many fail normality; spot check one
hyp.2 %>% filter(Key == "CUD") %>% 
  pull(DRD2) %>% shapiro.test()

# set working directory for graph output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/")

# create the graph & save
ggpairs(hyp.2, columns = 4:9, aes(alpha = 0.5), 
        diag = list(continuous = wrap("densityDiag", alpha = 0.5)), 
        lower = list(continuous = "smooth"), 
        upper = list(continuous = wrap("cor", size = 5.5,
                                       method = "spearman"))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 12),
        axis.text.y = element_text(hjust=1, size = 12),
        strip.text = element_text(size = 18))

ggsave("sf3_cud_scatterplot_correlations_ln_res.pdf")

# Individual correlations for significant hits---------------------------------

#' `COMT v Age onset`
# all samples:
cor.test(hyp.2$COMT, hyp.2$`Age onset`, alternative = "two.sided",
         method = "spearman") %>% report()
# S = 670446, p-value = 0.6231, rho = -0.0395219 (Spearman)
# technically these variables both pass normality so could use pearson
# t = -0.37714, df = 155, p-value = 0.7066, cor = -0.03027865  (Pearson)
