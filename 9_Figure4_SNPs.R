#' this script is used to generate F5, SF5, SF6, and accompanying stats
#' 
#' info about SNPs coming from DOI: 10.1016/j.humimm.2015.09.032
#' Lori Turner and Alex Hatton generated genotyping data. There are
#' more genotyped samples than there were RNA samples for qPCR.

# Set up ------------------------------------------------------------------

#' load libraries
library(tidyverse)
library(magrittr)
library(ggpubr)
library(rstatix)
library(broom)
library(chisq.posthoc.test)
library(report)


# Import & tidy data 
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
df <- read_csv("7_2023_Addiction_clean.csv") #349 observations

# make Key a factor:
df$Key <- factor(df$Key)

# select variables of interest (molecular and PC data + BIS11, SSSV):
target <- df %>% select(Key, Patient_ID, Study_ID, 
                        BIS11_total, SSSV_total,
                        42:47,68:71,80:83)


# pivot long:
long <- target %>% rename(BIS11 = BIS11_total,
                          SSSV = SSSV_total,
                          `Psych PC1` = psych_PC1,
                          `Psych PC2` = psych_PC2,
                          `PCR PC1` = pcr_PC1, 
                          `PCR PC2` = pcr_PC2,
                          DRD2 = ddDRD2.res.sq, 
                          DRD3 = ddDRD3.res.sq,
                          DRD4 = ddDRD4.res.sq, 
                          COMT = ddCOMT.res.sq) %>% 
  pivot_longer(-c(Key, Patient_ID, Study_ID, 
                  BIS11, SSSV,
                  `Psych PC1`, `Psych PC2`, 
                  `PCR PC1`, `PCR PC2`, 
                  DRD2, DRD3, DRD4, COMT),
               values_to = "genotype", 
               names_to = "SNP") %>% 
  tidyr::drop_na(genotype)


# create empty vectors for planned comparisons:
p.chi <- c(NA, NA, NA) # for chi-square results
p.pc1 <- c(NA, NA, NA) # for psych pc1 ~ genotype x group
p.pc2 <- c(NA, NA, NA) # for psych pc2 ~ genotype x group
p.bis <- c(NA, NA, NA) # for bis11 ~ genotype x group
p.rna <- c(NA, NA, NA) # for RNA ~ genotype x group
p.snp <- c(NA, NA, NA) # for RNA ~ genotype

# *rs4680 (COMT-related)--------------------------------------------------------
#' this is a COMT-related SNP. Not sure about associations with disease

# set output directory:
setwd("~/~r_projects/2023_kigar_cud_paper/results/f5/")

# filter to SNP of interest: 
rs4680 <- long %>% dplyr::filter(SNP == "rs4680") #278

rs4680$genotype <- factor(rs4680$genotype)

#### Chi square for frequency ####
table(rs4680$Key, rs4680$genotype)
#         A/A A/G G/G
# control  37  73  41
# CUD      33  58  36

rs4680 %>%  
  ggplot(aes(x = Key, fill = genotype)) + 
  geom_bar(position = "fill") +
  labs(y = "Proportion") + theme_pubclean() +
  ggpubr::font("xy.text", size = 18) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) + 
  scale_fill_manual(values = c("#cccccc", "#666666", "#000000")) +
  xlab("") + ylab("Frequency\n") 
ggsave("chi_square_genotype/rs4680_frequency_stackedbar.pdf")


# create count summary of group x genotype, stored in temp object M:
M <- rs4680 %>% select(Key, genotype) %>% 
  group_by(Key, genotype) %>% 
  summarise(n = n())

# extract counts:
col <- M$n

# convert counts into a 2 row x 2 column matrix:
snp <- matrix(col, nrow = 3, ncol = 2)

# name rows and columns:
dimnames(snp) = list(genotype = M$genotype[1:3], 
                     group = c("control", "CUD"))

# perform chi squared test:
chisq.test(snp) 
# X-squared = 0.20035, df = 2, p-value = 0.9047
p.chi[1] <- 0.9047

# get counts for gene expression x snp data
temp <- rs4680 %>% tidyr::drop_na(COMT)
table(temp$Key, temp$genotype)
#         A/A A/G G/G
# control  21  39  23
# CUD      20  33  21

# get counts for psychometric data x snp d
temp <- rs4680 %>% tidyr::drop_na(BIS11)
table(temp$Key, temp$genotype)
#         A/A A/G G/G
# control  36  71  40
# CUD      32  57  34

#### Visualize relationships ####

# graph genotype x Key --> Psych PC1
rs4680 %>% ggplot(aes(x = genotype, y = `Psych PC1`, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs4680") + ylab("Psych PC1\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("snp_v_psychPC1/f4_rs4680_psychPC1.pdf")

# graph genotype x Key --> Psych PC2
rs4680 %>% ggplot(aes(x = genotype, y = `Psych PC2`, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs4680") + ylab("Psych PC2\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("snp_v_psychPC2/f4_rs4680_psychPC2.pdf")

# graph genotype x Key --> COMT
rs4680 %>% ggplot(aes(x = genotype, y = COMT, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs4680") + ylab("COMT mRNA (batch corrected)\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("bespoke_snp_gene/f4_rs4680_COMT.pdf")

# graph genotype --> COMT
rs4680 %>% ggplot(aes(x = genotype, y = COMT, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs4680") + ylab("COMT mRNA (batch corrected)\n") +
  theme(legend.position="top") 
ggsave("bespoke_snp_gene/f5_rs4680_COMT.pdf")

# graph genotype x Key --> BIS11
rs4680 %>% ggplot(aes(x = genotype, y = BIS11, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs4680") + ylab("BIS11 total\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("snp_v_bis11/f4_rs4680_BIS11.pdf")

# graph PCR PC2 x genotype --> Psych PC1
rs4680 %>% 
  ggplot(aes (x = `PCR PC2`, y = `Psych PC1`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nPCR PC2") + ylab("Psych PC1\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("pcrPC2_v_psychPC1/f4_rs4680_pcrPC2_v_psychPC1.pdf")

# graph PCR PC2 x genotype --> Psych PC2
rs4680 %>% 
  ggplot(aes (x = `PCR PC2`, y = `Psych PC2`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nPCR PC2") + ylab("Psych PC2\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("pcrPC2_v_psychPC2/f4_rs4680_pcrPC2_v_psychPC2.pdf")

# graph COMT x genotype --> Psych PC1
rs4680 %>% 
  ggplot(aes (x = COMT, y = `Psych PC1`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nCOMT mRNA (batch corrected)") + ylab("Psych PC1\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("bespoke_v_psychPC1/f4_rs4680_COMT_v_psychPC1.pdf")


# graph COMT x genotype --> Psych PC2
rs4680 %>% 
  ggplot(aes (x = COMT, y = `Psych PC2`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nCOMT mRNA (batch corrected)") + ylab("Psych PC2\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("bespoke_v_psychPC2/f4_rs4680_COMT_v_psychPC2.pdf")

# graph COMT x genotype --> BIS11
rs4680 %>% 
  ggplot(aes (x = COMT, y = BIS11, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nCOMT mRNA (batch corrected)") + ylab("BIS11 total\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("bespoke_v_BIS11/f4_rs4680_COMT_v_BIS11.pdf")

#### Run Anova ####

aov(COMT ~ genotype, data = rs4680) %>% report()
p.snp[1] <- 0.127
  
aov(COMT ~ Key * genotype,
    data = rs4680) %>% report()
# The ANOVA (formula: COMT ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically not significant and very small (F(1, 151) =
#                                                                               1.24e-03, p = 0.972; Eta2 (partial) = 8.20e-06, 95% CI [0.00, 1.00])
# - The main effect of genotype is statistically not significant and small (F(2, 151) =
#                                                                             2.05, p = 0.132; Eta2 (partial) = 0.03, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and very
# small (F(2, 151) = 0.14, p = 0.870; Eta2 (partial) = 1.85e-03, 95% CI [0.00, 1.00])
# 
# Effect sizes were labelled following Field's (2013) recommendations.
p.rna[1] <- 0.870 # interaction


aov(`Psych PC1` ~ genotype, data = rs4680) %>% report()
# The ANOVA (formula: `Psych PC1` ~ genotype) suggests that:
#   
#   - The main effect of genotype is statistically not significant and very small
# (F(2, 266) = 0.64, p = 0.527; Eta2 = 4.80e-03, 95% CI [0.00, 1.00])

aov(`Psych PC1` ~ Key * genotype,
    data = rs4680) %>% report()
# The ANOVA (formula: `Psych PC1` ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and large (F(1, 263) = 187.32, p <
#                                                                      .001; Eta2 (partial) = 0.42, 95% CI [0.34, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 263)
#                                                                                = 1.20, p = 0.302; Eta2 (partial) = 9.06e-03, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and small
# (F(2, 263) = 1.72, p = 0.182; Eta2 (partial) = 0.01, 95% CI [0.00, 1.00])
# 
# Effect sizes were labelled following Field's (2013) recommendations.
p.pc1[1] <- 0.182 # interaction

aov(`Psych PC2` ~ genotype, data = rs4680) %>% report()
# The ANOVA (formula: `Psych PC2` ~ genotype) suggests that:
#   
#   - The main effect of genotype is statistically not significant and very small
# (F(2, 266) = 0.33, p = 0.717; Eta2 = 2.50e-03, 95% CI [0.00, 1.00])

aov(`Psych PC2` ~ Key * genotype,
            data = rs4680) %>% report()
# The ANOVA (formula: `Psych PC2` ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and small (F(1, 263) = 14.56, p <
#                                                                      .001; Eta2 (partial) = 0.05, 95% CI [0.02, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 263)
#                                                                                = 0.32, p = 0.730; Eta2 (partial) = 2.39e-03, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and very
# small (F(2, 263) = 0.51, p = 0.604; Eta2 (partial) = 3.83e-03, 95% CI [0.00, 1.00])
# 
# Effect sizes were labelled following Field's (2013) recommendations.
p.pc2[1] <- 0.604 # interaction

aov(BIS11 ~ genotype, data = rs4680) %>% report()
# The ANOVA (formula: BIS11 ~ genotype) suggests that:
#   
#   - The main effect of genotype is statistically not significant and very small
# (F(2, 267) = 0.47, p = 0.624; Eta2 = 3.52e-03, 95% CI [0.00, 1.00])

aov(BIS11 ~ Key * genotype,
    data = rs4680) %>% report()
# The ANOVA (formula: BIS11 ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and large (F(1, 264) = 145.19, p <
#                                                                      .001; Eta2 (partial) = 0.35, 95% CI [0.28, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 264)
#                                                                                = 0.90, p = 0.407; Eta2 (partial) = 6.79e-03, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and small
# (F(2, 264) = 1.67, p = 0.190; Eta2 (partial) = 0.01, 95% CI [0.00, 1.00])
# 
# Effect sizes were labelled following Field's (2013) recommendations.
p.bis[1] <- 0.190 # interaction

aov(SSSV ~ Key * genotype,
    data = rs4680) %>% report()
# The ANOVA (formula: SSSV ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and medium (F(1, 270) = 33.48, p <
#                                                                       .001; Eta2 (partial) = 0.11, 95% CI [0.06, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 270)
#                                                                                = 0.49, p = 0.614; Eta2 (partial) = 3.61e-03, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and very
# small (F(2, 270) = 0.77, p = 0.462; Eta2 (partial) = 5.71e-03, 95% CI [0.00, 1.00])


#### Test assumptions ####
# set output directory:
setwd("~/~r_projects/2023_kigar_cud_paper/results/f4/QC/")

#' Linearity & homogeneity check
#' looking for relatively straight lines and insignificant 
#' interaction terms!
rs4680 %>% tidyr::drop_na(`PCR PC2`) %>% 
  tidyr::drop_na(`Psych PC1`) %>% 
  ggscatter(x = "PCR PC2", y = "Psych PC1",
            title = "rs4680",
            facet.by  = c("genotype", "Key"), 
            short.panel.labs = FALSE) +
  stat_smooth(method = "loess", span = 0.9) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("linearitycheck_rs4680_psychpc1.pdf")

long %>% dplyr::filter(SNP == "rs4680") %>% 
  tidyr::drop_na(DRD2) %>% tidyr::drop_na(`Psych PC1`) %>% 
  unite(col = "group", genotype, Key) %>% 
  anova_test(`Psych PC1` ~ group*DRD2)

#       Effect DFn DFd      F        p p<.05   ges
# 1      group   5 138 20.349 3.54e-15     * 0.424
# 2       DRD2   1 138  3.127 7.90e-02       0.022
# 3 group:DRD2   5 138  0.617 6.87e-01       0.022

#' Normality & homogeneity check
rs4680 <- long %>% dplyr::filter(SNP == "rs4680") %>% 
  tidyr::drop_na(SNP) %>% tidyr::drop_na(DRD2) %>% tidyr::drop_na(`Psych PC1`)

# Fit the model (covariate goes first)
model <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs4680,
            na.action = na.exclude)

# Create model diagnostic metrics
model.metrics <- augment(model) %>% select(-.hat, -.sigma, -.fitted) 

# Assess normality of residuals using shapiro wilk test
shapiro_test(model.metrics$.resid) #non-normal. p = 0.00161

# Assess homogeneity of variance with Levene's test
levene_test(.resid ~ genotype*Key, data = model.metrics)
# non-homogeneous. p = 0.0385

# Assess outliers
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

# one possible outlier:
# Psych PC1       DRD2 genotype Key    .resid   .cooksd .std.resid
# 1 -3.143709 0.03701038      A/G CUD -2.581016 0.0578607  -3.483283


#### Run ANCOVA ####
ancova.4680 <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs4680,
                 na.action = na.exclude)
anova(ancova.4680)

# Response: Psych PC1
#               Df Sum Sq Mean Sq F value  Pr(>F)
# DRD2           1  0.184   0.184  0.3311 0.56590
# genotype       2  1.689   0.845  1.5174 0.22280
# Key            1 54.686  54.686 98.2499 < 2e-16 ***
# genotype:Key   2  3.676   1.838  3.3019 0.03964 *
# Residuals    143 79.593   0.557

par(mfrow=c(1,3))
plot(ancova.4680, add.smooth = FALSE, which = 1)
plot(ancova.4680, which = 2)
plot(ancova.4680, add.smooth = FALSE, which = 3)
#' QQ snakes at the ends. a little bending for scale-location but looks okay.
#'  manually saved as rs4680_ancova_psychPC1_drd2.pdf
dev.off()



# *rs6280 (DRD3-related) -----------------------------------------------------
#' this is a DRD3-related SNP. T is associated with higher dopamine
#' binding in vitro + alcohol and heroin dependence

# set output directory:
setwd("~/~r_projects/2023_kigar_cud_paper/results/f5/")

# filter to SNP of interest:
rs6280 <- long %>% dplyr::filter(SNP == "rs6280") %>% 
  dplyr::filter(genotype != "G/G") 
# one control G/G and no others. assume is error. n = 276

rs6280$genotype <- factor(rs6280$genotype)

#### *Chi square for frequency ####
table(rs6280$Key, rs6280$genotype)
#          C/C C/T G/G T/T
# control  13  70   1  66
# CUD      25  55   0  47

rs6280 %>%  
  ggplot(aes(x = Key, fill = genotype)) + 
  geom_bar(position = "fill") +
  labs(y = "Proportion") + theme_pubclean() +
  ggpubr::font("xy.text", size = 18) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) + 
  scale_fill_manual(values = c("#cccccc", "#666666", "#000000")) +
  xlab("") + ylab("Frequency\n") 
ggsave("chi_square_genotype/rs6280_frequency_stackedbar.pdf")


# create count summary of group x genotype, stored in temp object M:
M <- rs6280 %>% select(Key, genotype) %>% 
  group_by(Key, genotype) %>% 
  summarise(n = n())

# extract counts:
col <- M$n

# convert counts into a 2 row x 2 column matrix:
snp <- matrix(col, nrow = 3, ncol = 2)

# name rows and columns:
dimnames(snp) = list(genotype = M$genotype[1:3], 
                     group = c("control", "CUD"))

# perform chi squared test:
chisq.test(snp)
# X-squared = 7.0755, df = 2, p-value = 0.02908
p.chi[2] <- 0.02908

# run post-hoc pairwise comparisons using fdr correction method:
chisq.posthoc.test(snp, method = "fdr")
# Dimension     Value    control        CUD
# 1       C/C Residuals -2.6338317  2.6338317
# 2       C/C  p values  0.0253280  0.0253280 *
# 3       C/T Residuals  0.6109432 -0.6109432
# 4       C/T  p values  1.0000000  1.0000000
# 5       T/T Residuals  1.2271306 -1.2271306
# 6       T/T  p values  0.6593210  0.6593210

# get counts for gene x SNP data
temp <- rs6280 %>% tidyr::drop_na(DRD3)
table(temp$Key, temp$genotype)
#           C/C C/T T/T
# control   9  36  37
# CUD      16  35  23

# get counts for psychometric x SNP data
temp <- rs6280 %>% tidyr::drop_na(BIS11)
table(temp$Key, temp$genotype)
#           C/C C/T T/T
# control  12  68  65
# CUD      24  53  46

#### Visualize relationships ####

# graph genotype x Key --> Psych PC1
rs6280 %>% ggplot(aes(x = genotype, y = `Psych PC1`, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs6280") + ylab("Psych PC1\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("snp_v_psychPC1/f4_rs6280_psychPC1.pdf")

# graph genotype x Key --> Psych PC2
rs6280 %>% ggplot(aes(x = genotype, y = `Psych PC2`, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs6280") + ylab("Psych PC2\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("snp_v_psychPC2/f4_rs6280_psychPC2.pdf")

# graph genotype x Key --> DRD3
rs6280 %>% ggplot(aes(x = genotype, y = DRD3, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs6280") + ylab("DRD3 mRNA (batch corrected)\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("bespoke_snp_gene/f4_rs6280_DRD3.pdf")

# graph genotype --> DRD3
rs6280 %>% ggplot(aes(x = genotype, y = DRD3, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs6280") + ylab("DRD3 mRNA (batch corrected)\n") +
  theme(legend.position="top") 
ggsave("bespoke_snp_gene/f5_rs6280_DRD3.pdf")

# graph genotype x Key --> BIS11
rs6280 %>% ggplot(aes(x = genotype, y = BIS11, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs6280") + ylab("BIS11 total\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("snp_v_bis11/f4_rs6280_BIS11.pdf")

# graph PCR PC2 x genotype --> Psych PC1
rs6280 %>% 
  ggplot(aes (x = `PCR PC2`, y = `Psych PC1`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nPCR PC2") + ylab("Psych PC1\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("pcrPC2_v_psychPC1/f4_rs6280_pcrPC2_v_psychPC1.pdf")

# graph PCR PC2 x genotype --> Psych PC2
rs6280 %>% 
  ggplot(aes (x = `PCR PC2`, y = `Psych PC2`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nPCR PC2") + ylab("Psych PC2\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("pcrPC2_v_psychPC2/f4_rs6280_pcrPC2_v_psychPC2.pdf")

# graph DRD3 x genotype --> Psych PC1
rs6280 %>% 
  ggplot(aes (x = DRD3, y = `Psych PC1`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nDRD3 mRNA (batch corrected)") + ylab("Psych PC1\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("bespoke_v_psychPC1/f4_rs6280_DRD3_v_psychPC1.pdf")


# graph DRD3 x genotype --> Psych PC2
rs6280 %>% 
  ggplot(aes (x = DRD3, y = `Psych PC2`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nDRD3 mRNA (batch corrected)") + ylab("Psych PC2\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("bespoke_v_psychPC2/f4_rs6280_DRD3_v_psychPC2.pdf")

# graph DRD3 x genotype --> BIS11
rs6280 %>% 
  ggplot(aes (x = DRD3, y = BIS11, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nDRD3 mRNA (batch corrected)") + ylab("BIS11 total\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("bespoke_v_BIS11/f4_rs6280_DRD3_v_BIS11.pdf")

#### Run Anova ####
aov(DRD3 ~ genotype, data = rs6280) %>% report()
p.snp[2] <- 0.040
aov(DRD3 ~ genotype, data = rs6280) %>% tukey_hsd()
# term     group1 group2 null.value estimate conf.low conf.high  p.adj p.adj.signif
# * <chr>    <chr>  <chr>       <dbl>    <dbl>    <dbl>     <dbl>  <dbl> <chr>       
# 1 genotype C/C    C/T             0  -0.0120  -0.0530   0.0290  0.767  ns          
# 2 genotype C/C    T/T             0  -0.0392  -0.0811   0.00280 0.0729 ns          
# 3 genotype C/T    T/T             0  -0.0271  -0.0580   0.00378 0.098  ns 

aov(DRD3 ~ Key * genotype,
    data = rs6280) %>% report()
# The ANOVA (formula: DRD3 ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically not significant and very small (F(1, 150) =
#                                                                               0.51, p = 0.475; Eta2 (partial) = 3.41e-03, 95% CI [0.00, 1.00])
# - The main effect of genotype is statistically significant and small (F(2, 150) = 3.11,
#                                                                       p = 0.048; Eta2 (partial) = 0.04, 95% CI [2.27e-04, 1.00])
# - The interaction between Key and genotype is statistically not significant and small
# (F(2, 150) = 2.23, p = 0.111; Eta2 (partial) = 0.03, 95% CI [0.00, 1.00])
# 
# Effect sizes were labelled following Field's (2013) recommendations.
p.rna[2] <- 0.111 # interaction


aov(`Psych PC1` ~ Key * genotype,
    data = rs6280) %>% report()
# The ANOVA (formula: `Psych PC1` ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and large (F(1, 261) = 184.68, p <
#                                                                      .001; Eta2 (partial) = 0.41, 95% CI [0.34, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 261)
#                                                                                = 0.82, p = 0.441; Eta2 (partial) = 6.26e-03, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and small
# (F(2, 261) = 2.62, p = 0.075; Eta2 (partial) = 0.02, 95% CI [0.00, 1.00])
# 
# Effect sizes were labelled following Field's (2013) recommendations.
p.pc1[2] <- 0.075 # interaction


aov(`Psych PC2` ~ Key * genotype,
    data = rs6280) %>% report()
# The ANOVA (formula: `Psych PC2` ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and small (F(1, 261) = 14.86, p <
#                                                                      .001; Eta2 (partial) = 0.05, 95% CI [0.02, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 261)
#                                                                                = 0.13, p = 0.880; Eta2 (partial) = 9.79e-04, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and small
# (F(2, 261) = 2.20, p = 0.113; Eta2 (partial) = 0.02, 95% CI [0.00, 1.00])
# 
# Effect sizes were labelled following Field's (2013) recommendations.
p.pc2[2] <- 0.113 # interaction


aov(BIS11 ~ Key * genotype,
    data = rs6280) %>% report()
# The ANOVA (formula: BIS11 ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and large (F(1, 262) = 143.79, p <
#                                                                      .001; Eta2 (partial) = 0.35, 95% CI [0.28, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 262)
#                                                                                = 0.47, p = 0.625; Eta2 (partial) = 3.59e-03, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and small
# (F(2, 262) = 2.01, p = 0.136; Eta2 (partial) = 0.02, 95% CI [0.00, 1.00])
# 
# Effect sizes were labelled following Field's (2013) recommendations.
p.bis[2] <- 0.136 # interaction

aov(SSSV ~ genotype, data = rs6280) %>% report()
# The ANOVA (formula: SSSV ~ genotype) suggests that:
#   
#   - The main effect of genotype is statistically not significant and very small (F(2, 271)
#                                                                                  = 0.01, p = 0.988; Eta2 = 8.77e-05, 95% CI [0.00, 1.00])

aov(SSSV ~ Key * genotype,
    data = rs6280) %>% report()
# The ANOVA (formula: SSSV ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and medium (F(1, 268) = 33.64, p <
#                                                                       .001; Eta2 (partial) = 0.11, 95% CI [0.06, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 268)
#                                                                                = 0.38, p = 0.685; Eta2 (partial) = 2.82e-03, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and small
# (F(2, 268) = 2.16, p = 0.118; Eta2 (partial) = 0.02, 95% CI [0.00, 1.00])

#### Test assumptions ####
# set output directory:
setwd("~/~r_projects/2023_kigar_cud_paper/results/f4/QC/")

#' Linearity & homogeneity check
#' looking for relatively straight lines and insignificant 
#' interaction terms!
long %>% dplyr::filter(SNP == "rs6280") %>% tidyr::drop_na(DRD2) %>% 
  tidyr::drop_na(`Psych PC1`) %>%
  ggscatter(x = "DRD2", y = "Psych PC1",
            title = "rs6280",
            facet.by  = c("genotype", "Key"), 
            short.panel.labs = FALSE) +
  stat_smooth(method = "loess", span = 0.9) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("linearitycheck_rs6280_psychpc1.pdf") # very few (~8) C/C controls

long %>% dplyr::filter(SNP == "rs6280") %>%
  tidyr::drop_na(DRD2) %>% 
  tidyr::drop_na(`Psych PC1`) %>%
  unite(col = "group", genotype, Key) %>% 
  anova_test(`Psych PC1` ~ group*DRD2)

#       Effect DFn DFd      F        p p<.05   ges
# 1      group   5 137 22.301 2.60e-16     * 0.449
# 2       DRD2   1 137  3.418 6.70e-02       0.024
# 3 group:DRD2   5 137  1.217 3.04e-01       0.043

#' Normality & homogeneity check
rs6280 <- long %>% dplyr::filter(SNP == "rs6280") %>% 
  tidyr::drop_na(SNP) %>% tidyr::drop_na(DRD2) %>% tidyr::drop_na(`Psych PC1`)

# Fit the model (covariate goes first)
model <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs6280,
            na.action = na.exclude)

# Create model diagnostic metrics
model.metrics <- augment(model) %>% select(-.hat, -.sigma, -.fitted) 

# Assess normality of residuals using shapiro wilk test
shapiro_test(model.metrics$.resid) #non-normal. p = 0.00125

# Assess homogeneity of variance with Levene's test
levene_test(.resid ~ genotype*Key, data = model.metrics)
# homogeneous. p = 0.0826

# Assess outliers
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

# one possible outlier:
# Psych PC1       DRD2 genotype Key   .resid    .cooksd .std.resid
# 1 -3.143709 0.03701038      T/T CUD -2.26301 0.06465553  -3.120389


#### Run ANCOVA ####
ancova.6280 <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs6280,
                 na.action = na.exclude)
anova(ancova.6280)

# Response: Psych PC1
#               Df Sum Sq Mean Sq F value  Pr(>F)
# DRD2           1  0.184   0.184  0.3311 0.56590
# genotype       2  1.689   0.845  1.5174 0.22280
# Key            1 54.686  54.686 98.2499 < 2e-16 ***
# genotype:Key   2  3.676   1.838  3.3019 0.03964 *
# Residuals    143 79.593   0.557

par(mfrow=c(1,3))
plot(ancova.6280, add.smooth = FALSE, which = 1)
plot(ancova.6280, which = 2)
plot(ancova.6280, add.smooth = FALSE, which = 3)
#' QQ snakes at the ends. a little bending for scale-location but looks okay.
#'  manually saved as rs6280_ancova_psychPC1_drd2.pdf
dev.off()



# *rs6277 (DRD2-related) ------------------------------------------------------
#' this is a DRD2-related SNP. A is associated with decreased
#' DRD2 mRNA stability and translation

# set output directory:
setwd("~/~r_projects/2023_kigar_cud_paper/results/f5/")

# filter to SNP of interest: 
rs6277 <- long %>% dplyr::filter(SNP == "rs6277")  

rs6277$genotype <- factor(rs6277$genotype)

#### **Chi square for frequency ####
table(rs6277$Key, rs6277$genotype)
#          A/A A/G G/G
# control  47  69  34
# CUD      20  61  46

rs6277 %>%  
  ggplot(aes(x = Key, fill = genotype)) + 
  geom_bar(position = "fill") +
  labs(y = "Proportion") + theme_pubclean() +
  ggpubr::font("xy.text", size = 18) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) + 
  scale_fill_manual(values = c("#cccccc", "#666666", "#000000")) +
  xlab("") + ylab("Frequency\n") 
ggsave("chi_square_genotype/rs6277_frequency_stackedbar.pdf")


# create count summary of group x genotype, stored in temp object M:
M <- rs6277 %>% select(Key, genotype) %>% 
  group_by(Key, genotype) %>% 
  summarise(n = n())

# extract counts:
col <- M$n

# convert counts into a 2 row x 2 column matrix:
snp <- matrix(col, nrow = 3, ncol = 2)

# name rows and columns:
dimnames(snp) = list(genotype = M$genotype[1:3], 
                     group = c("control", "CUD"))

# perform chi squared test:
chisq.test(snp)
# X-squared = 11.341, df = 2, p-value = 0.003446
p.chi[3] <- 0.003446

# run post-hoc pairwise comparisons using fdr correction method:
chisq.posthoc.test(snp, method = "fdr")
# Dimension     Value    control        CUD
# 1       A/A Residuals  3.0182523 -3.0182523
# 2       A/A  p values  0.0076270  0.0076270 **
# 3       A/G Residuals -0.3375771  0.3375771
# 4       A/G  p values  1.0000000  1.0000000
# 5       G/G Residuals -2.4801076  2.4801076
# 6       G/G  p values  0.0394030  0.0394030 *

# get counts for gene x SNP data
temp <- rs6277 %>% tidyr::drop_na(DRD2)
table(temp$Key, temp$genotype)
#         A/A A/G G/G
# control  26  40  17
# CUD      15  34  25

# get counts for psychometric x SNP data
temp <- rs6277 %>% tidyr::drop_na(BIS11)
table(temp$Key, temp$genotype)
#         A/A A/G G/G
# control  45  67  34
# CUD      20  60  43

#### Visualize relationships ####

# graph genotype x Key --> Psych PC1
rs6277 %>% ggplot(aes(x = genotype, y = `Psych PC1`, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs6277") + ylab("Psych PC1\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("snp_v_psychPC1/f4_rs6277_psychPC1.pdf")

# graph genotype x Key --> Psych PC2
rs6277 %>% ggplot(aes(x = genotype, y = `Psych PC2`, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs6277") + ylab("Psych PC2\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("snp_v_psychPC2/f4_rs6277_psychPC2.pdf")

# graph genotype x Key --> DRD2
rs6277 %>% ggplot(aes(x = genotype, y = DRD2, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs6277") + ylab("DRD2 mRNA (batch corrected)\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("bespoke_snp_gene/f4_rs6277_DRD2.pdf")

# graph genotype  --> DRD2
rs6277 %>% ggplot(aes(x = genotype, y = DRD2, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs6277") + ylab("DRD2 mRNA (batch corrected)\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) 
ggsave("bespoke_snp_gene/f5_rs6277_DRD2.pdf")

# graph genotype x Key --> BIS11
rs6277 %>% ggplot(aes(x = genotype, y = BIS11, col = genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
  theme_pubclean() +
  ggpubr::font("xy.text", size = 14) +
  ggpubr::font("ylab", size = 18) +
  ggpubr::font("xlab", size = 18) +
  ggpubr::font("legend.title", size = 18) +
  ggpubr::font("legend.text", size = 18) +
  xlab("\nrs6277") + ylab("BIS11 total\n") +
  theme(legend.position="top",
        strip.text = element_text(size = rel(1.5))) +
  facet_wrap(. ~ Key, scales = "fixed") 
ggsave("snp_v_bis11/f4_rs6277_BIS11.pdf")

# graph PCR PC2 x genotype --> Psych PC1
rs6277 %>% 
  ggplot(aes (x = `PCR PC2`, y = `Psych PC1`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nPCR PC2") + ylab("Psych PC1\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("pcrPC2_v_psychPC1/f4_rs6277_pcrPC2_v_psychPC1.pdf")

# graph PCR PC2 x genotype --> Psych PC2
rs6277 %>% 
  ggplot(aes (x = `PCR PC2`, y = `Psych PC2`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nPCR PC2") + ylab("Psych PC2\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("pcrPC2_v_psychPC2/f4_rs6277_pcrPC2_v_psychPC2.pdf")

# graph DRD2 x genotype --> Psych PC1
rs6277 %>% 
  ggplot(aes (x = DRD2, y = `Psych PC1`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nDRD2 mRNA (batch corrected)") + ylab("Psych PC1\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("bespoke_v_psychPC1/f4_rs6277_DRD2_v_psychPC1.pdf")


# graph DRD2 x genotype --> Psych PC2
rs6277 %>% 
  ggplot(aes (x = DRD2, y = `Psych PC2`, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nDRD2 mRNA (batch corrected)") + ylab("Psych PC2\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("bespoke_v_psychPC2/f4_rs6277_DRD2_v_psychPC2.pdf")

# graph DRD2 x genotype --> BIS11
rs6277 %>% 
  ggplot(aes (x = DRD2, y = BIS11, color=genotype)) + 
  geom_point(size = 3) + 
  geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
  theme_pubclean() +
  font("xy.text", size = 14) +
  font("ylab", size = 18) +
  font("xlab", size = 18) +
  font("legend.text", size = 18) +
  font("legend.title", size = 18) +
  xlab("\nDRD2 mRNA (batch corrected)") + ylab("BIS11 total\n") + 
  facet_wrap(. ~ Key, scales = "fixed") + 
  theme(strip.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("bespoke_v_BIS11/f4_rs6277_DRD2_v_BIS11.pdf")


#### Run Anova ####
aov(DRD2 ~ genotype, data = rs6277) %>% report()
p.snp[3] <- 0.414

aov(DRD2 ~ Key * genotype,
    data = rs6277) %>% report()
# The ANOVA (formula: DRD2 ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and small (F(1, 151) = 6.28, p =
#                                                                      0.013; Eta2 (partial) = 0.04, 95% CI [4.61e-03, 1.00])
# - The main effect of genotype is statistically not significant and small (F(2, 151) =
#                                                                             1.24, p = 0.294; Eta2 (partial) = 0.02, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and very
# small (F(2, 151) = 0.53, p = 0.591; Eta2 (partial) = 6.94e-03, 95% CI [0.00, 1.00])
# 
# Effect sizes were labelled following Field's (2013) recommendations.
p.rna[3] <- 0.591 # interaction


aov(`Psych PC1` ~ Key * genotype,
    data = rs6277) %>% report()
# The ANOVA (formula: `Psych PC1` ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and large (F(1, 262) = 181.75, p <
#                                                                      .001; Eta2 (partial) = 0.41, 95% CI [0.34, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 262)
#                                                                                = 1.87e-03, p = 0.998; Eta2 (partial) = 1.43e-05, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and very
# small (F(2, 262) = 0.12, p = 0.889; Eta2 (partial) = 9.00e-04, 95% CI [0.00, 1.00])
p.pc1[3] <- 0.889 # interaction


aov(`Psych PC2` ~ Key * genotype,
    data = rs6277) %>% report()
# The ANOVA (formula: `Psych PC2` ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and small (F(1, 262) = 14.65, p <
#                                                                      .001; Eta2 (partial) = 0.05, 95% CI [0.02, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 262)
#                                                                                = 0.72, p = 0.489; Eta2 (partial) = 5.45e-03, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and small
# (F(2, 262) = 1.48, p = 0.229; Eta2 (partial) = 0.01, 95% CI [0.00, 1.00])
p.pc2[3] <- 0.229 # interaction


aov(BIS11 ~ Key * genotype,
    data = rs6277) %>% report()
# The ANOVA (formula: BIS11 ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and large (F(1, 263) = 141.58, p <
#                                                                      .001; Eta2 (partial) = 0.35, 95% CI [0.28, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 263)
#                                                                                = 3.35e-03, p = 0.997; Eta2 (partial) = 2.55e-05, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and very
# small (F(2, 263) = 0.11, p = 0.895; Eta2 (partial) = 8.43e-04, 95% CI [0.00, 1.00])
# 
# Effect sizes were labelled following Field's (2013) recommendations.
p.bis[3] <- 0.895 # interaction

aov(SSSV ~ genotype, data = rs6277) %>% report()
# The ANOVA (formula: SSSV ~ genotype) suggests that:
#   
#   - The main effect of genotype is statistically not significant and very small (F(2, 272)
#                                                                                  = 1.13, p = 0.323; Eta2 = 8.26e-03, 95% CI [0.00, 1.00])

aov(SSSV ~ Key * genotype,
    data = rs6277) %>% report()
# The ANOVA (formula: SSSV ~ Key * genotype) suggests that:
#   
#   - The main effect of Key is statistically significant and medium (F(1, 269) = 33.66, p <
#                                                                       .001; Eta2 (partial) = 0.11, 95% CI [0.06, 1.00])
# - The main effect of genotype is statistically not significant and very small (F(2, 269)
#                                                                                = 0.84, p = 0.432; Eta2 (partial) = 6.23e-03, 95% CI [0.00, 1.00])
# - The interaction between Key and genotype is statistically not significant and small
# (F(2, 269) = 1.99, p = 0.139; Eta2 (partial) = 0.01, 95% CI [0.00, 1.00])


#### Test assumptions ####
# set output directory:
setwd("~/~r_projects/2023_kigar_cud_paper/results/f4/QC/")

#' Linearity & homogeneity check
#' looking for relatively straight lines and insignificant 
#' interaction terms!
long %>% dplyr::filter(SNP == "rs6277") %>% tidyr::drop_na(DRD2) %>% 
  tidyr::drop_na(`Psych PC1`) %>%
  ggscatter(x = "DRD2", y = "Psych PC1",
            title = "rs6277",
            facet.by  = c("genotype", "Key"), 
            short.panel.labs = FALSE) +
  stat_smooth(method = "loess", span = 0.9) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("linearitycheck_rs6277_psychpc1.pdf")

long %>% dplyr::filter(SNP == "rs6277") %>%
  tidyr::drop_na(DRD2) %>% 
  tidyr::drop_na(`Psych PC1`) %>%
  unite(col = "group", genotype, Key) %>% 
  anova_test(`Psych PC1` ~ group*DRD2)

#        Effect DFn DFd      F        p p<.05   ges
# 1      group   5 138 18.548 4.60e-14     * 0.402
# 2       DRD2   1 138  2.600 1.09e-01       0.018
# 3 group:DRD2   5 138  0.643 6.67e-01       0.023

#' Normality & homogeneity check
rs6277 <- long %>% dplyr::filter(SNP == "rs6277") %>% 
  tidyr::drop_na(SNP) %>% tidyr::drop_na(DRD2) %>% tidyr::drop_na(`Psych PC1`)

# Fit the model (covariate goes first)
model <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs6277,
            na.action = na.exclude)

# Create model diagnostic metrics
model.metrics <- augment(model) %>% select(-.hat, -.sigma, -.fitted) 

# Assess normality of residuals using shapiro wilk test
shapiro_test(model.metrics$.resid) #non-normal. p = 0.00125


# Assess homogeneity of variance with Levene's test
levene_test(.resid ~ genotype*Key, data = model.metrics)
# homogeneous. p = 0.0826

# Assess outliers
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

# one possible outlier:
# Psych PC1       DRD2 genotype Key   .resid    .cooksd .std.resid
# 1 -3.143709 0.03701038      T/T CUD -2.26301 0.06465553  -3.120389


#### Run ANCOVA ####
ancova.6277 <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs6277,
                 na.action = na.exclude)
anova(ancova.6277)

# Response: Psych PC1
#               Df Sum Sq Mean Sq F value  Pr(>F)
# DRD2           1  0.184   0.184  0.3311 0.56590
# genotype       2  1.689   0.845  1.5174 0.22280
# Key            1 54.686  54.686 98.2499 < 2e-16 ***
# genotype:Key   2  3.676   1.838  3.3019 0.03964 *
# Residuals    143 79.593   0.557

par(mfrow=c(1,3))
plot(ancova.6277, add.smooth = FALSE, which = 1)
plot(ancova.6277, which = 2)
plot(ancova.6277, add.smooth = FALSE, which = 3)
#' QQ snakes at the ends. a little bending for scale-location but looks okay.
#'  manually saved as rs6277_ancova_psychPC1_drd2.pdf
dev.off()




# multiple comparisons corrections ----------------------------------------
ls()

# order of input is rs4680, rs6280, rs6277:

p.adjust(p.bis, method="BH")
# [1] 0.285 0.285 0.895
p.adjust(p.chi, method="BH")
# [1] 0.904700 0.043620 0.010338
# p.adjust(p.pc1, method="BH")
# p.adjust(p.pc2, method="BH")
# p.adjust(p.rna, method="BH")
p.adjust(p.snp, method="BH")
# [1] 0.1905 0.1200 0.4140



#' # rs1800497 (DRD2-related)-----------------------------------------------------
#' #' this is a DRD2-related SNP. A has relatively low (15%) frequency and is
#' #' associated with lower striatal DR D2 density in healthy subjects
#' 
#' # set output directory:
#' setwd("~/~r_projects/2023_kigar_cud_paper/results/f4/")
#' 
#' # filter to SNP of interest: 
#' rs1800497 <- long %>% dplyr::filter(SNP == "rs1800497") # 277 
#' 
#' #### Chi square for frequency ####
#' table(rs1800497$Key, rs1800497$genotype)
#' #          A/A A/G G/G
#' # control   5  42 103
#' # CUD      11  42  74
#' 
#' rs1800497 %>%  
#'   ggplot(aes(x = Key, fill = genotype)) + 
#'   geom_bar(position = "fill") +
#'   labs(y = "Proportion") + theme_pubclean() +
#'   ggpubr::font("xy.text", size = 18) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) + 
#'   scale_fill_manual(values = c("#cccccc", "#666666", "#000000")) +
#'   xlab("") + ylab("Frequency\n") 
#' ggsave("chi_square_genotype/rs1800497_frequency_stackedbar.pdf")
#' 
#' 
#' # create count summary of group x genotype, stored in temp object M:
#' M <- rs1800497 %>% select(Key, genotype) %>% 
#'   group_by(Key, genotype) %>% 
#'   summarise(n = n())
#' 
#' # extract counts:
#' col <- M$n
#' 
#' # convert counts into a 2 row x 2 column matrix:
#' snp <- matrix(col, nrow = 3, ncol = 2)
#' 
#' # name rows and columns:
#' dimnames(snp) = list(genotype = M$genotype[1:3], 
#'                      group = c("control", "CUD"))
#' 
#' # perform chi squared test:
#' chisq.test(snp)
#' # X-squared = 5.127, df = 2, p-value = 0.07703
#' p.chi[4] <- 0.07703
#' 
#' #' (what if you combine genotypes with A allele together?)
#' # convert counts into a 2 row x 2 column matrix:
#' snp2 <- matrix(c(47,103,53,74), nrow = 2, ncol = 2)
#' 
#' # name rows and columns:
#' dimnames(snp2) = list(genotype = c("A/_", "G/G"), 
#'                      group = c("control", "CUD"))
#' 
#' # perform chi squared test:
#' chisq.test(snp2)
#' # X-squared = 2.7889, df = 1, p-value = 0.09492
#' 
#' 
#' #### Visualize relationships ####
#' 
#' # graph genotype x Key --> Psych PC1
#' rs1800497 %>% ggplot(aes(x = genotype, y = `Psych PC1`, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs1800497") + ylab("Psych PC1\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("snp_v_psychPC1/f4_rs1800497_psychPC1.pdf")
#' 
#' # graph genotype x Key --> Psych PC2
#' rs1800497 %>% ggplot(aes(x = genotype, y = `Psych PC2`, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs1800497") + ylab("Psych PC2\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("snp_v_psychPC2/f4_rs1800497_psychPC2.pdf")
#' 
#' # graph genotype x Key --> DRD2
#' rs1800497 %>% ggplot(aes(x = genotype, y = DRD2, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs1800497") + ylab("DRD2 mRNA (batch corrected)\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("bespoke_snp_gene/f4_rs1800497_DRD2.pdf")
#' 
#' # graph genotype x Key --> BIS11
#' rs1800497 %>% ggplot(aes(x = genotype, y = BIS11, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs1800497") + ylab("BIS11 total\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("snp_v_bis11/f4_rs1800497_BIS11.pdf")
#' 
#' # graph PCR PC2 x genotype --> Psych PC1
#' rs1800497 %>% 
#'   ggplot(aes (x = `PCR PC2`, y = `Psych PC1`, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nPCR PC2") + ylab("Psych PC1\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("pcrPC2_v_psychPC1/f4_rs1800497_pcrPC2_v_psychPC1.pdf")
#' 
#' # graph PCR PC2 x genotype --> Psych PC2
#' rs1800497 %>% 
#'   ggplot(aes (x = `PCR PC2`, y = `Psych PC2`, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nPCR PC2") + ylab("Psych PC2\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("pcrPC2_v_psychPC2/f4_rs1800497_pcrPC2_v_psychPC2.pdf")
#' 
#' # graph DRD2 x genotype --> Psych PC1
#' rs1800497 %>% 
#'   ggplot(aes (x = DRD2, y = `Psych PC1`, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nDRD2 mRNA (batch corrected)") + ylab("Psych PC1\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("bespoke_v_psychPC1/f4_rs1800497_DRD2_v_psychPC1.pdf")
#' 
#' 
#' # graph DRD2 x genotype --> Psych PC2
#' rs1800497 %>% 
#'   ggplot(aes (x = DRD2, y = `Psych PC2`, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nDRD2 mRNA (batch corrected)") + ylab("Psych PC2\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("bespoke_v_psychPC2/f4_rs1800497_DRD2_v_psychPC2.pdf")
#' 
#' # graph DRD2 x genotype --> BIS11
#' rs1800497 %>% 
#'   ggplot(aes (x = DRD2, y = BIS11, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nDRD2 mRNA (batch corrected)") + ylab("BIS11 total\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("bespoke_v_BIS11/f4_rs1800497_DRD2_v_BIS11.pdf")
#' 
#' #### Test assumptions ####
#' # set output directory:
#' setwd("~/~r_projects/2023_kigar_cud_paper/results/f4/QC/")
#' 
#' #' Linearity & homogeneity check
#' #' looking for relatively straight lines and insignificant 
#' #' interaction terms!
#' long %>% dplyr::filter(SNP == "rs1800497") %>% tidyr::drop_na(DRD2) %>% 
#'   tidyr::drop_na(`Psych PC1`) %>%
#'   dplyr::filter(genotype != "A/A") %>% 
#'   ggscatter(x = "DRD2", y = "Psych PC1",
#'             title = "rs1800497",
#'             facet.by  = c("genotype", "Key"), 
#'             short.panel.labs = FALSE) +
#'   stat_smooth(method = "loess", span = 0.9) + 
#'   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#' ggsave("linearitycheck_rs1800497_psychpc1.pdf") 
#' # note: there aren't enough of the A/A genotype for LOESS; excluded
#' 
#' long %>% dplyr::filter(SNP == "rs1800497") %>%
#'   tidyr::drop_na(DRD2) %>% 
#'   tidyr::drop_na(`Psych PC1`) %>%
#'   unite(col = "group", genotype, Key) %>% 
#'   anova_test(`Psych PC1` ~ group*DRD2)
#' 
#' #      Effect DFn DFd      F        p p<.05   ges
#' # 1      group   5 138 20.806 1.87e-15     * 0.430
#' # 2       DRD2   1 138  2.751 9.90e-02       0.020
#' # 3 group:DRD2   5 138  1.588 1.67e-01       0.054
#' 
#' #' Normality & homogeneity check
#' rs1800497 <- long %>% dplyr::filter(SNP == "rs1800497") %>% 
#'     tidyr::drop_na(SNP) %>% tidyr::drop_na(DRD2) %>% tidyr::drop_na(`Psych PC1`)
#' 
#' # Fit the model (covariate goes first)
#' model <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs1800497,
#'             na.action = na.exclude)
#' 
#' # Create model diagnostic metrics
#' model.metrics <- augment(model) %>% select(-.hat, -.sigma, -.fitted) 
#' 
#' # Assess normality of residuals using shapiro wilk test
#' shapiro_test(model.metrics$.resid) #non-normal. p = 0.0192
#' 
#' # Assess homogeneity of variance with Levene's test
#' levene_test(.resid ~ genotype*Key, data = model.metrics)
#' # non-homogeneous. p = 0.0139
#' 
#' # Assess outliers
#' model.metrics %>% 
#'   filter(abs(.std.resid) > 3) %>%
#'   as.data.frame()
#' 
#' # no outliers
#' 
#' 
#' #### Run ANCOVA ####
#' aov.1800497 <- anova_test(`Psych PC1` ~ DRD2 + genotype*Key, data = rs1800497)
#' get_anova_table(aov.1800497)
#' # Effect DFn DFd      F                     p p<.05   ges
#' # 1         DRD2   1 143  2.696 0.1029999999999999943       0.019
#' # 2     genotype   2 143  1.424 0.2439999999999999947       0.020
#' # 3          Key   1 143 90.642 0.0000000000000000602     * 0.388
#' # 4 genotype:Key   2 143  1.250 0.2899999999999999800       0.017
#' 
#' par(mfrow=c(2,2))
#' plot(aov.1800497) 
#' #' QQ snakes at the ends. residuals are okay -some outlier cases, scale-location
#' #' has downward slope. manually saved as rs1800497_ancova_psychPC1_drd2.pdf
#' dev.off()
#' 
#' #' Normality & homogeneity check (no A/A genotype)
#' alt_rs1800497 <- long %>% dplyr::filter(SNP == "rs1800497") %>% 
#'   dplyr::filter(genotype != "A/A") %>% 
#'   tidyr::drop_na(SNP) %>% tidyr::drop_na(DRD2) %>% tidyr::drop_na(`Psych PC1`)
#' 
#' # Fit the model (covariate goes first)
#' model <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = alt_rs1800497,
#'             na.action = na.exclude)
#' 
#' # Create model diagnostic metrics
#' model.metrics <- augment(model) %>% select(-.hat, -.sigma, -.fitted) 
#' 
#' # Assess normality of residuals using shapiro wilk test
#' shapiro_test(model.metrics$.resid) #normal. p = 0.0511
#' 
#' # Assess homogeneity of variance with Levene's test
#' levene_test(.resid ~ genotype*Key, data = model.metrics)
#' # non-homogeneous. p = 0.0113
#' 
#' # Assess outliers
#' model.metrics %>% 
#'   filter(abs(.std.resid) > 3) %>%
#'   as.data.frame()
#' 
#' # potential outlier:
#' # Psych PC1        DRD2 genotype Key    .resid    .cooksd .std.resid
#' # 1 -2.804579 0.007381844      G/G CUD -2.224197 0.05012467  -3.099419
#' 
#' #' `run the model:`
#' ancova.1800497 <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs1800497,
#'                  na.action = na.exclude)
#' anova(ancova.1800497)
#' 
#' # Response: Psych PC1
#' #               Df Sum Sq Mean Sq F value  Pr(>F)
#' # DRD2           1  0.184   0.184  0.3311 0.56590
#' # genotype       2  1.689   0.845  1.5174 0.22280
#' # Key            1 54.686  54.686 98.2499 < 2e-16 ***
#' # genotype:Key   2  3.676   1.838  3.3019 0.03964 *
#' # Residuals    143 79.593   0.557
#' 
#' par(mfrow=c(1,3))
#' plot(ancova.1800497, add.smooth = FALSE, which = 1)
#' plot(ancova.1800497, which = 2)
#' plot(ancova.1800497, add.smooth = FALSE, which = 3)
#' #' QQ snakes at the ends. a little bending for scale-location but looks okay.
#' #'  manually saved as rs1800497_ancova_psychPC1_drd2.pdf
#' dev.off()
#' 
#' 
#' # rs12364283 (DRD2-related)---------------------------------------------------
#' # this is a DRD2-related SNP. Unclear what associations to expect.
#' 
#' # set output directory:
#' setwd("~/~r_projects/2023_kigar_cud_paper/results/f4/")
#' 
#' # filter to SNP of interest: 
#' rs12364283 <- long %>% dplyr::filter(SNP == "rs12364283") #277 
#' 
#' #### Chi square for frequency ####
#' table(rs12364283$Key, rs12364283$genotype)
#' #         A/A A/G G/G
#' # control 124  23   3
#' # CUD     108  17   2
#' 
#' rs12364283 %>%  
#'   ggplot(aes(x = Key, fill = genotype)) + 
#'   geom_bar(position = "fill") +
#'   labs(y = "Proportion") + theme_pubclean() +
#'   ggpubr::font("xy.text", size = 18) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) + 
#'   scale_fill_manual(values = c("#cccccc", "#666666", "#000000")) +
#'   xlab("") + ylab("Frequency\n") 
#' ggsave("chi_square_genotype/rs12364283_frequency_stackedbar.pdf")
#' 
#' 
#' # create count summary of group x genotype, stored in temp object M:
#' M <- rs12364283 %>% select(Key, genotype) %>% 
#'   group_by(Key, genotype) %>% 
#'   summarise(n = n())
#' 
#' # extract counts:
#' col <- M$n
#' 
#' # convert counts into a 2 row x 2 column matrix:
#' snp <- matrix(col, nrow = 3, ncol = 2)
#' 
#' # name rows and columns:
#' dimnames(snp) = list(genotype = M$genotype[1:3], 
#'                      group = c("control", "CUD"))
#' 
#' # perform chi squared test:
#' chisq.test(snp)
#' # X-squared = 0.29574, df = 2, p-value = 0.8625
#' p.chi[5] <- 0.8625
#' 
#' 
#' #### Visualize relationships ####
#' 
#' # graph genotype x Key --> Psych PC1
#' rs12364283 %>% ggplot(aes(x = genotype, y = `Psych PC1`, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs12364283") + ylab("Psych PC1\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("snp_v_psychPC1/f4_rs12364283_psychPC1.pdf")
#' 
#' # graph genotype x Key --> Psych PC2
#' rs12364283 %>% ggplot(aes(x = genotype, y = `Psych PC2`, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs12364283") + ylab("Psych PC2\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("snp_v_psychPC2/f4_rs12364283_psychPC2.pdf")
#' 
#' # graph genotype x Key --> DRD2
#' rs12364283 %>% ggplot(aes(x = genotype, y = DRD2, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs12364283") + ylab("DRD2 mRNA (batch corrected)\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("bespoke_snp_gene/f4_rs12364283_DRD2.pdf")
#' 
#' # graph genotype x Key --> BIS11
#' rs12364283 %>% ggplot(aes(x = genotype, y = BIS11, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs12364283") + ylab("BIS11 total\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("snp_v_bis11/f4_rs12364283_BIS11.pdf")
#' 
#' # graph PCR PC2 x genotype --> Psych PC1
#' rs12364283 %>% 
#'   ggplot(aes (x = `PCR PC2`, y = `Psych PC1`, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nPCR PC2") + ylab("Psych PC1\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("pcrPC2_v_psychPC1/f4_rs12364283_pcrPC2_v_psychPC1.pdf")
#' 
#' # graph PCR PC2 x genotype --> Psych PC2
#' rs12364283 %>% 
#'   ggplot(aes (x = `PCR PC2`, y = `Psych PC2`, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nPCR PC2") + ylab("Psych PC2\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("pcrPC2_v_psychPC2/f4_rs12364283_pcrPC2_v_psychPC2.pdf")
#' 
#' # graph DRD2 x genotype --> Psych PC1
#' rs12364283 %>% 
#'   ggplot(aes (x = DRD2, y = `Psych PC1`, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nDRD2 mRNA (batch corrected)") + ylab("Psych PC1\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("bespoke_v_psychPC1/f4_rs12364283_DRD2_v_psychPC1.pdf")
#' 
#' 
#' # graph DRD2 x genotype --> Psych PC2
#' rs12364283 %>% 
#'   ggplot(aes (x = DRD2, y = `Psych PC2`, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nDRD2 mRNA (batch corrected)") + ylab("Psych PC2\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("bespoke_v_psychPC2/f4_rs12364283_DRD2_v_psychPC2.pdf")
#' 
#' # graph DRD2 x genotype --> BIS11
#' rs12364283 %>% 
#'   ggplot(aes (x = DRD2, y = BIS11, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nDRD2 mRNA (batch corrected)") + ylab("BIS11 total\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("bespoke_v_BIS11/f4_rs12364283_DRD2_v_BIS11.pdf")
#' 
#' 
#' #### Test assumptions ####
#' # set output directory:
#' setwd("~/~r_projects/2023_kigar_cud_paper/results/f4/QC/")
#' 
#' #' Linearity & homogeneity check
#' #' looking for relatively straight lines and insignificant 
#' #' interaction terms!
#' long %>% dplyr::filter(SNP == "rs12364283") %>% tidyr::drop_na(DRD2) %>% 
#'   tidyr::drop_na(`Psych PC1`) %>%
#'   dplyr::filter(genotype != "G/G") %>% 
#'   ggscatter(x = "DRD2", y = "Psych PC1",
#'             title = "rs12364283",
#'             facet.by  = c("genotype", "Key"), 
#'             short.panel.labs = FALSE) +
#'   stat_smooth(method = "loess", span = 0.9)
#' ggsave("linearitycheck_rs12364283_psychpc1.pdf")
#' #' note: G/G genotype is too rare for computations (3 people total); 
#' #' G/A also pretty rare
#' 
#' long %>% dplyr::filter(SNP == "rs12364283") %>%
#'   tidyr::drop_na(DRD2) %>% 
#'   tidyr::drop_na(`Psych PC1`) %>%
#'   dplyr::filter(genotype != "G/G") %>% 
#'   unite(col = "group", genotype, Key) %>% 
#'   anova_test(`Psych PC1` ~ group*DRD2)
#' 
#' #      Effect DFn DFd      F        p p<.05      ges
#' # 1      group   3 139 29.995 5.14e-15     * 3.93e-01
#' # 2       DRD2   1 139  1.950 1.65e-01       1.40e-02
#' # 3 group:DRD2   3 139  0.003 1.00e+00       6.56e-05
#' 
#' #' Normality & homogeneity check
#' rs12364283 <- long %>% dplyr::filter(SNP == "rs12364283") %>% 
#'   tidyr::drop_na(SNP) %>% tidyr::drop_na(DRD2) %>% tidyr::drop_na(`Psych PC1`)
#' 
#' # Fit the model (covariate goes first)
#' model <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs12364283,
#'             na.action = na.exclude)
#' 
#' # Create model diagnostic metrics
#' model.metrics <- augment(model) %>% select(-.hat, -.sigma, -.fitted) 
#' 
#' # Assess normality of residuals using shapiro wilk test
#' shapiro_test(model.metrics$.resid) #non-normal. p = 0.000731
#' 
#' 
#' # Assess homogeneity of variance with Levene's test
#' levene_test(.resid ~ genotype*Key, data = model.metrics)
#' # non-homogeneous. p = 0.000663
#' 
#' 
#' # Assess outliers
#' model.metrics %>% 
#'   filter(abs(.std.resid) > 3) %>%
#'   as.data.frame()
#' 
#' # one possible outlier:
#' # Psych PC1       DRD2 genotype Key    .resid    .cooksd .std.resid
#' # 1 -3.143709 0.03701038      A/A CUD -2.542769 0.02688876  -3.357516
#' 
#' #' run the model:
#' aov.12364283 <- anova_test(`Psych PC1` ~ DRD2 + genotype*Key, data = rs12364283)
#' get_anova_table(aov.12364283)
#' # Effect DFn DFd      F                     p p<.05      ges
#' # 1         DRD2   1 143  2.022 0.1570000000000000007       0.014000
#' # 2     genotype   2 143  0.952 0.3880000000000000115       0.013000
#' # 3          Key   1 143 94.099 0.0000000000000000208     * 0.397000
#' # 4 genotype:Key   2 143  0.066 0.9360000000000000542       0.000926
#' 
#' par(mfrow=c(2,2))
#' plot(aov.12364283) #Warning: not plotting observations with leverage one
#' #' QQ snakes at the ends. residuals are okay -some outlier cases, scale-location
#' #' has downward slope. manually saved as rs12364283_ancova_psychPC1_drd2.pdf
#' dev.off()
#' 
#' #' Normality & homogeneity check (no A/A genotype)
#' alt_rs12364283 <- long %>% dplyr::filter(SNP == "rs12364283") %>% 
#'   dplyr::filter(genotype != "G/G") %>% 
#'   tidyr::drop_na(SNP) %>% tidyr::drop_na(DRD2) %>% tidyr::drop_na(`Psych PC1`)
#' 
#' # Fit the model (covariate goes first)
#' model <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = alt_rs12364283,
#'             na.action = na.exclude)
#' 
#' # Create model diagnostic metrics
#' model.metrics <- augment(model) %>% select(-.hat, -.sigma, -.fitted) 
#' 
#' # Assess normality of residuals using shapiro wilk test
#' shapiro_test(model.metrics$.resid) #normal. p = 0.00111
#' 
#' # Assess homogeneity of variance with Levene's test
#' levene_test(.resid ~ genotype*Key, data = model.metrics)
#' # non-homogeneous. p = 0.000428
#' 
#' # Assess outliers
#' model.metrics %>% 
#'   filter(abs(.std.resid) > 3) %>%
#'   as.data.frame()
#' 
#' # potential outlier:
#' # Psych PC1       DRD2 genotype Key   .resid    .cooksd .std.resid
#' # 1 -3.143709 0.03701038      A/A CUD -2.54278 0.03738208   -3.34579
#' 
#' 
#' #### Run ANCOVA ####
#' ancova.12364283 <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs12364283,
#'                  na.action = na.exclude)
#' anova(ancova.12364283)
#' 
#' # Response: Psych PC1
#' #               Df Sum Sq Mean Sq F value  Pr(>F)
#' # DRD2           1  0.184   0.184  0.3311 0.56590
#' # genotype       2  1.689   0.845  1.5174 0.22280
#' # Key            1 54.686  54.686 98.2499 < 2e-16 ***
#' # genotype:Key   2  3.676   1.838  3.3019 0.03964 *
#' # Residuals    143 79.593   0.557
#' 
#' par(mfrow=c(1,3))
#' plot(ancova.12364283, add.smooth = FALSE, which = 1)
#' plot(ancova.12364283, which = 2)
#' plot(ancova.12364283, add.smooth = FALSE, which = 3)
#' #' QQ snakes at the ends. a little bending for scale-location but looks okay.
#' #'  manually saved as rs12364283_ancova_psychPC1_drd2.pdf
#' dev.off()
#' 
#' 
#' # rs686 (DRD1-related) -------------------------------------------------------------------
#' #' this is a DRD1-related SNP. A should be related to higher 
#' #' DRD1 expression. Associated with nicotine dependence, 
#' #' alcohol dependence, tobacco smoking in SZ
#' 
#' # set output directory:
#' setwd("~/~r_projects/2023_kigar_cud_paper/results/f4/")
#' 
#' # filter to SNP of interest: 
#' rs686 <- long %>% dplyr::filter(SNP == "rs686") #276
#' 
#' #### Chi square for frequency ####
#' table(rs686$Key, rs686$genotype)
#' #         A/A A/G G/G
#' # control  47  80  23
#' # CUD      38  61  27
#' 
#' rs686 %>%  
#'   ggplot(aes(x = Key, fill = genotype)) + 
#'   geom_bar(position = "fill") +
#'   labs(y = "Proportion") + theme_pubclean() +
#'   ggpubr::font("xy.text", size = 18) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) + 
#'   scale_fill_manual(values = c("#cccccc", "#666666", "#000000")) +
#'   xlab("") + ylab("Frequency\n") 
#' ggsave("chi_square_genotype/rs686_frequency_stackedbar.pdf")
#' 
#' 
#' # create count summary of group x genotype, stored in temp object M:
#' M <- rs686 %>% select(Key, genotype) %>% 
#'   group_by(Key, genotype) %>% 
#'   summarise(n = n())
#' 
#' # extract counts:
#' col <- M$n
#' 
#' # convert counts into a 2 row x 2 column matrix:
#' snp <- matrix(col, nrow = 3, ncol = 2)
#' 
#' # name rows and columns:
#' dimnames(snp) = list(genotype = M$genotype[1:3], 
#'                      group = c("control", "CUD"))
#' 
#' # perform chi squared test:
#' chisq.test(snp)
#' # X-squared = 1.7596, df = 2, p-value = 0.4149
#' p.chi[6] <- 0.4149
#' 
#' 
#' #### Visualize relationships ####
#' 
#' # graph genotype x Key --> Psych PC1
#' rs686 %>% ggplot(aes(x = genotype, y = `Psych PC1`, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs686") + ylab("Psych PC1\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("snp_v_psychPC1/f4_rs686_psychPC1.pdf")
#' 
#' # graph genotype x Key --> Psych PC2
#' rs686 %>% ggplot(aes(x = genotype, y = `Psych PC2`, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs686") + ylab("Psych PC2\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("snp_v_psychPC2/f4_rs686_psychPC2.pdf")
#' 
#' 
#' # graph genotype x Key --> BIS11
#' rs686 %>% ggplot(aes(x = genotype, y = BIS11, col = genotype)) + 
#'   geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) + 
#'   theme_pubclean() +
#'   ggpubr::font("xy.text", size = 14) +
#'   ggpubr::font("ylab", size = 18) +
#'   ggpubr::font("xlab", size = 18) +
#'   ggpubr::font("legend.title", size = 18) +
#'   ggpubr::font("legend.text", size = 18) +
#'   xlab("\nrs686") + ylab("BIS11 total\n") +
#'   theme(legend.position="top",
#'         strip.text = element_text(size = rel(1.5))) +
#'   facet_wrap(. ~ Key, scales = "fixed") 
#' ggsave("snp_v_bis11/f4_rs686_BIS11.pdf")
#' 
#' # graph PCR PC2 x genotype --> Psych PC1
#' rs686 %>% 
#'   ggplot(aes (x = `PCR PC2`, y = `Psych PC1`, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nPCR PC2") + ylab("Psych PC1\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("pcrPC2_v_psychPC1/f4_rs686_pcrPC2_v_psychPC1.pdf")
#' 
#' # graph PCR PC2 x genotype --> Psych PC2
#' rs686 %>% 
#'   ggplot(aes (x = `PCR PC2`, y = `Psych PC2`, color=genotype)) + 
#'   geom_point(size = 3) + 
#'   geom_smooth(method=lm, se=T, linewidth = 2, aes(fill = genotype)) +
#'   theme_pubclean() +
#'   font("xy.text", size = 14) +
#'   font("ylab", size = 18) +
#'   font("xlab", size = 18) +
#'   font("legend.text", size = 18) +
#'   font("legend.title", size = 18) +
#'   xlab("\nPCR PC2") + ylab("Psych PC2\n") + 
#'   facet_wrap(. ~ Key, scales = "fixed") + 
#'   theme(strip.text = element_text(size = rel(1.5)),
#'         axis.text.x = element_text(angle = 45, hjust=1)) 
#' ggsave("pcrPC2_v_psychPC2/f4_rs686_pcrPC2_v_psychPC2.pdf")
#' 
#' 
#' #### Test assumptions ####
#' # set output directory:
#' setwd("~/~r_projects/2023_kigar_cud_paper/results/f4/QC/")
#' 
#' #' Linearity & homogeneity check
#' #' looking for relatively straight lines and insignificant 
#' #' interaction terms!
#' long %>% dplyr::filter(SNP == "rs686") %>% tidyr::drop_na(DRD2) %>% 
#'   tidyr::drop_na(`Psych PC1`) %>%
#'   ggscatter(x = "DRD2", y = "Psych PC1",
#'             title = "rs686",
#'             facet.by  = c("genotype", "Key"), 
#'             short.panel.labs = FALSE) +
#'   stat_smooth(method = "loess", span = 0.9)
#' ggsave("linearitycheck_rs686_psychpc1.pdf")
#' 
#' long %>% dplyr::filter(SNP == "rs686") %>%
#'   tidyr::drop_na(DRD2) %>% 
#'   tidyr::drop_na(`Psych PC1`) %>%
#'   unite(col = "group", genotype, Key) %>% 
#'   anova_test(`Psych PC1` ~ group*DRD2)
#' 
#' #       Effect DFn DFd      F        p p<.05   ges
#' # 1      group   5 138 20.928 1.59e-15     * 0.431
#' # 2       DRD2   1 138  2.629 1.07e-01       0.019
#' # 3 group:DRD2   5 138  0.139 9.83e-01       0.005
#' 
#' #' Normality & homogeneity check
#' rs686 <- long %>% dplyr::filter(SNP == "rs686") %>% 
#'   tidyr::drop_na(SNP) %>% tidyr::drop_na(DRD2) %>% tidyr::drop_na(`Psych PC1`)
#' 
#' # Fit the model (covariate goes first)
#' model <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs686,
#'             na.action = na.exclude)
#' 
#' # Create model diagnostic metrics
#' model.metrics <- augment(model) %>% select(-.hat, -.sigma, -.fitted) 
#' 
#' # Assess normality of residuals using shapiro wilk test
#' shapiro_test(model.metrics$.resid) #non-normal. p = 0.00106
#' 
#' # Assess homogeneity of variance with Levene's test
#' levene_test(.resid ~ genotype*Key, data = model.metrics)
#' # homogeneous. p = 0.0125
#' 
#' # Assess outliers
#' model.metrics %>% 
#'   filter(abs(.std.resid) > 3) %>%
#'   as.data.frame()
#' 
#' # one possible outlier:
#' # Psych PC1       DRD2 genotype Key   .resid    .cooksd .std.resid
#' # 1 -3.143709 0.03701038      A/A CUD -2.373264 0.07212519  -3.255955
#' 
#' 
#' #### Run ANCOVA ####
#' ancova.686 <- lm(`Psych PC1` ~ DRD2 + genotype*Key, data = rs686,
#'                  na.action = na.exclude)
#' anova(ancova.686)
#' 
#' # Response: Psych PC1
#' #               Df Sum Sq Mean Sq F value  Pr(>F)
#' # DRD2           1  0.184   0.184  0.3311 0.56590
#' # genotype       2  1.689   0.845  1.5174 0.22280
#' # Key            1 54.686  54.686 98.2499 < 2e-16 ***
#' # genotype:Key   2  3.676   1.838  3.3019 0.03964 *
#' # Residuals    143 79.593   0.557
#' 
#' par(mfrow=c(1,3))
#' plot(ancova.686, add.smooth = FALSE, which = 1)
#' plot(ancova.686, which = 2)
#' plot(ancova.686, add.smooth = FALSE, which = 3)
#' #' QQ snakes at the ends. a little bending for scale-location but looks okay.
#' #'  manually saved as rs686_ancova_psychPC1_drd2.pdf
#' dev.off()




