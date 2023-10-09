#' this script is used to generate Figure 3, Figure 4, and accompanying
#' stats
#' 
#' NB: different genes have different calibrator samples used in the 
#' ddCt calculation. These have been converted to NA data, so each 
#' gene may have varying NA sample data.

# Set-up ------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(psych)
library(factoextra)
library(pracma)
library(ggplot2)
library(ggpubr)

# Import & tidy data 
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
all <- read_csv("5_2023_Addiction_clean.csv")

# Set working directory for image output:
setwd("~/~r_projects/2023_kigar_cud_paper/results/f3/")

# Create psychometric PCA -------------------------------------------------

# Select psychometric data and drop nas
psych <- all %>% 
  drop_na(Key) %>% 
  select(BDI_total, STAI_S_total, STAI_T_total, 
         BIS11_attention, BIS11_motor, BIS11_nonplan, 
         SSSV_Tas, SSSV_Es, SSSV_Dis, SSSV_Bs) %>% 
  rename(`SSSV-TAS` = SSSV_Tas,
         `SSSV-ES` = SSSV_Es,
         `SSSV-DIS` = SSSV_Dis,
         `SSSV-BS` = SSSV_Bs,
         `BIS11-AI` = BIS11_attention,
         `BIS11-MI` = BIS11_motor,
         `BIS11-NPI` = BIS11_nonplan,
         `STAI-S` = STAI_S_total,
         `STAI-T` = STAI_T_total,
         BDI = BDI_total) %>%
  drop_na()
# 338 observations (from 349 starting), so 11 have incomplete data.

KMO_psych <- KMO(psych)
KMO_psych # MSA = 0.84, which is good.
bartlett.test(psych) ## p < 2.2e-16, which is good.

# generate the PCA
pca_psych <- prcomp(psych, center = T, scale. = T)

# examine eigenvalues/make scree plot:
eig_val_psych <- get_eigenvalue(pca_psych)
eig_val_psych # first 2 are >1

scree_psych <- fviz_eig(pca_psych,
                        barfill = "#738b9b",
                        barcolor = "#738b9b",
                        font.main = 36,
                        main = "",
                        font.x = 36, 
                        font.tickslab = 24,
                        xlab = "\nComponent number",
                        font.y = 36, 
                        ylab = "% explained variance\n",
                        ggtheme = theme_linedraw()) 
scree_psych
ggsave("f3_screeplot_psychometric_pca.pdf")

# make graph of how individual variables contribute to PC1&2:
circ.psych <- fviz_pca_var(pca_psych, axes = c(1,2),
                           col.var = "contrib",
                           repel = T,
                           col.circle = "#dde2e7",
                           title = "",
                           arrowsize = 1,
                           labelsize = 11,
                           ggtheme = theme_light()) +
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=9) +
  theme(text = element_text(size = 24),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 36)) 
circ.psych
ggsave("f3_factormap_psychometric_pca.pdf")

# pull out identical columns + key so that everything lines up
# with the psych data frame:
temp <- all %>% drop_na(Key) %>% 
  select(3:4, Key, BDI_total, STAI_S_total, STAI_T_total, 
         BIS11_attention, BIS11_motor, BIS11_nonplan, 
         SSSV_Tas, SSSV_Es, SSSV_Dis, SSSV_Bs) %>% drop_na()

# make graph with individual values
psych.graph <- fviz_pca_ind(pca_psych, axes = c(1,2), #axes may change depending on which components you choose to keep or want to plot
             label = "true",
             habillage = temp$Key,
             addEllipses = T, 
             repel = T, 
             geom = "point",
             pointsize = 5,
             title = "",
             ggtheme = theme_minimal() + 
               theme(text = element_text(size = 24),
                     axis.text = element_text(size = 36),
                     axis.title = element_text(size = 36))) + 
  scale_color_manual(values = c('#00BFC4','#F8766D')) +
  scale_fill_manual(values = c('#00BFC4','#F8766D')) 
psych.graph
ggsave("f3_individualmap_psychometric_pca.pdf")

# now you want to see if these new 'scores' correlate with gene expression
#rotate to get individual PK scores
ncomp_psych <- 2 # this is the # of components you choose to keep

raw_loadings_psych <- pca_psych$rotation[,1:ncomp_psych] %*% 
  diag(pca_psych$sdev, ncomp_psych, ncomp_psych)

# applies a function to raw loadings to make them 
# 'more biologically relevant'. this is the most widely used way of 'rotating'
# the data apparently.
rotated_loadings_psych <- varimax(raw_loadings_psych)$loadings

inv_loadings_psych <- t(pracma::pinv(rotated_loadings_psych))

scores_psych <- scale(psych) %*% inv_loadings_psych

# convert to a tibble; rename columns to something sensible:
scores_psych <- as_tibble(scores_psych)
scores_psych %<>% rename(psych_PC1 = V1,
                         psych_PC2 = V2)

# merge with temp database:
temp <- cbind(temp, scores_psych)

# drop extraneous columns in temp database to merge into complete dataframe:
temp %<>% select(1:2,14:15)

all <- full_join(all, temp, by = c("Patient_ID", "Study_ID"))
colnames(all)


# Create gene expression PCA -----------------------------------------------

# use 2 group/technical variables regressed out values; 
# filter out users and relatives: 
pcr <- all %>%
  drop_na(Key) %>% 
  select(68:71) %>% 
  rename(DRD2 = ddDRD2.res.sq,
         DRD3 = ddDRD3.res.sq,
         DRD4 = ddDRD4.res.sq,
         COMT = ddCOMT.res.sq ) %>% 
  drop_na() # n = 164

# test assumptions for running PCA:
KMO_pcr <- KMO(pcr) # 0.52 - bad, borderline unacceptable
KMO_pcr
bartlett.test(pcr) # p < 2.2e-16

# run the pca: 
pca_pcr <- prcomp(pcr, center = T, scale = T)

# examine eigenvalues/make scree plot:
eig_val_pcr <- get_eigenvalue(pca_pcr)
eig_val_pcr # first 2 are >1

scree_pcr <- fviz_eig(pca_pcr,
                        barfill = "#738b9b",
                        barcolor = "#738b9b",
                        font.main = 36,
                        main = "",
                        font.x = 36, 
                        font.tickslab = 24,
                        xlab = "\nComponent number",
                        font.y = 36, 
                        ylab = "% explained variance\n",
                        ggtheme = theme_linedraw()) 
scree_pcr
ggsave("f3_screeplot_pcr_pca.pdf")

# make graph of how individual variables contribute to PC1&2:
circ.pcr <- fviz_pca_var(pca_pcr, axes = c(1,2),
                           col.var = "contrib",
                           repel = T,
                           col.circle = "#dde2e7",
                           title = "",
                           arrowsize = 1,
                           labelsize = 11,
                           ggtheme = theme_light()) +
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=20) +
  theme(text = element_text(size = 24),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 36)) 
circ.pcr
ggsave("f3_factormap_pcr_pca.pdf")


# pull out identical columns + key so that everything lines up
# with the pcr data frame:
temp <- all %>% drop_na(Key) %>% 
  select(3:4, Key, 64:67) %>% drop_na()


# make graph with individual values
pcr.graph <- fviz_pca_ind(pca_pcr, axes = c(1,2), #axes may change depending on which components you choose to keep or want to plot
             label = "true",
             habillage = temp$Key,
             addEllipses = T, 
             repel = T, 
             geom = "point",
             pointsize = 5,
             title = "",
             ggtheme = theme_minimal() + 
               theme(text = element_text(size = 24),
                     axis.text = element_text(size = 36),
                     axis.title = element_text(size = 36))) + 
  scale_color_manual(values = c('#00BFC4','#F8766D')) +
  scale_fill_manual(values = c('#00BFC4','#F8766D')) 
pcr.graph
ggsave("f3_individualmap_pcr_pca.pdf")

# now you want to see if these new 'scores' correlate with gene expression
#rotate to get individual PK scores
ncomp_pcr <- 2 # this is the # of components you choose to keep

raw_loadings_pcr <- pca_pcr$rotation[,1:ncomp_pcr] %*% 
  diag(pca_pcr$sdev, ncomp_pcr, ncomp_pcr)

# applies a function to raw loadings to make them 
# 'more biologically relevant'. this is the most widely used way of 'rotating'
# the data apparently.
rotated_loadings_pcr <- varimax(raw_loadings_pcr)$loadings

inv_loadings_pcr <- t(pracma::pinv(rotated_loadings_pcr))

scores_pcr <- scale(pcr) %*% inv_loadings_pcr

# convert to a tibble; rename columns to something sensible:
scores_pcr <- as_tibble(scores_pcr)
scores_pcr %<>% rename(pcr_PC1 = V1,
                         pcr_PC2 = V2)

# merge with temp database:
temp <- cbind(temp, scores_pcr)

# drop extraneous columns in temp database to merge into complete dataframe:
temp %<>% select(1:2,8:9)

all <- full_join(all, temp, by = c("Patient_ID", "Study_ID"))
colnames(all)


# Group comparison: psychometric PCs --------------------------------------
# create empty vector for planned comparisons:
p.values <- c(NA, NA, NA, NA)


#' drop superfluous variables,
#' drop NAs for psychometric PC values,
#' rename to something more intelligible,
#' pivot long for graphing: 
long <- all %>% select(1,3:4,80:81) %>% 
  drop_na(psych_PC1) %>% 
  rename(`Psych PC1` = psych_PC1,
         `Psych PC2` = psych_PC2) %>% 
  drop_na() %>% 
  pivot_longer(-c(Key, Patient_ID, Study_ID),
               values_to = "PC", 
               names_to = "target") 

## Graph relative/control combined vs addicted
long %>%
  ggplot(aes (x = Key, y = PC)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Key)) +
  theme_pubclean() +
  font("xy.text", size = 12) +
  font("ylab", size = 18) +
  font("legend.title", size = 18) +
  font("legend.text", size = 18) +
  xlab("") + ylab("Variance\n") +
  facet_wrap(. ~ target, scales = "free_y", dir = "v") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))
ggsave("f3_psychometric_factor_ttests.pdf.pdf")

pc1 <- long %>% filter(target == "Psych PC1")
t.test(PC~Key, data=pc1, alternative="two.sided",
       var.equal=TRUE)
# t = 15.553, df = 336, p-value < 2.2e-16
p.values[1] <- 2.2e-16


pc2 <- long %>% filter(target == "Psych PC2")
t.test(PC~Key, data=pc2, alternative="two.sided",
       var.equal=TRUE)
# t = -4.5335, df = 336, p-value = 8.078e-06
p.values[2] <- 8.078e-06

# Group comparison: pcr PCs -----------------------------------------------

#' drop superfluous variables,
#' drop NAs for psychometric PC values,
#' rename to something more intelligible,
#' pivot long for graphing: 
long <- all %>% select(1,3:4,82:83) %>% 
  drop_na(pcr_PC1) %>% 
  rename(`PCR PC1` = pcr_PC1,
         `PCR PC2` = pcr_PC2) %>% 
  drop_na() %>% 
  pivot_longer(-c(Key, Patient_ID, Study_ID),
               values_to = "PC", 
               names_to = "target") 


## Graph relative/control combined vs addicted
long %>%
  ggplot(aes (x = Key, y = PC)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Key)) +
  theme_pubclean() +
  font("xy.text", size = 12) +
  font("ylab", size = 18) +
  font("legend.title", size = 18) +
  font("legend.text", size = 18) +
  xlab("") + ylab("Variance\n") +
  facet_wrap(. ~ target, scales = "free_y", dir = "v") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))
ggsave("f3_pcr_factor_ttests.pdf")

pc1 <- long %>% filter(target == "PCR PC1")
t.test(PC~Key, data=pc1, alternative="two.sided",
       var.equal=TRUE)
# t = -0.59698, df = 162, p-value = 0.5514
p.values[3] <- 0.5514

pc2 <- long %>% filter(target == "PCR PC2")
t.test(PC~Key, data=pc2, alternative="two.sided",
       var.equal=TRUE)
# t = 2.4804, df = 162, p-value = 0.01415
p.values[4] <- 0.01415


# make R stop doing scientific notation:
options(scipen=999)

# fdr adjustment on p values
p.adjust(p.values, method = "BH") %>% round(digits = 5)
# 0.00000**** 0.00002**** 0.55140 0.01887*

# Save and export ---------------------------------------------------------

# Export data w/PCs:
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
all %>% write_csv("7_2023_Addiction_clean.csv")
all %>% write_csv("7_2023_Addiction_clean.csv")