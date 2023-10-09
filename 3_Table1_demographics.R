#' this script is used to generate the numbers that go in Table 1 
#' it also cleans up the database a bit

# Set-up ------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(ggsignif)
library(ggpubr)
library(ggplot2)
library(broom)
library(sciplot)
library(chisq.posthoc.test)
library(vcd)

#' import data:
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
df <- read_csv(file = "Addiction_demograph_molecular_complete.csv")

#' drop 'users' from analysis (starting with 380 samples): 
df %<>% filter(!is.na(Key.combo)) #353 samples

#' remove 'Key'; will be doing analysis on controls+relatives (Key.combo)
df %<>% select(-Key)

#' rename Key.combo to Key:
df %<>% rename(Key = Key.combo)

#' modify Key variable to read 'CUD' instead of addicted:
df$Key[df$Key == "addicted"] <- "CUD"

#' Make 'Key' a factor:
df$Key <- factor(df$Key)

# Categorical variables - wrangle -------------------------------------------

# drop rows that have no encoding for cocaine dependence:
df %<>% filter(!is.na(Dep_Cocaine)) #349 samples
#' lose 4 samples, 2070 from study 377
#' 1053 and 1054 from study 231
#' 1012 from study 299

# select variables of interest:
df2 <- df %>% select(Key, Gender, Ethnic, Smoker, 
                     Dep_Alcohol, Dep_Opiates, Dep_Cocaine) 

# Make gender a factor and decode:
df2$Gender <- factor(df2$Gender)
levels(df2$Gender)[1:2] <- c("male", "female")

# Make ethnicity a factor, decode, combine all non-white
# ethnicities into 'other' category:
df2$Ethnic <- factor(df2$Ethnic)
levels(df2$Ethnic)[1:6] <- c("white", 
                             "non-white", 
                             "non-white", 
                             "non-white", 
                             "non-white", 
                             "non-white")

# smoker, 0 = non-smoker, 1 = smoker, 2 = former smoker
# BUT at Karen's request, simplify former smokers and non-smokers:
df2$Smoker <- factor(df2$Smoker)
levels(df2$Smoker)[1:3] <- c("non-smoker", 
                             "smoker",
                             "non-smoker")

# Categorical variables - cocaine dependence REDO -----------------------------

# create count summary of dep x group, stored in temp object M:
M <- df2 %>% select(Key, Dep_Cocaine) %>% 
  group_by(Key, Dep_Cocaine) %>% 
  summarise(n = n())
# 186 controls w/no dependence, 0 with dependence
# 163 CUD with dependence

# extract counts:
col.dep <- M$n

# there were no controls w/dependence. Add zero for table:
col.dep <- append(col.dep, 0, after = 1)

# there were no CUDs w/o dependence. Add zero for table:
col.dep <- append(col.dep, 0, after = 2)

# convert counts into a 2 row x 2 column matrix:
dep <- matrix(col.dep, ncol = 2, nrow = 2)

# name rows and columns:
dimnames(dep) = list(`cocaine status` = c("not dependent", "dependent"), 
                     group = c("control", "CUD"))

# perform chi squared test:
chisq.test(dep)
# X-squared = 344.99, df = 1, p-value < 2.2e-16

# run post-hoc pairwise comparisons using fdr correction method:
chisq.posthoc.test(dep, method = "fdr")

# Dimension     Value   control       CUD
# 1 not dependent Residuals  18.68154 -18.68154
# 2 not dependent  p values   0.00000   0.00000
# 3     dependent Residuals -18.68154  18.68154
# 4     dependent  p values   0.00000   0.00000

# visualize departure from expected results:
dep; addmargins(dep)
expect <- round(chisq.test(dep)$expected)
addmargins(expect)
round(chisq.test(dep)$residuals, 1)
mosaic(dep, shade=T, gp = shading_max)
# manually saved as Addiction_mosaic_cocainedep_v_group.pdf

# Categorical variables - gender ------------------------------------------

# create count summary of sex x group, stored in temp object M:
M <- df2 %>% select(Key, Gender) %>% 
  group_by(Key, Gender) %>% 
  summarise(n = n())

# extract counts:
col.sex <- M$n

# convert counts into a 2 row x 2 column matrix:
sex <- matrix(col.sex, ncol = 2, nrow = 2)

# name rows and columns:
dimnames(sex) = list(gender = c("male", "female"), 
                     group = c("control", "CUD"))

# perform chi squared test:
chisq.test(sex)
# X-squared = 16.19, df = 1, p-value = 5.73e-05

# run post-hoc pairwise comparisons using fdr correction method:
chisq.posthoc.test(sex, method = "fdr")

# Dimension     Value   control       CUD
# 1      male Residuals -4.166786  4.166786
# 2      male  p values  0.000062  0.000062
# 3    female Residuals  4.166786 -4.166786
# 4    female  p values  0.000062  0.000062

# visualize departure from expected results:
sex; addmargins(sex)
expect <- round(chisq.test(sex)$expected)
addmargins(expect)
round(chisq.test(sex)$residuals, 1)
mosaic(sex, shade=T, gp = shading_max)
# manually saved as Addiction_mosaic_sex_v_group.pdf

# Categorical variables - ethnicity ------------------------------------------

# create count summary of ethnicity x group, stored in temp object M:
M <- df2 %>% select(Key, Ethnic) %>% 
  group_by(Key, Ethnic) %>% 
  summarise(n = n())

#' drop NA sample for stats:
M <- drop_na(M)

# extract counts:
col.eth <- M$n

# convert counts into a 2 row x 2 column matrix:
eth <- matrix(col.eth, ncol = 2, nrow = 2)

# name rows and columns:
dimnames(eth) = list(ethnicity = c("white", "other"), 
                     group = c("control", "CUD"))

# perform chi squared test:
chisq.test(eth)
# X-squared = 18.112, df = 1, p-value = 2.083e-05

# run post-hoc pairwise comparisons using fdr correction method:
chisq.posthoc.test(eth, method = "fdr")

# Dimension     Value   control       CUD
# 1     white Residuals  4.398911 -4.398911
# 2     white  p values  0.000022  0.000022
# 3     other Residuals -4.398911  4.398911
# 4     other  p values  0.000022  0.000022

# visualize departure from expected results:
eth; addmargins(eth)
expect <- round(chisq.test(eth)$expected)
addmargins(expect)
round(chisq.test(eth)$residuals, 1)
mosaic(eth, shade=T, gp = shading_max)
# manually saved as Addiction_mosaic_ethnicity_v_group.pdf

# Categorical variables - smokers ------------------------------------------

# create count summary of Smoker x group, stored in temp object M:
M <- df2 %>% select(Key, Smoker) %>% 
  group_by(Key, Smoker) %>% 
  summarise(n = n())

# extract counts:
col.smo <- M$n

# convert counts into a 2 row x 2 column matrix:
smo <- matrix(col.smo, ncol = 2, nrow = 2)

# name rows and columns:
dimnames(smo) = list(`smoking status` = c("non-smoker", "smoker"), 
                     group = c("control", "CUD"))

# perform chi squared test:
chisq.test(smo)
# X-squared = 170.48, df = 1, p-value < 2.2e-16

# run post-hoc pairwise comparisons using fdr correction method:
chisq.posthoc.test(smo, method = "fdr")

# Dimension     Value   control       CUD
# 1 non-smoker Residuals  13.16465 -13.16465
# 2 non-smoker  p values   0.00000   0.00000
# 3     smoker Residuals -13.16465  13.16465
# 4     smoker  p values   0.00000   0.00000

# visualize departure from expected results:
smo; addmargins(smo)
expect <- round(chisq.test(smo)$expected)
addmargins(expect)
round(chisq.test(smo)$residuals, 1)
mosaic(smo, shade=T, gp = shading_max)
# manually saved as Addiction_mosaic_smoker_v_group.pdf

# Categorical variables - alcohol dependence ----------------------------------

# create count summary of dep x group, stored in temp object M:
M <- df2 %>% select(Key, Dep_Alcohol) %>% 
  group_by(Key, Dep_Alcohol) %>% 
  summarise(n = n())
# 186 controls w/no dependence, 0 with dependence
# 136 CUD w/no dependence, 27 with dependence

# extract counts:
col.dep <- M$n

# there were no controls w/dependence. Add zero for table:
col.dep <- append(col.dep, 0, after = 1)

# convert counts into a 2 row x 2 column matrix:
dep <- matrix(col.dep, ncol = 2, nrow = 2)

# name rows and columns:
dimnames(dep) = list(`alcohol status` = c("not dependent", "dependent"), 
                     group = c("control", "CUD"))

# perform chi squared test:
chisq.test(dep)
# X-squared = 31.113, df = 1, p-value = 2.434e-08

# run post-hoc pairwise comparisons using fdr correction method:
chisq.posthoc.test(dep, method = "fdr")

# Dimension     Value   control       CUD
# 1 not dependent Residuals  5.778689 -5.778689
# 2 not dependent  p values  0.000000  0.000000
# 3     dependent Residuals -5.778689  5.778689
# 4     dependent  p values  0.000000  0.000000

# visualize departure from expected results:
dep; addmargins(dep)
expect <- round(chisq.test(dep)$expected)
addmargins(expect)
round(chisq.test(dep)$residuals, 1)
mosaic(dep, shade=T, gp = shading_max)
# manually saved as Addiction_mosaic_alcoholdep_v_group.pdf

# Categorical variables - opiate dependence ----------------------------------

# create count summary of dep x group, stored in temp object M:
M <- df2 %>% select(Key, Dep_Opiates) %>% 
  group_by(Key, Dep_Opiates) %>% 
  summarise(n = n())
# 186 controls w/no dependence, 0 with dependence
# 68 CUD w/no dependence, 95 with dependence

# extract counts:
col.dep <- M$n

# there were no controls w/dependence. Add zero for table:
col.dep <- append(col.dep, 0, after = 1)

# convert counts into a 2 row x 2 column matrix:
dep <- matrix(col.dep, ncol = 2, nrow = 2)

# name rows and columns:
dimnames(dep) = list(`opiate status` = c("not dependent", "dependent"), 
                     group = c("control", "CUD"))

# perform chi squared test:
chisq.test(dep)
# X-squared = 146.02, df = 1, p-value < 2.2e-16

# run post-hoc pairwise comparisons using fdr correction method:
chisq.posthoc.test(dep, method = "fdr")

# Dimension     Value   control       CUD
# 1 not dependent Residuals  12.20451 -12.20451
# 2 not dependent  p values   0.00000   0.00000
# 3     dependent Residuals -12.20451  12.20451
# 4     dependent  p values   0.00000   0.00000

# visualize departure from expected results:
dep; addmargins(dep)
expect <- round(chisq.test(dep)$expected)
addmargins(expect)
round(chisq.test(dep)$residuals, 1)
mosaic(dep, shade=T, gp = shading_max)
# manually saved as Addiction_mosaic_opiatedep_v_group.pdf

rm(col.dep, col.eth, col.sex, col.smo, dep, eth, M, sex, smo)

# Continuous variables - descriptive  ---------------------------------------

# isolate the variables you want to run descriptive stats on,
# store names as a vector:
col.table1 <- df %>% select(Age, Yrs_Ed, NART, 
                            BDI_total, BIS11_total, 
                            AUDIT, DAST_20, Age_stimulants,
                            Stimulants_yrs, OCDUS_total, OCI_total,
                            STAI_S_total, STAI_T_total,
                            CTQabuse, SSSV_total) %>% colnames()

# write a function to extract count, average, SD, and SEM on a per group 
# basis: 
get.stats <- function(x, column){
  result <- x %>% dplyr::group_by(Key) %>% 
    summarise(count = n(),
              avg_column = mean(.data[[column]], na.rm = TRUE),
              SD_column = sd(.data[[column]], na.rm = TRUE),
              SE_column = se(.data[[column]], na.rm = TRUE)
    )
  colnames(result) = gsub("column", column, colnames(result))
  return(result)
}


# create table1 object, populating 
# it with the first two variables of interest (group names and Age); 
table1 <- get.stats(df, col.table1[1])

# create for loop to repeat the function over the remaining variables;
# note: Key and count do not need to be repeated. bind to existing
# table1 object: 
for(name in col.table1[2:length(col.table1)]){
  temp <- get.stats(df,name)
  table1 <- cbind(table1, temp %>% select(-Key, -count))
}

# save data 
setwd("~/~r_projects/2023_kigar_cud_paper/data/")
table1 %>% write_csv(file = "Addiction_Table1_contin_stats.csv")


# Continuous variables - graph --------------------------------------------

# Pivot table long for subsequent facet-wrapping of continuous data:
long <- df %>% select(Key, Patient_ID, Study_ID, Age, Yrs_Ed, 
                      AUDIT, DAST_20, Stimulants_yrs, OCDUS_total, 
                      NART, STAI_S_total, STAI_T_total,
                      BDI_total, BIS11_total, SSSV_total) %>% 
  pivot_longer(-c(Key, Study_ID, Patient_ID), 
               values_to = "value", 
               names_to = "variable")

# Make names look nicer for graph:
long$variable[long$variable == "BDI_total"] <- "BDI-II"
long$variable[long$variable == "BIS11_total"] <- "BIS11"
long$variable[long$variable == "DAST_20"] <- "DAST-20"
long$variable[long$variable == "OCDUS_total"] <- "OCDUS"
long$variable[long$variable == "SSSV_total"] <- "SSS-V"
long$variable[long$variable == "STAI_S_total"] <- "STAI-S"
long$variable[long$variable == "STAI_T_total"] <- "STAI-T"
long$variable[long$variable == "Stimulants_yrs"] <- "Years of use"
long$variable[long$variable == "Yrs_Ed"] <- "Years of education"

# set working directory to save graph:
setwd("~/~r_projects/2023_kigar_cud_paper/results/")

# Generate graph, for the moment assuming data is parametric:
long %>% drop_na() %>% 
  ggplot(aes (x = Key, y = value)) + 
  geom_violin() + geom_jitter(width = 0.2, aes(color=Key)) +
  theme_pubclean() +
  font("xy.text", size = 12) +
  font("ylab", size = 18) +
  font("legend.title", size = 18) +
  font("legend.text", size = 18) +
  stat_compare_means(method = "t.test",
                     hide.ns = T,
                     label = "p.signif",
                     label.x.npc = "middle",
                     label.y.npc = "top",
                     vjust = 1) +
  xlab("") + ylab("") +
  facet_wrap(. ~ variable, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(2))) +
  theme(strip.text = element_text(size = rel(1.5))) + 
  scale_color_manual(values=c('#00BFC4', '#F8766D'))

ggsave("Addiction_demographic_contin_stats.pdf")


# Continuous variables - paired t tests -----------------------------------

# Age:
bartlett.test(Age~Key, data = df,
              na.action = na.exclude) #p-value = 0.01425

t.test(Age ~ Key, 
       data = df,
       alternative = "two.sided", 
       var.equal = FALSE,
       na.action = na.exclude)
# t = -0.31827, df = 345.94, p-value = 0.7505

# AUDIT:
bartlett.test(AUDIT~Key, data = df,
              na.action = na.exclude) #p-value < 2.2e-16

t.test(AUDIT ~ Key, 
       data = df,
       alternative = "two.sided", 
       var.equal = FALSE,
       na.action = na.exclude)
# t = -5.1534, df = 194.1, p-value = 6.27e-07

# BDI:
bartlett.test(BDI_total~Key, data = df,
              na.action = na.exclude) #p-value < 2.2e-16

t.test(BDI_total ~ Key, 
       data = df,
       alternative = "two.sided", 
       var.equal = FALSE,
       na.action = na.exclude)
# t = -15.823, df = 202.74, p-value < 2.2e-16

# BIS11:
bartlett.test(BIS11_total~Key, data = df, 
              na.action = na.exclude) #p-value = 0.617

t.test(BIS11_total ~ Key, 
       data = df,
       alternative = "two.sided", 
       var.equal = TRUE,
       na.action = na.exclude)
# t = -14.033, df = 337, p-value < 2.2e-16

# NART:
bartlett.test(NART~Key, data = df, 
              na.action = na.exclude) #p-value = 0.1901

t.test(NART ~ Key, 
       data = df,
       alternative = "two.sided", 
       var.equal = TRUE,
       na.action = na.exclude)
# t = 8.1684, df = 329, p-value = 6.78e-15

# SSSV:
bartlett.test(SSSV_total~Key, data = df, 
              na.action = na.exclude) #p-value = 0.9061

t.test(SSSV_total ~ Key, 
       data = df,
       alternative = "two.sided", 
       var.equal = TRUE,
       na.action = na.exclude)
# t = -6.7067, df = 344, p-value = 8.191e-11

# STAI-S:
bartlett.test(STAI_S_total~Key, data = df, 
              na.action = na.exclude) #p-value = 5.241e-14

t.test(STAI_S_total ~ Key, 
       data = df,
       alternative = "two.sided", 
       var.equal = FALSE,
       na.action = na.exclude)
# t = -8.3133, df = 240.84, p-value = 6.867e-15

# STAI-T:
bartlett.test(STAI_T_total~Key, data = df, 
              na.action = na.exclude) #p-value = 1.649e-05

t.test(STAI_T_total ~ Key, 
       data = df,
       alternative = "two.sided", 
       var.equal = FALSE,
       na.action = na.exclude)
# t = -12.838, df = 283.95, p-value < 2.2e-16

# Years education:
bartlett.test(Yrs_Ed~Key, data = df, 
              na.action = na.exclude) #p-value = 2.836e-08

t.test(Yrs_Ed ~ Key, 
       data = df,
       alternative = "two.sided", 
       var.equal = FALSE,
       na.action = na.exclude)
# t = 7.9314, df = 320.76, p-value = 3.639e-14


# Save modified dataframe -------------------------------------------------
# changed 'addiction' to 'CUD' in Key
# dropped 'users'
# swapped to one Key column that does not include relative encoding
# dropped 4 samples that had no encoded variable for cocaine addiction 

# save data 
setwd("~/~r_projects/2023_kigar_cud_paper/clean/")
df %>% write_csv(file = "3_Addiction_clean.csv")
