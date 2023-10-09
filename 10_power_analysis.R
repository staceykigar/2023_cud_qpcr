

# Set up ------------------------------------------------------------------

#' load libraries
library(tidyverse)
library(magrittr)
library(pwr)


# Power calculations for F2 -----------------------------------------------
# p = 0.05


#' `DRD2 v BIS11 (controls only)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 82, r = -0.2400895, sig.level = 0.05,
           alternative = "two.sided")
# power = 0.5901909, 59% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.8, r = -0.2400895, sig.level = 0.05,
           alternative = "two.sided")
# n = 133.0283, 134 individuals
134-82 # need 52 more people



#' `DRD2 v BIS11 (CUD only)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 75, r = -0.2072227, sig.level = 0.05,
           alternative = "two.sided")
# power = 0.4342765, 43.4% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.8, r = -0.2072227, sig.level = 0.05,
           alternative = "two.sided")
# n = 179.6466, 180individuals
180-75 # need 105 more people



#' `DRD3 v BIS11 (CUD only)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 75, r = -0.2823405, sig.level = 0.05,
           alternative = "two.sided")
# power = 0.6974403, 69.7% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.8, r = -0.2823405, sig.level = 0.05,
           alternative = "two.sided")
# n = 95.32449, 96 individuals
96-75 # need 21 more people



#' `DRD3 v SSSV (overall)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 157, r = -0.1714575, sig.level = 0.05,
           alternative = "two.sided")
# power = 0.5773, 57.7% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.8, r = -0.1714575, sig.level = 0.05,
           alternative = "two.sided")
# n = 263.8545, 264 individuals
264-157 # need 107 more people



#' `COMT v SSSV (CUD only)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 75, r = 0.29, sig.level = 0.05,
           alternative = "two.sided")
# power = 0.7218275, 72.2% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.8, r = 0.29, sig.level = 0.05,
           alternative = "two.sided")
# n = 90.19199, 91 individuals
91-75 # need 16 more people



# Power calculations for FDR-adjusted SF3 (9 comparisons) -------------------------

# for 9 comparisons (as in figure 2/SF3):
.05/9 #p = 0.005555556


#' `DRD2 v BDI (overall)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 157, r = 0.153, sig.level = 0.005555556,
           alternative = "two.sided")
# power = 0.1948662, 19.5% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.8, r = 0.153, sig.level = 0.005555556,
           alternative = "two.sided")
# n = 551.9944, 552 individuals
552-157 # need 395 more people



#' `DRD2 v BIS11 (controls only)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 82, r = -0.2400895, sig.level = 0.005555556,
           alternative = "two.sided")
# power = 0.2753197, 27.5% power

# how many samples would I need to match alpha = 0.05?
pwr.r.test(power = 0.59, r = -0.2400895, sig.level = 0.005555556,
           alternative = "two.sided")
# n = 152.8937, 153 individuals
153-82 # need 71 more people



#' `DRD2 v BIS11 (CUD only)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 75, r = -0.2072227, sig.level = 0.005555556,
           alternative = "two.sided")
# power = 0.160687, 16.1% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.5, r = -0.2072227, sig.level = 0.005555556,
           alternative = "two.sided")
# n = 176.7142, 177 individuals
177-75 # need 102 more people



#' `DRD3 v BIS11 (CUD only)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 75, r = -0.2823405, sig.level = 0.005555556,
           alternative = "two.sided")
# power = 0.3787014, 37.9% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.697, r = -0.2823405, sig.level = 0.005555556,
           alternative = "two.sided")
# n = 131.1017, 132 individuals
132-75 # need 57 more people



#' `DRD3 v SSSV (overall)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 157, r = -0.1714575, sig.level = 0.005555556,
           alternative = "two.sided")
# power = 0.2662738, 26.6% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.5773, r = -0.1714575, sig.level = 0.005555556,
           alternative = "two.sided")
# n = 296.4867, 297 individuals
297-157 # need 140 more people



#' `DRD3 v SSSV (control only)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 82, r = -0.213, sig.level = 0.005555556,
           alternative = "two.sided")
# power = 0.1970772, 19.7% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.8, r = -0.213, sig.level = 0.005555556,
           alternative = "two.sided")
# n = 281.8392, 282 individuals
282-82 # need 200 more people



#' `DRD3 v SSSV (CUD only)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 75, r = -0.192, sig.level = 0.005555556,
           alternative = "two.sided")
# power = 0.1298644, 13.0% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.8, r = -0.192, sig.level = 0.005555556,
           alternative = "two.sided")
# n = 348.2816, 349 individuals
349-75 # need 274 more people



#' `COMT v SSSV (CUD only)`
#'how much power for DRD2 v controls?
pwr.r.test(n = 75, r = 0.29, sig.level = 0.005555556,
           alternative = "two.sided")
# power = 0.4060711, 40.6% power

# how many samples would I need for, say, 80% power?
pwr.r.test(power = 0.7218275, r = 0.29, sig.level = 0.005555556,
           alternative = "two.sided")
# n = 129.4209, 130 individuals
130-75 # need 55 more people


# Power calculations for FDR-adjusted SF3 - CUD only (7 comparisons) ---------

# for 9 comparisons (as in figure 2/SF3):
.05/7 #p = 0.007142857


#' `COMT v Age onset`
#'how much power do I have?
pwr.r.test(n = 75, r = -0.343, sig.level = 0.05,
           alternative = "two.sided")
# power = 0.8624256, 86.2%% power

pwr.r.test(power = 0.8624256, r = -0.343, sig.level = 0.007142857,
           alternative = "two.sided")
# n = 114.4412, 115 individuals
115 - 75 # 40 more samples

pwr.r.test(n = 75, r = -0.343, sig.level = 0.007142857,
           alternative = "two.sided")
# power = 0.6368295

# Power calculations for sf5 ----------------------------------------------


