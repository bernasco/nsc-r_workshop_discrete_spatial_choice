# nsc-r_workshop_discrete_spatial_choice.R
# Wim Bernasco, April 2022

# Use key shortcut       Alt+O to collapse all sections of this script
# Use key shortcut Shift+Alt+O to expand all sections of this script

# Load libraries --------------------------------------------------
library(tidyverse) # data handling (read_csv, mutate, full_join, ..)
library(broom)     # post-estimation (tidy, glance)
library(survival)  # clogit() function
library(mlogit)    # mlogit() function
library(car)       # coefficient tests (linearHypothesis)

# Set working directory -------------------------------------------
setwd("C:/Users/wimbe/KINGSTON/R/NSC-R_Workshops/Series_2/Workshop_Materials/nsc-r_workshop_discrete_spatial_choice")

# Read data ------------------------------------------------
# For details: see DataFileDescriptions.txt

# File in which each row represents a burglary offence
burglaries     <- read_csv(file = "TheHagueBurglars.csv",
                           show_col_types = FALSE)
# Peek into data interactively
#View(burglaries)

# File in which each row represents a neighborhood
neighborhoods  <- read_csv(file = "TheHagueNeighborhoods.csv",
                           show_col_types = FALSE)

# Peek into data interactively
#View(neighborhoods)

# Create data structure for discrete choice estimation --------
BXN <- 
  # Create the cross-product (Carthesion product)
  full_join(burglaries,neighborhoods, by=as.character()) %>% 
  # Create dependent variable: was neighborhood target or not
  mutate(CHOSEN = as.numeric(NHOODBUR==NHOODID),
         # Create independent variable: distance home - potential target
         # Pythagoras' Theorem, result in kilometres
         DISTANCE = sqrt((XRESID-X)^2 + (YRESID-Y)^2) / 1000,
         # Correction: Distance within the same neighborhood (Gosh)
         DISTANCE_ADJUST = if_else(DISTANCE==0,
                                   sqrt(SURFACE) / 2,
                                   DISTANCE),
         # Transform distance measure to proximity measure for convenience
         PROXIMITY =   -DISTANCE_ADJUST
)

# Peek into data interactively
#View(BXN)

# Estimate model 1 with clogit -----------------------------------
# (Bernasco & Nieuwbeerta (2005) Table 2 (p. 308)
Model1_clogit <- 
  clogit(formula = CHOSEN ~ PROPVAL + SINGFAM + RESMOBIL + 
           ETNHETERO + PROXIMITY + PROXCITY + RESUNITS + 
           strata(CASE),
         data=BXN, method="exact")

# Process output -------------------------------------------------

summary(Model1_clogit)

# Obtain model output in dataframes
estimates_Model1_clogit <- tidy(Model1_clogit)
stats_Model1_clogit     <- glance(Model1_clogit)

# Transform into odds ratio + confidence interval
estimates_Model1_clogit_ECI <- tidy(Model1_clogit) %>%
  mutate(odds_ratio = exp(estimate),
         ci95low    = exp(estimate - 2 * std.error),
         ci95high   = exp(estimate + 2 * std.error)
  ) %>%
  select(term, odds_ratio,ci95low, ci95high)

# Peek into data interactively
#View(estimates_Model1_clogit_ECI)

# Save to disk as CSV
write_csv(estimates_Model1_clogit_ECI, 
          file = "estimates_Model1_clogit_ECI.csv")

# Using tidy helps to post-process output, e.g. in graph 
tidy(Model1_clogit) %>% 
  ggplot() +
  geom_point(aes(x=estimate, y=term)) +
  geom_errorbarh(aes(xmin=estimate-(2*std.error),
                     xmax=estimate+(2*std.error),
                     y=term))

estimates_Model1_clogit_ECI %>% 
  ggplot() +
  geom_point(aes(x=odds_ratio, y=term)) +
  geom_errorbarh(aes(xmin=ci95low,
                     xmax=ci95high,
                     y=term)) 

# robust SE estimates
Model1_robust_clogit<-
  clogit(formula = CHOSEN ~ PROPVAL + SINGFAM + RESMOBIL +
           ETNHETERO + PROXIMITY + PROXCITY + RESUNITS + 
           strata(CASE),
         data=BXN, 
         method="efron", cluster=PERSONID)

# Compare standard error and robust standard error
tidy(Model1_robust_clogit)

# Estimate model 2 with clogit -----------------------------------
# Bernasco & Nieuwbeerta (2005), Table 3 (p. 309)
BXN_Extended <- 
  BXN %>%
  mutate(
    # interaction terms
    EH_NATIVE  = ETNHETERO * B_NATIVE,
    EH_FOREIGN = ETNHETERO * B_FOREIGN,
    PR_MINOR   = PROXIMITY * B_MINOR,
    PR_ADULT   = PROXIMITY * B_ADULT
)

Model2_clogit<-
  clogit(formula = CHOSEN ~ PROPVAL + SINGFAM + RESMOBIL + PROXCITY + 
           RESUNITS + PR_ADULT + PR_MINOR + EH_NATIVE + EH_FOREIGN + 
           strata(CASE),
         data=BXN_Extended)

tidy(Model2_clogit)

# test coefficient equality
linearHypothesis(Model2_clogit, "EH_NATIVE = EH_FOREIGN")
linearHypothesis(Model2_clogit, "PR_MINOR = PR_ADULT")

# Random sampling from alternatives -------------------------------

BXN_SA <-
  # select the non-chosen part of the choice set
  BXN %>% 
  filter(CHOSEN==0) %>%
  group_by(CASE) %>% 
  # random sample of 29 (out of 89-1=88), stratified by `CASE`
  slice_sample(n=29) %>%
  # and add the chosen part
  bind_rows(BXN %>% filter(CHOSEN==1))

# estimate model on subset of alternatives
Model1_SA_clogit <- 
  clogit(formula = CHOSEN ~ PROPVAL + SINGFAM + RESMOBIL + 
           ETNHETERO + PROXIMITY + PROXCITY + RESUNITS + 
           strata(CASE),
         data=BXN_SA, method="exact")
estimates_Model1_SA_clogit <- tidy(Model1_SA_clogit)

# Compare outcomes with and without sampling from alternatives
compare_models <-
   full_join(estimates_Model1_clogit, 
             estimates_Model1_SA_clogit, 
          by="term")

# Estimate Model 1 Using mlogit -----------------------------------

# sort by CASEID, NHOODID
BXN <- BXN %>%
  arrange(CASE, NHOODID)

# Create data frame for use with mlogit
BXN.mldata <- mlogit.data(data=BXN, 
                          shape="long",
                          choice = "CHOSEN", 
                          alt.var="NHOODID", 
                          chid.var="CASE",
                          id.var="PERSONID")

# estimate conditional logit model
Model1_mlogit <-
  mlogit(formula = 
           CHOSEN ~ PROPVAL + SINGFAM + RESMOBIL +
           ETNHETERO + PROXIMITY + PROXCITY + RESUNITS |0,
         data=BXN.mldata)
summary(Model1_mlogit)
tidy(Model1_mlogit)
glance(Model1_mlogit)




