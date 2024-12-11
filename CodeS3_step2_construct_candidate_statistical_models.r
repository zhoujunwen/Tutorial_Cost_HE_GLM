
# Start #####################################################################
# Title: Estimating Costs Associated with Disease Model States Using Generalized Linear Models: A Tutorial
# Step 2: Candidate statistical model
# > A: Construct candidate statistical models with initial set of covariate
# Author: Junwen Zhou
# Date: 1 August 2023

rm(list = ls())

# Set path
# setwd("C:\\XXX\\XXX") # If you don't set it, it will be the directory where you first open the r program
path_output <- getwd()

# Source packages
library(tidyverse)
select <- dplyr::select

# Import data ----

dat <- readRDS(file.path(path_output, "step1_ana.rds"))

# Two-part model ----

# Model two parts 
# > Part 1: probability of any costs incurring in the period
# > Part 2: costs conditional on any costs incurring in the period

## Part 1 ----

# Probability of any costs incurring in the period

# > Logistic regression ----

# Convert cost outcome to 1 or 0 for the part 1 model
ana <- dat %>% mutate(cost = ifelse(cost > 0, 1, 0))

# Define the formula: outcome ~ covariate
var_y <- "cost"
var_x <- c("male", "race", "townsend", "smoke", "pa", "unhealthy_diet", 
           "bmi", "ldl", "hdl", "creatinine", "sbp", "dbp", "atht", "db", "cancer", "mental", 
           "cur_age",
           "mi", "stroke", "vd", "nvd")
form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))

# Fit the logistic regression model to the data 
mod <- glm(data = ana, 
           formula = form, 
           family = binomial(link = "logit")) 

# Save the initial model for later use
saveRDS(mod, file = file.path(path_output, "step2_mod_2p1_logit.rds"))

## Part 2 ----

# Costs conditional on any costs incurring in the period

# Select the records with positive cost outcome for the part 2 model
ana <- dat %>% filter(cost > 0)

# Define the formula: outcome ~ covariate
var_y <- "cost"
var_x <- c("male", "race", "townsend", "smoke", "pa", "unhealthy_diet", 
           "bmi", "ldl", "hdl", "creatinine", "sbp", "dbp", "atht", "db", "cancer", "mental", 
           "cur_age",
           "mi", "stroke", "vd", "nvd")
form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))

# Define the candidate GLMs
list_test <- list(gau_id = gaussian("identity"),
                  gau_log = gaussian("log"),
                  poi_id = poisson("identity"),
                  poi_log = poisson("log"),
                  gam_id = Gamma("identity"),
                  gam_log = Gamma("log"))

# > Gaussian - Identity ----

name_test <- "gau_id"
mod <- glm(data = ana, 
           formula = form, 
           family = list_test[[name_test]])
saveRDS(mod, file = file.path(path_output, str_c("step2_mod_2p2_", name_test, ".rds")))

lm_coef <- coef(mod) # Keep the coefficient from linear regression for later use

# > Gaussian - LOG ----

name_test <- "gau_log"
mod <- glm(data = ana, 
           formula = form, 
           family = list_test[[name_test]])
saveRDS(mod, file = file.path(path_output, str_c("step2_mod_2p2_", name_test, ".rds")))

# > Poisson - Identity (LM start) ----

name_test <- "poi_id"
mod <- glm(data = ana %>% mutate(cost = round(cost, 0)), 
           # Poisson regression usually require outcome to be rounded to 1
           formula = form, 
           family = list_test[[name_test]],
           start = lm_coef
           # Sometimes poisson - identity GLM will need a good initial coefficients 
           # as identity link is not the default link for poisson regression
           )
saveRDS(mod, file = file.path(path_output, str_c("step2_mod_2p2_", name_test, ".rds")))

# > Poisson - Log ----

name_test <- "poi_log"
mod <- glm(data = ana %>% mutate(cost = round(cost, 0)), 
           # Poisson regression usually require outcome to be rounded to 1
           formula = form, 
           family = list_test[[name_test]])
saveRDS(mod, file = file.path(path_output, str_c("step2_mod_2p2_", name_test, ".rds")))

# > Gamma - Identity (LM start) ----

name_test <- "gam_id"
mod <- glm(data = ana, 
           formula = form, 
           family = list_test[[name_test]],
           start = lm_coef
           # Sometimes Gamma - identity GLM will need a good initial coefficients 
           # as identity link is not the default link for Gamma regression
           )
saveRDS(mod, file = file.path(path_output, str_c("step2_mod_2p2_", name_test, ".rds")))

# > Gamma - LOG ----

name_test <- "gam_log"
mod <- glm(data = ana, 
           formula = form, 
           family = list_test[[name_test]])
saveRDS(mod, file = file.path(path_output, str_c("step2_mod_2p2_", name_test, ".rds")))

# One-part model  ----

# Model the costs incurring in the period directly 

## linear regression ----

ana <- dat
var_y <- "cost"
var_x <- c("male", "race", "townsend", "smoke", "pa", "unhealthy_diet", 
           "bmi", "ldl", "hdl", "creatinine", "sbp", "dbp", "atht", "db", "cancer", "mental", 
           "cur_age",
           "mi", "stroke", "vd", "nvd")
form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))

name_test <- "gau_id"
mod <- glm(data = ana, 
           formula = form, 
           family = gaussian("identity"))
saveRDS(mod, file = file.path(path_output, str_c("step2_mod_1p_", name_test, ".rds")))

