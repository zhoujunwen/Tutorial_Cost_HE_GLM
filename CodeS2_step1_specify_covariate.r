
# Start #####################################################################
# Title: Estimating Costs Associated with Disease Model States Using Generalized Linear Models: A Tutorial
# Step 1: Preparation of dataset
# > Specify covariate
# Author: Junwen Zhou
# Date: 1 August 2023

rm(list = ls())

# Set path
# setwd("C:\\XXX\\XXX") # If you don't set it, it will be the directory where you first open the r program
path_output <- getwd()

# Source package
library(tidyverse)
select <- dplyr::select

# Import data ----

dat <- readRDS(file.path(path_output, "step0_ana.rds"))

# Specification ----

tp1 <- dat %>% 
  # Initial specification of covariates
  mutate(
    # Standardardize continuous covariate
    cur_age = (cur_age - 60) / 10,
    ldl = (ldl - 3.6) / 1,
    hdl = log(hdl),
    creatinine = ((log(creatinine) - 4.4) / 0.2),
    sbp = (sbp - 140) / 20,
    dbp = (dbp - 80) / 10,
    # Set reference level for discrete covariate
    male = factor(male, level = c("0", "1")),
    race = factor(race, level = c("white", "black", "asian", "other")),
    townsend = factor(townsend, level = str_c("q", c(3,1,2,4,5))),
    smoke = factor(smoke, level = c("none","former","current")),
    pa = factor(pa, level = c("moderate", "low", "high")),
    unhealthy_diet = factor(unhealthy_diet, level = c("0", "1")),
    bmi = factor(bmi, level = c("normal","underweight","overweight",
                                "obesity1", "obesity2","obesity3")),
    atht = factor(atht, level = c("0", "1")),
    db = factor(db, level = c("0", "1")),
    cancer = factor(cancer, level = c("0", "1")),
    mental = factor(mental, level = c("0", "1")),
    mi = factor(mi, level = as.character(0:10)),
    stroke = factor(stroke, level = as.character(0:10)),
    vd = factor(vd, level = c("0", "1")),
    nvd = factor(nvd, level = c("0", "1"))
    ) %>% 
  # Further specify some covariates
  mutate(
    # Combine categories for temporal history covarites
    mi = fct_collapse(mi, "4" = as.character(4:10)),
    stroke = fct_collapse(stroke, "4" = as.character(4:10))
    )
         
saveRDS(tp1, file = file.path(path_output, "step1_ana.rds"))

# Summary baseline characteristics ----

# Create a function to describe each covariate depending on whether it is continuous
f_des_cov <- function(cov, is_con = TRUE){
  if(is_con){
    # Describe mean and sd for continuous covariate
    tmp1 <- mean(cov)
    tmp2 <- sd(cov)
    output <- tibble(term = "Z", # There is no category for continuous covariate, set Z here.
                     value = str_c(round(tmp1, 1), " (", round(tmp2, 1), ")"))
  } else {
    # Describe total number of each category and their proportion for discrete covariate
    tmp1 <- table(cov)
    tmp2 <- prop.table(tmp1) * 100
    output <- tibble(term = names(tmp1),
                     value = str_c(tmp1, " (", round(tmp2, 1), ")"))
  }
  return(output)
}

# Name of baseline characteristics 
cov_dbl <- c("cur_age", "ldl", "hdl", "creatinine", "sbp", "dbp") # Continuous covariate
cov_cat <- c("male", "race", "townsend", "smoke", "pa", "unhealthy_diet", 
             "bmi", "atht","db", "cancer", "mental") # discrete covariate

# Describe original dataset
temp1 <- map_df(dat %>% filter(year == 1) %>% select(all_of(cov_dbl)), 
                ~f_des_cov(.x, is_con = TRUE), 
                .id = "cov")
temp2 <- map_df(dat %>% filter(year == 1) %>% select(all_of(cov_cat)),
                ~f_des_cov(.x, is_con = FALSE), 
                .id = "cov")
tmp1 <- bind_rows(temp1, temp2)

# Describe the dataset after specification
temp1 <- map_df(tp1 %>% filter(year == 1) %>% select(all_of(cov_dbl)), 
                ~f_des_cov(.x, is_con = TRUE), 
                .id = "cov")
temp2 <- map_df(tp1 %>% filter(year == 1) %>% select(all_of(cov_cat)), 
                ~f_des_cov(.x, is_con = FALSE), 
                .id = "cov")
tmp2 <- bind_rows(temp1, temp2)

# Covariate and term

# > Continuous covariate
temp1 <- tibble(cov = cov_dbl, # cov_dbl defined above
                term = "Z", 
                term2 = c("(Z - 60) / 10", 
                          "(Z - 3.6) / 1",
                          "Ln(Z)",
                          "(Ln(Z) - 4.4) / 0.2",
                          "(Z - 140) / 20",
                          "(Z - 80) / 10"))

# > Discrete covariate
temp2 <- tp1 %>% 
  select(all_of(cov_cat)) %>% 
  mutate(pa = fct_relevel(pa, "low"),
         bmi = fct_relevel(bmi, "underweight"),
         townsend = fct_relevel(townsend, "q1", "q2")) %>% 
  map_df(~tibble(term = levels(.x)), .id = "cov") %>% 
  mutate(term2 = term) # No difference in the term before and after specification except the reference level

# Order of covariate
temp3 <- c("cur_age",
           "male", "race", "townsend", "smoke", "pa", "unhealthy_diet", "bmi", 
           "ldl", "hdl", "creatinine", "sbp", "dbp",
           "atht","db", "cancer", "mental")

# Final term
tmp3 <- bind_rows(temp1, temp2) %>% 
  mutate(cov = factor(cov, levels = temp3)) %>% 
  arrange(cov)

# Final summary data
output <- left_join(tmp3, tmp1, by = c("cov", "term")) %>% 
  left_join(tmp2 %>% rename(value2 = value), by = c("cov", "term")) %>% 
  relocate(cov, term, value, term2, value2)

write.csv(output, file = file.path(path_output, "step1_tbl_cov_specification.csv"))



