
# Start #####################################################################
# Title: Estimating Costs Associated with Disease Model States Using Generalized Linear Models: A Tutorial
# Step 0: Generation of synthetic dataset
# Author: Junwen Zhou
# Date: 1 August 2023

rm(list = ls())

# Set path
# setwd("C:\\XXX\\XXX") # If you don't set it, it will be the directory where you first open the r program
path_output <- getwd()

# Source packages 
library(tidyverse)
select <- dplyr::select

# Prepare data ###############################################################

input <- list()

input$base <- list(
  age = c(56, 8.1),
  male = 0.439,
  race = list(white = 0.94, 
              black = 0.016, 
              asian = 0.016, 
              other = 0.022),
  townsend = list(q1 = 0.374, 
                  q2 = 0.201, 
                  q3 = 0.163, 
                  q4 = 0.145, 
                  q5 = 0.117),
  smoke = list(none = 0.559, 
               former = 0.332, 
               current = 0.103),
  pa = list(low = 0.148, 
            moderate = 0.329, 
            high = 0.327),
  unhealthy_diet = 0.335 / (0.335 + 0.643),
  bmi = list(underweight = 0.005,
             normal = 0.335,
             overweight = 0.423,
             obesity1 = 0.167,
             obesity2 = 0.046,
             obesity3 = 0.018),
  ldl = c(3.6, 0.8),
  hdl = c(1.5, 0.4),
  creatinine = c(71.5, 15.1),
  sbp = c(137.8, 18.6),
  dbp = c(82.4, 10.1),
  atht = 0.162,
  db = 0.049 ,
  cancer = 0.074,
  mental = 0.081
)

input$evt <- tibble(mi = 0.01,
                    stroke = 0.009,
                    vd = 0.004,
                    nvd = 0.02) %>% 
  mutate_all(~(-log( 1- .) / 7.1))

input$cost_p1 <- c(
  base = 0.13,
  male = 0.92,
  cur_age = 1.38,
  mi1 = 47.09,
  mi2 = 1.76,
  mi3 = 1.44,
  mi4 = 1.35,
  stroke1 = 47.08,
  stroke2 = 2.58,
  stroke3 = 1.78,
  stroke4 = 1.49,
  vd1 = 2.32,
  nvd1 = 11.4) %>% 
  log() # Perform natural log transformation easier for matrix calculation

input$cost_p2 <- c(
  base = 2102,
  male = -65,
  cur_age = 173,
  mi1 = 3054,
  mi2 = 670,
  mi3 = 304,
  mi4 = 304,
  stroke1 = 4485,
  stroke2 = 2192,
  stroke3 = 833,
  stroke4 = 833,
  vd1 = 4318,
  nvd1 = 4792)

# Simulate participant ----

set.seed(1234)

# > Simulate baseline ----

f_simulate_patient <- function(npat, input_base){
  
  para_n <- map_dbl(input_base, ~length(.x))
  para_type <- ifelse(para_n == 1, "Binary", ifelse(para_n == 2, "Continuous", "Categorical"))
  
  dat_con <- input_base[para_type == "Continuous"] 
  dat_bin <- input_base[para_type == "Binary"]
  dat_cat <- map(input_base[para_type == "Categorical"], 
                 ~bind_cols(.x) %>% 
                   gather(key = "cat", value = "prob") %>% 
                   mutate(prob = prob / sum(prob)))
  
  pat_con <- map(dat_con, ~rnorm(n = npat, mean = .x[[1]], sd = .x[[2]]))
  neg_con <- map_dbl(pat_con, ~sum(.x <0))
  pat_con[neg_con>0] <- map(pat_con[neg_con>0], ~ifelse(.x <0, min(.x[.x > 0]), .x))
  
  pat_bin <-  map(dat_bin, ~sample(c("1","0"), npat, replace=TRUE, prob=c(.x,1-.x)))
  pat_cat <- map(dat_cat, ~sample(.x$cat, npat, replace=TRUE, prob= .x$prob))
  
  output <- bind_cols(id = 1:npat, c(pat_con, pat_bin, pat_cat)[names(input_base)]) 
  return(output) 
}

pat_base <- f_simulate_patient(10000, input$base)

# > Simulate incident event ----

f_simulate_event <- function(npat = 10000, nyear = 10, input_evt = input$evt){
  
  evt_t <- map(input_evt, ~rexp(npat, rate = .x))
  evt_year <- map(evt_t, ~ tibble(id = 1:npat, evtyear = ifelse(.x > nyear, nyear + 1, ceiling(.x))))
  
  lastyear <- map2_dbl(evt_year$vd$evtyear, evt_year$nvd$evtyear, ~min(.x, .y, nyear))
  
  pat_year <- expand_grid(id = 1:npat, year = 1:nyear) %>% 
    left_join(tibble(id = 1:npat, yearlast = lastyear), by = "id")  %>% 
    filter(year <= yearlast) %>% 
    select(id, year)
  
  pat_evt_year <- map2(evt_year, names(evt_year),
                       ~left_join(pat_year, .x, by = "id") %>% 
                         mutate(evtyear2 = ifelse(year >= evtyear, 1, 0)) %>% 
                         group_by(id) %>% 
                         arrange(id, year) %>% 
                         mutate(evtyear3 = cumsum(evtyear2)) %>% 
                         ungroup() %>% 
                         arrange(id, year) %>% 
                         pull(evtyear3)) 
  output <- bind_cols(pat_year, pat_evt_year)
  
  return(output) 
}

pat_evt <-  f_simulate_event(npat = 10000, nyear = 10, input_evt = input$evt)


# > Simulate costs ----

f_simulate_cost <- function(sim_pat = pat_base %>% select(id, age, male),
                            sim_evt = pat_evt, 
                            input_p1 = input$cost_p1, 
                            input_p2 = input$cost_p2){
  
  mx_evt <- left_join(sim_evt, sim_pat, by = "id") %>% 
    transmute(male = male,
              cur_age = (age + year - 1 - 60) / 10,
              mi1 = mi == 1,
              mi2 = mi == 2,
              mi3 = mi == 3,
              mi4 = mi >= 4,
              stroke1 = stroke == 1,
              stroke2 = stroke == 2,
              stroke3 = stroke == 3,
              stroke4 = stroke >= 4,
              vd1 = vd == 1,
              nvd1 = nvd == 1) %>% 
    mutate_all(~as.numeric(.)) %>% 
    mutate(base = 1) %>% 
    relocate(base) %>% 
    as.matrix()
  
  pat_year_odd <- exp(mx_evt %*% input_p1) 
  pat_year_prob <- pat_year_odd / (pat_year_odd + 1)
  pat_year_cost <- mx_evt %*% input_p2
  
  npatyr <- nrow(mx_evt)
  
  pat_year_cost_p1 <- ifelse(runif(npatyr) < pat_year_prob, 1, 0)
  pat_year_cost_p2 <- rnorm(npatyr, pat_year_cost, pat_year_cost / 2)
  pat_year_cost_p2 <- ifelse(pat_year_cost_p2 < 0, 0, pat_year_cost_p2)
  output <- bind_cols(
    sim_evt %>% select(id, year),
    cost = (pat_year_cost_p1 * pat_year_cost_p2) %>% as.vector()
  )
  
  return(output) 
}

pat_cost <- f_simulate_cost(sim_pat = pat_base %>% select(id, age, male),
                            sim_evt = pat_evt, 
                            input_p1 = input$cost_p1, 
                            input_p2 = input$cost_p2)

output <- left_join(pat_base, pat_evt, by = "id") %>% 
  left_join(pat_cost, by = c("id", "year")) %>% 
  relocate(id, year) %>% 
  mutate(cur_age = age + year - 1)

saveRDS(output, file.path(path_output, "step0_ana.rds"))
