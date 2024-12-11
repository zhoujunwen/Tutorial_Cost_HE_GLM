
# Start #####################################################################
# Title: Estimating Costs Associated with Disease Model States Using Generalized Linear Models: A Tutorial
# Step 4: Use of developed model
# > B: Estimate marginal effect of a disease state
# Author: Junwen Zhou
# Date: 1 August 2023

rm(list = ls())

# Set path
# setwd("C:\\XXX\\XXX") # If you don't set it, it will be the directory where you first open the r program
path_output <- getwd()

# Source packages
library(tidyverse)
select <- dplyr::select

# Marginal effect estimation ----

# > Prepare function ----

f_pred_2pcost_byevt <- function(mod_p1, mod_p2, dat, evt, lv){
  
  # Set baseline and event to target level
  if(evt == "vd"){  
    dat <- dat %>% mutate(nvd = "0")
  } else if(evt == "nvd"){ 
    dat <- dat %>% mutate(vd = "0")
  }
  
  # Set event to target level
  dat <- dat %>% mutate_at(evt,~lv) 
  
  # Part1 
  rst_p1 <- predict(mod_p1, newdata = dat, type = "response")
  
  # Part 2
  rst_p2 <- predict(mod_p2, newdata = dat, type = "response")
  
  # Final
  rst <- rst_p1 * rst_p2
  return(rst)
}

# Calculate marginal effect ----

# > Perform estimation ----

evt_list <- c("mi", "stroke","vd", "nvd")

# Import data 
mod_p1 <- readRDS(file.path(path_output, "step3_mod_fin_2p1.rds"))
mod_p2 <- readRDS(file.path(path_output, "step3_mod_fin_2p2.rds"))
mod_data <- readRDS(file = file.path(path_output, "step1_ana.rds"))

# Set levels for marginal effect calculation

# > Main event
evt_list_lv <- map(evt_list %>% set_names(), ~levels(mod_data[[.x]]))

# > Calculate marginal effect

tp1 <- map_df(
  evt_list %>% set_names(), 
  function(evt) { 
    # Predict cost for each category for each event 
    evt_lv <- evt_list_lv[[evt]] 
    rst <- map(evt_lv %>% set_names(), ~f_pred_2pcost_byevt(mod_p1, mod_p2, mod_data, evt, .x))
    # Calculate marginal effect for each event 
    output <- map_df(rst[2:length(rst)], ~tibble(me.mean = round(mean(.x - rst[[1]]),0)), .id = "lv")
    return(output) },
  .id = "evt")
# We can directly output this results if we are just interested in the mean effect

# Bootstrapping standard error ----

# > Preparation ----

evt_list <- c("mi", "stroke","vd", "nvd")
evt_list_lv <- map(evt_list %>% set_names(), ~levels(mod_data[[.x]]))

# Import data 
mod_p1 <- readRDS(file.path(path_output, "step3_mod_fin_2p1.rds"))
mod_p2 <- readRDS(file.path(path_output, "step3_mod_fin_2p2.rds"))
mod_data <- readRDS(file = file.path(path_output, "step1_ana.rds"))

x_start_p1 <- coef(mod_p1)
x_start_p2 <- coef(mod_p2)
mod_id <- unique(mod_data$id)

# > Bootstrapping ----

# Set the number of bootstrap samples
# - Set 1000 to get the same results shown in Table 5 (1000 bootstrappings are usuallly required for the non-parametric standard error)
# - However, 1000 bootstrappings took ~40 minutes to run.  
n_bt <- 10

set.seed(1234)
tmp <- list()

for(i in 1:n_bt){
  
  sample_id <- sample(mod_id, length(mod_id), replace = T)
  sample_data <- left_join(tibble(id = sample_id),
                           mod_data %>% group_by(id) %>% nest(), by = "id") %>%
    unnest(data)
  
  sample_p1 <- glm(formula = mod_p1$formula, 
                   data = sample_data %>% mutate(cost = ifelse(cost>0, 1, 0)),
                   family = binomial(link = "logit"),
                   start = x_start_p1)
  
  n_vd <- sample_data %>% filter(cost > 0 & vd == 1) %>% nrow()
  if(n_vd > 0){
    # To fit the model to the data, the covariate needs to have at least two records with different value
    sample_p2 <- glm(formula = mod_p2$formula, 
                     data = sample_data %>% filter(cost > 0),
                     family = Gamma(link = "identity"),
                     start = x_start_p2)
  } else {
    tmp_form <- as.character(mod_p2$formula)
    tmp_cov <- str_remove(tmp_form[[3]], " \\+ vd")
    new_form <- as.formula(str_c(tmp_form[[2]], " ~ ", tmp_cov))
    sample_p2 <- glm(formula = new_form, 
                     data = sample_data %>% filter(cost > 0),
                     family = Gamma(link = "identity"),
                     start = x_start_p2[!str_detect(names(x_start_p2), "^vd1\\b")])
  }
  
  # > Calculate marginal effect
  
  tmp[[i]] <- map_df(
    evt_list %>% set_names(), 
    function(evt) { 
      # Predict cost for each category for each event 
      evt_lv <- evt_list_lv[[evt]] 
      rst <- map(evt_lv %>% set_names(), ~f_pred_2pcost_byevt(sample_p1, sample_p2, sample_data, evt, .x))
      # Calculate marginal effect for each event 
      output <- map_df(rst[2:length(rst)], ~tibble(me.mean = mean(.x - rst[[1]])), .id = "lv")
      return(output) },
    .id = "evt")
  if(i %in% (n_bt/10 * (1:10))){ print(i) }
  
}

tp2 <- bind_rows(tmp) %>% 
  group_by(evt, lv) %>% 
  summarize(m = mean(me.mean), 
            se = sd(me.mean)) %>% 
  mutate(evt = factor(evt, levels = evt_list),
         m = round(m, 0),
         se = round(se, 0)) %>% 
  arrange(evt)

# > Output results
output <- left_join(tp1, tp2, by = c("evt", "lv")) %>% 
  mutate(l = me.mean - 1.96 * se,
         h = me.mean + 1.96 * se,
         out = str_c(round(me.mean,0), " (", round(l,0), ", ", round(h,0), ")")) %>% 
  bind_rows(tibble(evt = "*To replicate results shown in Table 5: Set 'n_bt <- 1000'"))
write.csv(output, file = file.path(path_output, "step4_marginal_effect_bt.csv"))

