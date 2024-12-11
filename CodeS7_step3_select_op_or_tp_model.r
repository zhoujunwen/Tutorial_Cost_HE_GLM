
# Start #####################################################################
# Title: Estimating Costs Associated with Disease Model States Using Generalized Linear Models: A Tutorial
# Step 3: Model selection
# > C: Test for selection between one-part and two-part model
# Author: Junwen Zhou
# Date: 1 August 2023

rm(list = ls())

# Set path
# setwd("C:\\XXX\\XXX") # If you don't set it, it will be the directory where you first open the r program
path_output <- getwd()

# Source packages
library(tidyverse)
select <- dplyr::select

# Performance tests ----

## Prepare function ----

f_test_gof <- function(res){
  
  # ME: mean error
  # MAE: mean absolute error
  # RMSE: root mean squared error
  output <- tibble(me = round(mean(res)),
                   mae = round(mean(abs(res))), 
                   rmse = round(sqrt(mean(res^2))))
  return(output)
  
}

## Perform tests for overall ----

## > One-part ----

name_mod <- c("gau_id")
n <- length(name_mod)
tmp1 <- rep(list(NA), n)

for(i in 1:n){
  mod <- readRDS(file.path(path_output, str_c("step3_mod_1p_",name_mod[[i]], ".rds"))) 
  ana <- with(mod, tibble(y = data$cost,
                          y_hat = fitted.values)) %>% 
    mutate(res = y_hat - y)
  tmp1[[i]] <- f_test_gof(ana$res)
}

names(tmp1) <- str_c("op_", name_mod)

## > Two-part - GAM - ID ----

mod_data <- readRDS(file.path(path_output, "step1_ana.rds")) 
mod_p1 <- readRDS(file.path(path_output, "step3_mod_2p1_logit.rds")) 

name_mod <- c("gam_id")
n <- length(name_mod)
tmp2 <- rep(list(NA), n)

ana2_p1 <- mod_data %>% 
  select(id,  year, cost) %>% 
  bind_cols(cost_p1 = mod_p1$fitted.values)

for(i in 1:n){
  mod_p2 <- readRDS(file.path(path_output, str_c("step3_mod_2p2_",name_mod[[i]], ".rds"))) 
  ana2_p2 <- predict(mod_p2, newdata= mod_data, type = "response")
  ana2 <- bind_cols(ana2_p1, tibble(cost_p2 = ana2_p2)) %>% 
    mutate(y_hat = cost_p1 * cost_p2,
           res = y_hat - cost)
  tmp2[[i]] <- f_test_gof(ana2$res)
}

names(tmp2) <- str_c("tp_p2_", name_mod)

## > Output ----

output <- bind_rows(c(tmp1, tmp2), .id = "mod")
write.csv(output, file.path(path_output, "step3_tbl_gof_overall.csv"),
          row.names = FALSE)

# Performance by predicted results -----

## Prepare function ----

f_test_gof_by_y <- function(res, y = NULL, by_n = NULL){
  
  if(is.null(y)){
    ana <- tibble(res = res, y_n = 1) 
  } else {
    ana <- tibble(res = res, 
                  y_n = factor(dplyr::ntile(y, n = by_n), levels = 1:by_n))
  }
  output <- ana %>% 
    group_by(y_n) %>% 
    summarize(me = round(mean(res)),
              mae = round(mean(abs(res))),
              rmse = round(sqrt(mean(res^2))))
  return(output)
}

## Performance for overall ----

## > One-part ----

name_mod <- c("gau_id")
n <- length(name_mod)
tmp1 <- rep(list(NA), n)

for(i in 1:n){
  mod <- readRDS(file.path(path_output, str_c("step3_mod_1p_",name_mod[[i]], ".rds"))) 
  ana <- with(mod, tibble(y = data$cost,
                          y_hat = fitted.values)) %>% 
    mutate(res = y_hat - y)
  tmp1[[i]] <- f_test_gof_by_y(res= ana$res, y = ana$y_hat, by_n = 10)
}

names(tmp1) <- str_c("op_", name_mod)

## > Two-part - GAM - ID ----

mod_data <- readRDS(file.path(path_output, "step1_ana.rds")) 
mod_p1 <- readRDS(file.path(path_output, "step3_mod_2p1_logit.rds")) 
ana2_p1 <- mod_data %>% 
  select(id,  year, cost) %>% 
  bind_cols(cost_p1 = mod_p1$fitted.values)

name_mod <- c("gam_id")
n <- length(name_mod)
tmp2 <- rep(list(NA), n)

for(i in 1:n){
  
  mod_p2 <- readRDS(file.path(path_output, str_c("step3_mod_2p2_",name_mod[[i]], ".rds"))) 
  ana2_p2 <- predict(mod_p2, newdata= mod_data, type = "response")
  ana2 <- bind_cols(ana2_p1, tibble(cost_p2 = ana2_p2)) %>% 
    mutate(y_hat = cost_p1 * cost_p2,
           res = y_hat - cost)
  tmp2[[i]] <- f_test_gof_by_y(res= ana2$res, y = ana2$y_hat, by_n = 10)
}

names(tmp2) <- str_c("tp_p2_", name_mod)

## > Output ----

tp <- bind_rows(c(tmp1, tmp2), .id = "mod")
p <- ggplot(data = tp %>% mutate(mod = factor(mod, level = unique(tp$mod))),
            aes(x = y_n, y = me, group = mod)) + 
  geom_line(aes(color = mod), linetype = "dotted", linewidth =1) +
  geom_point(aes(color = mod),size = 2) + 
  geom_hline(yintercept = 0, color = "red") + 
  scale_color_manual(name = "GLM (Distribution - Link)", 
                     values = c("op_gau_id" = "#F8766D",
                                "tp_p2_gam_id" = "#00BFC4"),
                     labels = c("op_gau_id" ="One-part (Gaussian - Identity)",
                                "tp_p2_gam_id" = "Two-part (P1: Logistic; P2: Gamma - Identity)")) + 
  labs(x = "Deciles of predicted annual costs",
       y = "Mean difference of predicted vs. actual annual costs ") + 
  theme_bw()  +
  theme(legend.position = "bottom",
        legend.title = element_text(size=17),
        legend.text=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17,face="bold"))

ggsave(p, filename  = file.path(path_output, "step3_fig_me_by_decile_overall.png") , 
       height = 200, width = 300, units = "mm")

# Final model ----

# Prepare function 

f_clx <- function(mod,cluster){
  # library(sandwich)
  # library(lmtest)
  
  ## https://rdrr.io/cran/ivpack/src/R/clx.R
  dfcw=1
  M <- length(unique(cluster))
  N <- length(cluster)
  dfc <- (M/(M-1))*((N-1)/(N-mod$rank))
  u <- apply(estfun(mod),2,
             function(x) tapply(x, cluster, sum))
  vcovCL <- dfc*sandwich(mod, meat.=crossprod(u)/N)*dfcw
  test <- coeftest(mod, vcovCL)
  list(vcovCL = vcovCL, coeftest = test)
}

## Part 1 ----

mod <- readRDS(file = file.path(path_output, "step3_mod_2p1_logit.rds"))
saveRDS(mod, file= file.path(path_output, "step3_mod_fin_2p1.rds"))

# Estimate cluster robust standard error as we include multiple records from the same participant
tmp <- f_clx(mod, cluster = mod$data$id)$coeftest
tbl <- tibble(term = rownames(tmp),
              est = tmp[,1],
              se = tmp[,2]) %>% 
  mutate(l = est - 1.96 * se,
         h = est + 1.96 * se) %>% 
  mutate_at(c("est", "l", "h"), ~round(exp(.),2)) %>% 
  mutate(out = str_c(est, " (", l, ", ", h, ")")) %>% 
  select(term, out)

## Part 2 ----

mod <- readRDS(file = file.path(path_output, "step3_mod_2p2_gam_id.rds"))
saveRDS(mod, file= file.path(path_output, "step3_mod_fin_2p2.rds"))

# Estimate cluster robust standard error as we include multiple records from the same participant
tmp <- f_clx(mod, cluster = mod$data$id)$coeftest
tbl2 <- tibble(term = rownames(tmp),
              est = tmp[,1],
              se = tmp[,2]) %>% 
  mutate(l = est - 1.96 * se,
         h = est + 1.96 * se) %>% 
  mutate_at(c("est", "l", "h"), ~round(.,0)) %>% 
  mutate(out = str_c(est, " (", l, ", ", h, ")")) %>% 
  select(term, out)

# Output 
output <- bind_rows(tibble(mod = "mod_2p1"),
                    tbl, 
                    tibble(mod = "mod_2p2"),
                    tbl2)

write.csv(output, file= file.path(path_output, "step3_tbl_mod_fin.csv"),
          row.names = FALSE)
