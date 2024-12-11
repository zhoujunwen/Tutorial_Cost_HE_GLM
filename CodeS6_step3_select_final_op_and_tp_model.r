
# Start #####################################################################
# Title: Estimating Costs Associated with Disease Model States Using Generalized Linear Models: A Tutorial
# Step 3: Model selection
# > B: Test for selection within one-part and two-part model respectively
# Author: Junwen Zhou
# Date: 1 August 2023

rm(list = ls())

# Set path
# setwd("C:\\XXX\\XXX") # If you don't set it, it will be the directory where you first open the r program
path_output <- getwd()

# Source packages
library(tidyverse)
library(sandwich)
library(lmtest)
select <- dplyr::select

# Specification and performance tests ----

## Prepare function ----

## > Cluster robust standard error ----

f_clx <- function(mod,cluster){
  
  ## https://rdrr.io/cran/ivpack/src/R/clx.R
  ## https://www.ne.su.se/polopoly_fs/1.216115.1426234213!/menu/standard/file/clustering1.pdf
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


## > Specification test ----

f_test_spe_hl_clx <- function(id, y_hat, res){
  
  # Required package
  # library(sandwich) # f_clx
  # library(lmtest) # f_clx
  # library(aod)  # f_clx
  
  # Required function: f_clx
  
  # - checks fit on raw scale for systematic bias 
  ana <- tibble(id = id, res = res)
  ana$y_hat_decile <- factor(dplyr::ntile(y_hat, n = 10), levels = 1:10)
  hl_mod <- lm(res ~ y_hat_decile - 1, data = ana)
  hl_test <- aod::wald.test(b = coef(hl_mod), Sigma = f_clx(hl_mod, ana$id)$vcovCL, Terms = 1:10)
  output <- tibble(hl_p = hl_test$result[[1]][3])
  return(output)
}

f_test_spe_pl_clx <- function(family, id, y, y_link){
  
  # - determines linearity of response 
  ana <- tibble(y = y,
                y_link = y_link,
                id = id)
  
  pl_mod <- glm(formula = y ~ + y_link + I(y_link^2), data = ana, 
                family = family, start = c(0,1,0))
  pl_test <- f_clx(pl_mod, ana$id)$coeftest
  
  output <- tibble(pl_p = pl_test[3,4], pl_est = pl_test[3,1])
  return(output)  
  
}

f_test_spe_mp_clx <- function(id, y_hat, res){
  
  # Modified Park test
  # - determines family 
  ana <- tibble(y_hat = y_hat, res = res, id = id)
  n_le0 <- sum(ana$y_hat <= 0)
  
  output <- tibble(mp_slope = NA)
  
  if(n_le0 == 0){
    mp_mod <- glm(I(res^2) ~ I(log(y_hat)), data = ana %>% filter(res != 0), 
                  family = Gamma(link="log"), start = rep(2, 2))
    mp_test <- f_clx(mp_mod, ana %>% filter(res !=0) %>% pull(id))
    
    output$mp_slope <- mp_test$coeftest[2]
    # Slope
    # 0: Gaussian
    # 1: Poisson
    # 2: Gamma
    # 3: Inverse Gaussian or Wald
    
  }
  return(output)  
}

# > Performance test ----

f_test_gof <- function(res){
  
  # ME: mean error
  # MAE: mean absolute error
  # RMSE: root mean squared error
  
  output <- tibble(me = round(mean(res)),
                   mae = round(mean(abs(res))), 
                   rmse = round(sqrt(mean(res^2))))
  return(output)
  
}

## Perform tests ----

name_mod <- c("gam_id", "gam_log")
n <- length(name_mod)
tmp <- rep(list(NA), n)

for(i in 1:n){
  
  mod <- readRDS(file.path(path_output, str_c("step3_mod_2p2_",name_mod[[i]], ".rds"))) 
  ana <- with(mod, tibble(id = data$id,
                          y = data$cost,
                          y_hat = fitted.values,
                          y_link = linear.predictors)) %>% 
    mutate(res = y - y_hat)
  family <- mod$family
  test_hl <- with(ana, f_test_spe_hl_clx(id, y_hat, res))
  test_pl <- with(ana, f_test_spe_pl_clx(family = family, id, y, y_link))
  test_mp <- with(ana, f_test_spe_mp_clx(id, y_hat, res))
  test_gof <- f_test_gof(ana$res)
  tmp[[i]] <- bind_cols(test_hl, test_pl, test_mp, test_gof)
}

names(tmp) <- name_mod
tp <- bind_rows(tmp, .id = "mod")
output <- tp %>% 
  mutate_at(c("hl_p", "pl_p", "mp_slope"), 
            ~ifelse(. < 0.01, "<0.01", as.character(round(., 2)))) %>% 
  select(mod,mp_slope, hl_p, pl_p, me, mae, rmse)
write.csv(output, file.path(path_output, "step3_tbl_gof_test_mod_p2.csv"),
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

## Perform tests ----

name_mod <- c("gam_id", "gam_log")
n <- length(name_mod)
tmp1 <- rep(list(NA), n)

for(i in 1:n){
  
  mod <- readRDS(file.path(path_output, str_c("step2_mod_2p2_",name_mod[[i]], ".rds"))) 
  ana <- with(mod, tibble(y = data$cost,
                          y_hat = fitted.values)) %>% 
    mutate(res = y_hat - y)
  tmp1[[i]] <- f_test_gof_by_y(res= ana$res, y = ana$y_hat, by_n = 10)
}

names(tmp1) <- name_mod

## Output plot ----

tp <- bind_rows(tmp1, .id = "mod")
p <- ggplot(data = tp %>% mutate(mod = factor(mod, level = c(name_mod))),
            aes(x = y_n, y = me, group = mod)) + 
  geom_line(aes(color = mod), linetype = "dotted", linewidth =1) +
  geom_point(aes(color = mod),size = 2) + 
  geom_hline(yintercept = 0, color = "red") + 
  scale_color_manual(name = "GLM (Distribution - Link)", 
                       values = c("gam_id" = "#F8766D",
                                  "gam_log" = "#00BFC4"),
                       labels = c("gam_id" = "Gamma - Identity",
                                  "gam_log" = "Gamma - LOG")) + 
  labs(x = "Deciles of predicted annual costs",
       y = "Mean difference of predicted vs. actual annual costs ") + 
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size=17),
        legend.text=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17,face="bold"))

ggsave(p, filename  = file.path(path_output, "step3_fig_me_by_decile_mod_p2.png") , 
       height = 200, width = 300, units = "mm")



