
# Start #####################################################################
# Title: Estimating Costs Associated with Disease Model States Using Generalized Linear Models: A Tutorial
# Step 3: Model selection
# > A: Covariate selection for promising models 
# Author: Junwen Zhou
# Date: 1 August 2023

rm(list = ls())

# Set path
# setwd("C:\\XXX\\XXX") # If you don't set it, it will be the directory where you first open the r program
path_output <- getwd()

# Source packages
library(tidyverse)
library(sandwich) # f_clx
library(lmtest) # f_clx
library(aod)  # f_clx
library (miceadds) # LMtest
select <- dplyr::select

# Function ----

# > Cluster robust standard error ----

f_clx <- function(mod,cluster){
  # library(sandwich)
  # library(lmtest)
  
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

# > Stepwise backward selection in GLM ----

form <- function(y, x){ as.formula(str_c(y, "~", str_c(x, collapse = "+")))}

f_glm_select_cov <- function(mod, cluster, pval.in, pval.out, x.fix = NULL, opt_detail = FALSE){
  # 
  # mod = mod
  # cluster = "id"
  # pval.in = 0.05
  # pval.out = 0.05
  # x.fix = var_x_fix
  # opt_detail = FALSE
  
  mod_data <- mod$data
  mod_family <- mod$family
  mod_y <- names(mod$model[1])
  
  # Generate terms per covariate
  mod_x <- attr(mod$terms , "term.labels")
  mod_x.fct_lv <- map_df(mod$xlevels, ~tibble(lv = .x), .id = "x")
  mod_x_lv <- left_join(tibble(x = mod_x), mod_x.fct_lv, by = "x")
  mod_x_term <- mod_x_lv %>% mutate(term = if_else(is.na(lv), x, str_c(x, lv)))

  # Link each term to the term used for coef
  mod_coef <- coef(mod)
  mod_term <- names(mod_coef)
  mod_x.term <- map(mod_x %>% set_names(), ~mod_x_term %>% filter(x == .x  & term %in% mod_term) %>% pull(term))
  mod_x.coef <- tibble(term = mod_term, coef = mod_coef) %>% left_join(mod_x_term, by = "term") 
  
  # Initialize parameter for looping
  temp_mod <- mod
  test_rst_long <- tibble()
  test_rst_short <- tibble()
  
  x.droplist <- c()
  x.drop <- "init"
  x.add <- NA
  
  i <- 0
  
  # Start the loop 
  while(!is.na(x.drop)|!is.na(x.add)){
    
    i <- i + 1
    
    tmp_short <- tibble(step = i, var_drop = NA, p_drop = NA, var_add = NA, p_add = NA)
    
    # Drop covariate
    # - First update model info
    x <- attr(temp_mod$terms , "term.labels")
    coef <- coef(temp_mod)
    term <- names(coef)
    
    # - If cluster: take cluster SE
    if(is.null(cluster)) { 
      sigma <- vcov(temp_mod)
    } else {
      sigma <- f_clx(temp_mod, mod_data[[cluster]])$vcovCL
    }
    
    # - Test covariate one by one (except for newly added and fixed)
    x.testlist <- setdiff(x, c(x.fix, x.add)) 
    test <- map_df(
      x.testlist %>% set_names(), 
      function(x){
        if(sum(is.na(coef)) == 0){
          tibble(p = aod::wald.test(sigma, coef, which(term %in% mod_x.term[[x]]))$result$chi2[3])
        } else {
          new_coef <- coef[!is.na(coef)]
          new_term <- names(new_coef)
          na_term <- setdiff(term, new_term)
          tibble(p = aod::wald.test(sigma, new_coef, 
                                    which(new_term %in% setdiff(mod_x.term[[x]], na_term)))$result$chi2[3])
        }
      },
      .id = "x") %>% 
      arrange(desc(p))
    tmp_rst <- bind_cols(step = i, dir = "back", test, rst = "in")
    tmp_drop <- test %>% filter(p>= pval.out)
    x.drop <- tmp_drop$x[1]
    
    # If any need to drop
    if(!is.na(x.drop)) {
      
      # Print the dropped covariate
      x.drop_p <- tmp_drop$p[1]
      
      # x.drop_p <- test$p[1]
      print(paste0("drop:", x.drop, " (P=", round(x.drop_p, 2),")")) # 2 digits
      
      # Update selection results
      tmp_short$var_drop <- x.drop
      tmp_short$p_drop <- x.drop_p
      tmp_rst$rst[tmp_rst$x == x.drop] <- "out"
      
      # Update list of dropped covariates
      x.droplist <- c(x.droplist, x.drop)
      
      # Update temporal model
      x <- setdiff(x, x.drop)
      x_startval <- c(coef[1], coef[names(coef) %in% unlist(mod_x.term[x])])
      temp_mod <- glm(form(mod_y, x), data = mod_data, family = mod_family, start = x_startval)
    }
    
    # Update long selection results
    test_rst_long <- bind_rows(test_rst_long, tmp_rst)
    
    # Add back previously dropped covariate
    x.testlist <- setdiff(x.droplist, x.drop)
    
    if(!is.null(x.testlist)){
      # Test the added covariate one by one
      test <- tibble(x = x.testlist, p = NA)
      test.coef <- coef(temp_mod)
      for(x.test in x.testlist){ 
        # use for loop to avoid generating large list
        
        # Get start value for the new covariate from base case value
        x_startval_add <- mod_x.coef %>% filter(x == x.test) %>% pull(coef)
        x_startval <- c(test.coef[1], x_startval_add, test.coef[-1])
        temp_mod.add <- glm(
            form(y = mod_y, x = c(x.test, x)), 
            data = mod_data, family = mod_family, 
            start = x_startval)
        
        coef <- coef(temp_mod.add)
        term <- names(coef)
        
        # - If cluster: take cluster SE
        if(is.null(cluster)) { 
          sigma <- vcov(temp_mod.add)
        } else {
          sigma <- f_clx(temp_mod.add, mod_data[[cluster]])$vcovCL
        }
        
        if(sum(is.na(coef)) == 0){
          test$p[test$x == x.test] <- aod::wald.test(
            sigma, coef, Terms = which(term %in% mod_x.term[[x.test]]))$result$chi2[3]
        } else {
          new_coef <- coef[!is.na(coef)]
          new_term <- names(new_coef)
          na_term <- setdiff(term, new_term)
          test$p[test$x == x.test] <- aod::wald.test(
            sigma, new_coef, 
            which(new_term %in% setdiff(mod_x.term[[x]], na_term)))$result$chi2[3]
        }
        
      }
      test <- test %>% arrange(p)
      tmp_rst <- bind_cols(step = i, dir = "forward", test, rst = "out")
      tmp_add <- test %>% filter(p < pval.in)
      x.add <- tmp_add$x[1]
    }
    
    # if anything added back
    if(!is.na(x.add)){
      
      # Print add results
      x.add_p <- tmp_add$p[1]
      print(paste0("add:", x.add, " (P = ", round(x.add_p, 2),")")) # 2 digits
      
      # Update selection results
      tmp_short$var_add <- x.add
      tmp_short$p_add <- x.add_p
      
      tmp_rst$rst[tmp_rst$x == x.add] <- "in"
      
      # Update droplist
      x.droplist <- setdiff(x.droplist, x.add)
      
      # Get previous coef for starting value
      coef <- coef(temp_mod)
      x_startval_add <- mod_x.coef %>% filter(x == x.add) %>% pull(coef)
      
      # Update model
      x_startval <- c(coef[1], x_startval_add, coef[-1])
      temp_mod <- glm(form(y = mod_y, x = c(x.add, x)), 
                      data = mod_data, family = mod_family, 
                      start = x_startval)
    }
    
    if(!is.null(x.testlist)){ test_rst_long <- bind_rows(test_rst_long, tmp_rst) }
    test_rst_short <- bind_rows(test_rst_short, tmp_short)
    
  }
  
  if(opt_detail == TRUE){
    output <- list(long = test_rst_long, short = test_rst_short)
  } else {
    output <- test_rst_short
  }
  return(output)
}

f_glm_select_cov_bt <- function(mod, cluster, pval.in, pval.out, x.fix = NULL, n_bt = 10){
  
  # mod = mod
  # cluster = "id"
  # pval.in = 0.05
  # pval.out = 0.05
  # x.fix = NULL
  # n_bt = 10
  
  # Model info
  mod_family <- mod$family
  mod_y <- names(mod$model[1])
  mod_x <- attr(mod$terms , "term.labels")
  mod_data0 <- mod$data
  
  # Link term to the term used for coef
  
  # > Generate terms per covariate
  mod_x.fct_lv <- map_df(mod$xlevels, ~tibble(lv = .x), .id = "x")
  mod_x_lv <- left_join(tibble(x = mod_x), mod_x.fct_lv, by = "x")
  mod_x_term <- mod_x_lv %>% mutate(term = if_else(is.na(lv), x, str_c(x, lv)))
  
  # > Link each term to the term used for coef
  mod_coef <- coef(mod)
  mod_term <- names(mod_coef)
  mod_x.term <- map(mod_x %>% set_names(), ~mod_x_term %>% filter(x == .x  & term %in% mod_term) %>% pull(term))
  mod_x.coef <- tibble(term = mod_term, coef = mod_coef) %>% left_join(mod_x_term, by = "term") 
  
  # Create model data nested by id (for later sampling)
  if(!is.null(cluster)){
    colnames(mod_data0)[which(colnames(mod_data0)==cluster)] <- "id0"
  } else {
    mod_data0 <- mod_data0 %>% mutate(id0 = row_number())
  }
  id0_uni <- unique(mod_data0$id0)
  mod_data0_nest <- mod_data0 %>% group_by(id0) %>% nest()
  
  # Bootstrapping stepwise backward selection

  test_rst_bt <- vector(mode = "list", length = n_bt)
  for(j in 1:n_bt){
    
    if(j %in% ((1:10) * (n_bt/10))){ print(str_c("Completed bootstrapping: ", j-1, " out of ", n_bt)) }
    
    id_sample <- sample(x = id0_uni, size = length(id0_uni), replace = TRUE)
    mod_data <- tibble(id0 = id_sample) %>% 
      mutate(id = row_number()) %>% 
      left_join(mod_data0_nest, by = "id0") %>% 
      unnest(data) %>% 
      arrange(id, year) %>%
      select(-id0)
    
    tmp_val_n <- map_dbl(mod_data %>% select(-id), ~length(unique(.x)))
    if(sum(tmp_val_n == 1) >0){
      next
    }
    
    # Initialize parameter for looping
    temp_mod <- glm(form(mod_y, mod_x), data = mod_data, family = mod_family, start = mod_coef)
    
    test_rst_long <- tibble()
    test_rst_short <- tibble()
    
    x.droplist <- c()
    x.drop <- "init"
    x.add <- NA
    
    i <- 0
    # Start the loop 
    while(!is.na(x.drop)|!is.na(x.add)){
      
      i <- i + 1
      
      tmp_short <- tibble(step = i, var_drop = NA, p_drop = NA, var_add = NA, p_add = NA)
      
      # Drop covariate
      # - First update model info
      x <- attr(temp_mod$terms , "term.labels")
      coef <- coef(temp_mod)
      term <- names(coef)
      
      # - If cluster: take cluster SE
      if(is.null(cluster)) { 
        sigma <- vcov(temp_mod)
      } else {
        sigma <- f_clx(temp_mod, mod_data$id)$vcovCL
      }
      
      # - Test covariate one by one (except for newly added and fixed)
      x.testlist <- setdiff(x, c(x.fix, x.add)) 
      test <- map_df(
        x.testlist %>% set_names(), 
        function(x){
          if(sum(is.na(coef)) == 0){
            tibble(p = aod::wald.test(sigma, coef, which(term %in% mod_x.term[[x]]))$result$chi2[3])
          } else {
            new_coef <- coef[!is.na(coef)]
            new_term <- names(new_coef)
            na_term <- setdiff(term, new_term)
            tibble(p = aod::wald.test(sigma, new_coef, 
                                      which(new_term %in% setdiff(mod_x.term[[x]], na_term)))$result$chi2[3])
          }
        },
        .id = "x") %>% 
        arrange(desc(p))
      tmp_rst <- bind_cols(step = i, dir = "back", test, rst = "in")
      tmp_drop <- test %>% filter(p>= pval.out)
      x.drop <- tmp_drop$x[1]
      
      # If any need to drop
      if(!is.na(x.drop)) {
        
        # Print the dropped covariate
        x.drop_p <- tmp_drop$p[1]
        
        # x.drop_p <- test$p[1]
        # print(paste0("drop:", x.drop, " (P=", round(x.drop_p, 2),")")) # 2 digits
        
        # Update selection results
        tmp_short$var_drop <- x.drop
        tmp_short$p_drop <- x.drop_p
        tmp_rst$rst[tmp_rst$x == x.drop] <- "out"
        
        # Update list of dropped covariates
        x.droplist <- c(x.droplist, x.drop)
        
        # Update temporal model
        x <- setdiff(x, x.drop)
        x_startval <- c(coef[1], coef[names(coef) %in% unlist(mod_x.term[x])])
        temp_mod <- glm(form(mod_y, x), data = mod_data, family = mod_family, start = x_startval)
      }
      
      # Update long selection results
      test_rst_long <- bind_rows(test_rst_long, tmp_rst)
      
      # Add back previously dropped covariate
      x.testlist <- setdiff(x.droplist, x.drop)
      
      if(!is.null(x.testlist)){
        # Test the added covariate one by one
        test <- tibble(x = x.testlist, p = NA)
        test.coef <- coef(temp_mod)
        for(x.test in x.testlist){ 
          # use for loop to avoid generating large list
          
          # Get start value for the new covariate from base case value
          x_startval_add <- mod_x.coef %>% filter(x == x.test) %>% pull(coef)
          x_startval <- c(test.coef[1], x_startval_add, test.coef[-1])
          temp_mod.add <- glm(
            form(y = mod_y, x = c(x.test, x)), 
            data = mod_data, family = mod_family, 
            start = x_startval)
          
          coef <- coef(temp_mod.add)
          term <- names(coef)
          
          # - If cluster: take cluster SE
          if(is.null(cluster)) { 
            sigma <- vcov(temp_mod.add)
          } else {
            sigma <- f_clx(temp_mod.add, mod_data$id)$vcovCL
          }
          
          if(sum(is.na(coef)) == 0){
            test$p[test$x == x.test] <- aod::wald.test(
              sigma, coef, Terms = which(term %in% mod_x.term[[x.test]]))$result$chi2[3]
          } else {
            new_coef <- coef[!is.na(coef)]
            new_term <- names(new_coef)
            na_term <- setdiff(term, new_term)
            test$p[test$x == x.test] <- aod::wald.test(
              sigma, new_coef, 
              which(new_term %in% setdiff(mod_x.term[[x]], na_term)))$result$chi2[3]
          }
          
        }
        test <- test %>% arrange(p)
        tmp_rst <- bind_cols(step = i, dir = "forward", test, rst = "out")
        tmp_add <- test %>% filter(p < pval.in)
        x.add <- tmp_add$x[1]
      }
      
      # if anything added back
      if(!is.na(x.add)){
        
        # Print add results
        x.add_p <- tmp_add$p[1]
        # print(paste0("add:", x.add, " (P = ", round(x.add_p, 2),")")) # 2 digits
        
        # Update selection results
        tmp_short$var_add <- x.add
        tmp_short$p_add <- x.add_p
        
        tmp_rst$rst[tmp_rst$x == x.add] <- "in"
        
        # Update droplist
        x.droplist <- setdiff(x.droplist, x.add)
        
        # Get previous coef for starting value
        coef <- coef(temp_mod)
        x_startval_add <- mod_x.coef %>% filter(x == x.add) %>% pull(coef)
        
        # Update model
        x_startval <- c(coef[1], x_startval_add, coef[-1])
        temp_mod <- glm(form(y = mod_y, x = c(x.add, x)), 
                        data = mod_data, family = mod_family, 
                        start = x_startval)
      }
      
      if(!is.null(x.testlist)){ test_rst_long <- bind_rows(test_rst_long, tmp_rst) }
      test_rst_short <- bind_rows(test_rst_short, tmp_short)
    }
    
    # Results of the j bootstrapping
    test_rst_bt[[j]] <- mod_x %in% x.droplist
    
    # To inspect the results
    
    # > Long: showing P value of all the covariates at each step
    # test_rst_long 
    
    # > Short: showing the dropped or added covariate at each step and their p value
    # test_rst_short 
    
  }
  
  output <- map_df(test_rst_bt, function(z) {
    if(is.null(z)) {
      tibble(x = mod_x, include = NA)
    } else {
      tibble(x = mod_x, include = as.numeric(!z))
    }
  },
  .id = "test") 
  
  return(output)
}

# Preparation ----

# Create a list to store the covariate selection results
name_mod <- c("mod_2p1_logit", "mod_2p2_gam_id", "mod_2p2_gam_log", "mod_1p_gau_id") 
tp <- map(name_mod %>% set_names, ~list())
  
# Selection for two-part model ----

# > Part 1 - Logistic ----

name_mod_i <- "mod_2p1_logit"

# Model for covariate selection
mod <- readRDS(file = file.path(path_output, str_c("step2_", name_mod_i, ".rds")))

# Perform selection 
var_id <- "id" 
rst  <- f_glm_select_cov(mod, var_id, 0.05, 0.05)

# Selected covariate
var_x <- setdiff(attr(mod$terms , "term.labels"), unique(rst$var_drop))
  # This is a simple way to output selected covariate since no previously 
  # dropped covariate added back; 
  # We should do it differently when there is any dropped covariate added back

# Store selection results 
tp[[name_mod_i]] <- bind_rows(tibble(mod = name_mod_i, 
                                     var_fin = str_c(var_x, collapse = ",")),
                              rst)

# Update model with selected covariate 
var_y <- names(mod$model[1])
mod2 <- glm(form(var_y, var_x), 
            data = mod$data, 
            family = binomial("logit"))
saveRDS(mod2, file = file.path(path_output, str_c("step3_", name_mod_i, ".rds")))

# > Part 2 - Gamma-id ----

name_mod_i <- "mod_2p2_gam_id"

# Model for covariate selection
mod <- readRDS(file = file.path(path_output, str_c("step2_", name_mod_i, ".rds")))

# Perform selection
var_id <- "id"
rst  <- f_glm_select_cov(mod, var_id, 0.05, 0.05)

# Selected covariate
var_x <- setdiff(attr(mod$terms , "term.labels"), unique(rst$var_drop)) 

# Store selection results 
tp[[name_mod_i]] <- bind_rows(tibble(mod = name_mod_i, 
                                     var_fin = str_c(var_x, collapse = ",")),
                              rst)

# Update model with selected covariate 
var_y <- names(mod$model[1])

# Gamma identity link usually needs a good starting value so that it can work
tmp <- lm(form(var_y, var_x), data = mod$data)
tmp2 <- mod$coefficients[names(mod$coefficients) %in% names(coef(tmp))]
mod2 <- glm(form(var_y, var_x), 
            data = mod$data, 
            family = Gamma(link = "identity"), 
            start = tmp2) 
saveRDS(mod2, file = file.path(path_output, str_c("step3_", name_mod_i, ".rds")))

# > Part 2 - Gamma-LOG ----

name_mod_i <- "mod_2p2_gam_log"

# Model for covariate selection
mod <- readRDS(file = file.path(path_output, str_c("step2_", name_mod_i, ".rds")))

# Perform selection
var_id <- "id"
rst  <- f_glm_select_cov(mod, var_id, 0.05, 0.05)

# Selected covariate
var_x <- setdiff(attr(mod$terms , "term.labels"), unique(rst$var_drop))

# Store selection results 
tp[[name_mod_i]] <- bind_rows(tibble(mod = name_mod_i, 
                                     var_fin = str_c(var_x, collapse = ",")),
                              rst)

# Update model with selected covariate 
var_y <- names(mod$model[1])
mod2 <- glm(form(var_y, var_x), 
            data = mod$data, 
            family = Gamma(link = "log")) 
saveRDS(mod2, file = file.path(path_output, str_c("step3_", name_mod_i, ".rds")))

# Selection for one-part model ----

# > Gaussian - Identity ----

name_mod_i <- "mod_1p_gau_id"

# Model for covariate selection
mod <- readRDS(file = file.path(path_output, str_c("step2_", name_mod_i, ".rds")))

# Perform covariate selection
var_id <- "id"
rst  <- f_glm_select_cov(mod, var_id, 0.05, 0.05)

# Selected covariate
var_x <- setdiff(attr(mod$terms , "term.labels"), unique(rst$var_drop))

# Store selection results 
tp[[name_mod_i]] <- bind_rows(tibble(mod = name_mod_i, 
                                     var_fin = str_c(var_x, collapse = ",")),
                              rst)

# Update model with selected covariate 
var_y <- names(mod$model[1])
mod2 <- glm(form(var_y, var_x), 
            data = mod$data, 
            family = gaussian(link = "identity")) 
saveRDS(mod2, file = file.path(path_output, str_c("step3_", name_mod_i, ".rds")))

# Output covariate selection results ----

output <- map(tp, ~ .x %>% 
                filter(!(is.na(var_drop) & is.na(var_add) & is.na(var_fin))) %>%
                mutate(p_drop_chr = as.character(round(p_drop, 2)),
                       p_add_chr = as.character(round(p_add, 2)),
                       var_add = ifelse(is.na(var_add), "None", var_add)) %>%
                select(mod, var_fin, step, var_drop, p_drop_chr, var_add, p_add_chr)
              ) %>% 
  bind_rows()

write.csv(output, 
          file = file.path(path_output, "step3_tbl_cov_sel.csv"), 
          row.names = FALSE,
          na = "")

# Bootstrapping stepwise backward ----

# Set the number of bootstrap samples
# - Set 100 to get the same results shown in the supplement table 2
# - However, 100 bootstrapping takes ~1 hour to run
n_bt <- 10

# Part2 Gamma id
name_mod_i <- "mod_2p2_gam_id"

# Model for covariate selection
mod <- readRDS(file = file.path(path_output, str_c("step2_", name_mod_i, ".rds")))
  
# Perform selection
var_id <- "id"
set.seed(1234)
rst  <- f_glm_select_cov_bt(mod, var_id, 0.05, 0.05, x.fix = NULL, n_bt = n_bt)
  
output <- rst %>% 
  group_by(x) %>% 
  summarize(include = mean(include, na.rm = T)) %>% 
  arrange(desc(include)) %>% 
  bind_rows(tibble(x = "*To replicate results shown in Supplementary Table 2: Set 'n_bt <- 100'"))

write.csv(output, 
          file = file.path(path_output, "step3_tbl_cov_sel_bt.csv"), 
          row.names = FALSE,
          na = "")



