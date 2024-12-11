
# Start #####################################################################
# Title: Estimating Costs Associated with Disease Model States Using Generalized Linear Models: A Tutorial
# Step 4: Use of developed model
# > A: Predict cost for individual
# Author: Junwen Zhou
# Date: 1 August 2023

rm(list = ls())

# Set path
# setwd("C:\\XXX\\XXX") # If you don't set it, it will be the directory where you first open the r program
path_output <- getwd()

# Source packages
library(tidyverse)
select <- dplyr::select

# Individual profiles ----

dat <- tibble(age = 50,
              male = 0,
              sbp = 120,
              db = 1,
              mi = 1,
              stroke = 2,
              vd = 0,
              nvd = 0) 
# for the disease state descriptor (e.g. MI)
# > 1: same year of event
# > 2: one year after event
# > 3: two years after event 
# > 4 to more: same pattern as above

# Model coefficients ----

mod_p1 <- readRDS(file.path(path_output, "step3_mod_fin_2p1.rds"))
mod_p2 <- readRDS(file.path(path_output, "step3_mod_fin_2p2.rds"))
coef_p1 <- coef(mod_p1)
coef_p2 <- coef(mod_p2)

# Perform prediction ----

# > Specified individual profiles ----

ana <- dat %>% 
  transmute("(Intercept)" = 1,
            male1 = ifelse(male == 1, 1, 0),
            sbp = (sbp - 140) / 20,
            db1 = ifelse(db == 1, 1, 0),
            cur_age = (age - 60) / 10,
            mi1 = ifelse(mi == 1, 1, 0),
            mi2 = ifelse(mi == 2, 1, 0),
            mi3 = ifelse(mi == 3, 1, 0),
            mi4 = ifelse(mi >= 4, 1, 0),
            stroke1 = ifelse(stroke == 1, 1, 0),
            stroke2 = ifelse(stroke == 2, 1, 0),
            stroke3 = ifelse(stroke == 3, 1, 0),
            stroke4 = ifelse(stroke >= 4, 1, 0),
            vd1 = ifelse(vd == 1, 1, 0),
            nvd1 = ifelse(nvd == 1, 1, 0)) %>% 
  as.matrix()

# > Predict part 1 ----

rst_p1_odd <- exp(coef_p1 %*% ana[,names(coef_p1)])
rst_p1_prob <- rst_p1_odd / (rst_p1_odd + 1) 

# > Predict part 2 ----

rst_p2 <- coef_p2 %*% ana[, names(coef_p2)]

# > Final predicted costs ----

rst <- rst_p1_prob * rst_p2

print(round(rst_p1_prob, 2))
print(round(rst_p2, 0))
print(round(rst, 0))

# Comparison of crude costs and costs from model ---- 

# > Crude costs ----

dat <- readRDS("step1_ana.rds")
rst1 <- dat %>% 
  group_by(mi) %>% 
  summarize(cost = mean(cost))

# > Predicted costs ----

ana <- dat %>% 
  transmute("(Intercept)" = 1,
            male1 = ifelse(male == 1, 1, 0),
            sbp = sbp, # already transformed
            db1 = ifelse(db == 1, 1, 0),
            cur_age = cur_age,
            mi1 = ifelse(mi == "1", 1, 0),
            mi2 = ifelse(mi == "2", 1, 0),
            mi3 = ifelse(mi == "3", 1, 0),
            mi4 = ifelse(mi == "4", 1, 0),
            stroke1 = ifelse(stroke == "1", 1, 0),
            stroke2 = ifelse(stroke == "2", 1, 0),
            stroke3 = ifelse(stroke == "3", 1, 0),
            stroke4 = ifelse(stroke == "4", 1, 0),
            vd1 = ifelse(vd == 1, 1, 0),
            nvd1 = ifelse(nvd == 1, 1, 0)) %>% 
  as.matrix()

rst_p1_odd <- exp(ana[,names(coef_p1)] %*% coef_p1)
rst_p1_prob <- rst_p1_odd / (rst_p1_odd + 1) 
rst_p2 <- ana[, names(coef_p2)] %*% coef_p2
rst2_0 <- rst_p1_prob * rst_p2

rst2 <- bind_cols(dat %>% 
                      mutate(age2 = age + year - 1,
                             age2 = ifelse(age2 <30, 30, ifelse(age2>70, 70, age2)),
                             age2 = as.character(floor(age2))) %>% 
                      select(age2, mi), 
                    cost = as.vector(rst2_0)) %>% 
  group_by(mi) %>% 
  summarize(cost = mean(cost))

# > Predicted costs for 1 patient ----

# Patient profiles

pat <- expand_grid(age = 30:70,
                   male = 0:1,
                   sbp = 140,
                   db = 0,
                   mi = 0:4,
                   stroke = 0,
                   vd = 0,
                   nvd = 0) 

ana <- pat %>% 
  transmute("(Intercept)" = 1,
            male1 = ifelse(male == 1, 1, 0),
            sbp = (sbp - 140) / 20,
            db1 = ifelse(db == 1, 1, 0),
            cur_age = (age - 60) / 10,
            mi1 = ifelse(mi == "1", 1, 0),
            mi2 = ifelse(mi == "2", 1, 0),
            mi3 = ifelse(mi == "3", 1, 0),
            mi4 = ifelse(mi == "4", 1, 0),
            stroke1 = ifelse(stroke == "1", 1, 0),
            stroke2 = ifelse(stroke == "2", 1, 0),
            stroke3 = ifelse(stroke == "3", 1, 0),
            stroke4 = ifelse(stroke == "4", 1, 0),
            vd1 = ifelse(vd == 1, 1, 0),
            nvd1 = ifelse(nvd == 1, 1, 0)) %>% 
  as.matrix()

rst_p1_odd <- exp(ana[,names(coef_p1)] %*% coef_p1)
rst_p1_prob <- rst_p1_odd / (rst_p1_odd + 1) 
rst_p2 <- ana[, names(coef_p2)] %*% coef_p2
rst3_0 <- rst_p1_prob * rst_p2

rst3 <- bind_cols(pat %>% select(age, male, mi) %>% mutate(mi = factor(mi, levels = c(0:4))),
                  tibble(cost = as.vector(rst3_0))) 


# > Plot ----

dat_plot <- bind_rows(rst1 %>% mutate(type = "crude"),
                      rst2 %>% mutate(type = "model"),
                      rst3 %>% mutate(type = ifelse(male == 0, "model_f", "model_m")) %>% select(-male)
                      ) %>% 
  mutate(mi = factor(mi, levels = as.character(0:4),
                     labels = c("No MI", "Same year of MI", "1 year after MI", "2 years after MI", "\u22653 years after MI"))
         )

p <- ggplot() +
  geom_line(data = dat_plot %>% filter(type %in% c("model_f", "model_m")), 
            aes(x = age, y = cost, color = type),
            linetype = 2, linewidth = 1, alpha = 0.5) + 
  geom_hline(data = dat_plot %>% filter(!type %in% c("model_f", "model_m")), 
             aes(yintercept = cost, color = type),
            linetype = 1, linewidth = 1, alpha = 0.5) + 
  facet_wrap(~mi, nrow = 1) + 
  scale_color_manual(name = "Estimation", 
                       values = c("crude" = "black", 
                                  "model" = "#00BFC4",
                                  "model_f" = "#F8766D", 
                                  "model_m" = "#619CFF"),
                       labels = c("crude" = "Crude population mean",
                                  "model" = "Predicted population mean",
                                  "model_f" = "Predicted for a hypothetical female",
                                  "model_m" = "Predicted for a hypothetical male")) + 
  labs(x = "Age, years",
       y = "Annual costs (£) of health state") + 
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size=17),
        legend.text=element_text(size=15),
        strip.text = element_text(size = 15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17,face="bold")) + 
  guides(color=guide_legend(nrow=2,byrow=TRUE))
  
ggsave(p, filename  = file.path(path_output, "step4_fig_esimate_mod_vs_raw.png") , 
       height = 150, width = 300, units = "mm")

