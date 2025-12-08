#Propensity score matching survival analysis for vedollizumab/anti-TNF/ustekinumab

library(haven)
library(MatchIt)
library(dplyr)
library(cobalt)
library(survival)
library(lubridate)
library(survminer)

df <- read.csv("analytic_extract_antitnf.csv")

#PSM formula
ps_formula <- treat ~
  Age_at_tnf +
  Gender + Race +
  diabetes + obesity + nicotine + psc + ra + psoriasis +
  chronic_resp + ihd + hf + ckd + stroke +
  pred_prior_any + mesalamine_prior_any + uste_prior_any + vedo_prior_any + iv_steroid_prior_any

#PSM using MatchIt
m.out <- matchit(
  ps_formula,
  data    = df,
  method  = "nearest",        
  distance = "glm",           
  link    = "logit",   
  caliper = 0.1,              
  std.caliper = TRUE,         
  ratio   = 1,                
  replace = FALSE,             
  m.order = 'random'
)

#summary(m.out)
#love.plot(m.out, threshold = 0.1)
matched_df <- match.data(m.out)

#cox model
fit_cox <- coxph(
  Surv(time_col, event_col) ~ treat,
  data    = matched_df,
  weights = weights,    
  robust  = TRUE,
  cluster = subclass    
)
summary(fit_cox)

#Kaplan-Meier curve
fit_km <- survfit(
  Surv(time_col, event_col) ~ treat,
  data = matched_df,
  weights = weights        
)
ggsurvplot(
  fit_km,
  data = matched_df,
  conf.int = TRUE,                 
  pval = TRUE,                     
  pval.coord = c(5, 0.825),
  risk.table = TRUE,              
  risk.table.col = "strata",
  ggtheme = theme_minimal(),
  palette = c("deepskyblue", "hotpink"),
  xlab = "Time (days)",
  ylab = "Colectomy-free survival",
  ylim = c(0.8, 1),
  legend.labs = c("No concomitant ACEi/ARB", "Concomitant ACEi/ARB"),
  legend.title = "Group"
)
