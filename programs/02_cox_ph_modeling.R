# ==============================================================================
# Script: 02_cox_ph_modeling.R
# Purpose: Multivariable Cox PH model evaluating continuous CXCL17P on OS,
#          linear contrasts, adjusted survival curves, and model diagnostics.
# Dataset: hepatoCellular (from asaur package)
# ==============================================================================

# 1. Load Required Libraries
library(survival)
library(survminer)
library(asaur)
library(dplyr)
library(ggplot2)
library(doBy)
library(Hmisc)

# ==============================================================================
# 2. Data Preparation and Cleaning (Complete-Case Analysis)
# ==============================================================================

data("hepatoCellular")

hcc_clean = hepatoCellular |> 
  select(OS, Death, CXCL17P, Tumormultiplicity, AFP, AST, 
         Tumorsize, Vascularinvasion, TNM, BCLC) |> 
  na.omit() |> 
  # default reference cell / dummy coding in R)
  mutate(
    Tumormultiplicity = factor(Tumormultiplicity), # 1 = Solitary (Ref), 2 = Multiple
    AFP               = factor(AFP),               # 1 = <=25 (Ref), 2 = >25
    AST               = factor(AST),               # 1 = <=42 (Ref), 2 = >42
    Tumorsize         = factor(Tumorsize),         # 1 = <=5cm (Ref), 2 = >5cm
    Vascularinvasion  = factor(Vascularinvasion),  # 0 = Absent (Ref), 1 = Present
    TNM               = factor(TNM),               # 1 = I+II (Ref), 2 = III+IV
    BCLC              = factor(BCLC)               # 1 = 0-A (Ref), 2 = B-C
  )
# Create the Survival Object
surv_obj_hcc = Surv(hcc_clean$OS, hcc_clean$Death)

# ==============================================================================
# 3. Fit the Multivariable Cox Proportional Hazards Model
# ==============================================================================

cox_model = coxph(
  surv_obj_hcc ~ CXCL17P + Tumormultiplicity + AFP + AST + Tumorsize + 
    Vascularinvasion + TNM + BCLC, 
  data = hcc_clean,
  x = TRUE # Stores the data matrix; required for post-hoc testing (like esticon/cox.zph)
)

summary(cox_model)

# Output Harrell's C Index (Concordance)
concordance(cox_model)

# ==============================================================================
# 4. Contrast: 75th vs. 25th Percentile of Continuous CXCL17P
# ==============================================================================

cxcl17p_quantiles = quantile(hcc_clean$CXCL17P, probs = c(0.25, 0.75))
p25 = as.numeric(cxcl17p_quantiles[1])
p75 = as.numeric(cxcl17p_quantiles[2])
delta = p75 - p25

# Create the L matrix for linear combination l^T * beta
# CXCL17P is the 1st coefficient (weighted by delta); the 7 confounders are set to 0
L = rbind(
  c(delta, 0, 0, 0, 0, 0, 0, 0)
)

# Run the contrast using esticon()
cxcl17p_75_v_25 = esticon(cox_model, L, conf.int = TRUE)

# Output raw contrast table (provides estimate, std.error, statistic, and p.value)
cxcl17p_75_v_25

# Output log-hazard 95% confidence interval bounds
c(cxcl17p_75_v_25$lwr, cxcl17p_75_v_25$upr)

# Exponentiate to get the Hazard Ratio and 95% CI for the 75th vs 25th percentile
paste0("Hazard Ratio for CXCL17P (75th vs 25th Percentile): ",
       round(exp(cxcl17p_75_v_25$estimate), 2), " (",
       round(exp(cxcl17p_75_v_25$lwr), 2), ", ",
       round(exp(cxcl17p_75_v_25$upr), 2), "), p-value = ",
       signif(cxcl17p_75_v_25$p.value, 3))

# ==============================================================================
# 5. Adjusted Survival Curves for the Two CXCL17P Contrast Levels
# ==============================================================================

# Direct Adjusted Survival (Average Survival across Cohort)

hcc_trick = hcc_clean %>% 
  mutate(CXCL17P = c(p25, rep(p75, nrow(hcc_clean) - 1)))

ggadjustedcurves(
  cox_model, 
  variable = "CXCL17P", 
  data = hcc_trick
) + 
  labs(title = "Direct Adjusted Survival: 75th vs 25th Percentile of CXCL17P",
       x = "Time (Months)",
       y = "Survival Rate",
       color = "CXCL17P Percentile")
