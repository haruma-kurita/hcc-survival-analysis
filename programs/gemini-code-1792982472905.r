# ==============================================================================
# Script: 02_cox_ph_modeling.R
# Purpose: Multivariable Cox PH model evaluating continuous CXCL17P on OS,
#          linear contrasts, adjusted survival curves, and model diagnostics.
# Dataset: hepatoCellular (from asaur package)
# ==============================================================================

# 1. Load Required Libraries
# Source: Module 4, Slide 4 (survival, survminer)
# Source: Module 8, Slide 24 (asaur)
# Source: Module 9, Slide 22 (doBy) & Slide 30 (Hmisc)
# Source: Module 11, Slide 13 (ggplot2) & Module 12, Slide 12 (dplyr)
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
# Source: Module 8, Slide 24 & Module 12, Slide 12 (Loading & selecting from asaur)
data("hepatoCellular")

hcc_clean = hepatoCellular %>%
  select(OS, Death, CXCL17P, Tumormultiplicity, AFP, AST, 
         Tumorsize, Vascularinvasion, TNM, BCLC) %>%
  na.omit() %>%
  # Source: Module 6, Slide 19 & Module 9, Slides 9 & 11 (Factor conversion for 
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
# Source: Module 4, Slides 5-6; Module 8, Slide 18; Module 12, Slide 12
surv_obj_hcc = Surv(hcc_clean$OS, hcc_clean$Death)

# ==============================================================================
# 3. Fit the Multivariable Cox Proportional Hazards Model
# ==============================================================================
# Source: Module 8, Slide 19 & Module 9, Slides 8 & 11
# Note: x = TRUE is included per Module 10, Slide 24 to retain the covariate 
# matrix in the output object if needed for adjusted curve packages.
cox_model = coxph(
  surv_obj_hcc ~ CXCL17P + Tumormultiplicity + AFP + AST + 
                 Tumorsize + Vascularinvasion + TNM + BCLC, 
  data = hcc_clean,
  x = TRUE
)

# Output log-hazard coefficients, Hazard Ratios exp(coef), 95% CIs, Wald/LRT/Score tests
# Source: Module 8, Slides 15-19 & Module 9, Slide 8
summary(cox_model)

# Output Harrell's C Index (Concordance)
# Source: Module 10, Slide 27
concordance(cox_model)

# ==============================================================================
# 4. Contrast: 75th vs. 25th Percentile of Continuous CXCL17P
# ==============================================================================
# Pre-specified contrast levels (25th and 75th percentiles of CXCL17P)
cxcl17p_quantiles = quantile(hcc_clean$CXCL17P, probs = c(0.25, 0.75))
p25 = as.numeric(cxcl17p_quantiles[1])
p75 = as.numeric(cxcl17p_quantiles[2])
delta = p75 - p25

# Create the L matrix for linear combination l^T * beta
# Source: Module 9, Slides 20-22 & Slide 32
# CXCL17P is the 1st coefficient (weighted by delta); the 7 confounders are set to 0
L = rbind(
  c(delta, 0, 0, 0, 0, 0, 0, 0)
)

# Run the contrast using esticon()
# Source: Module 9, Slide 22
cxcl17p_75_v_25 = esticon(cox_model, L, conf.int = TRUE)

# Output raw contrast table (provides estimate, std.error, statistic, and p.value)
# Source: Module 9, Slide 22
cxcl17p_75_v_25

# Output log-hazard 95% confidence interval bounds
# Source: Module 9, Slide 22
c(cxcl17p_75_v_25$lwr, cxcl17p_75_v_25$upr)

# Exponentiate to get the Hazard Ratio and 95% CI for the 75th vs 25th percentile
# Source: Module 9, Slide 22 & Slide 32
paste0("Hazard Ratio for CXCL17P (75th vs 25th Percentile): ",
       round(exp(cxcl17p_75_v_25$estimate), 2), " (",
       round(exp(cxcl17p_75_v_25$lwr), 2), ", ",
       round(exp(cxcl17p_75_v_25$upr), 2), "), p-value = ",
       signif(cxcl17p_75_v_25$p.value, 3))

# ==============================================================================
# 5. Adjusted Survival Curves for the Two CXCL17P Contrast Levels
# ==============================================================================

# --- Method A: Direct Adjusted Survival (Average Survival across Cohort) ---
# Source: Module 10, Slides 15, 21, & 23 ("Trick ggadjustedcurves into thinking 
# you have 2 values only")
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

# --- Method B: Conditional Adjusted Survival (Reference Patient Profile) ---
# Source: Module 10, Slides 13, 19-20, & 22 (Using survfit with newdata and ggsurvplot)
surv_ref_hcc = survfit(
  cox_model,
  newdata = data.frame(
    CXCL17P           = c(p25, p75),
    Tumormultiplicity = factor(c(1, 1), levels = levels(hcc_clean$Tumormultiplicity)),
    AFP               = factor(c(1, 1), levels = levels(hcc_clean$AFP)),
    AST               = factor(c(1, 1), levels = levels(hcc_clean$AST)),
    Tumorsize         = factor(c(1, 1), levels = levels(hcc_clean$Tumorsize)),
    Vascularinvasion  = factor(c(0, 0), levels = levels(hcc_clean$Vascularinvasion)),
    TNM               = factor(c(1, 1), levels = levels(hcc_clean$TNM)),
    BCLC              = factor(c(1, 1), levels = levels(hcc_clean$BCLC))
  )
)

ggsurvplot(
  surv_ref_hcc, 
  data = hcc_clean,
  legend.labs = c("25th Percentile (Low)", "75th Percentile (High)"),
  xlab = "Time (Months)",
  ylab = "Survival probability"
)

# ==============================================================================
# 6. Model Assumptions and Diagnostics (Modules 11 & 12)
# ==============================================================================

# ------------------------------------------------------------------------------
# A. Proportional Hazards (PH) Assumption: Scaled Schoenfeld Residuals
# Source: Module 11, Slides 18-20 & 23; Module 12, Slides 5 & 7
# ------------------------------------------------------------------------------
# Formal hypothesis test and Schoenfeld residual plots for all predictors
cox.zph(cox_model)
ggcoxzph(cox.zph(cox_model))

# Complementary Log-Log Plot (Visual check for categorical covariates, e.g., AFP)
# Source: Module 11, Slides 10, 13-14, & 22
surv_plot_afp = ggsurvplot(survfit(surv_obj_hcc ~ AFP, data = hcc_clean), 
                           data = hcc_clean)
ggplot(data = surv_plot_afp$data.survplot, 
       aes(log(time), log(-log(surv)), color = AFP)) +
  geom_point() +
  geom_path() +
  labs(title = "Complementary Log-Log Plot for AFP")

# ------------------------------------------------------------------------------
# B. Linearity Assumption (Functional Form of Continuous CXCL17P): Martingale Residuals
# Source: Module 12, Slides 8-23
# ------------------------------------------------------------------------------

# 1. Null Martingale Residuals using ggcoxfunctional()
# Source: Module 12, Slides 11, 13, & 14
ggcoxfunctional(surv_obj_hcc ~ CXCL17P, ylim = c(-1, 1), data = hcc_clean)

# 2. Martingale Residuals without CXCL17P (Reduced Model + Spline Overlays)
# Source: Module 12, Slides 11 & 17-18
cox_reduced <- coxph(
  surv_obj_hcc ~ Tumormultiplicity + AFP + AST + 
                 Tumorsize + Vascularinvasion + TNM + BCLC, 
  data = hcc_clean
)
martingale_resid_red <- residuals(cox_reduced, type = "martingale")

spline5 <- rcspline.eval(hcc_clean$CXCL17P, nk = 5, inclx = TRUE)
spline3 <- rcspline.eval(hcc_clean$CXCL17P, nk = 3, inclx = TRUE)

plot_data_red <- data.frame(
  CXCL17P = hcc_clean$CXCL17P,
  resid   = martingale_resid_red,
  spline5 = spline5,
  spline3 = spline3
)

plot_reduced = ggplot(plot_data_red, aes(x = CXCL17P, y = resid)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y ~ rcspline.eval(x, knots = attributes(spline5)$knots),
              se = FALSE, color = "blue") +
  geom_smooth(method = "lm",
              formula = y ~ rcspline.eval(x, knots = attributes(spline3)$knots),
              se = FALSE, color = "darkgreen") +
  theme_minimal() + 
  geom_smooth(color = "red") +
  geom_smooth(method = "lm", color = "orange", se = FALSE)
plot_reduced

# 3. Full Model Martingale Residuals (Checking for randomness about y = 0)
# Source: Module 12, Slides 11 & 19-20
martingale_resid_full <- residuals(cox_model, type = "martingale")

plot_data_full <- data.frame(
  CXCL17P = hcc_clean$CXCL17P,
  resid   = martingale_resid_full
)

plot_full = ggplot(plot_data_full, aes(x = CXCL17P, y = resid)) +
  geom_point() + 
  geom_hline(yintercept = 0) +
  theme_minimal() + 
  geom_smooth(color = "red")
plot_full

# ------------------------------------------------------------------------------
# C. Outliers Check: Deviance Residuals
# Source: Module 12, Slides 24-27
# Rule of thumb: absolute deviance residuals > 2 or 3 warrant investigation
# ------------------------------------------------------------------------------
ggcoxdiagnostics(
  cox_model,
  type = "deviance",
  ox.scale = "linear.predictions"
)

# ------------------------------------------------------------------------------
# D. Influential Points Check: Scaled DFBETAS
# Source: Module 12, Slides 28-31
# Rule of thumb: compare against 2 / sqrt(n) using scaled "dfbetas"
# ------------------------------------------------------------------------------
ggcoxdiagnostics(
  cox_model,
  type = "dfbetas",
  ox.scale = "observation.id"
)