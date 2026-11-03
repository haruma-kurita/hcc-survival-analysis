# 1. SETUP & DATA PREPARATION
library(survival)
library(survminer)
library(doBy) 
library(asaur)
library(dplyr)
library(Hmisc) # Added for spline adjustments

data("hepatoCellular")

# Data Cleaning
df_clean <- hepatoCellular |> 
  select(OS, Death, CXCL17P, Tumormultiplicity, AFP, AST, Tumorsize, Vascularinvasion, TNM, BCLC) |> 
  na.omit()

# Convert clinical covariates to factors (reference cell coding)
df_clean <- df_clean %>%
  mutate(across(c(Tumormultiplicity, AFP, AST, Tumorsize, Vascularinvasion, TNM, BCLC), as.factor))

summary(df_clean)


# 2. AIM 1: INITIAL MULTIVARIABLE COX MODEL
# Fit the initial unadjusted Cox proportional hazards model
cox_model <- coxph(Surv(df_clean$OS, df_clean$Death) ~ CXCL17P + Tumormultiplicity + AFP + AST + Tumorsize + 
                     Vascularinvasion + TNM + BCLC, data = df_clean)

summary(cox_model)


# 3. INITIAL CONTRAST: 75th vs 25th Percentile of CXCL17P
# Get IQR specifically for CXCL17P
q1_cxcl <- quantile(df_clean$CXCL17P, 0.25)
q3_cxcl <- quantile(df_clean$CXCL17P, 0.75)
L_matrix <- rbind(c(q3_cxcl - q1_cxcl, 0, 0, 0, 0, 0, 0, 0))

# Calculate and print the contrast estimate, CI, and p-value
cxcl_contrast <- esticon(cox_model, L_matrix, conf.int = TRUE)
print(cxcl_contrast)

# Print exponentiated Hazard Ratio and CI for interpretation
paste0("Initial Hazard Ratio for 75th vs 25th Percentile: ",
       round(exp(cxcl_contrast$estimate), 2), " (",
       round(exp(cxcl_contrast$lwr), 2), ", ",
       round(exp(cxcl_contrast$upr), 2), ")")


# 4. INITIAL VISUAL: ADJUSTED SURVIVAL CURVES (Average Patient Approach)
new_data_plot <- data.frame(
  CXCL17P = c(q1_cxcl, q3_cxcl),
  Tumormultiplicity = factor(c(1, 1), levels = levels(df_clean$Tumormultiplicity)),
  AFP = factor(c(1, 1), levels = levels(df_clean$AFP)),
  AST = factor(c(1, 1), levels = levels(df_clean$AST)),
  Tumorsize = factor(c(1, 1), levels = levels(df_clean$Tumorsize)),
  Vascularinvasion = factor(c(0, 0), levels = levels(df_clean$Vascularinvasion)),
  TNM = factor(c(1, 1), levels = levels(df_clean$TNM)),
  BCLC = factor(c(1, 1), levels = levels(df_clean$BCLC))
)

surv_avg <- survfit(cox_model, newdata = new_data_plot)

ggsurvplot(surv_avg, data = new_data_plot,
           legend.title = "CXCL17P Percentile",
           legend.labs = c("25th Percentile", "75th Percentile"))


# 5. DIAGNOSTICS & ASSUMPTIONS
# 5a. Proportional Hazards: Scaled Schoenfeld Residuals
ph_test <- cox.zph(cox_model)
print(ph_test)
ggcoxzph(ph_test)

# 5b. Linearity of Continuous Predictor: Null Martingale Residuals
ggcoxfunctional(Surv(OS, Death) ~ CXCL17P, data = df_clean)

# 5c. Outliers: Deviance Residuals
ggcoxdiagnostics(cox_model, 
                 type = "deviance", 
                 ox.scale = "linear.predictions")

# 5d. Influential Points: DFBETAS
ggcoxdiagnostics(cox_model, 
                 type = "dfbetas", 
                 ox.scale = "observation.id")


# ---------------------------------------------------------
# 6. UPDATED MODEL: FIXING LINEARITY & PROPORTIONAL HAZARDS
# ---------------------------------------------------------

# 6a. Fix Linearity: Generate Restricted Cubic Splines for CXCL17P (3 knots)
# using norm = 0 and extracting knots 
spline_cov = rcspline.eval(df_clean$CXCL17P, nk = 3, norm = 0)
knots = attr(spline_cov, "knots")

# 6b. Fix Proportional Hazards: Stratify violators and refit the Cox model
# Evaluating the spline directly inside the formula alongside the main effect
cox_model_fixed = coxph(Surv(OS, Death) ~ CXCL17P + rcspline.eval(CXCL17P, knots, norm = 0) + 
                          AFP + AST + 
                          strata(Tumormultiplicity) + 
                          strata(Tumorsize) + 
                          strata(Vascularinvasion) + 
                          strata(TNM) + 
                          strata(BCLC), 
                        data = df_clean)

summary(cox_model_fixed)

# 6c. Updated Contrast: 75th vs 25th Percentile (Spline Adjusted)
cxcl_q3vq1 = rbind(
  c(q3_cxcl - q1_cxcl, 
    rcspline.eval(q3_cxcl, knots, norm = 0) - rcspline.eval(q1_cxcl, knots, norm = 0), 
    0, 0)
)

cxcl_contrasts = esticon(cox_model_fixed, cxcl_q3vq1, conf.int = TRUE)

print("--- ADJUSTED FINAL HAZARD RATIO ---")
paste0("Adjusted HR for 75th vs 25th Percentile: ",
       round(exp(cxcl_contrasts$estimate), 2), " (",
       round(exp(cxcl_contrasts$lwr), 2), ", ",
       round(exp(cxcl_contrasts$upr), 2), ")")

# 6d. Updated Visuals: Adjusted Survival Curves
new_data_plot_fixed = data.frame(
  CXCL17P = c(q1_cxcl, q3_cxcl),
  AFP = factor(c(1, 1), levels = levels(df_clean$AFP)),
  AST = factor(c(1, 1), levels = levels(df_clean$AST)),
  Tumormultiplicity = factor(c(1, 1), levels = levels(df_clean$Tumormultiplicity)),
  Tumorsize = factor(c(1, 1), levels = levels(df_clean$Tumorsize)),
  Vascularinvasion = factor(c(0, 0), levels = levels(df_clean$Vascularinvasion)),
  TNM = factor(c(1, 1), levels = levels(df_clean$TNM)),
  BCLC = factor(c(1, 1), levels = levels(df_clean$BCLC))
)

surv_avg_fixed = survfit(cox_model_fixed, newdata = new_data_plot_fixed)

ggsurvplot(surv_avg_fixed, data = new_data_plot_fixed,
           legend.title = "CXCL17P Percentile",
           legend.labs = c("25th Percentile", "75th Percentile"))
