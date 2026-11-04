# ==============================================================================
# Program Name:   01_Project_Script.R
# Project:        Prognostic Value of CXCL17 in Hepatocellular Carcinoma (SAP Execution)
# Author:         Haruma Kurita
# Institution:    Duke University, Department of Biostatistics & Bioinformatics
# Date:           November 3, 2026
#
# Description:    Evaluates the association between continuous peritumoral CXCL17 
#                 expression and overall survival (OS) in HCC complete cases (n=227). 
#                 Implements a multivariable Cox proportional hazards model utilizing 
#                 restricted cubic splines (3 knots) for the continuous exposure and 
#                 stratification for clinical covariates to satisfy PH assumptions.
#
# Input Data:     hepatoCellular (via 'asaur' package)
# Outputs:        - Table 1: Baseline Characteristics (gtsummary)
#                 - Figure 1: Adjusted KM Survival Curves (75th vs 25th percentile)
#                 - Appendix: Model Diagnostics (Schoenfeld & Martingale residuals)
#
# Disclosure:     Code structuring, troubleshooting, and formatting assisted by 
#                 Google Gemini.
# ==============================================================================

# 1. SETUP & DATA PREPARATION
library(gtsummary)
library(survival)
library(survminer)
library(doBy) 
library(asaur)
library(dplyr)
library(Hmisc)

data("hepatoCellular")

# Data Cleaning
df_clean <- hepatoCellular |> 
  select(OS, Death, CXCL17P, Tumormultiplicity, AFP, AST, Tumorsize, Vascularinvasion, TNM, BCLC) |> 
  na.omit()

# Convert clinical covariates to factors (reference cell coding)
df_clean <- df_clean %>%
  mutate(across(c(Tumormultiplicity, AFP, AST, Tumorsize, Vascularinvasion, TNM, BCLC), as.factor))

summary(df_clean)

# Generate Table 1: Baseline Characteristics
table_1 <- df_clean |>
  select(CXCL17P, Tumormultiplicity, AFP, AST, Tumorsize, Vascularinvasion, TNM, BCLC) |>
  tbl_summary(
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      CXCL17P ~ "Peritumoral CXCL17 Expression",
      Tumormultiplicity ~ "Tumor Multiplicity",
      Tumorsize ~ "Tumor Size",
      Vascularinvasion ~ "Vascular Invasion",
      TNM ~ "TNM Stage",
      BCLC ~ "BCLC Stage"
    )
  ) |>
  add_n() |> # Adds a column showing total valid N for each variable
  modify_header(label = "**Baseline Characteristic**") |>
  bold_labels()

# Print Table 1 to the Viewer pane
table_1

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

# 5. DIAGNOSTICS & ASSUMPTIONS (INITIAL MODEL)
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
spline_cov = rcspline.eval(df_clean$CXCL17P, nk = 3, norm = 0)
knots = attr(spline_cov, "knots")

# 6b. Fix Proportional Hazards: Stratify violators and refit the Cox model
cox_model_fixed = coxph(Surv(OS, Death) ~ CXCL17P + rcspline.eval(CXCL17P, knots, norm = 0) + 
                          AFP + AST + 
                          strata(Tumormultiplicity) + 
                          strata(Tumorsize) + 
                          strata(Vascularinvasion) + 
                          strata(TNM) + 
                          strata(BCLC), 
                        data = df_clean)

summary(cox_model_fixed)

tbl_regression(cox_model_fixed, exponentiate = TRUE)

# 6c. Updated Contrast: 75th vs 25th Percentile (Spline Adjusted)
cxcl_q3vq1 = rbind(
  c(q3_cxcl - q1_cxcl, 
    rcspline.eval(q3_cxcl, knots, norm = 0) - rcspline.eval(q1_cxcl, knots, norm = 0), 
    0, 0)
)

cxcl_contrasts = esticon(cox_model_fixed, cxcl_q3vq1, conf.int = TRUE)
cxcl_contrasts

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


# ---------------------------------------------------------
# 7. APPENDIX DIAGNOSTICS (FINAL MODEL)
# ---------------------------------------------------------

# 7a. Final Proportional Hazards Assumption (Schoenfeld Residuals)
# Note: This will only plot the unstratified variables (CXCL17, AFP, AST)
ph_check_fixed <- cox.zph(cox_model_fixed)
print("--- FINAL SCHOENFELD RESIDUALS ---")
print(ph_check_fixed)

# Print the plot to the Viewer pane for manual export
ggcoxzph(ph_check_fixed)

# 7b. Final Linearity Assumption (Martingale Residuals)
# Plotting the continuous variable against a null model's residuals to show necessity of splines
ggcoxfunctional(Surv(OS, Death) ~ CXCL17P, data = df_clean)
