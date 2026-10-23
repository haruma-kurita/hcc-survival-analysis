# ==============================================================================
# Script: 02_cox_ph_modeling.R
# Purpose: Multivariable Cox PH model evaluating continuous CXCL17P on overall survival
# Dataset: hepatoCellular (from asaur package)
# ==============================================================================

# 1. Load Required Libraries
library(survival)
library(survminer)
library(asaur)
library(doBy)      # For esticon() contrasts
library(dplyr)     # For data manipulation

# 2. Data Preparation and Cleaning (Complete-Case Analysis)
data("hepatoCellular")

hcc_clean <- hepatoCellular |> 
  # Select only the variables required for Aim 1
  select(OS, Death, CXCL17P, Tumormultiplicity, AFP, AST, 
         Tumorsize, Vascularinvasion, TNM, BCLC) |> 
  # Handle missing data via complete-case analysis
  na.omit() |> 
  # Convert categorical clinical covariates to factors for dummy coding
  mutate(
    Tumormultiplicity = factor(Tumormultiplicity),
    AFP               = factor(AFP),
    AST               = factor(AST),
    Tumorsize         = factor(Tumorsize),
    Vascularinvasion  = factor(Vascularinvasion),
    TNM               = factor(TNM),
    BCLC              = factor(BCLC)
  )

# 3. Descriptive Summaries for CXCL17P Contrast Levels
# Calculate the 25th and 75th percentiles (IQR) for CXCL17P
cxcl17p_quantiles
