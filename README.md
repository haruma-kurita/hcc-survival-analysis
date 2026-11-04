# Prognostic Value of CXCL17 in Hepatocellular Carcinoma

## Overview
This repository contains a clinical statistical analysis evaluating the association between continuous peritumoral chemokine CXC ligand 17 (CXCL17) expression and overall survival (OS) in patients with hepatocellular carcinoma (HCC). 

The project demonstrates advanced survival analysis techniques in R, specifically addressing non-linear dose-response relationships using restricted cubic splines and resolving proportional hazards (PH) assumption violations through model stratification.

## Data Source
The analysis utilizes the `hepatoCellular` dataset from the **`asaur`** R package. The final analytic dataset consists of a complete-case cohort (n=227), which includes 97 observed mortality events.

## Statistical Analysis Plan (SAP)
A multivariable Cox proportional hazards regression model was developed to quantify the independent effect of peritumoral CXCL17 expression on mortality risk. 
* **Covariate Adjustment:** The model adjusts for baseline clinical covariates including tumor multiplicity, alpha-fetoprotein (AFP), aspartate aminotransferase (AST), tumor size, vascular invasion, TNM stage, and BCLC stage.
* **Assumption Diagnostics & Corrections:**
  * **Linearity:** Evaluated using Martingale residuals. Non-linearity in the continuous CXCL17 exposure was corrected using restricted cubic splines (3 knots).
  * **Proportional Hazards:** Evaluated using scaled Schoenfeld residuals. PH violations among clinical covariates were resolved by stratifying the baseline hazard by tumor multiplicity, tumor size, vascular invasion, TNM stage, and BCLC stage.

## Key Findings
* **Primary Contrast:** When comparing patients at the 75th percentile of CXCL17 expression to those at the 25th percentile, the spline-adjusted model yielded a hazard ratio of 1.42 (95% CI: 0.85, 2.38; p = 0.186). 
* **Conclusion:** Because the confidence interval crosses 1.0 and the p-value is greater than 0.05, continuous peritumoral CXCL17 expression did not demonstrate statistically significant independent prognostic value when rigorously adjusted for standard clinical confounders. 
* **Secondary Findings:** Elevated alpha-fetoprotein (AFP Level 2) emerged as a highly significant independent predictor of poorer overall survival (HR = 2.63, p = 0.00045).

## Repository Structure
* `HCC_CXCL17_Survival_Analysis.R`: The primary execution script containing data preparation, descriptive statistics (Table 1), modeling, contrast estimation, and diagnostic testing.
* `output/`: Directory containing generated clinical and diagnostic visuals.
  * `Table1_Baseline_Characteristics`
  * `Fig1_Adjusted_Survival_Curves.png`
  * `Fig2_Schoenfeld_Residuals.png` (PH Assumption Check)
  * `Fig3_Martingale_Residuals.png` (Linearity Assumption Check)

## Technologies & Packages
* **Language:** R
* **Packages:** `survival`, `survminer` (Modeling and Kaplan-Meier visuals) | `gtsummary` (Clinical table generation) | `Hmisc` (Restricted cubic splines) | `doBy` (Contrast matrices/esticon) | `dplyr` (Data wrangling)

## Author
**Haruma Kurita**  
Master of Biostatistics Candidate  
Duke University School of Medicine, Department of Biostatistics & Bioinformatics
