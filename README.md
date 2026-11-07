# Prognostic Value of Peritumoral CXCL17 in Hepatocellular Carcinoma

## Overview
This repository contains an end-to-end clinical survival analysis evaluating the independent association between continuous peritumoral chemokine CXC ligand 17 (CXCL17) expression and overall survival (OS) in patients with hepatocellular carcinoma (HCC) following curative surgical resection.

Whereas the original published study ([Li et al., 2014](https://doi.org/10.1371/journal.pone.0110064)) dichotomized peritumoral CXCL17 expression at the sample median, this project preserves the continuous dose-response relationship. The workflow demonstrates advanced survival modeling techniques in R, specifically addressing non-linearity via restricted cubic splines, resolving proportional hazards (PH) violations through multivariable stratification, and verifying model stability via deviance and DFBETAS residual diagnostics.

## Data Source
The analysis utilizes the `hepatoCellular` dataset from the **`asaur`** R package, derived from 227 patients who underwent curative resection for HCC at Sun Yat-sen University Cancer Center. The final complete-case analytic cohort ($n = 227$) captures 97 observed mortality events over a median follow-up of 36 months (range: 2–83 months).

## Statistical Methodology
A multivariable Cox proportional hazards regression model was developed according to a pre-specified Statistical Analysis Plan (SAP) to quantify the independent effect of continuous peritumoral CXCL17 density (cells/megapixel) on overall survival:

* **Covariate Adjustment:** Controlled for seven pre-specified clinical and pathological confounders: tumor multiplicity (solitary vs. multiple), serum alpha-fetoprotein (AFP: $\le 25$ vs. $> 25\text{ ng/mL}$), aspartate aminotransferase (AST: $\le 42$ vs. $> 42\text{ U/L}$), tumor size ($\le 5$ vs. $> 5\text{ cm}$), vascular invasion (absent vs. present), TNM stage (I–II vs. III–IV), and BCLC stage (0–A vs. B–C).
* **Assumption Diagnostics & Remediations:**
  * **Linearity (Martingale Residuals):** Null Martingale residual plots with LOWESS smoothing revealed a non-linear functional form for continuous peritumoral CXCL17 density. This was remediated by modeling CXCL17 with a **3-knot restricted cubic spline**.
  * **Proportional Hazards (Scaled Schoenfeld Residuals):** Initial `cox.zph()` testing indicated a significant global PH violation ($\chi^2 = 17.61, p = 0.024$) driven by five categorical staging and pathology covariates ($p < 0.05$). Stratifying the baseline hazard across tumor multiplicity, tumor size, vascular invasion, TNM stage, and BCLC stage fully resolved the violation in the final model (global $\chi^2 = 3.21, \text{df} = 4, p = 0.52$; all individual $p > 0.10$).
  * **Outliers & Influential Observations (Deviance & DFBETAS Residuals):** Deviance residuals were symmetrically bounded between $-1.79$ and $2.72$ (within $\pm 3.0$), and maximum absolute DFBETAS across all parameters were $< 0.25$, confirming no extreme outliers or unduly influential observations.

## Key Findings
* **Primary Biomarker Contrast (75th vs. 25th Percentile):** Comparing patients at the 75th percentile of peritumoral CXCL17 density ($141.5\text{ cells/megapixel}$) to those at the 25th percentile ($15.4\text{ cells/megapixel}$) yielded a spline-adjusted hazard ratio of **1.42 (95% CI: 0.85, 2.38; $p = 0.186$)**.
* **Clinical Conclusion:** Because the 95% confidence interval spans 1.0, continuous peritumoral CXCL17 expression does not demonstrate statistically significant independent prognostic value for overall survival once non-linearity and non-proportional clinical confounders are properly accounted for. This suggests that the significant association reported in the original dichotomized analysis ($\text{HR} = 2.07, p = 0.002$) may have been sensitive to median cut-point categorization and unaddressed PH violations.
* **Secondary Clinical Predictors:** Among unstratified covariates, elevated preoperative serum **AFP ($> 25\text{ ng/mL}$)** remained a strong, statistically significant independent predictor of mortality (**HR = 2.63; 95% CI: 1.53, 4.52; $p < 0.001$**), whereas elevated **AST ($> 42\text{ U/L}$)** was not significantly associated with overall survival (**HR = 1.21; 95% CI: 0.76, 1.92; $p = 0.422$**).

## Repository Structure
```text
├── 01_Project_Script.R              # Primary R script: data prep, Table 1, Cox modeling, contrasts, & plots
├── Statistical_Analysis_Plan.pdf    # Pre-specified SAP detailing study aims, hypotheses, and methods
├── Analysis_Report.pdf              # Final clinical report with executive summary, results, and critique
├── output/                          # Generated publication-quality tables and diagnostic figures
│   ├── Table1_Baseline_Characteristics.png
│   ├── Fig1_Adjusted_Survival_Curves.png
│   ├── AppFig1_Martingale_Residuals.png
│   ├── AppFig2_Schoenfeld_Residuals.png
│   ├── AppFig3_Deviance_Residuals.png
│   └── AppFig4_DFBETAS_Residuals.png
└── README.md
