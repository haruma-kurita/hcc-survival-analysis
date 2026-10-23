# src/01_data_setup.R

# Load required packages
library(asaur)
library(dplyr)

# 1. Load the dataset
data("hepatoCellular")

# 2. Save a local copy to the data/ folder for your GitHub portfolio
# This ensures reproducibility for reviewers who don't have the asaur package
write.csv(hepatoCellular, "data/hepatoCellular.csv", row.names = FALSE)

# 3. Inspect the continuous primary predictor to pick our contrast levels
summary(hepatoCellular$CXCL17P)
quantile(hepatoCellular$CXCL17P, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

cleaned_hepato <- hepatoCellular |> 
  select(OS, Death, CXCL17P, Tumormultiplicity, AFP, AST, Tumorsize, Vascularinvasion, TNM, BCLC)
nrow(cleaned_hepato)
