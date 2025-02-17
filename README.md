# PHACS Project

This repository contains the code and documentation for the PHACS project. The project focuses on processing and analyzing biomarker data, calculating pace of aging, and assessing HIV viral load.

## Table of Contents

- [Main Workflow](#main-workflow)
  - [Helper Functions](#helper-functions)
  - [Initial Biomarker Merge](#initial-biomarker-merge)
  - [Biomarker Cleaning](#biomarker-cleaning)
  - [Pace of Aging Calculations](#pace-of-aging-calculations)
  - [HIV Viral Load Assessments](#hiv-viral-load-assessments)
    - [Outcomes](#outcomes)
- [Regression Templates](#regression-templates)
- [Contributing](#contributing)
- [License](#license)

## Main Workflow

### Helper Functions

- **Helper_Functions.R**  
  Contains all custom helper functions.  
  **Note:** Run this file at `root.R` first.

### Initial Biomarker Merge

- **Initial_Biomarker.Rmd**  
  Merges all the necessary biomarker files and handles NaN values and unnecessary columns.

  **Input Files:**
  - `evw0173.sas7bdat` (blood pressure)
  - *(Routine Chemistries)*
  - `lbw0071.sas7bdat` (HOMA-Insulin Resistance)
  - `lbw0073.sas7bdat` (Routine Hematologies)
  - `evw0169.sas7bdat` (Height, Weight, Body measurements)
  - `dmw0080.sas7bdat` (HIV Status)
  - `master.sas7bdat` (Master Sheet for Participant Demographics)

  **Output Files:**
  - `m4_dem.csv` (Participant demographics and all merged biomarkers)
  - `m4_fin.csv` (Merged biomarkers only)

### Biomarker Cleaning

- **Biomarker_Cleaning.Rmd**  
  Reads in the merged biomarker file and cleans the data for modeling. This file also cleans the viral load data and produces three different cleaned versions based on feature selection.

  **Input Files:**
  - `m4_dem.csv`
  - `m4_fin.csv`
  - `f3109.sas7bdat` (viral load)
  - `hxw0101.sas7bdat`

  **Output Files:**
  - `ivs_stand_flipped.csv` (Participant and Biomarker Data ready for modeling)
  - `ivs_stand_flipped_x.csv` (Alternate version for feature selection)

### Pace of Aging Calculations

Each script in this section outputs its own Pace of Aging file:

- **PoA_V1.Rmd**  
  Calculates the Pace of Aging for all three methods.

  **Input Files:**
  - `ivs_stand_flipped.csv` (from Biomarker_Cleaning.Rmd)
  - `hxw0101.sas7bdat` (Viral Load)

  **Output:**  
  Data frames with Pace of Aging for each method.

- **PoA_V2.Rmd**  
  Uses biomarkers chosen purely by correlation with age. Only biomarkers with a correlation above the median are used.

- **PoA_V3.Rmd** *(Current Working Model)*  
  Biomarkers chosen based on stepAIC.

- **PoA_V4.Rmd**  
  Biomarkers chosen based on `bm~weeks…` that are not singular or have convergence issues.

### HIV Viral Load Assessments

#### Outcomes

- **Outcome Measures – Copy.Rmd**  
  Organizes and cleans outcome measures.  
  *(Originally done by Kevin)*

- **PoA_Outcome_Merge.Rmd**  
  Reads in the cleaned outcomes and merges the Pace of Aging and outcome files, preparing the data for analyses.

## Regression Templates

- **run_regression.R**  
  Reads in merged files from `PoA_Outcome_Merge.Rmd` and runs the regression templates.

- **Regression_Template.qmd**  
  Template for running regressions and outputting comparison metrics.  
  **Usage:** Run from `run_regression.R`.

- **ZeroInflation_Regression_Temp.qmd**  
  Template for running zero-inflated regressions and outputting comparison metrics.  
  **Usage:** Run from `run_regression.R`.

- **BinLog_Reg_Template.qmd**  
  Template for running binomial logistic regression.  
  **Usage:** Run from `run_regression.R`.

## Contributing

Contributions are welcome! Please follow the standard GitHub guidelines for contributing to this repository.

## License

[Specify your license here]

---

*This project is part of the PHACS initiative.*

