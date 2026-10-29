# Colorectal Cancer Survival Analysis (Hospital Cancer Registry, Mato Grosso, Brazil)

R analysis of overall survival (OS) and cancer-specific survival (CSS) in patients with colorectal cancer
(ICD-O-3 C18-C20), using Kaplan-Meier curves, log-rank tests and Cox proportional hazards models.

## Contents
- `colorectal_cancer_survival.R`: full pipeline (cohort definition, recoding, survival times, KM, log-rank, Cox, PH checks)

## Methods
1. **Cohort:** sites C18, C19, C20; rare morphologies and death-certificate-only cases excluded; age >= 18
2. **Time and events:** time from diagnosis to death or administrative censoring (2021-12-24), truncated at 10 years;
   OS = death from any cause, CSS = death from cancer (other deaths censored)
3. **Bivariate:** frequencies, event proportions and log-rank tests
4. **Cox models:** univariable, then multivariable (age, race, sex, education, city, morphology, stage, topography, marital status)
5. **Diagnostics:** proportional-hazards test (Schoenfeld residuals)

## Data (not included)
Patient-level data are confidential and are excluded by `.gitignore`. Place the Excel file at `data/banco.xlsx`
(or edit `DATA_FILE`). Expected columns: `Sexo, Idade, RacaCor, Grau_de_Instrucao, Estado_Civil, Cidade_Endereco,
Codigo_da_Topografia, Descricao_da_Morfologia, Codigo_da_Morfologia, Meio_de_Diagnostico, Extensao, Status_Vital,
Tipo_do_Obito, Data_do_obito, Data_de_Diagnostico, Ano_de_Diagnostico` (plus a few unused ones listed in the script).

## Run
Open the script in RStudio and run it section by section. Packages are installed automatically
(`tidyverse, janitor, readxl, survival, survminer, broom`). Figures and tables go to `output/`.
