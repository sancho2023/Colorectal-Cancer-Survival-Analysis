# =============================================================================
# Survival Analysis of Colorectal Cancer (ICD-O-3 C18-C20)
# Hospital-based cancer registry, Mato Grosso, Brazil
# =============================================================================
#
# Goal:     Estimate overall survival (OS) and cancer-specific survival (CSS)
#           and identify prognostic factors with Kaplan-Meier curves, log-rank
#           tests and Cox proportional hazards models.
#
# Input:    An Excel file with one row per patient (see README for the columns).
#           Patient-level data are NOT included in this repository.
#
# Sections:
#   1. Setup (packages, paths, parameters, helper functions)
#   2. Data import and cohort definition
#   3. Recoding of covariates
#   4. Survival time and event definitions
#   5. Kaplan-Meier: overall vs. cancer-specific survival
#   6. Bivariate analysis: descriptives and log-rank tests
#   7. Univariable Cox models
#   8. Multivariable Cox model: overall survival
#   9. Multivariable Cox model: cancer-specific survival
#  10. Time-varying effects (preparation)
#
# Run the script top to bottom, section by section, in RStudio.
# =============================================================================


# =============================================================================
# 1. SETUP ----
# =============================================================================

# --- 1.1 Packages ------------------------------------------------------------
pkgs <- c("tidyverse", "janitor", "readxl", "survival", "survminer", "broom")
missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)  # data wrangling (dplyr, forcats, purrr) and ggplot2
  library(janitor)    # tabyl() frequency tables
  library(readxl)     # read Excel files
  library(survival)   # Surv(), survfit(), survdiff(), coxph(), cox.zph()
  library(survminer)  # ggsurvplot(), ggcoxzph(), surv_median()
  library(broom)      # tidy model output
})

# --- 1.2 Paths ---------------------------------------------------------------
# Edit these lines for your machine. Do NOT commit the data file.
DATA_FILE <- "data/banco.xlsx"   # patient-level Excel file (not in the repository)
OUT_DIR   <- "output"            # figures and tables are saved here
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- 1.3 Analysis parameters -------------------------------------------------
CUTOFF_DATE <- as.Date("2021-12-24")  # administrative end of follow-up
HORIZON     <- 10                     # follow-up truncated at 10 years

# --- 1.4 Helper: save a plot as a TIFF file ----------------------------------
save_tiff <- function(plot, file, width, height, units = "mm", res = 600) {
  tiff(file.path(OUT_DIR, file), width = width, height = height,
       units = units, res = res, compression = "lzw")
  print(plot)
  dev.off()
}


# =============================================================================
# 2. DATA IMPORT AND COHORT DEFINITION ----
# =============================================================================

# --- 2.1 Import and keep the variables of interest ---------------------------
raw <- read_excel(DATA_FILE)

cohort <- raw |>
  select(Sexo, DTNASC, Idade, Idade_Categorizada, RacaCor, Naturalidade_Estado,
         Naturalidade, Grau_de_Instrucao, Estado_Civil, Cidade_Endereco,
         Descricao_da_Topografia, Codigo_da_Topografia, Descricao_da_Morfologia,
         Descricao_da_Doenca, Codigo_da_Doenca_categorizado, Meio_de_Diagnostico,
         Extensao, Status_Vital, Tipo_do_Obito, Data_do_obito,
         Data_de_Diagnostico, Ano_de_Diagnostico, Codigo_da_Morfologia)

# --- 2.2 Select the cancer sites ---------------------------------------------
# C18 = colon (C180-C189), C19 = rectosigmoid junction, C20 = rectum.
# C21 (anus) is excluded.
table(cohort$Codigo_da_Topografia)

cohort <- cohort |>
  filter(Codigo_da_Topografia %in% c("C180", "C181", "C182", "C183", "C184",
                                     "C185", "C186", "C187", "C188", "C189",
                                     "C199", "C209"))

# --- 2.3 Exclude rare / non-informative morphology codes ---------------------
cohort <- cohort |>
  filter(!Codigo_da_Morfologia %in% c(80001, 80102, 80702, 81402,
                                      82102, 82612, 82632))

# --- 2.4 Exclude cases diagnosed only by death certificate (SDO) -------------
# Survival time is not reliable when the diagnosis comes from the death certificate.
cohort <- cohort |> filter(Meio_de_Diagnostico != "SDO")


# =============================================================================
# 3. RECODING OF COVARIATES ----
# =============================================================================
# Labels are recoded to English. Source values (Portuguese) are kept in the
# case_when() conditions because they come from the registry.

cohort <- cohort |>
  mutate(
    # --- Sex ---
    sex = case_when(Sexo == "FEMININO"  ~ "Female",
                    Sexo == "MASCULINO" ~ "Male"),

    # --- Race / skin color ---
    race = case_when(RacaCor %in% c("AMARELA", "INDÍGENA", "PARDA") ~ "Other",
                     RacaCor == "BRANCO"         ~ "White",
                     RacaCor == "PRETA"          ~ "Black",
                     RacaCor == "SEM INFORMAÇÃO" ~ "No_information"),

    # --- City of residence ---
    city = case_when(Cidade_Endereco == "CUIABÁ"        ~ "Cuiaba",
                     Cidade_Endereco == "VÁRZEA GRANDE" ~ "Varzea_Grande",
                     Cidade_Endereco == "RONDONÓPOLIS"  ~ "Rondonopolis",
                     Cidade_Endereco == "SINOP"         ~ "Sinop",
                     TRUE                               ~ "Other"),

    # --- Education ---
    education = case_when(
      Grau_de_Instrucao %in% c("FUNDAMENTAL I (1ª A 4ª SÉRIE)",
                               "FUNDAMENTAL II (5ª A 8ª SÉRIE)") ~ "Elementary",
      Grau_de_Instrucao == "MÉDIO (ANTIGO SEGUNDO GRAU)"          ~ "High_school",
      Grau_de_Instrucao %in% c("SUPERIOR COMPLETO",
                               "SUPERIOR INCOMPLETO")             ~ "College",
      Grau_de_Instrucao == "SEM ESCOLARIDADE"                     ~ "No_schooling",
      Grau_de_Instrucao == "SEM INFORMAÇÃO"                       ~ "No_information"),

    # --- Marital status ---
    marital = case_when(
      Estado_Civil %in% c("CASADO", "UNIÃO CONSENSUAL")                           ~ "Married",
      Estado_Civil %in% c("SOLTEIRO", "SEPARADO JUDICIALMENTE", "VIÚVO")          ~ "Not_married",
      Estado_Civil == "SEM INFORMAÇÃO"                                            ~ "No_information"),

    # --- Topography (ICD-O site) ---
    topography = case_when(
      Codigo_da_Topografia %in% c("C180", "C181", "C182", "C183", "C184",
                                  "C185", "C186", "C187", "C188", "C189") ~ "C18",
      Codigo_da_Topografia == "C199" ~ "C19",
      Codigo_da_Topografia == "C209" ~ "C20"),

    # --- Histological class ---
    # AC = conventional adenocarcinoma; MC = mucinous (incl. signet-ring cell);
    # NOS = malignant neoplasm, not otherwise specified; OTHER = rare variants.
    morphology = case_when(
      Descricao_da_Morfologia %in% c(
        "ADENOCARCINOMA- SOE", "ADENOCARCINOMA TUBULAR", "ADENOCARCINOMA VILOSO",
        "ADENOCARCINOMA EM ADENOMA TUBULOVILOSO", "ADENOCARCINOMA EM ADENOMA VILOSO",
        "ADENOCARCINOMA EM POLIPO ADENOMATOSO", "ADENOCARCINOMA PAPILAR- SOE",
        "ADENOCARCINOMA TRABECULAR", "ADENOCARCINOMA DE GLANDULAS DA REGIAO ANAL",
        "ADENOCARCINOMA DE PROPAGACAO SUPERFICIAL") ~ "AC",
      Descricao_da_Morfologia %in% c(
        "ADENOCARCINOMA MUCINOSO", "ADENOCARCINOFIBROMA MUCINOSO",
        "ADENOCARCINOMA PRODUTOR DE MUCINA", "CISTADENOCARCINOMA MUCINOSO- SOE",
        "CARCINOMA DE CELULAS EM ANEL DE SINETE") ~ "MC",
      Descricao_da_Morfologia %in% c("NEOPLASIA MALIGNA",
                                     "NEOPLASMA MALIGNO") ~ "NOS",
      TRUE ~ "OTHER"),

    # --- Disease extent (stage at diagnosis) ---
    stage = case_when(Extensao == "IN SITU"        ~ "In_situ",
                      Extensao == "LOCALIZADO"     ~ "Localized",
                      Extensao == "METÁSTASE"      ~ "Metastasis",
                      Extensao == "SEM INFORMAÇÃO" ~ "No_information"),

    # --- Age categories ---
    age_cat = case_when(Idade < 50              ~ "Under50",
                        Idade >= 50 & Idade < 60 ~ "50_59",
                        Idade >= 60 & Idade < 70 ~ "60_69",
                        Idade >= 70              ~ "70_plus"),

    # --- Year of diagnosis (descriptive only) ---
    year_diag = case_when(Ano_de_Diagnostico <= 2005                          ~ "Up_to_2005",
                          Ano_de_Diagnostico > 2005 & Ano_de_Diagnostico <= 2011 ~ "2006_2011",
                          Ano_de_Diagnostico > 2011                           ~ "2012_2018")
  )

# --- Remaining exclusions: invalid stage, missing age, minors ----------------
# Stage values outside the four categories above (e.g. "not applicable") become NA.
cohort <- cohort |>
  filter(!is.na(stage), !is.na(Idade), Idade >= 18)

# Quick checks of the recoded variables
lapply(cohort[c("sex", "race", "city", "education", "marital",
                "topography", "morphology", "stage", "age_cat", "year_diag")], table)


# =============================================================================
# 4. SURVIVAL TIME AND EVENT DEFINITIONS ----
# =============================================================================

# --- 4.1 Follow-up time ------------------------------------------------------
# Deaths: death date - diagnosis date. Survivors are censored at CUTOFF_DATE.
cohort <- cohort |>
  mutate(
    time_days = if_else(
      Status_Vital == "MORTO",
      as.numeric(difftime(Data_do_obito, Data_de_Diagnostico, units = "days")),
      as.numeric(difftime(CUTOFF_DATE,   Data_de_Diagnostico, units = "days"))
    )
  ) |>
  filter(time_days >= 0)   # drops negative times (date errors) and missing times

# Time in years (and months, for reference)
cohort <- cohort |>
  mutate(time_months = round(time_days / 30.44, 2),
         time_years  = round(time_days / 365.25, 2))

# Patients with time = 0 (<1% of the sample) get 1 day, because Cox / KM
# procedures require positive times.
cohort <- cohort |>
  mutate(time_years = if_else(time_years == 0, 1 / 365.25, time_years),
         time_years = round(time_years, 3))

summary(cohort$time_years)

# --- 4.2 Event indicators ----------------------------------------------------
# death_any    = 1 if died from any cause   (overall survival, OS)
# death_cancer = 1 if died from cancer      (cancer-specific survival, CSS);
#                deaths from other causes are censored.
cohort <- cohort |>
  mutate(
    death_any    = case_when(Status_Vital == "MORTO" ~ 1,
                             Status_Vital == "VIVO"  ~ 0),
    death_cancer = case_when(Status_Vital == "VIVO"                           ~ 0,
                             Status_Vital == "MORTO" & Tipo_do_Obito == "CÂNCER" ~ 1,
                             Status_Vital == "MORTO"                          ~ 0)
  )

table(cohort$death_any, cohort$death_cancer, dnn = c("death_any", "death_cancer"))

# --- 4.3 Administrative censoring at the prediction horizon (10 years) -------
# Follow-up beyond HORIZON is truncated and events after it are censored.
cohort <- cohort |>
  mutate(
    time_years_10 = pmin(time_years, HORIZON),
    death_any     = if_else(time_years > HORIZON, 0, death_any),
    death_cancer  = if_else(time_years > HORIZON, 0, death_cancer)
  )

# Sanity checks
summary(cohort$time_years_10)
table(death_any = cohort$death_any, death_cancer = cohort$death_cancer)

# --- 4.4 Reference categories for the regression models ----------------------
as_ref_factor <- function(x, ref) relevel(factor(x), ref = ref)

cohort <- cohort |>
  mutate(
    age_cat    = as_ref_factor(age_cat,    "Under50"),
    race       = as_ref_factor(race,       "White"),
    sex        = as_ref_factor(sex,        "Female"),
    education  = as_ref_factor(education,  "College"),
    city       = as_ref_factor(city,       "Other"),
    marital    = as_ref_factor(marital,    "Married"),
    morphology = as_ref_factor(morphology, "AC"),
    stage      = as_ref_factor(stage,      "Localized"),
    topography = factor(topography)          # C18 (alphabetical first) is the reference
  )

cat("Final cohort:", nrow(cohort), "patients |",
    sum(cohort$death_any), "deaths (OS) |",
    sum(cohort$death_cancer), "cancer deaths (CSS)\n")

# Covariates used in the bivariate and Cox analyses
COVARS <- c("age_cat", "race", "sex", "education", "city",
            "marital", "morphology", "stage", "topography")


# =============================================================================
# 5. KAPLAN-MEIER: OVERALL vs. CANCER-SPECIFIC SURVIVAL ----
# =============================================================================

# --- 5.1 Overall survival ----------------------------------------------------
km_os <- survfit(Surv(time_years_10, death_any) ~ 1, data = cohort)
km_os                              # n, events, median
summary(km_os, times = 0:10)       # survival probability by year
surv_median(km_os)                 # median overall survival

# --- 5.2 Cancer-specific survival --------------------------------------------
km_css <- survfit(Surv(time_years_10, death_cancer) ~ 1, data = cohort)
km_css
summary(km_css, times = 0:10)
surv_median(km_css)

# --- 5.3 Combined plot -------------------------------------------------------
# ggsurvplot() needs a distinct formula in each fit's call to label the curves.
km_os$call$formula  <- Surv(time_years_10, death_any)    ~ "Overall"
km_css$call$formula <- Surv(time_years_10, death_cancer) ~ "Cancer-specific"

p_km <- ggsurvplot(
  fit = list(Overall = km_os, Cancer_specific = km_css),
  data = cohort,
  combine = TRUE,
  conf.int = TRUE, conf.int.style = "ribbon", conf.int.alpha = 0.2,
  legend.title = "Survival type",
  legend.labs = c("Overall survival", "Cancer-specific survival"),
  legend = "top",
  xlab = "Time (years)", ylab = "Probability of survival",
  palette = c("#0072B2", "#D55E00"),
  break.time.by = 2, ylim = c(0, 1),
  surv.median.line = "hv",
  risk.table = TRUE, risk.table.height = 0.30,
  risk.table.title = "Number at risk",
  risk.table.fontsize = 4.5, risk.table.y.text = FALSE,
  censor = TRUE, censor.shape = 124, censor.size = 3,
  ggtheme = theme_classic(base_size = 16) +
    theme(axis.title   = element_text(face = "bold"),
          axis.text    = element_text(color = "black"),
          legend.title = element_text(face = "bold"))
)
p_km

save_tiff(p_km, "KM_overall_vs_cancer_specific.tiff", width = 12, height = 8, units = "in")


# =============================================================================
# 6. BIVARIATE ANALYSIS: DESCRIPTIVES AND LOG-RANK TESTS ----
# =============================================================================

# For one covariate and one event, print: frequency table, events by category
# (n and row %), and the log-rank test. Returns the log-rank p-value.
describe_and_logrank <- function(data, var, event,
                                 time = "time_years_10") {
  cat("\n=====", var, "|", event, "=====\n")

  # Frequency of each category
  print(adorn_pct_formatting(tabyl(data, !!sym(var))))

  # Events by category (counts with margins, then row percentages)
  tab <- table(data[[var]], data[[event]], dnn = c(var, event))
  print(addmargins(tab))
  print(round(prop.table(tab, margin = 1) * 100, 2))

  # Log-rank test of survival differences between categories
  lr <- survdiff(as.formula(sprintf("Surv(%s, %s) ~ %s", time, event, var)),
                 data = data)
  print(lr)
  pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
}

# --- 6.1 Overall survival ----------------------------------------------------
logrank_os <- map_dbl(COVARS, ~ describe_and_logrank(cohort, .x, "death_any"))

# --- 6.2 Cancer-specific survival --------------------------------------------
logrank_css <- map_dbl(COVARS, ~ describe_and_logrank(cohort, .x, "death_cancer"))

# --- 6.3 Summary table of log-rank p-values ----------------------------------
logrank_table <- tibble(variable = COVARS,
                        p_OS  = signif(logrank_os, 3),
                        p_CSS = signif(logrank_css, 3))
logrank_table
write_csv(logrank_table, file.path(OUT_DIR, "logrank_pvalues.csv"))


# =============================================================================
# 7. UNIVARIABLE COX MODELS ----
# =============================================================================

# Fit one Cox model per covariate; return hazard ratios (HR) with 95% CIs.
univariable_cox <- function(data, vars, event, time = "time_years_10") {
  map_dfr(vars, function(v) {
    fit <- coxph(as.formula(sprintf("Surv(%s, %s) ~ %s", time, event, v)), data = data)
    tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
      mutate(variable = v, .before = 1)
  }) |>
    select(variable, term, HR = estimate, conf.low, conf.high, p.value) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

# --- 7.1 Overall survival ----------------------------------------------------
cox_uni_os <- univariable_cox(cohort, COVARS, "death_any")
print(cox_uni_os, n = Inf)
write_csv(cox_uni_os, file.path(OUT_DIR, "cox_univariable_OS.csv"))

# --- 7.2 Cancer-specific survival --------------------------------------------
cox_uni_css <- univariable_cox(cohort, COVARS, "death_cancer")
print(cox_uni_css, n = Inf)
write_csv(cox_uni_css, file.path(OUT_DIR, "cox_univariable_CSS.csv"))


# =============================================================================
# 8. MULTIVARIABLE COX MODEL: OVERALL SURVIVAL ----
# =============================================================================

# --- 8.1 Fully adjusted model ------------------------------------------------
model_os <- coxph(
  Surv(time_years_10, death_any) ~ age_cat + race + sex + education +
    city + morphology + stage + topography + marital,
  data = cohort
)
summary(model_os)

# Backward/forward stepwise selection by AIC (for reference only)
step(model_os)

# --- 8.2 Adjusted survival curve (at the mean covariate profile) -------------
fit_os <- survfit(model_os)
surv_median(fit_os)
summary(fit_os, times = c(1, 3, 5, 10))

p_os <- ggsurvplot(
  fit_os, data = cohort,
  conf.int = TRUE, conf.int.style = "ribbon", conf.int.alpha = 0.15,
  censor = TRUE, censor.shape = 124, censor.size = 2,
  palette = "#C00000", size = 1.3,
  xlab = "Time since diagnosis (years)", ylab = "Overall survival probability",
  break.time.by = 2, surv.scale = "percent",
  ggtheme = theme_classic(base_size = 15)
)
p_os$plot <- p_os$plot +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "none")
p_os

save_tiff(p_os, "Adjusted_survival_OS.tiff", width = 180, height = 150)

# --- 8.3 Proportional hazards assumption (Schoenfeld residuals) --------------
ph_os <- cox.zph(model_os, transform = "km", global = TRUE)
ph_os
p_ph_os <- ggcoxzph(ph_os)
p_ph_os

tiff(file.path(OUT_DIR, "Schoenfeld_OS.tiff"), width = 4000, height = 2300, res = 300)
print(p_ph_os)
dev.off()


# =============================================================================
# 9. MULTIVARIABLE COX MODEL: CANCER-SPECIFIC SURVIVAL ----
# =============================================================================

# --- 9.1 Fully adjusted model ------------------------------------------------
model_css <- coxph(
  Surv(time_years_10, death_cancer) ~ age_cat + race + sex + education +
    city + morphology + stage + topography + marital,
  data = cohort
)
summary(model_css)

step(model_css)

# --- 9.2 Adjusted survival curve with number at risk -------------------------
fit_css <- survfit(model_css)
surv_median(fit_css)
summary(fit_css, times = c(1, 3, 5, 10))

p_css <- ggsurvplot(
  fit_css, data = cohort,
  color = "#C00000", size = 1.4,
  conf.int = TRUE, conf.int.style = "ribbon", conf.int.alpha = 0.10,
  censor = TRUE, censor.shape = 3, censor.size = 1.8,
  risk.table = TRUE, risk.table.height = 0.22,
  risk.table.y.text = FALSE, risk.table.title = "Number at risk",
  cumevents = FALSE, cumcensor = FALSE,
  surv.scale = "percent", break.time.by = 1,
  xlab = "Time since diagnosis (years)", ylab = "Cancer-specific survival probability",
  ggtheme = theme_classic(base_size = 15)
)
p_css$plot <- p_css$plot +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 12, colour = "black"),
        legend.position = "none")
p_css

save_tiff(p_css, "Adjusted_survival_CSS.tiff", width = 180, height = 150, res = 300)

# --- 9.3 Proportional hazards assumption -------------------------------------
ph_css <- cox.zph(model_css, transform = "km", global = TRUE)
ph_css
p_ph_css <- ggcoxzph(ph_css)
p_ph_css

tiff(file.path(OUT_DIR, "Schoenfeld_CSS.tiff"), width = 4000, height = 2500, res = 300)
print(p_ph_css)
dev.off()


# =============================================================================
# 10. TIME-VARYING EFFECTS (PREPARATION) ----
# =============================================================================
# When the Schoenfeld test shows non-proportional hazards, the effect of a
# covariate can be allowed to change over time (covariate x time interaction).
# This step builds a dummy-variable matrix (sex is omitted here) to be used in
# those interaction models. The interaction models themselves are the next step.

dummies_os <- model.matrix(
  ~ age_cat + race + education + morphology + stage + topography + marital - 1,
  data = cohort
) |>
  as.data.frame()

cohort_dummies <- cbind(cohort, dummies_os)
names(cohort_dummies)

# TODO: fit coxph() with tt() terms or split follow-up with survSplit() using
# the dummy variables of the covariates that violate proportional hazards.
