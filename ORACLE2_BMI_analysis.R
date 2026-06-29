# =============================================================================
# ORACLE2 BMI Analysis – Publication Script
#
# Figures produced:
#   Figure 1  : Figure_3_v3                    (T2 biomarker 3×3 grid by BMI, p-badges)
#   Figure 2  : ACQ_fig1_adj_FVC_imputed       (ACQ means + item heatmap, FVC-adjusted)
#   Figure 3  : LungFunction_BMI_violin_trend  (violin + ANOVA trend p)
#   Figure 4  : RCS_ASAAR_FVC_final            (RCS spline + nadir + strip, FVC-adjusted)
#
# Supplement figures:
#   S1  : BMI ~ FeNO / BEC correlations
#   S2  : BEC density by FeNO stratum  (density plot)
#   S3  : Figure 3 variant – adults only (≥18 years)
#   S5  : Figure 3 variant – continuous BMI only
#   S6  : Delta ACQ vs Normal-weight reference (forest plot)
#   S7  : ACQ by sex × BMI class
#   S8  : Lung function violin (without trend annotation)
#   S9  : RCS ASAAR stratified by T2 group
#   S10 : 3×3 heatmaps (obesity × T2) – rate ratios from imputed data
#   S11 : Splines by obesity × FeNO (BEC modelled continuously)
#
# Table 1: Baseline characteristics by BMI group (4 categories)
#
# Inputs expected in the global environment
#   data1   – non-imputed baseline dataset  (one row per patient)
#   demo1   – baseline dataset with raw (non-log) FeNO/BEC columns:
#               FeNO_notlog, BEC_notlog
#   demo2   – cleaned baseline dataset (demo1 filtered, 9999 → NA, BMI_group4 added)
#   demo    – dataset with time-to-first-attack variables
#   imp_data_ORACLE_final_COMP     – long-format imputed data (.imp column)
#   imp_data_ORACLE_final_COMP_NR  – as above, with linear_predictors column
#   imp_data_ORACLE_final_COMP_upd – as above, additionally containing
#                                    FVC_postBD_PCT_Baseline_IMPUTED
#
# =============================================================================


# ── Libraries ─────────────────────────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS)
library(mice)
library(rms)
library(stringr)
library(tibble)
library(broom)
library(scales)
library(forcats)
library(patchwork)
library(emmeans)
library(ggpubr)
library(ggtext)
library(tableone)
library(gt)
library(gghalves)   # geom_half_violin  (for Figure S8)
library(car)        # vif()


# ── Global variable names (edit here if column names differ) ──────────────────
BMI_VAR   <- "BMI"
ATTACKS_Y <- "Number_severe_asthma_attacks_during_followup"
FU_DAYS   <- "Follow_up_duration_days_notlogged"
SEX_VAR   <- "Gender_0Female_1Male"
AGE_VAR   <- "Age"
BEC_VAR   <- "Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced"
FENO_VAR  <- "FeNO_baseline_ppb"
SEV_VAR   <- "Treatment_step"
FEV1PCT   <- "FEV1_preBD_PCT_Baseline"
PREVATT   <- "Any_severe_attack_previous_12m_0no_1yes"
TRIAL_VAR <- "Enrolled_Trial_name"
ACQ_Y     <- "ACQ_baseline_score_mean"
FVC_PCT   <- "FVC_postBD_PCT_Baseline_IMPUTED"   # used in FVC-adjusted models

OBESE_CUTOFF <- 30

# colour constant used throughout
COL_MAIN <- "#084081"


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Add WHO-6 BMI categories to any data frame
add_bmi_categories <- function(df) {
  df %>%
    mutate(
      bmi_value = .data[[BMI_VAR]],
      bmi_who6 = case_when(
        bmi_value < 18.5                  ~ "Underweight (<18.5)",
        bmi_value >= 18.5 & bmi_value < 25 ~ "Normal (18.5\u2013<25)",
        bmi_value >= 25   & bmi_value < 30 ~ "Overweight (25\u2013<30)",
        bmi_value >= 30   & bmi_value < 35 ~ "Obese I (30\u2013<35)",
        bmi_value >= 35   & bmi_value < 40 ~ "Obese II (35\u2013<40)",
        bmi_value >= 40                    ~ "Obese III (\u226540)",
        TRUE ~ NA_character_
      ),
      bmi_who6 = factor(
        bmi_who6,
        levels = c(
          "Underweight (<18.5)", "Normal (18.5\u2013<25)",
          "Overweight (25\u2013<30)", "Obese I (30\u2013<35)",
          "Obese II (35\u2013<40)", "Obese III (\u226540)"
        )
      )
    )
}

# Pool negative-binomial models across imputations
fit_pool_nb_by_imp <- function(long_df, formula_rhs, exposure_days_var = FU_DAYS) {
  stopifnot(".imp" %in% names(long_df))
  imps <- sort(unique(long_df$.imp))
  fits <- vector("list", length(imps))
  for (k in seq_along(imps)) {
    d <- long_df %>% filter(.imp == imps[k])
    fml <- as.formula(paste0(
      ATTACKS_Y, " ~ ", formula_rhs,
      " + offset(log(", exposure_days_var, "))"
    ))
    fits[[k]] <- MASS::glm.nb(fml, data = d)
  }
  mice::pool(fits)
}

# RCS knot positions (10th / 50th / 90th percentile)
get_rcs3_knots <- function(x) {
  as.numeric(stats::quantile(x, probs = c(0.10, 0.50, 0.90), na.rm = TRUE))
}

# Most-common value (used in prediction grids)
get_mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}

# Fit RCS NB model for one imputation slice (without FVC)
fit_rcs3_nb_oneimp <- function(d) {
  bmi <- d[[BMI_VAR]]
  kn  <- get_rcs3_knots(bmi)
  fml <- as.formula(paste0(
    ATTACKS_Y, " ~ rms::rcs(", BMI_VAR, ", knots = c(",
    paste(sprintf("%.6f", kn), collapse = ","), "))",
    " + factor(", SEX_VAR, ")",
    " + ", AGE_VAR, " + ", BEC_VAR, " + ", FENO_VAR,
    " + ", SEV_VAR, " + ", FEV1PCT, " + ", PREVATT,
    " + factor(", TRIAL_VAR, ")",
    " + offset(log(", FU_DAYS, "))"
  ))
  MASS::glm.nb(fml, data = d)
}

# Fit RCS NB model with imputed FVC adjustment
fit_rcs3_nb_oneimp_fvc <- function(d) {
  bmi <- d[[BMI_VAR]]
  kn  <- get_rcs3_knots(bmi)
  fml <- as.formula(paste0(
    ATTACKS_Y, " ~ rms::rcs(", BMI_VAR, ", knots = c(",
    paste(sprintf("%.6f", kn), collapse = ","), "))",
    " + factor(", SEX_VAR, ")",
    " + ", AGE_VAR, " + ", BEC_VAR, " + ", FENO_VAR,
    " + ", SEV_VAR, " + ", FEV1PCT, " + ", FVC_PCT,
    " + ", PREVATT,
    " + factor(", TRIAL_VAR, ")",
    " + offset(log(", FU_DAYS, "))"
  ))
  MASS::glm.nb(fml, data = d)
}

# Safe prediction: pins TRIAL_VAR to a level the fitted model actually saw
predict_rate_curve_safe <- function(fit, d_ref, bmi_grid) {
  coef_names   <- names(coef(fit))
  trial_prefix <- paste0("factor(", TRIAL_VAR, ")")
  known_trials <- sub(trial_prefix, "",
                      coef_names[startsWith(coef_names, trial_prefix)],
                      fixed = TRUE)
  trial_counts <- sort(table(d_ref[[TRIAL_VAR]]), decreasing = TRUE)
  safe_trial   <- names(trial_counts)[names(trial_counts) %in% known_trials]
  if (length(safe_trial) == 0) safe_trial <- known_trials[1] else safe_trial <- safe_trial[1]

  ref <- d_ref %>%
    summarise(
      across(all_of(AGE_VAR),   ~ mean(.x, na.rm = TRUE)),
      across(all_of(BEC_VAR),   ~ mean(.x, na.rm = TRUE)),
      across(all_of(FENO_VAR),  ~ mean(.x, na.rm = TRUE)),
      across(all_of(SEV_VAR),   ~ median(.x, na.rm = TRUE)),
      across(all_of(FEV1PCT),   ~ mean(.x, na.rm = TRUE)),
      across(all_of(PREVATT),   ~ median(.x, na.rm = TRUE))
    )
  ref[[SEX_VAR]]   <- get_mode(d_ref[[SEX_VAR]])
  ref[[TRIAL_VAR]] <- safe_trial

  newdat             <- ref[rep(1, length(bmi_grid)), , drop = FALSE]
  newdat[[BMI_VAR]]  <- bmi_grid
  newdat[[FU_DAYS]]  <- 365.25

  pr <- predict(fit, newdata = newdat, type = "link", se.fit = TRUE)
  tibble(bmi = bmi_grid,
         rate = exp(pr$fit),
         lo   = exp(pr$fit - 1.96 * pr$se.fit),
         hi   = exp(pr$fit + 1.96 * pr$se.fit))
}

# Safe prediction – with FVC covariate added to reference row
predict_rate_curve_fvc_safe <- function(fit, d_ref, bmi_grid) {
  coef_names   <- names(coef(fit))
  trial_prefix <- paste0("factor(", TRIAL_VAR, ")")
  known_trials <- sub(trial_prefix, "",
                      coef_names[startsWith(coef_names, trial_prefix)],
                      fixed = TRUE)
  trial_counts <- sort(table(d_ref[[TRIAL_VAR]]), decreasing = TRUE)
  safe_trial   <- names(trial_counts)[names(trial_counts) %in% known_trials]
  if (length(safe_trial) == 0) safe_trial <- known_trials[1] else safe_trial <- safe_trial[1]

  ref <- d_ref %>%
    summarise(
      across(all_of(AGE_VAR),   ~ mean(.x, na.rm = TRUE)),
      across(all_of(BEC_VAR),   ~ mean(.x, na.rm = TRUE)),
      across(all_of(FENO_VAR),  ~ mean(.x, na.rm = TRUE)),
      across(all_of(SEV_VAR),   ~ median(.x, na.rm = TRUE)),
      across(all_of(FEV1PCT),   ~ mean(.x, na.rm = TRUE)),
      across(all_of(FVC_PCT),   ~ mean(.x, na.rm = TRUE)),
      across(all_of(PREVATT),   ~ median(.x, na.rm = TRUE))
    )
  ref[[SEX_VAR]]   <- get_mode(d_ref[[SEX_VAR]])
  ref[[TRIAL_VAR]] <- safe_trial

  newdat             <- ref[rep(1, length(bmi_grid)), , drop = FALSE]
  newdat[[BMI_VAR]]  <- bmi_grid
  newdat[[FU_DAYS]]  <- 365.25

  pr <- predict(fit, newdata = newdat, type = "link", se.fit = TRUE)
  tibble(bmi = bmi_grid,
         rate = exp(pr$fit),
         lo   = exp(pr$fit - 1.96 * pr$se.fit),
         hi   = exp(pr$fit + 1.96 * pr$se.fit))
}

# p-value badge helper for Figure 3 family
format_p <- function(p) {
  if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
}

# T2 cell colour (green = T2-low, red = T2-high, yellow = intermediate)
cell_colour <- function(feno, eos) {
  case_when(
    feno == "<20 ppb" & eos == "<0.15"  ~ "#1a9641",
    feno == "\u226535 ppb" & eos == "\u22650.3" ~ "#d7191c",
    TRUE                                ~ "#f4a700"
  )
}

# Figure 3 badge layer (one call per badge colour)
make_badge_layer <- function(pval_df, col) {
  geom_label(
    data          = pval_df %>% filter(badge_col == col),
    aes(x = Eos_cat, y = FeNO_cat, label = pval_lab),
    inherit.aes   = FALSE,
    vjust         = 2.4,
    size          = 2.4,
    fontface      = "bold",
    fill          = "white",
    colour        = col,
    label.size    = 0.4,
    label.r       = unit(0.08, "lines"),
    label.padding = unit(0.12, "lines")
  )
}


# =============================================================================
# TABLE 1  –  Baseline characteristics by BMI group (4-category)
# =============================================================================

# ── Prepare working dataset ──────────────────────────────────────────────────
data2 <- data1 %>%
  filter(!is.na(BMI)) %>%
  mutate(across(everything(), ~ ifelse(. == 9999, NA, .))) %>%
  mutate(
    BMI_group4 = case_when(
      BMI < 18.5 ~ "Underweight",
      BMI < 25   ~ "Normal",
      BMI < 30   ~ "Overweight",
      BMI >= 30  ~ "Obese",
      TRUE       ~ NA_character_
    ),
    BMI_group4 = factor(BMI_group4,
                        levels = c("Underweight", "Normal", "Overweight", "Obese")),
    Smoking = case_when(
      Smoking_0never_1ex_2current == 0 ~ "Never",
      Smoking_0never_1ex_2current == 1 ~ "Ex-smoker",
      Smoking_0never_1ex_2current == 2 ~ "Current",
      TRUE ~ NA_character_
    ),
    Smoking = factor(Smoking, levels = c("Never", "Ex-smoker", "Current")),
    Treatment_step_f = factor(Treatment_step, levels = 1:5,
                              labels = paste("Step", 1:5)),
    Atopy           = ifelse(Atopy_history_0no_1yes_9999notknown == 1, "Yes", "No"),
    AllergyTest     = ifelse(Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown == 1, "Yes", "No"),
    Eczema          = ifelse(Eczema_0no_1yes_9999notknown == 1, "Yes", "No"),
    AllergicRhin    = ifelse(AllergicRhinitis__0no_1yes_9999notknown == 1, "Yes", "No"),
    ChronRhinosinus = ifelse(Chronic_Rhinosinusitis_0no_1yes_9999notknown == 1, "Yes", "No"),
    NasalPolyposis  = ifelse(Nasal_polyposis_0no_1yes_9999notknown == 1, "Yes", "No"),
    PsychDisease    = ifelse(Psychiatric_disease_0no_1yes_9999notknown_NOTIMPUTED == 1, "Yes", "No"),
    SevereAttack_12m = ifelse(Any_severe_attack_previous_12m_0no_1yes == 1, "Yes", "No"),
    mOCS             = ifelse(maintenance_OCS_prescribed__0no_1yes == 1, "Yes", "No"),
    BEC_cellsperµL   = Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced * 1000
  ) %>%
  filter(!is.na(BMI_group4))

vars_continuous_median <- intersect(
  c("Age", "BMI", "ACQ_baseline_score_mean",
    "FEV1_preBD_PCT_Baseline", "FEV1_FVC_ratio",
    "FEV1_PCT_reversibility_postBD", "FeNO_baseline_ppb",
    "BEC_cellsperµL", "Total_IgE"),
  names(data2)
)

vars_categorical <- intersect(
  c("Smoking", "Ethnicity", "Treatment_step_f", "mOCS",
    "Atopy", "AllergyTest", "Eczema", "AllergicRhin",
    "ChronRhinosinus", "NasalPolyposis", "PsychDisease",
    "SevereAttack_12m", "Number_severe_attack_previous_12m_con"),
  names(data2)
)
all_vars <- c(vars_continuous_median, vars_categorical)

t1_bmi4 <- CreateTableOne(
  vars       = all_vars,
  strata     = "BMI_group4",
  data       = data2,
  factorVars = vars_categorical,
  includeNA  = FALSE,
  addOverall = TRUE
)

print_bmi4 <- print(
  t1_bmi4,
  nonnormal     = vars_continuous_median,
  showAllLevels = TRUE,
  contDigits    = 1,
  catDigits     = 1,
  quote         = FALSE,
  noSpaces      = TRUE,
  printToggle   = FALSE,
  test          = TRUE
)

df_table1 <- as.data.frame(print_bmi4) %>%
  rownames_to_column("Variable") %>%
  mutate(Variable = gsub("_", " ", Variable)) %>%
  rename(`p-value` = p, Test = test)

write.csv(df_table1, "Table1_BMI4_groups.csv", row.names = FALSE)

gt_tbl <- df_table1 %>%
  dplyr::select(-any_of("Test")) %>%
  gt() %>%
  tab_header(
    title    = md("**Table 1. Baseline characteristics by BMI group**"),
    subtitle = md("Percentages calculated from participants with available data (non-missing denominator)")
  ) %>%
  tab_spanner(label = "BMI group",
              columns = c("Underweight", "Normal", "Overweight", "Obese")) %>%
  cols_label(Variable = "Characteristic", Overall = "Overall") %>%
  tab_options(table.font.size = 12, table.width = pct(100),
              row.striping.include_table_body = TRUE)

gt_tbl
gt::gtsave(gt_tbl, "Table1_BMI4.html")


# =============================================================================
# DATA PREPARATION FOR FIGURES 1, 3, 4, AND SUPPLEMENTS
# =============================================================================

# ── Imputed dataset: add BMI categories ──────────────────────────────────────
long_df2 <- add_bmi_categories(imp_data_ORACLE_final_COMP) %>%
  mutate(bmi_who6 = relevel(bmi_who6, ref = "Normal (18.5\u2013<25)"))

# Plot-order levels (top-to-bottom in coord_flip figures)
bmi_levels_plot <- c(
  "Obese III (\u226540)", "Obese II (35\u2013<40)", "Obese I (30\u2013<35)",
  "Overweight (25\u2013<30)", "Normal (18.5\u2013<25)", "Underweight (<18.5)"
)

# Imputed dataset WITH FVC (imp_data_ORACLE_final_COMP_upd must exist)
long_df2_fvc <- add_bmi_categories(imp_data_ORACLE_final_COMP_upd) %>%
  mutate(bmi_who6 = relevel(bmi_who6, ref = "Normal (18.5\u2013<25)"))

imps     <- sort(unique(long_df2$.imp))
imps_fvc <- sort(unique(long_df2_fvc$.imp))

# BMI grid for continuous modelling
bmi_grid <- seq(
  min(long_df2[[BMI_VAR]], na.rm = TRUE),
  max(long_df2[[BMI_VAR]], na.rm = TRUE),
  length.out = 200
)
bmi_grid_fvc <- seq(
  min(long_df2_fvc[[BMI_VAR]], na.rm = TRUE),
  max(long_df2_fvc[[BMI_VAR]], na.rm = TRUE),
  length.out = 200
)


# =============================================================================
# FIGURE 2  –  LungFunction_BMI_violin_trend
# Lung function (FEV1%, FVC%, FEV1/FVC) by 4-category BMI;
# full violin + boxplot + ANOVA linear-trend p-value per panel.
# =============================================================================

data1 <- data1 %>%
  mutate(
    bmi_cat = case_when(
      BMI < 18.5               ~ "Underweight",
      BMI >= 18.5 & BMI < 25  ~ "Normal",
      BMI >= 25   & BMI < 30  ~ "Overweight",
      BMI >= 30                ~ "Obese"
    ),
    bmi_cat = factor(bmi_cat,
                     levels = c("Underweight", "Normal", "Overweight", "Obese"))
  )

data_lung <- data1 %>%
  dplyr::select(
    bmi_cat,
    `FEV1 % predicted`  = FEV1_preBD_PCT_Baseline,
    `FVC % predicted`   = FVC_postBD_PCT_Baseline_NOTIMPUTED,
    `FEV1/FVC ratio`    = FEV1_FVC_ratio
  ) %>%
  drop_na(bmi_cat) %>%
  mutate(
    `FEV1/FVC ratio` = ifelse(`FEV1/FVC ratio` <= 1,
                              `FEV1/FVC ratio` * 100,
                              `FEV1/FVC ratio`)
  ) %>%
  pivot_longer(
    cols      = c(`FEV1 % predicted`, `FVC % predicted`, `FEV1/FVC ratio`),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(Metric = factor(Metric,
                         levels = c("FEV1 % predicted", "FVC % predicted", "FEV1/FVC ratio")))

mean_labels <- data_lung %>%
  filter(!is.na(Value)) %>%
  group_by(bmi_cat, Metric) %>%
  summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop")

bmi_order_lf <- c("Underweight", "Normal", "Overweight", "Obese")

trend_pvals <- data_lung %>%
  mutate(bmi_num = as.numeric(factor(bmi_cat, levels = bmi_order_lf))) %>%
  group_by(Metric) %>%
  summarise(
    p_trend = summary(lm(Value ~ bmi_num))$coefficients["bmi_num", "Pr(>|t|)"],
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      p_trend < 0.001 ~ "Trend p<0.001",
      TRUE            ~ sprintf("Trend p=%.3f", p_trend)
    )
  )

trend_annot <- data_lung %>%
  group_by(Metric) %>%
  summarise(y_pos = max(Value, na.rm = TRUE) * 1.02, .groups = "drop") %>%
  left_join(trend_pvals, by = "Metric")

p_lung_v3 <- ggplot(data_lung,
                    aes(x = bmi_cat, y = Value, fill = bmi_cat, colour = bmi_cat)) +
  geom_violin(alpha = 0.45, linewidth = 0.3, trim = TRUE, scale = "width") +
  geom_boxplot(aes(fill = bmi_cat),
               width = 0.12, outlier.alpha = 0.08, outlier.size = 0.6,
               linewidth = 0.35, colour = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.0, colour = "black") +
  geom_label(
    data          = mean_labels,
    aes(x = bmi_cat, y = mean_val, label = sprintf("%.1f", mean_val)),
    inherit.aes   = FALSE,
    vjust         = -0.5, size = 2.5, fontface = "bold",
    colour        = "#1a3a5c", fill = "white",
    label.size    = 0.25, label.r = unit(0.12, "lines"),
    label.padding = unit(0.18, "lines")
  ) +
  geom_text(
    data        = trend_annot,
    aes(x = Inf, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 0, size = 3.0, fontface = "italic", colour = "grey25"
  ) +
  facet_wrap(~ Metric, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "Underweight" = "#dce8f0", "Normal" = "#b4ccde",
    "Overweight"  = "#8cb0cc", "Obese"  = "#6494ba"
  ), guide = "none") +
  scale_colour_manual(values = c(
    "Underweight" = "#dce8f0", "Normal" = "#b4ccde",
    "Overweight"  = "#8cb0cc", "Obese"  = "#6494ba"
  ), guide = "none") +
  scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.10))) +
  labs(x = "", y = "Value (%)") +
  theme_classic(base_size = 10) +
  theme(
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", colour = NA),
    axis.text.x      = element_text(size = 8),
    panel.spacing    = unit(1.0, "lines")
  )

ggsave("LungFunction_BMI_violin_trend.jpg", p_lung_v3,
       width = 7, height = 9, dpi = 600)


# =============================================================================
# FIGURE 1  –  ACQ_fig1_adj_FVC_imputed
# Adjusted mean ACQ (panel A) + item heatmap (panel B), additionally
# adjusted for imputed FVC% predicted.
# =============================================================================

acq_item_vars <- c(
  "ACQ_baseline_score_item1_sleepawakenings_NOTIMPUTED",
  "ACQ_baseline_score_item2_morningsymptoms_NOTIMPUTED",
  "ACQ_baseline_score_item3_activitylimitation_NOTIMPUTED",
  "ACQ_baseline_score_item4_dyspnea_NOTIMPUTED",
  "ACQ_baseline_score_item5_wheezing_NOTIMPUTED"
)

item_labels <- c(
  "ACQ_baseline_score_item1_sleepawakenings_NOTIMPUTED"    = "1: Sleep awakenings",
  "ACQ_baseline_score_item2_morningsymptoms_NOTIMPUTED"    = "2: Morning symptoms",
  "ACQ_baseline_score_item3_activitylimitation_NOTIMPUTED" = "3: Activity limitation",
  "ACQ_baseline_score_item4_dyspnea_NOTIMPUTED"            = "4: Dyspnea",
  "ACQ_baseline_score_item5_wheezing_NOTIMPUTED"           = "5: Wheezing"
)

# Fit ACQ model (FVC-adjusted) per imputation
fits_acq_fvc <- lapply(imps_fvc, function(i) {
  d <- long_df2_fvc %>% filter(.imp == i) %>% droplevels()
  lm(
    as.formula(paste0(
      ACQ_Y, " ~ bmi_who6",
      " + factor(", SEX_VAR, ")",
      " + ", AGE_VAR, " + ", BEC_VAR, " + ", FENO_VAR,
      " + ", SEV_VAR, " + ", FEV1PCT, " + ", FVC_PCT,
      " + ", PREVATT, " + factor(", TRIAL_VAR, ")"
    )),
    data = d
  )
})

# Pool emmeans – total ACQ
emm_list_fvc <- lapply(seq_along(fits_acq_fvc), function(k) {
  out      <- as.data.frame(emmeans(fits_acq_fvc[[k]], ~ bmi_who6))
  out$.imp <- imps_fvc[k]
  out
}) %>% bind_rows()

emm_pooled_fvc <- emm_list_fvc %>%
  group_by(bmi_who6) %>%
  summarise(
    adj_mean = mean(emmean,   na.rm = TRUE),
    lo       = mean(lower.CL, na.rm = TRUE),
    hi       = mean(upper.CL, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    `BMI class` = factor(bmi_who6, levels = bmi_levels_plot),
    label       = sprintf("%.2f [%.2f\u2013%.2f]", adj_mean, lo, hi),
    label_pos   = hi + 0.05
  )

# Pool emmeans – per ACQ item
item_emm_list_fvc <- lapply(seq_along(fits_acq_fvc), function(k) {
  lapply(acq_item_vars, function(itv) {
    fml_item <- as.formula(paste0(
      itv, " ~ bmi_who6",
      " + factor(", SEX_VAR, ")",
      " + ", AGE_VAR, " + ", BEC_VAR, " + ", FENO_VAR,
      " + ", SEV_VAR, " + ", FEV1PCT, " + ", FVC_PCT,
      " + ", PREVATT, " + factor(", TRIAL_VAR, ")"
    ))
    d <- long_df2_fvc %>% filter(.imp == imps_fvc[k]) %>% droplevels()
    fit_item <- lm(fml_item, data = d)
    out <- as.data.frame(emmeans(fit_item, ~ bmi_who6))
    out$item <- itv
    out$.imp  <- imps_fvc[k]
    out
  }) %>% bind_rows()
}) %>% bind_rows()

item_emm_pooled_fvc <- item_emm_list_fvc %>%
  group_by(item, bmi_who6) %>%
  summarise(adj_mean = mean(emmean, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    bmi_who6   = factor(bmi_who6,   levels = rev(bmi_levels_plot)),
    item_label = factor(item_labels[item], levels = rev(item_labels))
  )

# Panel A: adjusted mean ACQ
p_A_fvc <- ggplot(emm_pooled_fvc, aes(x = `BMI class`, y = adj_mean)) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.18, colour = COL_MAIN, linewidth = 0.6) +
  geom_point(size = 2.8, colour = COL_MAIN) +
  geom_text(aes(y = label_pos, label = label),
            hjust = 0, size = 3.0, colour = COL_MAIN) +
  coord_flip(clip = "off") +
  labs(x = "", y = "Adjusted mean ACQ score", title = "A") +
  theme_classic() +
  theme(axis.text.y  = element_text(size = 9, face = "bold"),
        plot.title   = element_text(face = "bold", size = 10),
        plot.margin  = margin(r = 95))

# Panel B: item heatmap
p_B_fvc <- ggplot(item_emm_pooled_fvc,
                  aes(x = bmi_who6, y = item_label, fill = adj_mean)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", adj_mean)),
            size = 2.9, colour = "white", fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = COL_MAIN,
    midpoint = 1.5, limits = c(0, 3.3),
    name = "Adj. mean\nitem score"
  ) +
  scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
  labs(x = "BMI class", y = "", title = "B") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x  = element_text(size = 8, angle = 30, hjust = 1),
    axis.text.y  = element_text(size = 9, face = "bold"),
    plot.title   = element_text(face = "bold", size = 10),
    panel.grid   = element_blank()
  )

p_acq_fig1_fvc <- (p_A_fvc / p_B_fvc) + plot_layout(heights = c(1.1, 0.9))

ggsave("ACQ_fig1_adj_FVC_imputed.jpg", p_acq_fig1_fvc,
       width = 9, height = 8, dpi = 600)


# =============================================================================
# FIGURE 4  –  RCS_ASAAR_FVC_final
# RCS spline (FVC-adjusted), nadir annotation, and categorical-RR strip.
# =============================================================================

# Categorical NB model (FVC-adjusted) for strip RRs
rhs_adj_fvc <- paste(
  "bmi_who6",
  " + factor(", SEX_VAR, ")",
  " + ", AGE_VAR, " + ", BEC_VAR, " + ", FENO_VAR,
  " + ", SEV_VAR, " + ", FEV1PCT, " + ", FVC_PCT,
  " + ", PREVATT, " + factor(", TRIAL_VAR, ")"
)

pooled_cat_who6_fvc <- fit_pool_nb_by_imp(long_df2_fvc, rhs_adj_fvc)

df_forest_fvc <- summary(pooled_cat_who6_fvc, conf.int = TRUE, exponentiate = TRUE) %>%
  as_tibble() %>%
  rename(RR = estimate, CI_lower = conf.low, CI_upper = conf.high) %>%
  filter(str_detect(term, "^bmi_who6")) %>%
  mutate(
    Variable_clean = str_trim(str_remove(term, "^bmi_who6")),
    Variable_clean = factor(Variable_clean, levels = rev(
      levels(long_df2_fvc$bmi_who6)[levels(long_df2_fvc$bmi_who6) != "Normal (18.5\u2013<25)"]
    ))
  )

# RCS fits (FVC-adjusted)
fits_rcs_fvc <- lapply(imps_fvc, function(i) {
  d <- long_df2_fvc %>% filter(.imp == i) %>% droplevels()
  fit_rcs3_nb_oneimp_fvc(d)
})

pred_list_fvc <- lapply(seq_along(fits_rcs_fvc), function(k) {
  d <- long_df2_fvc %>% filter(.imp == imps_fvc[k])
  predict_rate_curve_fvc_safe(fits_rcs_fvc[[k]], d_ref = d,
                              bmi_grid = bmi_grid_fvc) %>%
    mutate(.imp = imps_fvc[k])
})

pred_pool_fvc <- bind_rows(pred_list_fvc) %>%
  group_by(bmi) %>%
  summarise(rate = mean(rate), lo = mean(lo), hi = mean(hi), .groups = "drop")

pred_pool_plot_fvc <- pred_pool_fvc %>% filter(bmi >= 15, bmi <= 50)

nadir_row_fvc  <- pred_pool_plot_fvc %>% slice_min(rate, n = 1)
nadir_bmi_fvc  <- round(nadir_row_fvc$bmi,  1)
nadir_rate_fvc <- round(nadir_row_fvc$rate, 3)

# N per category for strip labels (from non-imputed demo2)
n_bmi6 <- demo2 %>%
  mutate(bmi_who6 = case_when(
    BMI < 18.5                  ~ "Underweight (<18.5)",
    BMI >= 18.5 & BMI < 25     ~ "Normal (18.5\u2013<25)",
    BMI >= 25   & BMI < 30     ~ "Overweight (25\u2013<30)",
    BMI >= 30   & BMI < 35     ~ "Obese I (30\u2013<35)",
    BMI >= 35   & BMI < 40     ~ "Obese II (35\u2013<40)",
    BMI >= 40                   ~ "Obese III (\u226540)"
  )) %>%
  filter(!is.na(bmi_who6)) %>%
  count(bmi_who6) %>%
  deframe()

# Strip data
strip_fvc <- df_forest_fvc %>%
  mutate(Variable_clean = as.character(Variable_clean)) %>%
  left_join(
    tibble(
      Variable_clean = c("Underweight (<18.5)", "Overweight (25\u2013<30)",
                         "Obese I (30\u2013<35)", "Obese II (35\u2013<40)", "Obese III (\u226540)"),
      xmin = c(15, 25, 30, 35, 40),
      xmax = c(18.5, 30, 35, 40, 50)
    ),
    by = "Variable_clean"
  ) %>%
  mutate(
    xmid    = (xmin + xmax) / 2,
    label   = paste0(
      gsub("\\s*\\([^)]*\\)", "", Variable_clean), "\n",
      sprintf("%.2f", RR), "\n",
      sprintf("[%.2f\u2013%.2f]", CI_lower, CI_upper)
    ),
    n_label = paste0("n=", format(n_bmi6[Variable_clean], big.mark = ","))
  )

normal_box <- tibble(
  Variable_clean = "Normal (18.5\u2013<25)",
  xmin = 18.5, xmax = 25, xmid = 21.75,
  label   = "Normal\n(reference)",
  n_label = paste0("n=", format(n_bmi6["Normal (18.5\u2013<25)"], big.mark = ","))
)

y_lo <- 0.5
y_hi <- 1.6
strip_y0 <- y_lo + 0.01
strip_y1 <- strip_y0 + 0.09
strip_yt <- (strip_y0 + strip_y1) / 2
strip_yn <- strip_y0 - 0.025

p_rcs_lim_ASAAR_nad_fvc_v3 <- ggplot(pred_pool_plot_fvc, aes(x = bmi, y = rate)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#2C7FB8", alpha = 0.25) +
  geom_line(colour = "#084081", linewidth = 1.1) +
  annotate("segment",
           x = nadir_bmi_fvc, xend = nadir_bmi_fvc,
           y = nadir_rate_fvc, yend = 0.45,
           linetype = "dashed", colour = "#d7191c", linewidth = 0.5, alpha = 0.6) +
  annotate("point",
           x = nadir_bmi_fvc, y = nadir_rate_fvc,
           colour = "#d7191c", size = 3.5) +
  annotate("label",
           x = nadir_bmi_fvc, y = nadir_rate_fvc,
           label = paste0("BMI = ", nadir_bmi_fvc),
           vjust = -0.95, hjust = 0.15, size = 3.2, fontface = "bold",
           colour = "#d7191c", fill = "white",
           label.size = 0.3, label.r = unit(0.1, "lines"),
           label.padding = unit(0.2, "lines")) +
  # Non-Normal strips
  geom_rect(data = strip_fvc,
            aes(xmin = xmin, xmax = xmax, ymin = strip_y0, ymax = strip_y1),
            inherit.aes = FALSE, fill = "grey92", colour = "grey55", linewidth = 0.4) +
  geom_text(data = strip_fvc,
            aes(x = xmid, y = strip_yt, label = label),
            inherit.aes = FALSE, size = 2.8, colour = "#084081",
            fontface = "bold", lineheight = 0.9) +
  # Normal reference strip
  geom_rect(data = normal_box,
            aes(xmin = xmin, xmax = xmax, ymin = strip_y0, ymax = strip_y1),
            inherit.aes = FALSE, fill = "grey92", colour = "grey55", linewidth = 0.4) +
  geom_text(data = normal_box,
            aes(x = xmid, y = strip_yt, label = label),
            inherit.aes = FALSE, size = 2.8, colour = "#084081",
            fontface = "bold", lineheight = 0.9) +
  # n= labels
  geom_text(data = strip_fvc, aes(x = xmid, y = strip_yn, label = n_label),
            inherit.aes = FALSE, size = 2.6, colour = "grey30") +
  geom_text(data = normal_box, aes(x = xmid, y = strip_yn, label = n_label),
            inherit.aes = FALSE, size = 2.6, colour = "grey30") +
  scale_x_continuous(limits = c(15, 50), breaks = seq(15, 50, 5)) +
  scale_y_continuous(breaks = seq(0.6, 1.6, 0.2)) +
  labs(x = "BMI (kg/m\u00b2)",
       y = "Estimated annual severe asthma attack rate (ASAAR)") +
  coord_cartesian(ylim = c(y_lo, y_hi), clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    axis.title  = element_text(size = 11, face = "bold"),
    axis.text   = element_text(size = 10),
    panel.grid  = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 20, l = 10)
  )

ggsave("RCS_ASAAR_FVC_final.jpg", p_rcs_lim_ASAAR_nad_fvc_v3,
       width = 10, height = 6, dpi = 600)


# =============================================================================
# FIGURE 3  –  Figure_3_v3
# 3×3 T2 biomarker grid by 4-category BMI, p-value badges vs Normal.
# =============================================================================

fig3_data <- demo1 %>%
  transmute(
    BMI       = BMI,
    FeNO_ppb  = FeNO_notlog,
    Eos_x109L = BEC_notlog,
    bmi4 = case_when(
      BMI < 18.5 ~ "Underweight (<18.5)",
      BMI < 25   ~ "Normal (18.5\u2013<25)",
      BMI < 30   ~ "Overweight (25\u2013<30)",
      BMI >= 30  ~ "Obese (\u226530)"
    ),
    FeNO_cat = case_when(
      FeNO_ppb < 20 ~ "<20 ppb",
      FeNO_ppb < 35 ~ "20\u2013<35 ppb",
      TRUE          ~ "\u226535 ppb"
    ),
    Eos_cat = case_when(
      Eos_x109L < 0.15 ~ "<0.15",
      Eos_x109L < 0.30 ~ "0.15\u2013<0.3",
      TRUE             ~ "\u22650.3"
    )
  ) %>%
  filter(!is.na(bmi4), !is.na(FeNO_cat), !is.na(Eos_cat)) %>%
  mutate(
    bmi4 = factor(bmi4, levels = c(
      "Underweight (<18.5)", "Normal (18.5\u2013<25)",
      "Overweight (25\u2013<30)", "Obese (\u226530)"
    )),
    FeNO_cat = factor(FeNO_cat, levels = c("<20 ppb", "20\u2013<35 ppb", "\u226535 ppb")),
    Eos_cat  = factor(Eos_cat,  levels = c("<0.15", "0.15\u2013<0.3", "\u22650.3"))
  )

fig3_plot_data <- fig3_data %>%
  count(bmi4, FeNO_cat, Eos_cat, name = "Count") %>%
  group_by(bmi4) %>%
  mutate(
    Total = sum(Count),
    Pct   = 100 * Count / Total,
    lab   = paste0(Count, "/", Total, "\n(", sprintf("%.1f", Pct), "%)")
  ) %>%
  ungroup()

normal_cells_fig3 <- fig3_plot_data %>%
  filter(bmi4 == "Normal (18.5\u2013<25)") %>%
  dplyr::select(FeNO_cat, Eos_cat, Count_ref = Count, Total_ref = Total)

pval_data_fig3 <- fig3_plot_data %>%
  filter(bmi4 %in% c("Underweight (<18.5)", "Overweight (25\u2013<30)", "Obese (\u226530)")) %>%
  left_join(normal_cells_fig3, by = c("FeNO_cat", "Eos_cat")) %>%
  rowwise() %>%
  mutate(
    pval     = prop.test(x = c(Count, Count_ref), n = c(Total, Total_ref),
                         correct = FALSE)$p.value,
    pval_lab = format_p(pval)
  ) %>%
  ungroup() %>%
  mutate(badge_col = cell_colour(as.character(FeNO_cat), as.character(Eos_cat)))

p_bmi4_final_demo1_v3 <- ggplot(fig3_plot_data,
                                aes(x = Eos_cat, y = FeNO_cat, fill = Pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = lab), fontface = "bold", size = 3.0, vjust = 0.3) +
  make_badge_layer(pval_data_fig3, "#f4a700") +
  make_badge_layer(pval_data_fig3, "#1a9641") +
  make_badge_layer(pval_data_fig3, "#d7191c") +
  facet_wrap(~ bmi4, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#084081", name = "% of group") +
  labs(x = "Blood eosinophils (\u00d710\u2079 cells/L)", y = "FeNO (ppb)") +
  coord_fixed() +
  theme_classic() +
  theme(
    strip.text    = element_text(face = "bold", size = 10),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    axis.title    = element_text(face = "bold", size = 12),
    legend.title  = element_text(face = "bold"),
    plot.margin   = margin(t = 10, r = 10, b = 10, l = 10)
  )

ggsave("Figure_3_v3.jpg", p_bmi4_final_demo1_v3,
       dpi = 600, width = 14, height = 5)


# =============================================================================
# SUPPLEMENTARY FIGURES
# =============================================================================

# ── SUPPLEMENT: Figure S1  –  BMI ~ FeNO / BEC correlations ─────────────────
df_corr <- demo1 %>%
  mutate(
    FeNO_ppb  = FeNO_notlog,
    BEC_x109L = 10^Blood_Eos_Log10baseline_x10_9_cells_per_L_zeroreplaced
  ) %>%
  filter(!is.na(BMI))

p_feno_corr <- ggplot(df_corr %>% filter(!is.na(FeNO_ppb)),
                      aes(x = BMI, y = FeNO_ppb)) +
  geom_point(alpha = 0.25, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  stat_cor(method = "spearman", label.x = 15,
           label.y = max(df_corr$FeNO_ppb, na.rm = TRUE) * 0.95, size = 5) +
  labs(title = "BMI vs FeNO", x = "BMI (kg/m\u00b2)", y = "FeNO (ppb)") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 12))

p_bec_corr <- ggplot(df_corr %>% filter(!is.na(BEC_x109L)),
                     aes(x = BMI, y = BEC_x109L)) +
  geom_point(alpha = 0.25, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  stat_cor(method = "spearman", label.x = 15,
           label.y = max(df_corr$BEC_x109L, na.rm = TRUE) * 0.95, size = 5) +
  scale_y_log10() +
  labs(title = "BMI vs Blood eosinophils",
       x = "BMI (kg/m\u00b2)",
       y = expression("Blood eosinophils ("*10^9*" cells/L)")) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 12))

ggsave("BMI_T2biomarker_correlations.png", p_feno_corr | p_bec_corr,
       width = 11, height = 5, dpi = 600)


# ── SUPPLEMENT: Figure S2  –  BEC density by FeNO stratum (3 BMI groups) ────
df_density <- data1 %>%
  mutate(
    bmi_group3 = case_when(
      BMI < 25               ~ "Lean",
      BMI >= 25 & BMI < 30  ~ "Overweight",
      BMI >= 30              ~ "Obese"
    ),
    bmi_group3 = factor(bmi_group3, levels = c("Lean", "Overweight", "Obese")),
    FeNO_Category = case_when(
      FeNO_baseline_ppb <  log10(20) ~ "<20 ppb",
      FeNO_baseline_ppb >= log10(20) & FeNO_baseline_ppb < log10(35) ~ "20\u2013<35 ppb",
      FeNO_baseline_ppb >= log10(35) ~ "\u226535 ppb"
    ),
    FeNO_Category = factor(FeNO_Category, levels = c("<20 ppb", "20\u2013<35 ppb", "\u226535 ppb")),
    BEC_raw = 10^Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced * 1000
  ) %>%
  filter(!is.na(bmi_group3), !is.na(FeNO_Category), !is.na(BEC_raw), BEC_raw > 0)

p_density_bec <- ggplot(df_density,
                        aes(x = BEC_raw, colour = FeNO_Category, fill = FeNO_Category)) +
  geom_density(alpha = 0.35, linewidth = 0.8, adjust = 1) +
  scale_x_continuous(trans = "log10",
                     breaks = c(50, 100, 150, 300, 500, 1000, 2000),
                     labels = as.character(c(50, 100, 150, 300, 500, 1000, 2000))) +
  scale_colour_manual(values = c("green3", "orange2", "red3"), name = "FeNO (ppb)") +
  scale_fill_manual(values = c("green3", "orange2", "red3"), name = "FeNO (ppb)") +
  facet_wrap(~ bmi_group3, nrow = 1) +
  labs(x = "Blood eosinophils (cells/\u03bcL, log scale)", y = "Density") +
  theme_bw() +
  theme(axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "top")

ggsave("Density_BloodEosinophils_byFeNO_byBMIgroup.png", p_density_bec,
       width = 14, height = 5, dpi = 600)


# ── SUPPLEMENT: Figure S3  –  Figure 3 (adults ≥18 years only) ──────────────
fig3_data_adults <- demo1 %>%
  filter(Age >= 18) %>%
  transmute(
    BMI = BMI, FeNO_ppb = FeNO_notlog, Eos_x109L = BEC_notlog,
    bmi4 = case_when(
      BMI < 18.5 ~ "Underweight (<18.5)", BMI < 25 ~ "Normal (18.5\u2013<25)",
      BMI < 30   ~ "Overweight (25\u2013<30)", BMI >= 30 ~ "Obese (\u226530)"
    ),
    FeNO_cat = case_when(
      FeNO_ppb < 20 ~ "<20 ppb", FeNO_ppb < 35 ~ "20\u2013<35 ppb", TRUE ~ "\u226535 ppb"
    ),
    Eos_cat = case_when(
      Eos_x109L < 0.15 ~ "<0.15", Eos_x109L < 0.30 ~ "0.15\u2013<0.3", TRUE ~ "\u22650.3"
    )
  ) %>%
  filter(!is.na(bmi4), !is.na(FeNO_cat), !is.na(Eos_cat)) %>%
  mutate(
    bmi4     = factor(bmi4, levels = c("Underweight (<18.5)", "Normal (18.5\u2013<25)",
                                       "Overweight (25\u2013<30)", "Obese (\u226530)")),
    FeNO_cat = factor(FeNO_cat, levels = c("<20 ppb", "20\u2013<35 ppb", "\u226535 ppb")),
    Eos_cat  = factor(Eos_cat,  levels = c("<0.15", "0.15\u2013<0.3", "\u22650.3"))
  )

fig3_plot_data_adults <- fig3_data_adults %>%
  count(bmi4, FeNO_cat, Eos_cat, name = "Count") %>%
  group_by(bmi4) %>%
  mutate(Total = sum(Count), Pct = 100 * Count / Total,
         lab = paste0(Count, "/", Total, "\n(", sprintf("%.1f", Pct), "%)")) %>%
  ungroup()

normal_adults <- fig3_plot_data_adults %>%
  filter(bmi4 == "Normal (18.5\u2013<25)") %>%
  dplyr::select(FeNO_cat, Eos_cat, Count_ref = Count, Total_ref = Total)

pval_data_adults <- fig3_plot_data_adults %>%
  filter(!bmi4 %in% "Normal (18.5\u2013<25)") %>%
  left_join(normal_adults, by = c("FeNO_cat", "Eos_cat")) %>%
  rowwise() %>%
  mutate(pval = prop.test(x = c(Count, Count_ref), n = c(Total, Total_ref),
                          correct = FALSE)$p.value,
         pval_lab = format_p(pval)) %>%
  ungroup() %>%
  mutate(badge_col = cell_colour(as.character(FeNO_cat), as.character(Eos_cat)))

p_fig3_adults <- ggplot(fig3_plot_data_adults,
                        aes(x = Eos_cat, y = FeNO_cat, fill = Pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = lab), fontface = "bold", size = 3.0, vjust = 0.3) +
  make_badge_layer(pval_data_adults, "#f4a700") +
  make_badge_layer(pval_data_adults, "#1a9641") +
  make_badge_layer(pval_data_adults, "#d7191c") +
  facet_wrap(~ bmi4, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#084081", name = "% of group") +
  labs(x = "Blood eosinophils (\u00d710\u2079 cells/L)", y = "FeNO (ppb)") +
  coord_fixed() + theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold"))

ggsave("Figure_3_v3_adults.png", p_fig3_adults, dpi = 600, width = 14, height = 5)


# ── SUPPLEMENT: Figure S5  –  Figure 3 (continuous BMI only) ─────────────────
# Uses demo_clean (all ages, continuous BMI values only)
fig3_data_clean <- demo_clean %>%
  transmute(
    BMI = BMI,
    FeNO_ppb  = FeNO_baseline_ppb,
    Eos_x109L = Blood_Eos_baseline_x10_9_cells_per_L,
    bmi4 = case_when(
      BMI < 18.5 ~ "Underweight (<18.5)", BMI < 25 ~ "Normal (18.5\u2013<25)",
      BMI < 30   ~ "Overweight (25\u2013<30)", BMI >= 30 ~ "Obese (\u226530)"
    ),
    FeNO_cat = case_when(
      FeNO_ppb < 20 ~ "<20 ppb", FeNO_ppb < 35 ~ "20\u2013<35 ppb", TRUE ~ "\u226535 ppb"
    ),
    Eos_cat = case_when(
      Eos_x109L < 0.15 ~ "<0.15", Eos_x109L < 0.30 ~ "0.15\u2013<0.3", TRUE ~ "\u22650.3"
    )
  ) %>%
  filter(!is.na(bmi4), !is.na(FeNO_cat), !is.na(Eos_cat)) %>%
  mutate(
    bmi4     = factor(bmi4, levels = c("Underweight (<18.5)", "Normal (18.5\u2013<25)",
                                       "Overweight (25\u2013<30)", "Obese (\u226530)")),
    FeNO_cat = factor(FeNO_cat, levels = c("<20 ppb", "20\u2013<35 ppb", "\u226535 ppb")),
    Eos_cat  = factor(Eos_cat,  levels = c("<0.15", "0.15\u2013<0.3", "\u22650.3"))
  )

fig3_plot_data_clean <- fig3_data_clean %>%
  count(bmi4, FeNO_cat, Eos_cat, name = "Count") %>%
  group_by(bmi4) %>%
  mutate(Total = sum(Count), Pct = 100 * Count / Total,
         lab = paste0(Count, "/", Total, "\n(", sprintf("%.1f", Pct), "%)")) %>%
  ungroup()

normal_clean <- fig3_plot_data_clean %>%
  filter(bmi4 == "Normal (18.5\u2013<25)") %>%
  dplyr::select(FeNO_cat, Eos_cat, Count_ref = Count, Total_ref = Total)

pval_data_clean <- fig3_plot_data_clean %>%
  filter(!bmi4 %in% "Normal (18.5\u2013<25)") %>%
  left_join(normal_clean, by = c("FeNO_cat", "Eos_cat")) %>%
  rowwise() %>%
  mutate(pval = prop.test(x = c(Count, Count_ref), n = c(Total, Total_ref),
                          correct = FALSE)$p.value,
         pval_lab = format_p(pval)) %>%
  ungroup() %>%
  mutate(badge_col = cell_colour(as.character(FeNO_cat), as.character(Eos_cat)))

p_fig3_clean <- ggplot(fig3_plot_data_clean,
                       aes(x = Eos_cat, y = FeNO_cat, fill = Pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = lab), fontface = "bold", size = 3.0, vjust = 0.3) +
  make_badge_layer(pval_data_clean, "#f4a700") +
  make_badge_layer(pval_data_clean, "#1a9641") +
  make_badge_layer(pval_data_clean, "#d7191c") +
  facet_wrap(~ bmi4, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#084081", name = "% of group") +
  labs(x = "Blood eosinophils (\u00d710\u2079 cells/L)", y = "FeNO (ppb)") +
  coord_fixed() + theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold"))

ggsave("Figure_3_v3_clean_allages.png", p_fig3_clean, dpi = 600, width = 14, height = 5)


# ── SUPPLEMENT: Figure S6  –  Delta ACQ vs Normal BMI ───────────────────────
# Requires fits_acq (ACQ models without FVC, pooled via pooled_acq)

fits_acq <- lapply(imps, function(i) {
  d <- long_df2 %>% filter(.imp == i)
  lm(
    as.formula(paste0(
      ACQ_Y, " ~ bmi_who6 + factor(", SEX_VAR, ") + ", AGE_VAR, " + ", BEC_VAR,
      " + ", FENO_VAR, " + ", SEV_VAR, " + ", FEV1PCT, " + ", PREVATT,
      " + factor(", TRIAL_VAR, ")"
    )),
    data = d
  )
})
pooled_acq <- mice::pool(fits_acq)

reg_bmi <- summary(pooled_acq, conf.int = TRUE) %>%
  as_tibble() %>%
  filter(grepl("^bmi_who6", term)) %>%
  mutate(
    BMI_class = factor(
      trimws(sub("^bmi_who6", "", term)),
      levels = rev(bmi_levels_plot[bmi_levels_plot != "Normal (18.5\u2013<25)"])
    ),
    label_pos = pmax(conf.high, 0) + 0.04,
    label     = sprintf("%.2f [%.2f\u2013%.2f]", estimate, conf.low, conf.high),
    sig       = p.value < 0.05
  )

p_B_deltaACQ <- ggplot(reg_bmi, aes(x = BMI_class, y = estimate)) +
  geom_hline(yintercept = 0,    linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept =  0.5, linetype = "dotted", colour = "firebrick", linewidth = 0.7) +
  geom_hline(yintercept = -0.5, linetype = "dotted", colour = "firebrick", linewidth = 0.7) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.18, colour = COL_MAIN, linewidth = 0.6) +
  geom_point(aes(shape = sig), size = 2.8, colour = COL_MAIN) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  geom_text(aes(y = label_pos, label = label), hjust = 0, size = 3.0, colour = COL_MAIN) +
  annotate("text", x = 0.55, y =  0.5, label = "MCID +0.5",
           hjust = 0.5, vjust = -0.2, size = 2.7, colour = "grey40", fontface = "italic") +
  annotate("text", x = 0.55, y = -0.5, label = "MCID \u22120.5",
           hjust = 0.5, vjust = -0.2, size = 2.7, colour = "grey40", fontface = "italic") +
  coord_flip(clip = "off") +
  labs(x = "", y = "\u0394 ACQ vs normal BMI (adjusted)",
       caption = "Filled = p<0.05; open = n.s.") +
  theme_classic() +
  theme(axis.text.y  = element_text(size = 9, face = "bold"),
        plot.caption = element_text(size = 7, colour = "grey50"),
        plot.margin  = margin(r = 95))

ggsave("ACQ_fig2_deltaACQ.png", p_B_deltaACQ, width = 9, height = 3.5, dpi = 600)


# ── SUPPLEMENT: Figure S7  –  ACQ stratified by sex ─────────────────────────
sex_levels <- list(Female = 0, Male = 1)
COL_FEMALE <- "#B2182B"
COL_MALE   <- "#2166AC"

fits_acq_sex <- lapply(imps, function(i) {
  d <- long_df2 %>% filter(.imp == i)
  lapply(names(sex_levels), function(sx) {
    ds <- d %>% filter(.data[[SEX_VAR]] == sex_levels[[sx]]) %>% droplevels()
    lm(
      as.formula(paste0(
        ACQ_Y, " ~ bmi_who6 + ", AGE_VAR, " + ", BEC_VAR, " + ", FENO_VAR,
        " + ", SEV_VAR, " + ", FEV1PCT, " + ", PREVATT,
        " + factor(", TRIAL_VAR, ")"
      )),
      data = ds
    )
  }) %>% setNames(names(sex_levels))
})

emm_sex_list <- lapply(seq_along(fits_acq_sex), function(k) {
  lapply(names(sex_levels), function(sx) {
    out <- as.data.frame(emmeans(fits_acq_sex[[k]][[sx]], ~ bmi_who6))
    out$sex  <- sx
    out$.imp <- imps[k]
    out
  }) %>% bind_rows()
}) %>% bind_rows()

emm_pooled_sex <- emm_sex_list %>%
  group_by(sex, bmi_who6) %>%
  summarise(adj_mean = mean(emmean, na.rm = TRUE),
            lo = mean(lower.CL, na.rm = TRUE),
            hi = mean(upper.CL, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    `BMI class` = factor(bmi_who6, levels = rev(c(
      "Underweight (<18.5)", "Normal (18.5\u2013<25)", "Overweight (25\u2013<30)",
      "Obese I (30\u2013<35)", "Obese II (35\u2013<40)", "Obese III (\u226540)"
    ))),
    label     = sprintf("%.2f [%.2f\u2013%.2f]", adj_mean, lo, hi),
    label_pos = hi + 0.05,
    sex       = factor(sex, levels = c("Female", "Male"))
  )

p_acq_sex <- ggplot(emm_pooled_sex,
                    aes(x = `BMI class`, y = adj_mean, colour = sex, shape = sex)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.6,
                position = position_dodge(width = 0.55)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.55)) +
  geom_text(aes(y = label_pos, label = label), hjust = 0, size = 2.8,
            position = position_dodge(width = 0.55), show.legend = FALSE) +
  scale_colour_manual(name = "Sex", values = c(Female = COL_FEMALE, Male = COL_MALE)) +
  scale_shape_manual(name = "Sex", values = c(Female = 16, Male = 16)) +
  coord_flip(clip = "off") +
  labs(x = "", y = "Adjusted mean ACQ score") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9, face = "bold"),
        plot.margin = margin(r = 115))

ggsave("ACQ_sex_stratified.png", p_acq_sex, width = 10, height = 5, dpi = 600)


# ── SUPPLEMENT: Figure S8  –  Lung function violin (no trend annotation) ─────
p_lung_v2 <- ggplot(data_lung,
                    aes(x = bmi_cat, y = Value, fill = bmi_cat, colour = bmi_cat)) +
  geom_violin(alpha = 0.45, linewidth = 0.3, trim = TRUE, scale = "width") +
  geom_boxplot(aes(fill = bmi_cat), width = 0.12, outlier.alpha = 0.08,
               outlier.size = 0.6, linewidth = 0.35, colour = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.0, colour = "black") +
  geom_label(data = mean_labels, aes(x = bmi_cat, y = mean_val,
                                     label = sprintf("%.1f", mean_val)),
             inherit.aes = FALSE, vjust = -0.5, size = 2.5, fontface = "bold",
             colour = "#1a3a5c", fill = "white",
             label.size = 0.25, label.r = unit(0.12, "lines"),
             label.padding = unit(0.18, "lines")) +
  stat_compare_means(
    comparisons = list(c("Normal", "Underweight"), c("Normal", "Overweight"),
                       c("Normal", "Obese")),
    method = "t.test", p.adjust.method = "BH",
    label = "p.signif", tip.length = 0.01, size = 2.8, colour = "grey30"
  ) +
  facet_wrap(~ Metric, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "Underweight" = "#dce8f0", "Normal" = "#b4ccde",
    "Overweight"  = "#8cb0cc", "Obese"  = "#6494ba"
  ), guide = "none") +
  scale_colour_manual(values = c(
    "Underweight" = "#dce8f0", "Normal" = "#b4ccde",
    "Overweight"  = "#8cb0cc", "Obese"  = "#6494ba"
  ), guide = "none") +
  scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  labs(x = "", y = "Value (%)",
       caption = paste(
         "FEV1 pre-BD; FVC post-BD (pre-BD unavailable); FEV1/FVC pre-BD",
         "Adjacent-group comparisons: Benjamini-Hochberg correction;",
         "* p<0.05  ** p<0.01  *** p<0.001",
         sep = "\n"
       )) +
  theme_classic(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey95", colour = NA),
        axis.text.x = element_text(size = 8),
        plot.caption = element_text(size = 7, colour = "grey50"),
        panel.spacing = unit(1.0, "lines"))

ggsave("LungFunction_BMI_violin.png", p_lung_v2, width = 7, height = 9, dpi = 600)


# ── SUPPLEMENT: Figure S9  –  RCS ASAAR by T2 group ─────────────────────────
LOG_FENO_20 <- log10(20); LOG_FENO_35 <- log10(35)
LOG_BEC_015 <- log10(0.15); LOG_BEC_030 <- log10(0.30)

long_df2_t2 <- long_df2 %>%
  mutate(
    t2_group = case_when(
      .data[[FENO_VAR]] <  LOG_FENO_20 & .data[[BEC_VAR]] <  LOG_BEC_015 ~ "T2-low",
      .data[[FENO_VAR]] >= LOG_FENO_35 & .data[[BEC_VAR]] >= LOG_BEC_030 ~ "T2-high",
      !is.na(.data[[FENO_VAR]]) & !is.na(.data[[BEC_VAR]])               ~ "T2-intermediate",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(t2_group))

run_rcs_for_group <- function(df_subset, bmi_range = c(15, 50)) {
  df_subset <- droplevels(df_subset)
  imps_sub  <- sort(unique(df_subset$.imp))
  fits <- lapply(imps_sub, function(i) {
    d <- df_subset %>% filter(.imp == i)
    fit_rcs3_nb_oneimp(d)
  })
  bmi_grid_sub <- seq(
    max(bmi_range[1], min(df_subset[[BMI_VAR]], na.rm = TRUE)),
    min(bmi_range[2], max(df_subset[[BMI_VAR]], na.rm = TRUE)),
    length.out = 200
  )
  bind_rows(lapply(seq_along(fits), function(k) {
    d <- df_subset %>% filter(.imp == imps_sub[k])
    predict_rate_curve_safe(fits[[k]], d_ref = d, bmi_grid = bmi_grid_sub) %>%
      mutate(.imp = imps_sub[k])
  })) %>%
    group_by(bmi) %>%
    summarise(rate = mean(rate), lo = mean(lo), hi = mean(hi), .groups = "drop")
}

message("Fitting T2-low ...")
df_low  <- long_df2_t2 %>% filter(t2_group == "T2-low")
pred_t2low  <- run_rcs_for_group(df_low);  n_t2low  <- df_low  %>% filter(.imp == 1) %>% nrow()

message("Fitting T2-intermediate ...")
df_int  <- long_df2_t2 %>% filter(t2_group == "T2-intermediate")
pred_t2int  <- run_rcs_for_group(df_int);  n_t2int  <- df_int  %>% filter(.imp == 1) %>% nrow()

message("Fitting T2-high ...")
df_high <- long_df2_t2 %>% filter(t2_group == "T2-high")
pred_t2high <- run_rcs_for_group(df_high); n_t2high <- df_high %>% filter(.imp == 1) %>% nrow()

all_preds_t2   <- bind_rows(pred_t2low, pred_t2int, pred_t2high)
y_shared_min   <- floor(min(all_preds_t2$lo, na.rm = TRUE) * 10) / 10
y_shared_max   <- 2.5
strip_y0_sh    <- y_shared_min + 0.004 * (y_shared_max - y_shared_min)
strip_y1_sh    <- y_shared_min + 0.12  * (y_shared_max - y_shared_min)
strip_yt_sh    <- (strip_y0_sh + strip_y1_sh) / 2

make_spline_simple <- function(pred_df, title_label, n,
                               line_col, ribbon_col, show_y = FALSE) {
  ggplot(pred_df, aes(x = bmi, y = rate)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = ribbon_col, alpha = 0.25) +
    geom_line(colour = line_col, linewidth = 1.1) +
    scale_x_continuous(limits = c(15, 50), breaks = seq(15, 50, 5)) +
    coord_cartesian(ylim = c(y_shared_min, y_shared_max), clip = "off") +
    labs(title = paste0(title_label, "  (n=", n, ")"),
         x = "BMI (kg/m\u00b2)",
         y = if (show_y) "Estimated annual severe asthma\nattack rate (ASAAR)" else "") +
    theme_classic(base_size = 13) +
    theme(plot.title  = element_text(face = "bold", size = 13),
          axis.title  = element_text(size = 13, face = "bold"),
          axis.text   = element_text(size = 12),
          panel.grid  = element_blank(),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
}

p_t2low_s  <- make_spline_simple(pred_t2low,  "T2-low  (FeNO <20 ppb & BEC <0.15)",  n_t2low,
                                 "#1a9641", "#78c679", show_y = TRUE)
p_t2int_s  <- make_spline_simple(pred_t2int,  "T2-intermediate  (all others)",        n_t2int,
                                 "#f4a700", "#fed976")
p_t2high_s <- make_spline_simple(pred_t2high, "T2-high  (FeNO \u226535 ppb & BEC \u22650.30)", n_t2high,
                                 "#d7191c", "#fc8d59")

ggsave("RCS_ASAAR_by_T2group_simple.png",
       p_t2low_s | p_t2int_s | p_t2high_s,
       width = 16, height = 5, dpi = 600)


# ── SUPPLEMENT: Figure S10  –  3×3 heatmaps by obesity × T2 (imputed RRs) ───
# Uses imp_data_ORACLE_final_COMP_NR (contains linear_predictors)

imp_data_ORACLE_final_COMP_NR_obsex <- imp_data_ORACLE_final_COMP_NR %>%
  mutate(
    obesity = if_else(.data[[BMI_VAR]] >= OBESE_CUTOFF, "Obese", "Non-obese"),
    obesity = factor(obesity, levels = c("Non-obese", "Obese"))
  )

demo_type2_obsex <- imp_data_ORACLE_final_COMP_NR_obsex %>%
  filter(!is.na(FeNO_baseline_ppb), !is.na(Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced)) %>%
  mutate(
    FeNO_raw = 10^FeNO_baseline_ppb,
    Eos_raw  = 10^Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced,
    FeNO_Category_20_35 = case_when(
      FeNO_raw < 20  ~ "<20 ppb",
      FeNO_raw < 35  ~ "20-<35 ppb",
      TRUE           ~ "\u226535 ppb"
    ),
    Eos_Category_015_03 = case_when(
      Eos_raw < 0.15 ~ "<0.15",
      Eos_raw < 0.30 ~ "0.15-<0.3",
      TRUE           ~ "\u22650.3"
    )
  ) %>%
  mutate(
    FeNO_Category_20_35 = factor(FeNO_Category_20_35, levels = c("<20 ppb", "20-<35 ppb", "\u226535 ppb")),
    Eos_Category_015_03 = factor(Eos_Category_015_03, levels = c("<0.15", "0.15-<0.3", "\u22650.3")),
    sex = if_else(Gender_0Female_1Male == 1, "Male", "Female"),
    sex = factor(sex, levels = c("Female", "Male"))
  )

# T2-cell binary flags
cell_flags <- list(
  FeNO_low_Eos_low   = c("<20 ppb",   "<0.15"),
  FeNO_low_Eos_mid   = c("<20 ppb",   "0.15-<0.3"),
  FeNO_low_Eos_high  = c("<20 ppb",   "\u22650.3"),
  FeNO_mid_Eos_low   = c("20-<35 ppb","<0.15"),
  FeNO_mid_Eos_mid   = c("20-<35 ppb","0.15-<0.3"),
  FeNO_mid_Eos_high  = c("20-<35 ppb","\u22650.3"),
  FeNO_high_Eos_low  = c("\u226535 ppb","<0.15"),
  FeNO_high_Eos_mid  = c("\u226535 ppb","0.15-<0.3"),
  FeNO_high_Eos_high = c("\u226535 ppb","\u22650.3")
)

obesity_totals_obsex <- demo_type2_obsex %>%
  group_by(obesity) %>% summarise(Total = n(), .groups = "drop")

feno_eos_prev_obsex <- demo_type2_obsex %>%
  group_by(obesity, FeNO_Category_20_35, Eos_Category_015_03) %>%
  summarise(Count = n(), .groups = "drop") %>%
  left_join(obesity_totals_obsex, by = "obesity") %>%
  mutate(Percentage = 100 * Count / Total)

# Fit NB models per T2 cell × obesity × imputation
results_list_obsex <- list()
for (ob in c("Non-obese", "Obese")) {
  for (cat in names(cell_flags)) {
    imp_data_ORACLE_final_COMP_NR_obsex <- imp_data_ORACLE_final_COMP_NR_obsex %>%
      mutate(!!cat := as.factor(
        as.integer(FeNO_Category_20_35 == cell_flags[[cat]][1] &
                     Eos_Category_015_03 == cell_flags[[cat]][2])
      ))
    res_comb <- lapply(1:10, function(i) {
      d <- subset(imp_data_ORACLE_final_COMP_NR_obsex, obesity == ob & .imp == i)
      MASS::glm.nb(
        Number_severe_asthma_attacks_during_followup ~ get(cat) +
          offset(d$linear_predictors),
        data = d
      )
    })
    results_list_obsex[[paste(ob, cat, sep = "_")]] <-
      summary(mice::pool(res_comb), conf.int = TRUE, exp = TRUE)
  }
}

extract_rr <- function(res_list, obesity_label) {
  bind_rows(lapply(names(res_list), function(nm) {
    if (!startsWith(nm, paste0(obesity_label, "_"))) return(NULL)
    cat  <- sub(paste0("^", obesity_label, "_"), "", nm)
    rows <- res_list[[nm]]
    tr   <- which(grepl("get\\(cat\\)", rows$term))[1]
    if (is.na(tr)) return(NULL)
    data.frame(
      obesity             = obesity_label,
      FeNO_Category_20_35 = cell_flags[[cat]][1],
      Eos_Category_015_03 = cell_flags[[cat]][2],
      Estimate            = rows$estimate[tr],
      Lower_CI            = rows$conf.low[tr],
      Upper_CI            = rows$conf.high[tr]
    )
  }))
}

results_all_obsex <- bind_rows(
  extract_rr(results_list_obsex, "Non-obese"),
  extract_rr(results_list_obsex, "Obese")
) %>%
  mutate(
    FeNO_Category_20_35 = factor(FeNO_Category_20_35, levels = c("<20 ppb","20-<35 ppb","\u226535 ppb")),
    Eos_Category_015_03 = factor(Eos_Category_015_03, levels = c("<0.15","0.15-<0.3","\u22650.3"))
  )

feno_eos_prev_obsex_rr <- feno_eos_prev_obsex %>%
  left_join(results_all_obsex, by = c("obesity","FeNO_Category_20_35","Eos_Category_015_03"))

mid_rr <- (min(feno_eos_prev_obsex_rr$Estimate, na.rm=TRUE) +
             max(feno_eos_prev_obsex_rr$Estimate, na.rm=TRUE)) / 2

plot_heatmap_obsex <- function(data, title_txt) {
  ggplot(data, aes(x = Eos_Category_015_03, y = FeNO_Category_20_35, fill = Estimate)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "lightyellow", mid = "orange", high = "firebrick",
                         midpoint = mid_rr,
                         limits   = c(min(feno_eos_prev_obsex_rr$Estimate, na.rm=TRUE),
                                      max(feno_eos_prev_obsex_rr$Estimate, na.rm=TRUE)),
                         name = "Rate ratio") +
    ggtext::geom_richtext(
      aes(label = paste0("<b><span style='font-size:12pt;'>", round(Estimate, 2), "</span></b><br>",
                         round(Lower_CI, 2), "\u2013", round(Upper_CI, 2), "<br>",
                         Count, "/", Total, " (", round(Percentage, 1), "%)")),
      fill = NA, label.color = NA, size = 3
    ) +
    labs(title = title_txt, x = "Blood eosinophils (\u00d710\u2079 cells/L)", y = "FeNO (ppb)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(size = 14, face = "bold")) +
    coord_fixed()
}

plot_4panel_obsex <-
  (plot_heatmap_obsex(feno_eos_prev_obsex_rr %>% filter(obesity == "Non-obese"), "Non-obese") |
   plot_heatmap_obsex(feno_eos_prev_obsex_rr %>% filter(obesity == "Obese"),     "Obese")) +
  plot_layout(guides = "collect")

ggsave("FENO_BEC_3x3_obesityonly_RR_imputed.png", plot_4panel_obsex,
       width = 12, height = 6, dpi = 600)


# ── SUPPLEMENT: Figure S11  –  Splines by FeNO × obesity ─────────────────────
df_obonly_25 <- imp_data_ORACLE_final_COMP %>%
  mutate(
    FeNO_actual = 10^FeNO_baseline_ppb,
    FeNO_category_25 = factor(if_else(FeNO_actual < 25, "<25", "\u226525"),
                              levels = c("<25", "\u226525")),
    obesity = factor(if_else(.data[[BMI_VAR]] >= OBESE_CUTOFF, "Obese", "Non-obese"),
                     levels = c("Non-obese", "Obese"))
  )

fit_and_predict_bec_by_feno_obonly_25 <- function(dat) {
  m <- MASS::glm.nb(
    Number_severe_asthma_attacks_during_followup ~
      rms::rcs(Blood_Eos_baseline_x10_9_cells_per_L_NOTIMPUTED, 4) * FeNO_category_25 +
      offset(log(Follow_up_duration_days_notlogged)) +
      ACQ_baseline_score_mean + Any_severe_attack_previous_12m_0no_1yes +
      FEV1_preBD_PCT_Baseline + Treatment_step + as.factor(Enrolled_Trial_name),
    data = dat
  )
  pg <- expand.grid(
    Blood_Eos_baseline_x10_9_cells_per_L_NOTIMPUTED = seq(0.1, 1.5, length.out = 100),
    FeNO_category_25 = levels(dat$FeNO_category_25)
  )
  pg$ACQ_baseline_score_mean              <- mean(dat$ACQ_baseline_score_mean,              na.rm = TRUE)
  pg$Any_severe_attack_previous_12m_0no_1yes <- median(dat$Any_severe_attack_previous_12m_0no_1yes, na.rm = TRUE)
  pg$FEV1_preBD_PCT_Baseline              <- mean(dat$FEV1_preBD_PCT_Baseline,              na.rm = TRUE)
  pg$Treatment_step                       <- median(dat$Treatment_step,                      na.rm = TRUE)
  pg$Enrolled_Trial_name                  <- names(which.max(table(dat$Enrolled_Trial_name)))
  pg$Follow_up_duration_days_notlogged    <- 365.25
  pr <- predict(m, newdata = pg, type = "link", se.fit = TRUE)
  pg %>% mutate(predicted = exp(pr$fit),
                lower     = exp(pr$fit - 1.96 * pr$se.fit),
                upper     = exp(pr$fit + 1.96 * pr$se.fit))
}

pred_all_obonly_25 <- df_obonly_25 %>%
  filter(!is.na(obesity), !is.na(FeNO_category_25)) %>%
  group_split(obesity) %>%
  setNames(levels(df_obonly_25$obesity)) %>%
  purrr::map(~ {
    dat <- .x
    pg  <- fit_and_predict_bec_by_feno_obonly_25(dat)
    pg$obesity <- unique(dat$obesity)
    pg
  }) %>%
  bind_rows()

p_spline_obonly_25 <- ggplot(
  pred_all_obonly_25,
  aes(x = Blood_Eos_baseline_x10_9_cells_per_L_NOTIMPUTED,
      y = predicted, colour = FeNO_category_25, fill = FeNO_category_25)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linewidth = 0) +
  geom_line(linewidth = 1) +
  scale_x_continuous(trans = "log", limits = c(0.1, 1.5),
                     breaks = c(0.1, 0.15, 0.3, 0.6, 1, 1.5)) +
  scale_colour_manual(values = c("green3", "red3"), name = "FeNO (ppb)") +
  scale_fill_manual(values = c("green3", "red3"), name = "FeNO (ppb)") +
  facet_wrap(~ obesity, ncol = 2) +
  coord_cartesian(ylim = c(0, 3)) +
  labs(x = "Blood eosinophils (\u00d710\u2079 cells/L)",
       y = "Estimated annual rate of severe asthma attacks") +
  theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "top")

ggsave("ObesityOnly_splines_BEC_by_FeNOcut25.png", p_spline_obonly_25,
       width = 12, height = 6, dpi = 600)
