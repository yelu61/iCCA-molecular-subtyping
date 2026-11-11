#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE, warn = 1)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "6-ClinicalSurvival/publication/tables")
figure_dir <- file.path(root, "6-ClinicalSurvival/publication/figures")
report_dir <- file.path(root, "6-ClinicalSurvival/publication/reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

required_tables <- c(
  "step6_survival_model_contract.csv", "step6_survival_population_by_study_and_state.csv",
  "step6_clinical_variable_missingness.csv", "step6_clinical_baseline_by_main_state.csv",
  "step6_clinical_association_tests.csv", "step6_km_logrank_tests.csv",
  "step6_km_median_survival.csv", "step6_cox_main_states_core.csv",
  "step6_cox_main_states_pairwise.csv", "step6_cox_main_states_extended_sensitivity.csv",
  "step6_cox_all_states_exploratory.csv", "step6_cox_proportional_hazards_diagnostics.csv",
  "step6_cox_time_varying_subtype_effects.csv", "step6_c3_vs_c1_leave_one_study_out.csv"
)
required_figures <- c(
  "step6_km_main_states.pdf", "step6_km_main_states.png",
  "step6_cox_main_states_forest.pdf", "step6_cox_main_states_forest.png",
  "step6_cox_time_varying_subtype_effects.pdf", "step6_cox_time_varying_subtype_effects.png",
  "step6_c3_vs_c1_leave_one_study_out.pdf", "step6_c3_vs_c1_leave_one_study_out.png"
)
nonempty <- function(path) file.exists(path) && file.info(path)$size > 0

contract <- fread(file.path(table_dir, "step6_survival_model_contract.csv"))
core <- fread(file.path(table_dir, "step6_cox_main_states_core.csv"))
pairwise <- fread(file.path(table_dir, "step6_cox_main_states_pairwise.csv"))
extended <- fread(file.path(table_dir, "step6_cox_main_states_extended_sensitivity.csv"))
ph <- fread(file.path(table_dir, "step6_cox_proportional_hazards_diagnostics.csv"))
tv <- fread(file.path(table_dir, "step6_cox_time_varying_subtype_effects.csv"))
km <- fread(file.path(table_dir, "step6_km_logrank_tests.csv"))
medians <- fread(file.path(table_dir, "step6_km_median_survival.csv"))
loso <- fread(file.path(table_dir, "step6_c3_vs_c1_leave_one_study_out.csv"))
baseline <- fread(file.path(table_dir, "step6_clinical_baseline_by_main_state.csv"))

checks <- data.table(
  Item = c("Required clinical source tables", "Required clinical figures", "Core denominator and event count",
           "Study-stratified Cox specification", "Proportional-hazards violation adjudicated",
           "Main and all-state scopes separated", "Observed duct remains observed-only",
           "Canonical subtype naming"),
  Status = c(
    ifelse(all(vapply(file.path(table_dir, required_tables), nonempty, logical(1))), "PASS", "BLOCKED"),
    ifelse(all(vapply(file.path(figure_dir, required_figures), nonempty, logical(1))), "PASS", "BLOCKED"),
    ifelse(contract[Analysis == "Main states" & Scope == "Average-association summary", N] == 331L &&
             contract[Analysis == "Main states" & Scope == "Average-association summary", Events] == 153L, "PASS", "BLOCKED"),
    ifelse(all(grepl("strata\\(Study\\)", contract$Formula)), "PASS", "BLOCKED"),
    ifelse(ph[Model == "Core study-stratified Cox" & Term == "GLOBAL", PValue] < 0.05 &&
             nrow(tv) == 18L && all(c(3, 6, 12, 24, 36, 60) %in% tv$TimeMonths), "PASS", "BLOCKED"),
    ifelse(any(contract$Scope == "Exploratory") && any(contract$Scope == "Non-PH adjudication"), "PASS", "BLOCKED"),
    ifelse(any(baseline$Variable == "Observed duct type") && !any(grepl("pred", baseline$Variable, ignore.case = TRUE)), "PASS", "BLOCKED"),
    ifelse(!any(grepl("C2-Immune|Immune-sensitive|Immune-enriched", unlist(list(core, pairwise, extended, baseline)), ignore.case = TRUE)), "PASS", "BLOCKED")
  ),
  Evidence = c(paste(required_tables, collapse = ";"), paste(required_figures, collapse = ";"),
               "step6_survival_model_contract.csv", "step6_survival_model_contract.csv",
               "step6_cox_proportional_hazards_diagnostics.csv;step6_cox_time_varying_subtype_effects.csv",
               "step6_survival_model_contract.csv", "step6_clinical_baseline_by_main_state.csv",
               "publication-facing Step6 tables")
)
fwrite(checks, file.path(table_dir, "step6_clinical_survival_release_gate_status.csv"))

get_pair <- function(target, reference) pairwise[Target == target & Reference == reference]
get_tv <- function(contrast, time) tv[Contrast == contrast & TimeMonths == time]
c3c1 <- get_pair("C3", "C1")
c3c2 <- get_pair("C3", "C2")
c2c1 <- get_pair("C2", "C1")
c3_6 <- get_tv("C3 vs C1", 6)
c3_12 <- get_tv("C3 vs C1", 12)
c3_24 <- get_tv("C3 vs C1", 24)

methods <- c(
  "# Clinical and survival methods-ready evidence", "",
  "Overall-survival analyses used the canonical subtype master table and treated each tumor as one analytical unit. Valid survival required a positive overall-survival time and a binary vital-status value. Principal-state analyses were restricted to C1, C2, and C3. The headline complete-case population contained 331 tumors and 153 deaths from 2019Ahn, 2022Dong, 2023Deng, and TCGA-CHOL; cohorts without survival endpoints were excluded only from survival inference.", "",
  "Kaplan-Meier curves were displayed for descriptive available-case survival, while the global comparison used a study-stratified log-rank test. The core Cox model included subtype, continuous age, sex, and stage group (I/II versus III/IV), with study-specific baseline hazards through strata(Study). C1, female sex, and stage I/II were references. The extended complete-case sensitivity additionally included tumor size, HBV status, and vascular invasion. Six-state models were explicitly exploratory.", "",
  "Proportional hazards were evaluated with scaled Schoenfeld residuals using cox.zph. Because subtype and age violated the proportional-hazards assumption in the core model, constant Cox hazard ratios were retained only as average-association summaries. A prespecified time-varying Cox model added log-time interactions for C2, C3, and age and estimated adjusted subtype contrasts at 3, 6, 12, 24, 36, and 60 months. Leave-one-study-out models evaluated the direction and precision of the C3-versus-C1 association. Clinical feature comparisons were exploratory, unadjusted, available-case tests and were not interpreted as causal or independent of study composition.", "",
  "Two-sided P values were used. Benjamini-Hochberg correction was applied within the pairwise subtype, time-varying-effect, Cox-coefficient, and clinical-association families. Missingness, model populations, events, formulas, reference levels, and PH diagnostics are saved as source tables."
)
writeLines(methods, file.path(report_dir, "step6_methods_ready.md"))

results <- c(
  "# Clinical and survival results-ready evidence", "",
  sprintf("Among 331 principal-state tumors with valid survival data, 153 deaths were observed. Median overall survival was %.1f months for C1, %.1f months for C2, and %.1f months for C3. The study-stratified log-rank P value was %s.",
          medians[FinalSubtype == "C1", MedianMonths], medians[FinalSubtype == "C2", MedianMonths],
          medians[FinalSubtype == "C3", MedianMonths], format.pval(km[Test == "Study-stratified log-rank", PValue], digits = 3, eps = 0.001)), "",
  sprintf("The constant-coefficient core Cox model summarized an average C3-versus-C1 association of HR %.2f (95%% CI %.2f–%.2f; P=%s) and C3-versus-C2 HR %.2f (95%% CI %.2f–%.2f; P=%s). C2 did not differ significantly from C1 (HR %.2f, 95%% CI %.2f–%.2f; P=%s). These constant HRs are not interpreted as time-invariant because the global PH test was significant (P=%s).",
          c3c1$HR, c3c1$Lower95, c3c1$Upper95, format.pval(c3c1$PValue, digits = 3, eps = 0.001),
          c3c2$HR, c3c2$Lower95, c3c2$Upper95, format.pval(c3c2$PValue, digits = 3, eps = 0.001),
          c2c1$HR, c2c1$Lower95, c2c1$Upper95, format.pval(c2c1$PValue, digits = 3, eps = 0.001),
          format.pval(ph[Model == "Core study-stratified Cox" & Term == "GLOBAL", PValue], digits = 3, eps = 0.001)), "",
  sprintf("The time-varying model localized the excess C3 risk to earlier follow-up: C3 versus C1 HR %.2f (95%% CI %.2f–%.2f) at 6 months and %.2f (%.2f–%.2f) at 12 months, declining to %.2f (%.2f–%.2f) at 24 months. Thus the supported conclusion is early adverse survival associated with C3, not a constant lifetime hazard ratio.",
          c3_6$HR, c3_6$Lower95, c3_6$Upper95, c3_12$HR, c3_12$Lower95, c3_12$Upper95,
          c3_24$HR, c3_24$Lower95, c3_24$Upper95), "",
  sprintf("In the extended model (n=%d, %d deaths), the average C3-versus-C1 association remained directionally similar (HR %.2f, 95%% CI %.2f–%.2f). Leave-one-study-out C3-versus-C1 HRs ranged from %.2f to %.2f; omission of 2022Dong left only 59 C1/C3 tumors and produced a non-significant but directionally adverse estimate, so cross-study precision is limited.",
          extended$N[[1]], extended$Events[[1]], extended[Term == "FinalSubtypeC3", HR],
          extended[Term == "FinalSubtypeC3", Lower95], extended[Term == "FinalSubtypeC3", Upper95],
          min(loso$HR), max(loso$HR)), "",
  "Clinical covariate distributions are descriptive available-case evidence. Observed duct morphology is kept separate from Step3 predictions, and therapeutic-response predictions are not part of this survival gate."
)
writeLines(results, file.path(report_dir, "step6_results_ready.md"))

captions <- data.table(
  Figure = c("step6_km_main_states", "step6_cox_time_varying_subtype_effects",
             "step6_cox_main_states_forest", "step6_c3_vs_c1_leave_one_study_out"),
  Scope = c("Main", "Main", "Supplementary", "Supplementary"),
  Caption = c(
    "Available-case Kaplan-Meier curves for 331 principal-state tumors (153 deaths); P value from a study-stratified log-rank test; numbers at risk shown below the plot.",
    "Adjusted subtype hazard ratios at prespecified follow-up times from a study-stratified Cox model with log-time interactions for subtype and age; ribbons show 95% confidence intervals.",
    "Average-association Cox estimates adjusted for age, sex, and stage with study-specific baseline hazards; the proportional-hazards violation is addressed in the time-varying panel.",
    "C3-versus-C1 adjusted Cox sensitivity after omitting each survival-contributing study; estimates retain study-stratified baseline hazards where applicable."
  ),
  SourceTable = c("step6_km_logrank_tests.csv;step6_km_median_survival.csv",
                  "step6_cox_time_varying_subtype_effects.csv", "step6_cox_main_states_core.csv",
                  "step6_c3_vs_c1_leave_one_study_out.csv")
)
fwrite(captions, file.path(table_dir, "step6_figure_caption_drafts.csv"))

manifest <- data.table(
  Artifact = c(required_tables, "step6_clinical_survival_release_gate_status.csv", "step6_figure_caption_drafts.csv",
               "step6_methods_ready.md", "step6_results_ready.md", required_figures),
  Path = c(file.path(table_dir, required_tables), file.path(table_dir, "step6_clinical_survival_release_gate_status.csv"),
           file.path(table_dir, "step6_figure_caption_drafts.csv"), file.path(report_dir, "step6_methods_ready.md"),
           file.path(report_dir, "step6_results_ready.md"), file.path(figure_dir, required_figures))
)
manifest[, Exists := vapply(Path, nonempty, logical(1))]
fwrite(manifest, file.path(table_dir, "step6_output_manifest.csv"))

cat("Clinical/survival release gate:", ifelse(all(checks$Status == "PASS"), "PASS", "BLOCKED"), "\n")
if (!all(checks$Status == "PASS")) quit(status = 1L)
