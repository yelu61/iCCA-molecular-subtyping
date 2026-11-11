#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260715)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "6-ClinicalSurvival/publication")
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
report_dir <- file.path(out_dir, "reports")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

master_path <- file.path(root, "governance/data/subtype_master_table.csv")
dict_path <- file.path(root, "governance/data/subtype_label_color_dictionary.csv")
master <- fread(master_path, showProgress = FALSE)
dictionary <- fread(dict_path, showProgress = FALSE)
main_states <- c("C1", "C2", "C3")
all_states <- c(main_states, "C1_C2", "C1_C3", "C2_C3")
subtype_colors <- setNames(dictionary$Color, dictionary$FinalSubtype)

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
clean_binary <- function(x) {
  z <- num(x)
  fifelse(z %in% c(0, 1), z, NA_real_)
}

prepare <- function(x, state_levels) {
  z <- copy(x)
  z[, Sample := as.character(Sample)]
  z[, os_months := num(`Overall survival (months)`)]
  z[, os_event := clean_binary(`Vital state (1=Dead)`)]
  z[, age_years := num(Age_num)]
  z[, tumor_size_cm := num(`Tumor size (cm)`)]
  z[, hbv_positive := clean_binary(`HBV infection (1=Positive)`)]
  z[, vascular_invasion := clean_binary(`Vascular invasion (1=Positive)`)]
  z[, sex_clean := fifelse(Sex %chin% c("Female", "Male"), Sex, NA_character_)]
  z[, stage_group := fifelse(Stage %chin% c("I", "II"), "I/II",
                             fifelse(Stage %chin% c("III", "IV"), "III/IV", NA_character_))]
  z[, observed_duct := fifelse(`iCCA classification` %chin% c("small duct type", "large duct type"),
                               `iCCA classification`, NA_character_)]
  z[, FinalSubtype := factor(FinalSubtype, levels = state_levels)]
  z[, sex_clean := factor(sex_clean, levels = c("Female", "Male"))]
  z[, stage_group := factor(stage_group, levels = c("I/II", "III/IV"))]
  z
}

main <- prepare(master[SubtypeClass == "Main"], main_states)
all_state <- prepare(master, all_states)

survival_valid <- function(x) x[!is.na(os_months) & os_months > 0 & os_event %in% c(0, 1)]
main_survival <- survival_valid(main)
all_survival <- survival_valid(all_state)
core_vars <- c("os_months", "os_event", "FinalSubtype", "age_years", "sex_clean", "stage_group", "Study")
extended_vars <- c(core_vars, "tumor_size_cm", "hbv_positive", "vascular_invasion")
main_core <- main_survival[complete.cases(main_survival[, ..core_vars])]
main_extended <- main_survival[complete.cases(main_survival[, ..extended_vars])]
all_core <- all_survival[complete.cases(all_survival[, ..core_vars])]

if (nrow(main_core) != 331L || sum(main_core$os_event) != 153L) {
  stop("Unexpected core survival denominator: N=", nrow(main_core), " events=", sum(main_core$os_event))
}
if (anyDuplicated(main_core$Sample)) stop("Survival core contains duplicate samples")

fit_core <- coxph(Surv(os_months, os_event) ~ FinalSubtype + age_years + sex_clean + stage_group + strata(Study),
                  data = main_core, ties = "efron", x = TRUE)
fit_extended <- coxph(Surv(os_months, os_event) ~ FinalSubtype + age_years + sex_clean + stage_group +
                          tumor_size_cm + hbv_positive + vascular_invasion + strata(Study),
                        data = main_extended, ties = "efron", x = TRUE)
fit_all <- coxph(Surv(os_months, os_event) ~ FinalSubtype + age_years + sex_clean + stage_group + strata(Study),
                  data = all_core, ties = "efron", x = TRUE)

tidy_cox <- function(fit, analysis, model, n, events, reference, scope) {
  s <- summary(fit)
  data.table(
    Analysis = analysis, Model = model, Scope = scope, Reference = reference,
    Term = rownames(s$coefficients), HR = unname(s$conf.int[, "exp(coef)"]),
    Lower95 = unname(s$conf.int[, "lower .95"]), Upper95 = unname(s$conf.int[, "upper .95"]),
    PValue = unname(s$coefficients[, "Pr(>|z|)"]), N = n, Events = events
  )[, FDR := p.adjust(PValue, method = "BH")]
}

core_tidy <- tidy_cox(fit_core, "Principal states", "Core study-stratified Cox",
                      nrow(main_core), sum(main_core$os_event), "C1; Female; stage I/II", "Headline")
extended_tidy <- tidy_cox(fit_extended, "Principal states", "Extended study-stratified sensitivity",
                          nrow(main_extended), sum(main_extended$os_event), "C1; Female; stage I/II", "Sensitivity")
all_tidy <- tidy_cox(fit_all, "All six states", "Core study-stratified Cox",
                     nrow(all_core), sum(all_core$os_event), "C1; Female; stage I/II", "Exploratory")
fwrite(core_tidy, file.path(table_dir, "step6_cox_main_states_core.csv"))
fwrite(extended_tidy, file.path(table_dir, "step6_cox_main_states_extended_sensitivity.csv"))
fwrite(all_tidy, file.path(table_dir, "step6_cox_all_states_exploratory.csv"))

# Complete adjusted pairwise subtype contrasts use the same core population.
pairwise_rows <- list()
for (ref in main_states) {
  d <- copy(main_core)
  d[, SubtypeRefit := relevel(factor(as.character(FinalSubtype), levels = main_states), ref = ref)]
  fit <- coxph(Surv(os_months, os_event) ~ SubtypeRefit + age_years + sex_clean + stage_group + strata(Study),
               data = d, ties = "efron")
  s <- summary(fit)
  terms <- grep("^SubtypeRefit", rownames(s$coefficients), value = TRUE)
  for (term in terms) {
    target <- sub("^SubtypeRefit", "", term)
    pairwise_rows[[length(pairwise_rows) + 1L]] <- data.table(
      Target = target, Reference = ref,
      HR = unname(s$conf.int[term, "exp(coef)"]), Lower95 = unname(s$conf.int[term, "lower .95"]),
      Upper95 = unname(s$conf.int[term, "upper .95"]), PValue = unname(s$coefficients[term, "Pr(>|z|)"]),
      N = nrow(d), Events = sum(d$os_event),
      Formula = "Surv(OS months, event) ~ subtype + age + sex + stage group + strata(Study)"
    )
  }
}
pairwise <- unique(rbindlist(pairwise_rows), by = c("Target", "Reference"))
pairwise <- pairwise[match(c("C2|C1", "C3|C1", "C3|C2"), paste(Target, Reference, sep = "|"))]
pairwise[, FDR := p.adjust(PValue, method = "BH")]
fwrite(pairwise, file.path(table_dir, "step6_cox_main_states_pairwise.csv"))

ph_table <- function(fit, model, scope) {
  z <- cox.zph(fit, transform = "km")
  out <- as.data.table(z$table, keep.rownames = "Term")
  setnames(out, c("Term", "chisq", "df", "p"), c("Term", "ChiSquare", "DF", "PValue"))
  out[, `:=`(Model = model, Scope = scope)]
  setcolorder(out, c("Model", "Scope", "Term", "ChiSquare", "DF", "PValue"))
  out
}
ph <- rbindlist(list(
  ph_table(fit_core, "Core study-stratified Cox", "Headline"),
  ph_table(fit_extended, "Extended study-stratified sensitivity", "Sensitivity"),
  ph_table(fit_all, "All-state core study-stratified Cox", "Exploratory")
))
fwrite(ph, file.path(table_dir, "step6_cox_proportional_hazards_diagnostics.csv"))

# The core model violates proportional hazards for subtype and age. A
# prespecified log-time interaction model therefore estimates clinically
# interpretable subtype contrasts at fixed follow-up times rather than treating
# one constant HR as the full time course.
tv_data <- copy(main_core)
tv_data[, `:=`(C2_indicator = as.integer(FinalSubtype == "C2"),
               C3_indicator = as.integer(FinalSubtype == "C3"))]
fit_time_varying <- coxph(
  Surv(os_months, os_event) ~ C2_indicator + C3_indicator + age_years + sex_clean + stage_group +
    tt(C2_indicator) + tt(C3_indicator) + tt(age_years) + strata(Study),
  data = tv_data, ties = "efron", x = TRUE,
  tt = function(x, t, ...) x * log(pmax(t, 0.1))
)
tv_coef <- coef(fit_time_varying)
tv_vcov <- vcov(fit_time_varying)
tv_times <- c(3, 6, 12, 24, 36, 60)
contrast_vector <- function(contrast, time) {
  v <- setNames(rep(0, length(tv_coef)), names(tv_coef))
  lt <- log(time)
  if (contrast == "C2 vs C1") {
    v["C2_indicator"] <- 1; v["tt(C2_indicator)"] <- lt
  } else if (contrast == "C3 vs C1") {
    v["C3_indicator"] <- 1; v["tt(C3_indicator)"] <- lt
  } else if (contrast == "C3 vs C2") {
    v["C3_indicator"] <- 1; v["C2_indicator"] <- -1
    v["tt(C3_indicator)"] <- lt; v["tt(C2_indicator)"] <- -lt
  }
  v
}
tv_effects <- rbindlist(lapply(c("C2 vs C1", "C3 vs C1", "C3 vs C2"), function(contrast) {
  rbindlist(lapply(tv_times, function(time) {
    v <- contrast_vector(contrast, time)
    eta <- sum(v * tv_coef)
    se <- sqrt(as.numeric(t(v) %*% tv_vcov %*% v))
    data.table(Contrast = contrast, TimeMonths = time, HR = exp(eta),
               Lower95 = exp(eta - 1.96 * se), Upper95 = exp(eta + 1.96 * se),
               PValue = 2 * pnorm(-abs(eta / se)), N = nrow(tv_data), Events = sum(tv_data$os_event))
  }))
}))
tv_effects[, FDR := p.adjust(PValue, method = "BH")]
fwrite(tv_effects, file.path(table_dir, "step6_cox_time_varying_subtype_effects.csv"))

p_tv <- ggplot(tv_effects[Contrast %chin% c("C2 vs C1", "C3 vs C1")],
               aes(TimeMonths, HR, color = Contrast, fill = Contrast)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey55") +
  geom_ribbon(aes(ymin = Lower95, ymax = Upper95), alpha = 0.16, color = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16), limits = c(0.08, 25)) +
  scale_color_manual(values = c("C2 vs C1" = subtype_colors[["C2"]], "C3 vs C1" = subtype_colors[["C3"]])) +
  scale_fill_manual(values = c("C2 vs C1" = subtype_colors[["C2"]], "C3 vs C1" = subtype_colors[["C3"]])) +
  labs(title = "Subtype-associated mortality risk changes over follow-up",
       subtitle = "Study-stratified Cox with log-time interactions for subtype and age; n=331, events=153",
       x = "Follow-up time (months)", y = "Adjusted hazard ratio (log scale)", color = NULL, fill = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(figure_dir, "step6_cox_time_varying_subtype_effects.pdf"), p_tv, width = 7.0, height = 4.8, device = cairo_pdf)
ggsave(file.path(figure_dir, "step6_cox_time_varying_subtype_effects.png"), p_tv, width = 7.0, height = 4.8, dpi = 400)

# Study-stratified log-rank is reported alongside unadjusted KM curves.
km_fit <- survfit(Surv(os_months, os_event) ~ FinalSubtype, data = main_survival)
lr_unstrat <- survdiff(Surv(os_months, os_event) ~ FinalSubtype, data = main_survival)
lr_strat <- survdiff(Surv(os_months, os_event) ~ FinalSubtype + strata(Study), data = main_survival)
logrank <- data.table(
  Test = c("Unstratified log-rank", "Study-stratified log-rank"),
  N = nrow(main_survival), Events = sum(main_survival$os_event),
  ChiSquare = c(unname(lr_unstrat$chisq), unname(lr_strat$chisq)),
  DF = c(length(main_states) - 1L, length(main_states) - 1L),
  PValue = c(pchisq(lr_unstrat$chisq, df = length(lr_unstrat$n) - 1L, lower.tail = FALSE),
             pchisq(lr_strat$chisq, df = length(main_states) - 1L, lower.tail = FALSE))
)
fwrite(logrank, file.path(table_dir, "step6_km_logrank_tests.csv"))

km_summary <- summary(km_fit)$table
median_table <- data.table(FinalSubtype = sub("FinalSubtype=", "", rownames(km_summary)),
                           N = km_summary[, "records"], Events = km_summary[, "events"],
                           MedianMonths = km_summary[, "median"],
                           Lower95 = km_summary[, "0.95LCL"], Upper95 = km_summary[, "0.95UCL"])
fwrite(median_table, file.path(table_dir, "step6_km_median_survival.csv"))

stratified_p <- format.pval(logrank[Test == "Study-stratified log-rank", PValue], digits = 2, eps = 0.001)
stratified_p_label <- ifelse(startsWith(stratified_p, "<"), paste0("P", stratified_p), paste0("P=", stratified_p))
km_plot <- ggsurvplot(
  km_fit, data = main_survival, risk.table = TRUE, conf.int = FALSE,
  palette = unname(subtype_colors[main_states]), legend.title = "",
  legend.labs = main_states, xlab = "Overall survival (months)", ylab = "Survival probability",
  break.time.by = 12, risk.table.height = 0.25, censor.shape = 124, censor.size = 3,
  ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable(base_size = 9),
  title = "Overall survival across principal molecular states",
  subtitle = sprintf("Available-case KM; n=%d, events=%d; study-stratified log-rank %s",
                     nrow(main_survival), sum(main_survival$os_event), stratified_p_label)
)
pdf(file.path(figure_dir, "step6_km_main_states.pdf"), width = 7.2, height = 6.2, useDingbats = FALSE)
print(km_plot)
dev.off()
png(file.path(figure_dir, "step6_km_main_states.png"), width = 2880, height = 2480, res = 400)
print(km_plot)
dev.off()

term_label <- c("FinalSubtypeC2" = "C2 vs C1", "FinalSubtypeC3" = "C3 vs C1",
                "age_years" = "Age, per year", "sex_cleanMale" = "Male vs female",
                "stage_groupIII/IV" = "Stage III/IV vs I/II")
forest <- copy(core_tidy[Term %chin% names(term_label)])
forest[, Label := factor(term_label[Term], levels = rev(term_label))]
forest[, Text := sprintf("%.2f (%.2f–%.2f)", HR, Lower95, Upper95)]
p_forest <- ggplot(forest, aes(HR, Label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = Lower95, xmax = Upper95), orientation = "y", width = 0.18, color = "#444444") +
  geom_point(size = 2.5, color = "#0072B2") +
  scale_x_log10(breaks = c(0.5, 1, 2, 4)) +
  labs(title = "Study-stratified multivariable overall-survival model",
       subtitle = "Average Cox associations; PH violation is addressed in the time-varying effect panel",
       x = "Adjusted hazard ratio (log scale)", y = NULL) +
  theme_classic(base_size = 11)
ggsave(file.path(figure_dir, "step6_cox_main_states_forest.pdf"), p_forest, width = 7.0, height = 4.6, device = cairo_pdf)
ggsave(file.path(figure_dir, "step6_cox_main_states_forest.png"), p_forest, width = 7.0, height = 4.6, dpi = 400)

# Leave-one-study-out sensitivity for the C3 versus C1 adjusted contrast.
loso <- rbindlist(lapply(unique(main_core$Study), function(held) {
  d <- main_core[Study != held & FinalSubtype %chin% c("C1", "C3")]
  d[, FinalSubtype := relevel(droplevels(FinalSubtype), ref = "C1")]
  fit <- coxph(Surv(os_months, os_event) ~ FinalSubtype + age_years + sex_clean + stage_group + strata(Study),
               data = d, ties = "efron")
  s <- summary(fit)
  term <- "FinalSubtypeC3"
  data.table(HeldOutStudy = held, N = nrow(d), Events = sum(d$os_event),
             HR = unname(s$conf.int[term, "exp(coef)"]), Lower95 = unname(s$conf.int[term, "lower .95"]),
             Upper95 = unname(s$conf.int[term, "upper .95"]), PValue = unname(s$coefficients[term, "Pr(>|z|)"]))
}))
fwrite(loso, file.path(table_dir, "step6_c3_vs_c1_leave_one_study_out.csv"))
p_loso <- ggplot(loso, aes(HR, reorder(HeldOutStudy, HR))) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = Lower95, xmax = Upper95), orientation = "y", width = 0.18) +
  geom_point(size = 2.4, color = "#CC79A7") +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8)) +
  labs(title = "C3 versus C1 survival sensitivity",
       subtitle = "Core covariates with one survival-contributing study omitted at a time",
       x = "Adjusted hazard ratio (log scale)", y = "Held-out study") +
  theme_classic(base_size = 11)
ggsave(file.path(figure_dir, "step6_c3_vs_c1_leave_one_study_out.pdf"), p_loso, width = 6.6, height = 4.4, device = cairo_pdf)
ggsave(file.path(figure_dir, "step6_c3_vs_c1_leave_one_study_out.png"), p_loso, width = 6.6, height = 4.4, dpi = 400)

# Population, missingness, and clinical descriptive evidence.
population <- master[, .(TotalN = .N,
                         ValidOSN = sum(!is.na(num(`Overall survival (months)`)) & num(`Overall survival (months)`) > 0 &
                                          clean_binary(`Vital state (1=Dead)`) %in% c(0, 1)),
                         Events = sum(clean_binary(`Vital state (1=Dead)`) == 1 &
                                        !is.na(num(`Overall survival (months)`)) & num(`Overall survival (months)`) > 0, na.rm = TRUE)),
                     by = .(Study, FinalSubtype, SubtypeClass)]
fwrite(population, file.path(table_dir, "step6_survival_population_by_study_and_state.csv"))

missing_variables <- c(os_months = "Overall survival time", os_event = "Vital status", age_years = "Age",
                       sex_clean = "Sex", stage_group = "Stage group", tumor_size_cm = "Tumor size",
                       hbv_positive = "HBV", vascular_invasion = "Vascular invasion", observed_duct = "Observed duct type")
missingness <- rbindlist(lapply(c("All studies", unique(main$Study)), function(study) {
  d <- if (study == "All studies") main else main[Study == study]
  rbindlist(lapply(names(missing_variables), function(v) data.table(
    Study = study, Variable = missing_variables[[v]], EligibleN = nrow(d),
    AvailableN = sum(!is.na(d[[v]])), MissingN = sum(is.na(d[[v]])), MissingFraction = mean(is.na(d[[v]]))
  )))
}))
fwrite(missingness, file.path(table_dir, "step6_clinical_variable_missingness.csv"))

continuous_vars <- c(age_years = "Age (years)", tumor_size_cm = "Tumor size (cm)")
categorical_vars <- c(sex_clean = "Sex", stage_group = "Stage group", hbv_positive = "HBV positive",
                      vascular_invasion = "Vascular invasion", observed_duct = "Observed duct type")
baseline <- list()
tests <- list()
for (v in names(continuous_vars)) {
  baseline[[length(baseline) + 1L]] <- main[!is.na(get(v)), .(
    Variable = continuous_vars[[v]], Level = "Median (IQR)", AvailableN = .N,
    Summary = sprintf("%.1f (%.1f–%.1f)", median(get(v)), quantile(get(v), .25), quantile(get(v), .75))
  ), by = FinalSubtype]
  d <- main[!is.na(get(v))]
  kw <- kruskal.test(d[[v]] ~ d$FinalSubtype)
  tests[[length(tests) + 1L]] <- data.table(Variable = continuous_vars[[v]], Test = "Kruskal-Wallis",
                                             N = nrow(d), PValue = kw$p.value,
                                             Interpretation = "Exploratory unadjusted clinical association")
}
for (v in names(categorical_vars)) {
  d <- main[!is.na(get(v))]
  baseline[[length(baseline) + 1L]] <- d[, .N, by = .(FinalSubtype, Level = as.character(get(v)))][,
    `:=`(Variable = categorical_vars[[v]], AvailableN = sum(N), Summary = sprintf("%d (%.1f%%)", N, 100 * N / sum(N))),
    by = FinalSubtype]
  tab <- table(d$FinalSubtype, d[[v]])
  ft <- fisher.test(tab, simulate.p.value = any(dim(tab) > 2L), B = 20000)
  tests[[length(tests) + 1L]] <- data.table(Variable = categorical_vars[[v]], Test = "Fisher exact",
                                             N = nrow(d), PValue = ft$p.value,
                                             Interpretation = "Exploratory unadjusted clinical association")
}
baseline <- rbindlist(baseline, fill = TRUE)
setcolorder(baseline, c("Variable", "Level", "FinalSubtype", "AvailableN", "Summary", "N"), skip_absent = TRUE)
tests <- rbindlist(tests)
tests[, FDR := p.adjust(PValue, method = "BH")]
fwrite(baseline, file.path(table_dir, "step6_clinical_baseline_by_main_state.csv"))
fwrite(tests, file.path(table_dir, "step6_clinical_association_tests.csv"))

model_summary <- data.table(
  Analysis = c("Main states", "Main states", "Main states", "All six states"),
  Scope = c("Average-association summary", "Non-PH adjudication", "Sensitivity", "Exploratory"),
  N = c(nrow(main_core), nrow(tv_data), nrow(main_extended), nrow(all_core)),
  Events = c(sum(main_core$os_event), sum(tv_data$os_event), sum(main_extended$os_event), sum(all_core$os_event)),
  Formula = c(
    "Surv(OS months, event) ~ subtype + age + sex + stage group + strata(Study)",
    "Surv(OS months, event) ~ subtype indicators + age + sex + stage group + log(time) interactions for subtype and age + strata(Study)",
    "Surv(OS months, event) ~ subtype + age + sex + stage group + tumor size + HBV + vascular invasion + strata(Study)",
    "Surv(OS months, event) ~ six-state subtype + age + sex + stage group + strata(Study)"
  ),
  Reference = "C1; Female; stage I/II"
)
fwrite(model_summary, file.path(table_dir, "step6_survival_model_contract.csv"))

cat("Clinical/survival analysis complete: core N=", nrow(main_core), " events=", sum(main_core$os_event), "\n", sep = "")
