#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
})

set.seed(20260723)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "6-ClinicalSurvival/publication/tables")
figure_dir <- file.path(root, "6-ClinicalSurvival/publication/figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
clean_binary <- function(x) {
  z <- num(x)
  fifelse(z %in% c(0, 1), z, NA_real_)
}

master <- fread(file.path(root, "governance/data/subtype_master_table.csv"), showProgress = FALSE)
d <- copy(master[SubtypeClass == "Main"])
d[, `:=`(
  Sample = as.character(Sample),
  FinalSubtype = factor(FinalSubtype, levels = c("C1", "C2", "C3")),
  age = num(Age_num),
  sex = factor(fifelse(Sex %chin% c("Female", "Male"), Sex, NA_character_), levels = c("Female", "Male")),
  stage = factor(fifelse(Stage %chin% c("I", "II"), "I/II",
                         fifelse(Stage %chin% c("III", "IV"), "III/IV", NA_character_)),
                 levels = c("I/II", "III/IV")),
  stage_advanced = fifelse(Stage %chin% c("III", "IV"), 1,
                           fifelse(Stage %chin% c("I", "II"), 0, NA_real_)),
  tumor_size = num(`Tumor size (cm)`),
  vascular_invasion = clean_binary(`Vascular invasion (1=Positive)`),
  observed_large_duct = fifelse(`iCCA classification` == "large duct type", 1,
                                fifelse(`iCCA classification` == "small duct type", 0, NA_real_)),
  os_time = num(`Overall survival (months)`),
  os_event = clean_binary(`Vital state (1=Dead)`)
)]

contrasts <- data.table(Target = c("C2", "C3", "C3"), Reference = c("C1", "C1", "C2"))

random_effects_dl <- function(yi, sei) {
  keep <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[keep]; sei <- sei[keep]; k <- length(yi)
  if (k < 2) return(NULL)
  w <- 1 / sei^2; mu_fixed <- sum(w * yi) / sum(w)
  q <- sum(w * (yi - mu_fixed)^2); df <- k - 1
  cval <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (q - df) / cval)
  wr <- 1 / (sei^2 + tau2); mu <- sum(wr * yi) / sum(wr); se <- sqrt(1 / sum(wr))
  list(k = k, estimate = mu, se = se, ci.lb = mu - 1.96 * se, ci.ub = mu + 1.96 * se,
       pval = 2 * pnorm(-abs(mu / se)), tau2 = tau2,
       I2 = ifelse(q > 0, max(0, (q - df) / q) * 100, 0), QP = pchisq(q, df, lower.tail = FALSE))
}
endpoint_contract <- data.table(
  Endpoint = c("Tumor size", "Advanced stage", "Vascular invasion", "Observed large-duct morphology"),
  Variable = c("tumor_size", "stage_advanced", "vascular_invasion", "observed_large_duct"),
  Type = c("continuous", "binary", "binary", "binary"),
  PositiveLevel = c("larger size in cm", "stage III/IV", "vascular invasion present", "observed large-duct type")
)

hc3_contrast <- function(dd, variable, type, target, reference) {
  dd <- copy(dd[!is.na(get(variable))])
  dd[, SubtypeRefit := relevel(droplevels(FinalSubtype), ref = reference)]
  form <- as.formula(sprintf("%s ~ SubtypeRefit + factor(Study)", variable))
  fit <- if (type == "continuous") lm(form, data = dd) else glm(form, data = dd, family = binomial())
  robust <- coeftest(fit, vcov. = vcovHC(fit, type = "HC3"))
  term <- paste0("SubtypeRefit", target)
  if (!term %in% rownames(robust)) return(NULL)
  beta <- robust[term, 1]; se <- robust[term, 2]
  data.table(Target = target, Reference = reference, N = nrow(dd), MissingN = nrow(d) - nrow(dd),
             Studies = uniqueN(dd$Study), Estimate = if (type == "binary") exp(beta) else beta,
             Lower95 = if (type == "binary") exp(beta - 1.96 * se) else beta - 1.96 * se,
             Upper95 = if (type == "binary") exp(beta + 1.96 * se) else beta + 1.96 * se,
             PValue = robust[term, 4], EffectScale = ifelse(type == "binary", "Odds ratio", "Mean difference (cm)"),
             Model = "Outcome ~ subtype + Study; HC3 robust covariance")
}

phenotype_rows <- list()
for (i in seq_len(nrow(endpoint_contract))) {
  ep <- endpoint_contract[i]
  for (j in seq_len(nrow(contrasts))) {
    z <- hc3_contrast(d, ep$Variable, ep$Type, contrasts$Target[j], contrasts$Reference[j])
    if (!is.null(z)) {
      z[, `:=`(Endpoint = ep$Endpoint, PositiveLevel = ep$PositiveLevel)]
      phenotype_rows[[length(phenotype_rows) + 1L]] <- z
    }
  }
}
phenotype <- rbindlist(phenotype_rows, fill = TRUE)
phenotype[, FDR := p.adjust(PValue, method = "BH")]
setcolorder(phenotype, c("Endpoint", "PositiveLevel", "Target", "Reference", "N", "MissingN", "Studies",
                         "EffectScale", "Estimate", "Lower95", "Upper95", "PValue", "FDR", "Model"))
fwrite(phenotype, file.path(table_dir, "step6_study_adjusted_clinical_phenotype_contrasts.csv"))

availability <- rbindlist(lapply(seq_len(nrow(endpoint_contract)), function(i) {
  ep <- endpoint_contract[i]
  rbindlist(lapply(c("All studies", sort(unique(d$Study))), function(st) {
    dd <- if (st == "All studies") d else d[Study == st]
    data.table(Endpoint = ep$Endpoint, Study = st, EligibleN = nrow(dd),
               AvailableN = sum(!is.na(dd[[ep$Variable]])), MissingN = sum(is.na(dd[[ep$Variable]])),
               ContributingSubtypes = paste(levels(droplevels(dd[!is.na(get(ep$Variable))]$FinalSubtype)), collapse = ";"))
  }))
}))
fwrite(availability, file.path(table_dir, "step6_clinical_phenotype_availability.csv"))

study_effects <- list()
for (i in seq_len(nrow(endpoint_contract))) {
  ep <- endpoint_contract[i]
  for (j in seq_len(nrow(contrasts))) {
    target <- contrasts$Target[j]; reference <- contrasts$Reference[j]
    for (st in sort(unique(d$Study))) {
      dd <- d[Study == st & FinalSubtype %chin% c(target, reference) & !is.na(get(ep$Variable))]
      nt <- sum(dd$FinalSubtype == target); nr <- sum(dd$FinalSubtype == reference)
      if (nt < 3 || nr < 3) next
      if (ep$Type == "continuous") {
        xt <- dd[FinalSubtype == target][[ep$Variable]]; xr <- dd[FinalSubtype == reference][[ep$Variable]]
        yi <- mean(xt) - mean(xr); sei <- sqrt(var(xt) / length(xt) + var(xr) / length(xr)); method <- "Welch mean difference"
      } else {
        a <- sum(dd$FinalSubtype == target & dd[[ep$Variable]] == 1)
        b <- sum(dd$FinalSubtype == target & dd[[ep$Variable]] == 0)
        c0 <- sum(dd$FinalSubtype == reference & dd[[ep$Variable]] == 1)
        e0 <- sum(dd$FinalSubtype == reference & dd[[ep$Variable]] == 0)
        corrected <- any(c(a, b, c0, e0) == 0)
        if (corrected) { a <- a + .5; b <- b + .5; c0 <- c0 + .5; e0 <- e0 + .5 }
        yi <- log((a / b) / (c0 / e0)); sei <- sqrt(1 / a + 1 / b + 1 / c0 + 1 / e0)
        method <- ifelse(corrected, "Log odds ratio; Haldane-Anscombe correction", "Log odds ratio")
      }
      if (is.finite(yi) && is.finite(sei) && sei > 0) study_effects[[length(study_effects) + 1L]] <- data.table(
        Endpoint = ep$Endpoint, Type = ep$Type, Target = target, Reference = reference, Study = st,
        N = nrow(dd), NTarget = nt, NReference = nr, LogEffect = yi, SE = sei,
        Effect = ifelse(ep$Type == "binary", exp(yi), yi), Lower95 = ifelse(ep$Type == "binary", exp(yi - 1.96 * sei), yi - 1.96 * sei),
        Upper95 = ifelse(ep$Type == "binary", exp(yi + 1.96 * sei), yi + 1.96 * sei), Method = method)
    }
  }
}
study_effects <- rbindlist(study_effects, fill = TRUE)
fwrite(study_effects, file.path(table_dir, "step6_clinical_phenotype_study_specific_effects.csv"))

meta_rows <- list()
for (key in unique(study_effects[, paste(Endpoint, Target, Reference, sep = "|")])) {
  z <- study_effects[paste(Endpoint, Target, Reference, sep = "|") == key]
  if (nrow(z) < 2) next
  fit <- random_effects_dl(z$LogEffect, z$SE)
  if (is.null(fit)) next
  binary <- z$Type[1] == "binary"
  meta_rows[[length(meta_rows) + 1L]] <- data.table(
    Endpoint = z$Endpoint[1], Target = z$Target[1], Reference = z$Reference[1], Studies = nrow(z),
    EffectScale = ifelse(binary, "Random-effects odds ratio", "Random-effects mean difference (cm)"),
    Estimate = ifelse(binary, exp(fit$estimate), fit$estimate),
    Lower95 = ifelse(binary, exp(fit$ci.lb), fit$ci.lb), Upper95 = ifelse(binary, exp(fit$ci.ub), fit$ci.ub),
    PValue = fit$pval, Tau2 = fit$tau2, I2Percent = fit$I2, QPValue = fit$QP,
    MetaMethod = "DerSimonian-Laird random effects")
}
phenotype_meta <- rbindlist(meta_rows, fill = TRUE)
if (nrow(phenotype_meta)) phenotype_meta[, FDR := p.adjust(PValue, method = "BH")]
fwrite(phenotype_meta, file.path(table_dir, "step6_clinical_phenotype_random_effects_meta.csv"))

# Same-population repeated-CV comparison of clinical baseline versus clinical + subtype.
surv_dat <- d[!is.na(os_time) & os_time > 0 & os_event %in% c(0, 1) & complete.cases(age, sex, stage, Study)]
stopifnot(nrow(surv_dat) == 331L, sum(surv_dat$os_event) == 153L)
eval_times <- seq(3, 36, by = 3)
headline_times <- c(6, 12, 24)

make_folds <- function(dd, k, seed) {
  set.seed(seed); out <- integer(nrow(dd)); strata <- interaction(dd$Study, dd$os_event, drop = TRUE)
  for (lev in levels(strata)) { idx <- which(strata == lev); out[idx] <- sample(rep(seq_len(k), length.out = length(idx))) }
  out
}

predict_cox_risk <- function(fit, newdata, times) {
  bh <- basehaz(fit, centered = FALSE)
  lp <- as.numeric(predict(fit, newdata = newdata, type = "lp", reference = "zero"))
  out <- matrix(NA_real_, nrow(newdata), length(times))
  for (i in seq_len(nrow(newdata))) {
    st <- as.character(newdata$Study[i])
    b <- bh[as.character(bh$strata) == st | sub("^Study=", "", as.character(bh$strata)) == st, , drop = FALSE]
    if (!nrow(b)) stop("No baseline hazard for Study ", st)
    hfun <- stepfun(b$time, c(0, b$hazard), right = TRUE)
    out[i, ] <- 1 - exp(-hfun(times) * exp(lp[i]))
  }
  pmin(pmax(out, 0), 1)
}

pred_rows <- list(); fold_rows <- list()
for (repeat_id in 1:5) {
  fold_id <- make_folds(surv_dat, 5, 20260723 + repeat_id)
  fold_rows[[repeat_id]] <- data.table(Repeat = repeat_id, Sample = surv_dat$Sample, Study = surv_dat$Study,
                                      Event = surv_dat$os_event, Fold = fold_id)
  for (fold in 1:5) {
    tr <- surv_dat[fold_id != fold]; te <- surv_dat[fold_id == fold]
    fits <- list(
      Clinical = coxph(Surv(os_time, os_event) ~ age + sex + stage + strata(Study), data = tr, ties = "efron", x = TRUE, y = TRUE),
      ClinicalPlusSubtype = coxph(Surv(os_time, os_event) ~ age + sex + stage + FinalSubtype + strata(Study), data = tr, ties = "efron", x = TRUE, y = TRUE)
    )
    for (model_name in names(fits)) {
      risk <- predict_cox_risk(fits[[model_name]], newdata = te, times = eval_times)
      row <- data.table(Repeat = repeat_id, Fold = fold, Model = model_name, Sample = te$Sample, Study = te$Study,
                        Time = te$os_time, Event = te$os_event)
      for (ii in seq_along(eval_times)) row[[paste0("Risk", eval_times[ii])]] <- risk[, ii]
      pred_rows[[length(pred_rows) + 1L]] <- row
    }
  }
}
pred_long <- rbindlist(pred_rows, fill = TRUE)
fold_def <- rbindlist(fold_rows)
fwrite(pred_long, file.path(table_dir, "step6_incremental_prognostic_predictions_by_repeat.csv"))
fwrite(fold_def, file.path(table_dir, "step6_incremental_prognostic_fold_definitions.csv"))

risk_cols <- paste0("Risk", eval_times)
oof <- pred_long[, c(list(Study = Study[1], Time = Time[1], Event = Event[1]),
                          lapply(.SD, mean)), by = .(Model, Sample), .SDcols = risk_cols]
fwrite(oof, file.path(table_dir, "step6_incremental_prognostic_oof_predictions.csv"))

gfun <- function(time, event) {
  fit <- survfit(Surv(time, 1 - event) ~ 1)
  stepfun(fit$time, c(1, fit$surv), right = TRUE)
}
ipcw_brier <- function(time, event, risk, horizon) {
  G <- gfun(time, event); gt <- pmax(G(horizon), 1e-4); gi <- pmax(G(pmax(time - 1e-7, 0)), 1e-4)
  wcase <- as.numeric(time <= horizon & event == 1) / gi
  wcontrol <- as.numeric(time > horizon) / gt
  mean(wcase * (1 - risk)^2 + wcontrol * risk^2)
}
metric_set <- function(z) {
  risk_mat <- as.matrix(z[, ..risk_cols]); colnames(risk_mat) <- eval_times
  weighted_auc <- function(horizon, risk) {
    G <- gfun(z$Time, z$Event); cases <- which(z$Time <= horizon & z$Event == 1); controls <- which(z$Time > horizon)
    if (!length(cases) || !length(controls)) return(NA_real_)
    wc <- 1 / pmax(G(pmax(z$Time[cases] - 1e-7, 0)), 1e-4)
    wn <- rep(1 / pmax(G(horizon), 1e-4), length(controls))
    comp <- outer(risk[cases], risk[controls], ">") + 0.5 * outer(risk[cases], risk[controls], "==")
    sum(comp * outer(wc, wn)) / sum(outer(wc, wn))
  }
  auc <- vapply(headline_times, function(tt) weighted_auc(tt, risk_mat[, as.character(tt)]), numeric(1))
  brier <- vapply(seq_along(eval_times), function(i) ipcw_brier(z$Time, z$Event, risk_mat[, i], eval_times[i]), numeric(1))
  c(CIndex = concordance(Surv(Time, Event) ~ risk_mat[, which(eval_times == 12)], data = z, reverse = TRUE)$concordance,
    AUC6 = auc[1], AUC12 = auc[2], AUC24 = auc[3],
    Brier6 = brier[eval_times == 6], Brier12 = brier[eval_times == 12], Brier24 = brier[eval_times == 24],
    IBS3to36 = sum(diff(eval_times) * (head(brier, -1) + tail(brier, -1)) / 2) / diff(range(eval_times)))
}

models <- c("Clinical", "ClinicalPlusSubtype")
point <- setNames(lapply(models, function(m) metric_set(oof[Model == m])), models)
set.seed(20260723)
boot_n <- 500L
boot <- replicate(boot_n, {
  base <- oof[Model == "Clinical"]
  idx <- unlist(lapply(split(seq_len(nrow(base)), base$Study), function(ii) sample(ii, length(ii), replace = TRUE)))
  a <- metric_set(oof[Model == models[1]][idx])
  b <- metric_set(oof[Model == models[2]][idx])
  c(a, b, Delta = b - a)
})

perf_rows <- list(); metric_names <- names(point[[1]])
for (mi in seq_along(metric_names)) {
  met <- metric_names[mi]
  for (mm in seq_along(models)) {
    pos <- mi + (mm - 1) * length(metric_names)
    perf_rows[[length(perf_rows) + 1L]] <- data.table(Comparison = models[mm], Metric = met,
      Estimate = point[[models[mm]]][mi], Lower95 = quantile(boot[pos, ], .025, na.rm = TRUE),
      Upper95 = quantile(boot[pos, ], .975, na.rm = TRUE), N = nrow(surv_dat), Events = sum(surv_dat$os_event))
  }
  pos <- 2 * length(metric_names) + mi
  perf_rows[[length(perf_rows) + 1L]] <- data.table(Comparison = "Delta: clinical + subtype minus clinical", Metric = met,
    Estimate = point[[models[2]]][mi] - point[[models[1]]][mi], Lower95 = quantile(boot[pos, ], .025, na.rm = TRUE),
    Upper95 = quantile(boot[pos, ], .975, na.rm = TRUE), N = nrow(surv_dat), Events = sum(surv_dat$os_event))
}
perf <- rbindlist(perf_rows)
perf[, Interpretation := fifelse(grepl("Brier|IBS", Metric), "Negative delta favors adding subtype", "Positive delta favors adding subtype")]
fwrite(perf, file.path(table_dir, "step6_incremental_prognostic_performance.csv"))

calibration_rows <- list()
for (m in models) {
  z <- oof[Model == m]
  for (tt in headline_times) {
    p <- pmin(pmax(z[[paste0("Risk", tt)]], 1e-6), 1 - 1e-6)
    G <- gfun(z$Time, z$Event); gi <- pmax(G(pmax(z$Time - 1e-7, 0)), 1e-4); gt <- pmax(G(tt), 1e-4)
    y <- as.numeric(z$Time <= tt & z$Event == 1)
    w <- fifelse(z$Time <= tt & z$Event == 1, 1 / gi, fifelse(z$Time > tt, 1 / gt, 0))
    keep <- w > 0
    fit <- glm(y[keep] ~ qlogis(p[keep]), family = quasibinomial(), weights = w[keep])
    calibration_rows[[length(calibration_rows) + 1L]] <- data.table(Model = m, HorizonMonths = tt,
      CalibrationIntercept = coef(fit)[1], CalibrationSlope = coef(fit)[2], N = nrow(z),
      Method = "IPCW-weighted logistic calibration of repeated-CV OOF risk")
  }
}
calibration <- rbindlist(calibration_rows)
fwrite(calibration, file.path(table_dir, "step6_incremental_prognostic_calibration.csv"))

fit0 <- coxph(Surv(os_time, os_event) ~ age + sex + stage + strata(Study), data = surv_dat, ties = "efron")
fit1 <- coxph(Surv(os_time, os_event) ~ age + sex + stage + FinalSubtype + strata(Study), data = surv_dat, ties = "efron")
early <- copy(surv_dat); early[, `:=`(early_time = pmin(os_time, 24), early_event = as.integer(os_event == 1 & os_time <= 24))]
fit0e <- coxph(Surv(early_time, early_event) ~ age + sex + stage + strata(Study), data = early, ties = "efron")
fit1e <- coxph(Surv(early_time, early_event) ~ age + sex + stage + FinalSubtype + strata(Study), data = early, ties = "efron")
lrt_row <- function(a, b, scope, events) {
  lr <- 2 * (as.numeric(logLik(b)) - as.numeric(logLik(a))); df <- attr(logLik(b), "df") - attr(logLik(a), "df")
  data.table(Scope = scope, N = nrow(surv_dat), Events = events, ClinicalLogLik = as.numeric(logLik(a)),
             ClinicalPlusSubtypeLogLik = as.numeric(logLik(b)), LRChiSquare = lr, DF = df,
             PValue = pchisq(lr, df, lower.tail = FALSE),
             Interpretation = ifelse(scope == "Early 24-month restricted follow-up", "Addresses early C3 risk; not a lifetime effect", "Average-association benchmark; PH violation requires time-specific interpretation"))
}
lrt <- rbindlist(list(lrt_row(fit0, fit1, "Full follow-up", sum(surv_dat$os_event)),
                      lrt_row(fit0e, fit1e, "Early 24-month restricted follow-up", sum(early$early_event))))
fwrite(lrt, file.path(table_dir, "step6_incremental_prognostic_likelihood_ratio.csv"))

# Study-specific adjusted Cox effects and random-effects summaries.
cox_study <- list()
for (scope in c("Full follow-up", "Early 24-month restricted follow-up")) {
  for (st in sort(unique(surv_dat$Study))) {
    for (j in seq_len(nrow(contrasts))) {
      target <- contrasts$Target[j]; reference <- contrasts$Reference[j]
      dd <- copy(surv_dat[Study == st & FinalSubtype %chin% c(target, reference)])
      if (nrow(dd) < 15 || sum(dd$os_event) < 5 || min(table(droplevels(dd$FinalSubtype))) < 3) next
      dd[, SubtypeBinary := relevel(droplevels(FinalSubtype), ref = reference)]
      if (scope == "Early 24-month restricted follow-up") dd[, `:=`(analysis_time = pmin(os_time, 24), analysis_event = as.integer(os_event == 1 & os_time <= 24))]
      else dd[, `:=`(analysis_time = os_time, analysis_event = os_event)]
      fit <- tryCatch(coxph(Surv(analysis_time, analysis_event) ~ SubtypeBinary + age + sex + stage,
                             data = dd, ties = "efron"), error = function(e) NULL)
      if (is.null(fit)) next
      term <- paste0("SubtypeBinary", target); s <- summary(fit)
      if (!term %in% rownames(s$coefficients)) next
      beta <- unname(s$coefficients[term, "coef"]); se <- unname(s$coefficients[term, "se(coef)"])
      if (!is.finite(beta) || !is.finite(se) || se <= 0 || se > 4) next
      cox_study[[length(cox_study) + 1L]] <- data.table(Scope = scope, Study = st, Target = target, Reference = reference,
        N = nrow(dd), Events = sum(dd$analysis_event), LogHR = beta, SE = se, HR = exp(beta),
        Lower95 = exp(beta - 1.96 * se), Upper95 = exp(beta + 1.96 * se), PValue = 2 * pnorm(-abs(beta / se)),
        Formula = "Within-study Cox: subtype + age + sex + stage")
    }
  }
}
cox_study <- rbindlist(cox_study, fill = TRUE)
fwrite(cox_study, file.path(table_dir, "step6_subtype_cox_study_specific.csv"))

cox_meta <- list(); cox_influence <- list()
for (key in unique(cox_study[, paste(Scope, Target, Reference, sep = "|")])) {
  z <- cox_study[paste(Scope, Target, Reference, sep = "|") == key]
  if (nrow(z) < 2) next
  fit <- random_effects_dl(z$LogHR, z$SE)
  if (is.null(fit)) next
  cox_meta[[length(cox_meta) + 1L]] <- data.table(Scope = z$Scope[1], Target = z$Target[1], Reference = z$Reference[1],
    Studies = nrow(z), HR = exp(fit$estimate), Lower95 = exp(fit$ci.lb), Upper95 = exp(fit$ci.ub),
    PValue = fit$pval, Tau2 = fit$tau2, I2Percent = fit$I2, QPValue = fit$QP,
    MetaMethod = "DerSimonian-Laird", Interpretation = "Random-effects sensitivity; one-stage Study-stratified Cox remains primary")
  if (nrow(z) >= 3) for (held in z$Study) {
    zz <- z[Study != held]; ff <- random_effects_dl(zz$LogHR, zz$SE)
    if (!is.null(ff)) cox_influence[[length(cox_influence) + 1L]] <- data.table(Scope = z$Scope[1], Target = z$Target[1], Reference = z$Reference[1],
      HeldOutStudy = held, StudiesRemaining = nrow(zz), HR = exp(ff$estimate), Lower95 = exp(ff$ci.lb), Upper95 = exp(ff$ci.ub), PValue = ff$pval)
  }
}
cox_meta <- rbindlist(cox_meta, fill = TRUE); if (nrow(cox_meta)) cox_meta[, FDR := p.adjust(PValue, method = "BH")]
cox_influence <- rbindlist(cox_influence, fill = TRUE)
fwrite(cox_meta, file.path(table_dir, "step6_subtype_cox_random_effects_meta.csv"))
fwrite(cox_influence, file.path(table_dir, "step6_subtype_cox_meta_leave_one_study_out.csv"))

plot_perf <- perf[Comparison == "Delta: clinical + subtype minus clinical" & Metric %chin% c("CIndex", "AUC6", "AUC12", "AUC24", "IBS3to36")]
plot_perf[, Label := factor(Metric, levels = rev(c("CIndex", "AUC6", "AUC12", "AUC24", "IBS3to36")))]
p <- ggplot(plot_perf, aes(Estimate, Label)) + geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = Lower95, xmax = Upper95), orientation = "y", width = .18) +
  geom_point(size = 2.5, color = "#CC79A7") + theme_classic(base_size = 11) +
  labs(title = "Incremental prognostic value of molecular subtype", subtitle = "Repeated 5x5 CV; paired patient bootstrap; n=331, events=153",
       x = "Clinical + subtype minus clinical-only", y = NULL)
ggsave(file.path(figure_dir, "step6_incremental_prognostic_value.pdf"), p, width = 7, height = 4.6, device = cairo_pdf)
ggsave(file.path(figure_dir, "step6_incremental_prognostic_value.png"), p, width = 7, height = 4.6, dpi = 400)

forest <- cox_study[Scope == "Full follow-up" & Target == "C3" & Reference == "C1"]
if (nrow(forest)) {
  pooled <- cox_meta[Scope == "Full follow-up" & Target == "C3" & Reference == "C1"]
  fp <- ggplot(forest, aes(HR, reorder(Study, HR))) + geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
    geom_errorbar(aes(xmin = Lower95, xmax = Upper95), orientation = "y", width = .18) + geom_point(size = 2.4, color = "#CC79A7") +
    scale_x_log10() + theme_classic(base_size = 11) +
    labs(title = "Study-specific adjusted C3 versus C1 survival effects",
         subtitle = sprintf("Random-effects sensitivity HR %.2f (95%% CI %.2f–%.2f); one-stage model remains primary", pooled$HR, pooled$Lower95, pooled$Upper95),
         x = "Adjusted hazard ratio (log scale)", y = NULL)
  ggsave(file.path(figure_dir, "step6_subtype_cox_meta_forest.pdf"), fp, width = 7, height = 4.6, device = cairo_pdf)
  ggsave(file.path(figure_dir, "step6_subtype_cox_meta_forest.png"), fp, width = 7, height = 4.6, dpi = 400)
}

required <- c("step6_study_adjusted_clinical_phenotype_contrasts.csv", "step6_clinical_phenotype_availability.csv",
              "step6_clinical_phenotype_study_specific_effects.csv", "step6_clinical_phenotype_random_effects_meta.csv",
              "step6_incremental_prognostic_oof_predictions.csv", "step6_incremental_prognostic_performance.csv",
              "step6_incremental_prognostic_calibration.csv", "step6_incremental_prognostic_likelihood_ratio.csv",
              "step6_subtype_cox_study_specific.csv", "step6_subtype_cox_random_effects_meta.csv",
              "step6_subtype_cox_meta_leave_one_study_out.csv")
gate <- data.table(
  Check = c("all_targeted_tables_nonempty", "four_clinical_endpoints", "same_population_incremental_models",
            "repeated_cv_and_bootstrap", "time_specific_metrics_and_calibration", "study_specific_cox_and_meta",
            "early_24_month_non_ph_sensitivity", "one_stage_model_remains_primary"),
  Pass = c(all(file.exists(file.path(table_dir, required)) & file.info(file.path(table_dir, required))$size > 0),
           uniqueN(phenotype$Endpoint) == 4, all(oof[, uniqueN(Sample), by = Model]$V1 == 331),
           nrow(fold_def) == 5 * 331 && boot_n >= 500, all(c("AUC6", "AUC12", "AUC24", "IBS3to36") %chin% perf$Metric) && nrow(calibration) == 6,
           nrow(cox_study) > 0 && nrow(cox_meta) > 0, any(lrt$Scope == "Early 24-month restricted follow-up"), TRUE),
  Detail = c(paste(required, collapse = ";"), paste(unique(phenotype$Endpoint), collapse = ";"), "Clinical and ClinicalPlusSubtype use identical 331 patients",
             "5x5 CV; 500 paired patient bootstraps", "6/12/24-month AUC and IPCW calibration; 3-36-month IBS",
             sprintf("%d study-specific rows; %d random-effects summaries", nrow(cox_study), nrow(cox_meta)),
             "24-month restricted follow-up likelihood-ratio sensitivity", "Primary one-stage Study-stratified Cox retained in core Step6")
)
fwrite(gate, file.path(table_dir, "step6_targeted_upgrade_release_gate_status.csv"))
if (!all(gate$Pass)) stop("Step6 targeted upgrade release gate failed")
cat("Step6 targeted upgrades complete.\n")
