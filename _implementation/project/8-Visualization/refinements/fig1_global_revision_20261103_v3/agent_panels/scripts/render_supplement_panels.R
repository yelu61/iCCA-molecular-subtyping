#!/usr/bin/env Rscript
# R-only, display-only reconstruction from locked source tables.
# Contract: quantitative grid; source values are preserved, no models are fitted.
# Run in the project root with the supplied R 4.5 executable and --vanilla.
options(stringsAsFactors = FALSE, warn = 2)
library(grid)
stopifnot(requireNamespace("svglite", quietly = TRUE),
          requireNamespace("ragg", quietly = TRUE),
          requireNamespace("systemfonts", quietly = TRUE))
root <- normalizePath(getwd())
out <- file.path(root, "8-Visualization/refinements/fig1_global_revision_20261103_v3/agent_panels")
paths <- c(
  feature = "1-Clustering/publication/tables/step1_kmeans_feature_number_sensitivity.csv",
  seeds = "1-Clustering/publication/tables/step1_kmeans_k3_seed_sensitivity.csv",
  loso = "1-Clustering/publication/tables/step1_kmeans_leave_one_study_out.csv",
  refit = "1-Clustering/publication/tables/step1_locked_vs_refit_summary.csv",
  clinical = "6-ClinicalSurvival/publication/tables/step6_study_adjusted_clinical_phenotype_contrasts.csv",
  study_colors = "0-Data/qc/study_color_dictionary.csv",
  subtype_colors = "docs/subtype_color_dictionary.csv")
dat <- lapply(paths, function(p) read.csv(file.path(root, p), check.names = FALSE))
for (nm in names(paths)) file.copy(file.path(root, paths[[nm]]),
  file.path(out, "source_data", basename(paths[[nm]])), overwrite = TRUE)
source_manifest <- data.frame(Role = names(paths), Path = unname(paths),
  MD5 = unname(tools::md5sum(file.path(root, paths))),
  Rows = vapply(dat, nrow, integer(1)))
write.csv(source_manifest, file.path(out, "review/input_manifest.csv"), row.names = FALSE)

# Direct study labels identify cohorts. Neutral mark types identify analysis roles,
# so the centrally loaded study colors need not be assigned redundantly to marks.
stopifnot(all(dat$loso$HeldOutStudy %in% dat$study_colors$Study))
stopifnot(nrow(dat$feature) == 3L, identical(dat$feature$Features, c(1000L, 2000L, 5000L)),
          nrow(dat$seeds) == 100L, identical(dat$seeds$Seed, 1:100),
          nrow(dat$loso) == 7L, all(dat$loso$TrainingN + dat$loso$HeldOutN == 586L),
          nrow(dat$refit) == 3L, all(dat$refit$N == 586L))
expected_endpoints <- c("Tumor size", "Advanced stage", "Vascular invasion", "Observed large-duct morphology")
clinical_all <- dat$clinical
clinical <- clinical_all[clinical_all$Endpoint %in% expected_endpoints, ]
stopifnot(nrow(clinical) == 12L, all(table(clinical$Endpoint) == 3L),
          all(clinical$Lower95 < clinical$Estimate), all(clinical$Estimate < clinical$Upper95),
          all(clinical$N + clinical$MissingN == 497L),
          max(abs(p.adjust(clinical_all$PValue, method = "BH") - clinical_all$FDR)) < 1e-12)
write.csv(clinical, file.path(out, "source_data/clinical_display_rows.csv"), row.names = FALSE)

family <- "Arial"
font_file <- systemfonts::match_fonts(family)$path[[1]]
stopifnot(grepl("Arial", font_file, ignore.case = TRUE))
ink <- "#30343A"; soft <- "#71777E"; rule <- "#D8DCE0"; neutral <- "#525A64"
mm <- function(x) unit(x, "mm")
txt <- function(label, x, y, size = 7, face = "plain", col = ink, just = "left") {
  grid.text(label, x = mm(x), y = mm(y), just = just,
    gp = gpar(fontfamily = family, fontsize = size, fontface = face, col = col, lineheight = 1.08))
}
seg <- function(x0, y0, x1, y1, col = rule, lwd = .6, lty = 1) {
  grid.segments(mm(x0), mm(y0), mm(x1), mm(y1), gp = gpar(col = col, lwd = lwd, lty = lty))
}
point <- function(x, y, pch = 16, size = 1.35, col = ink, fill = ink) {
  grid.points(mm(x), mm(y), pch = pch, size = mm(size), gp = gpar(col = col, fill = fill, lwd = .7))
}
start_page <- function(w, h) {
  grid.newpage()
  pushViewport(viewport(x = 0, y = 0, width = mm(w), height = mm(h), just = c("left", "bottom")))
}
end_page <- function() popViewport()
export <- function(draw, name, w, h) {
  stem <- file.path(out, "figures", name)
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w/25.4, height = h/25.4, family = family)
  draw(); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w/25.4, height = h/25.4, bg = "white")
  draw(); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "mm", res = 600, background = "white")
  draw(); dev.off()
  ragg::agg_png(file.path(out, "review", paste0(name, "_preview.png")), width = w, height = h,
                units = "mm", res = 300, background = "white")
  draw(); dev.off()
}

draw_stability <- function() {
  start_page(174, 65)
  txt("Clustering sensitivity relative to stored k-means assignments", 1.3, 62, 8.1, "bold")
  txt("k = 3; n = 586 tumors. The corrected expression matrix is retained.", 1.3, 57.8, 7, col = soft)
  txt("Feature-number sensitivity", 1.3, 52.7, 7.5, "bold")
  txt("100 random seeds", 47.5, 52.7, 7.5, "bold")
  txt("Fixed-representation leave-study-out sensitivity", 89.7, 52.7, 7.4, "bold")
  sy <- function(z) 16 + (z - .65)/(.37) * 28
  ticks <- c(.7, .8, .9, 1)
  for (z in ticks) {
    seg(9, sy(z), 41, sy(z), col = rule, lwd = .45)
    txt(sprintf("%.1f", z), 6.7, sy(z), 7, just = "right")
    seg(50, sy(z), 83, sy(z), col = rule, lwd = .45)
  }
  txt("ARI", 1.3, 47.4, 7, col = soft)
  fx <- c(12, 25, 38)
  fy <- sy(dat$feature$ARIvsLockedKMeans)
  grid.lines(mm(fx), mm(fy), gp = gpar(col = neutral, lwd = .8))
  point(fx, fy, size = 1.8)
  for (i in seq_along(fx)) {
    txt(sprintf("%.3f", dat$feature$ARIvsLockedKMeans[i]), fx[i], fy[i] + 2.7, 7, just = "center")
    txt(format(dat$feature$Features[i], big.mark = ","), fx[i], 11.6, 7, just = "center")
    txt(paste0("\u0394n = ", dat$feature$ChangedSamples[i]), fx[i], 7.8, 6.8, col = soft, just = "center")
  }
  txt("Number of variable genes", 25, 3.1, 7, just = "center")
  sx <- 50 + (dat$seeds$Seed - 1)/99 * 33
  point(sx, sy(dat$seeds$ARIvsHistoricalKMeans), size = .7)
  for (s in c(1, 50, 100)) {
    xx <- 50 + (s - 1)/99 * 33
    txt(s, xx, 11.6, 7, just = "center")
  }
  txt("100/100 at ARI = 1.000", 66.5, 36.2, 7, just = "center")
  txt("25 starts per seed", 66.5, 31.9, 7, col = soft, just = "center")
  txt("Random seed", 66.5, 7.8, 7, just = "center")
  lx <- function(z) 119 + (z - .65)/.37 * 53
  for (z in ticks) seg(lx(z), 16, lx(z), 45, col = rule, lwd = .45)
  short <- c("2019Ahn" = "Ahn", "2022Dong" = "Dong", "2023Deng" = "Deng",
             "2023Jeon" = "Jeon", "2023Lee" = "Lee", "2024Fan" = "Fan", "TCGA-CHOL" = "TCGA-CHOL")
  yy <- seq(44, 17.3, length.out = 7)
  for (i in seq_len(nrow(dat$loso))) {
    d <- dat$loso[i, ]
    txt(paste0(short[[d$HeldOutStudy]], " (", d$HeldOutN, ")"), 89.7, yy[i], 7)
    point(lx(d$ARIAllVsHistoricalKMeans), yy[i] + .9, pch = 16, size = 1.4)
    point(lx(d$ARITrainingVsHistoricalKMeans), yy[i], pch = 1, size = 1.5, fill = "white")
    point(lx(d$ARIHeldOutVsHistoricalKMeans), yy[i] - .9, pch = 17, size = 1.7, col = neutral)
  }
  txt("Omitted study (n)", 89.7, 48.4, 6.8, col = soft)
  for (z in ticks) txt(sprintf("%.1f", z), lx(z), 12.4, 7, just = "center")
  txt("ARI in each comparison set", 146, 8.6, 7, just = "center")
  point(100, 3.3, 16, 1.4); txt("All", 102.6, 3.3, 6.9)
  point(119.5, 3.3, 1, 1.5, fill = "white"); txt("Retained", 122.1, 3.3, 6.9)
  point(149, 3.3, 17, 1.7, col = neutral); txt("Omitted", 151.6, 3.3, 6.9)
  end_page()
}

draw_refit <- function() {
  start_page(174, 26)
  txt("Exact refit of stored assignments", 1.3, 23, 8.1, "bold")
  cols <- c(1.3, 91, 112, 137, 165)
  txt("Assignment", cols[1], 17.9, 7.2, "bold")
  for (i in 2:5) txt(c("N", "ARI", "Agreement", "Changed")[i-1], cols[i], 17.9, 7.2, "bold", just = "center")
  seg(1.3, 15.5, 172, 15.5, col = rule, lwd = .7)
  labs <- c("NMF, k = 3 (50 runs)", "k-means, k = 3", "Governed six-state assignment")
  yy <- c(12.2, 7.5, 2.8)
  for (i in 1:3) {
    d <- dat$refit[i, ]
    txt(labs[i], cols[1], yy[i], 7.4)
    txt(d$N, cols[2], yy[i], 7.4, just = "center")
    txt(sprintf("%.3f", d$ARI), cols[3], yy[i], 7.4, just = "center")
    txt(sprintf("%.0f%%", 100*d$ExactAgreementAfterAlignment), cols[4], yy[i], 7.4, just = "center")
    txt(d$ChangedSamples, cols[5], yy[i], 7.4, just = "center")
  }
  end_page()
}

fmt_fdr <- function(z) if (z < .001) formatC(z, format = "e", digits = 1) else sprintf("%.3f", z)
draw_endpoint <- function(d, x0, y0, title) {
  pushViewport(viewport(x = mm(x0), y = mm(y0), width = mm(85), height = mm(32), just = c("left", "bottom")))
  stopifnot(length(unique(d$N)) == 1L, length(unique(d$Studies)) == 1L)
  txt(paste0(title, " (n = ", d$N[1], "; ", d$Studies[1], " studies)"), 0, 30.4, 7.7, "bold")
  is_or <- identical(unique(d$EffectScale), "Odds ratio")
  xlim <- if (is_or) c(-1, 2) else c(-1, 2.5)
  trans <- function(v) if (is_or) log10(v) else v
  xx <- function(v) 18 + (trans(v) - xlim[1])/diff(xlim)*25
  if (is_or) stopifnot(all(d$Lower95 > 0))
  vals <- if (is_or) c(.1, 1, 10, 100) else c(-1, 0, 1, 2)
  ref <- if (is_or) 1 else 0
  txt("Contrast", 0, 25.1, 7, "bold")
  txt(if (is_or) "Odds ratio" else "Mean difference (cm)", 30.5, 25.1, if (is_or) 7 else 6.8, "bold", just = "center")
  txt("Estimate [95% CI]", 46.2, 25.1, 7, "bold")
  txt("FDR", 84.5, 25.1, 7, "bold", just = "right")
  seg(xx(ref), 7, xx(ref), 22.4, col = soft, lwd = .6, lty = 2)
  yy <- c(20.2, 15.2, 10.2)
  for (i in 1:3) {
    r <- d[i, ]
    txt(paste(r$Target, "vs", r$Reference), 0, yy[i], 7.1)
    seg(xx(r$Lower95), yy[i], xx(r$Upper95), yy[i], col = ink, lwd = .8)
    point(xx(r$Estimate), yy[i], pch = 16, size = 1.5)
    txt(sprintf("%.2f [%.2f, %.2f]", r$Estimate, r$Lower95, r$Upper95), 46.2, yy[i], 7)
    txt(fmt_fdr(r$FDR), 84.5, yy[i], 7, just = "right")
  }
  seg(18, 6.6, 43, 6.6, col = ink, lwd = .6)
  for (v in vals) {
    seg(xx(v), 6.6, xx(v), 5.7, col = ink, lwd = .6)
    txt(format(v, trim = TRUE), xx(v), 2.7, 7, just = "center")
  }
  if (is_or) txt("log scale", 59, 2.7, 6.8, col = soft)
  popViewport()
}
draw_clinical <- function() {
  start_page(174, 75)
  txt("Study-adjusted clinical phenotype contrasts", 1.3, 72, 8.1, "bold")
  title <- c("Tumor size", "Stage III/IV", "Vascular invasion", "Observed large duct")
  origins <- list(c(1.3, 36.5), c(88, 36.5), c(1.3, 2), c(88, 2))
  for (i in seq_along(expected_endpoints)) {
    d <- clinical[clinical$Endpoint == expected_endpoints[i], ]
    d <- d[match(c("C2/C1", "C3/C1", "C3/C2"), paste(d$Target, d$Reference, sep = "/")), ]
    draw_endpoint(d, origins[[i]][1], origins[[i]][2], title[i])
  }
  end_page()
}

export(draw_stability, "SFig1_clustering_sensitivity", 174, 65)
export(draw_refit, "SFig1_exact_refit", 174, 26)
export(draw_clinical, "SFig1_adjusted_clinical", 174, 75)

manifest <- data.frame(Panel = c("Clustering sensitivity", "Exact refit", "Adjusted clinical phenotypes"),
  Stem = c("SFig1_clustering_sensitivity", "SFig1_exact_refit", "SFig1_adjusted_clinical"),
  Width_mm = 174, Height_mm = c(65, 26, 75), PNG_dpi = 600, Font = family,
  Font_path = font_file, Rows_displayed = c("3 feature settings + 100 seeds + 21 study/set ARIs", "3", "12"))
write.csv(manifest, file.path(out, "review/output_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "review/sessionInfo.txt"))
cat("Exported three panels in PDF, SVG and 600-dpi PNG. No upstream data changed.\n")
