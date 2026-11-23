#!/usr/bin/env Rscript
# ==============================================================================
# Script: plot_figure6.R
# Purpose: Generate Figure 6 panels (A-G) and full composite matching publication layout.
# Data Inputs: 8-Visualization/data/figure6/
# Outputs: 8-Visualization/figures/figure6/
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(png)
})

# Pinned font handling for macOS / sandboxed environment
if (!"Arial" %in% names(grDevices::pdfFonts())) grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])

# Paths & CLI Arguments
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default_val) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) args[idx + 1] else default_val
}

data_cands <- c(get_arg("--data-dir", ""), "figures/data/figure6", "8-Visualization/data/figure6", "data/figure6", "../data/figure6")
data_dir <- data_cands[dir.exists(data_cands) & data_cands != ""][1]
if (is.na(data_dir)) data_dir <- "8-Visualization/data/figure6"
data_dir <- normalizePath(data_dir)

out_cands <- c(get_arg("--out-dir", ""), if (dir.exists("8-Visualization/figures")) "8-Visualization/figures/figure6" else "", "figures/outputs/figure6", "outputs/figure6")
out_dir <- out_cands[out_cands != ""][1]
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Color palettes
dict_cands <- c("governance/data/subtype_label_color_dictionary.csv", "data/subtype_label_color_dictionary.csv", "../data/subtype_label_color_dictionary.csv")
dict_file <- dict_cands[file.exists(dict_cands)][1]
if (!is.na(dict_file)) {
  subtype_pal_df <- read.csv(dict_file, stringsAsFactors = FALSE)
  subtype_colors <- setNames(subtype_pal_df$Color, subtype_pal_df$FinalSubtype)
} else {
  subtype_colors <- c(C1 = "#0072B2", C2 = "#E69F00", C3 = "#CC79A7", C1_C2 = "#738859", C1_C3 = "#6675AC", C2_C3 = "#D98C53")
}
pal <- subtype_colors[c("C1", "C2", "C3")]

ink <- "#263E4B"; muted <- "#7C8992"; light <- "#E5EBEF"

mm <- function(x) unit(x, "mm")

txt <- function(s, x, y, size = 6.0, bold = FALSE, col = ink, just = "left", rot = 0, lineheight = 1) {
  grid.text(s, mm(x), mm(y), just = just, rot = rot,
            gp = gpar(fontfamily = "Arial", fontsize = size,
                      fontface = if (bold) "bold" else "plain", col = col, lineheight = lineheight))
}

rect_box <- function(x, y, w, h, fill, col = NA, lwd = 0.3, just = c("center", "center")) {
  grid.rect(mm(x), mm(y), mm(w), mm(h), just = just, gp = gpar(fill = fill, col = col, lwd = lwd))
}

line_seg <- function(x, y, col = "#A8B1B8", lwd = 0.5, lty = 1) {
  grid.lines(mm(x), mm(y), gp = gpar(col = col, lwd = lwd, lty = lty))
}

point_dot <- function(x, y, fill = ink, col = fill, size = 1.4, pch = 16, lwd = 0.65) {
  grid.points(mm(x), mm(y), pch = pch, size = mm(size), gp = gpar(col = col, fill = fill, lwd = lwd))
}

W <- 85; H <- 76
start <- function(w, h) { W <<- w; H <<- h; grid.newpage(); pushViewport(viewport(xscale = c(0, w), yscale = c(h, 0))) }
txt_grid <- function(s, x, y, size = 6.3, bold = FALSE, col = ink, just = 'left', italic = FALSE) {
  g <- textGrob(s, unit(x, 'native'), unit(y, 'native'), just = c(just, 'center'),
                gp = gpar(fontfamily = 'Arial', fontsize = size, fontface = if (italic) 'italic' else if (bold) 'bold' else 'plain', col = col, lineheight = 1.03))
  grid.draw(g)
}
ln_grid <- function(x, y, col = muted, lwd = .5, lty = 1) grid.lines(unit(x, 'native'), unit(y, 'native'), gp = gpar(col = col, lwd = lwd, lty = lty))
dot_grid <- function(x, y, col = ink, pch = 16, size = 1.25) grid.points(unit(x, 'native'), unit(y, 'native'), pch = pch, size = unit(size, 'mm'), gp = gpar(col = col, fill = col, lwd = .6))
axisx_grid <- function(xx, ticks, y, labels = as.character(ticks), size = 5.7) {
  ln_grid(xx(range(ticks)), c(y, y))
  for (i in seq_along(ticks)) {
    ln_grid(rep(xx(ticks[i]), 2), c(y, y + .6))
    txt_grid(labels[i], xx(ticks[i]), y + 2, size, just = 'center')
  }
}
interval_grid <- function(r, xx, y, col = ink, pch = 16, hr = FALSE) {
  e <- if (hr) r$HR else r$Estimate
  v <- xx(c(r$Lower95, e, r$Upper95))
  ln_grid(v[c(1, 3)], c(y, y), col, .75)
  dot_grid(v[2], y, col, pch)
}

# ------------------------------------------------------------------------------
# Panel A: Prognostic Information (331 Patients)
# ------------------------------------------------------------------------------
cat("Rendering Panel A...
")
draw_panelA <- function() {
  ah <- fread(file.path(data_dir, 'Figure6A_clinical_information_1.csv'))
  a <- fread(file.path(data_dir, 'Figure6A_clinical_information_2.csv'))
  
  start(85, 76); on.exit(popViewport())
  txt_grid('Prognostic information', 2, 3.2, 7.3, TRUE)
  txt_grid('331 patients', 83, 3.2, 6, just = 'right')
  txt_grid('24-month OOF prediction', 2, 8.5, 6.5, TRUE)
  txt_grid('AUC', 43, 13, 6.2, TRUE, just = 'center')
  txt_grid('Increment', 70, 13, 6.2, TRUE, just = 'center')
  
  ax <- function(v) 32 + (v - .5) / .35 * 18
  dx <- function(v) 60 + (v + .02) / .20 * 15
  ln_grid(rep(ax(.5), 2), c(16, 32), light, .6, 2)
  ln_grid(rep(dx(0), 2), c(16, 37.5), light, .6, 2)
  
  for (i in 1:3) {
    y <- c(18.5, 24.5, 30.5)[i]
    co <- c(muted, pal['C3'], ink)[i]
    cmp <- c('M0', 'MC3', 'MAll')[i]
    r <- a[Comparison == cmp]
    txt_grid(c('Clinical', '+ C3', '+ C1/C2/C3')[i], 2, y)
    interval_grid(r, ax, y, co, c(1, 16, 15)[i])
    txt_grid(sprintf('%.3f', r$Estimate), 57, y, 6.2, just = 'right')
    if (i == 1) {
      txt_grid('Ref.', 83, y, 6, just = 'right')
    } else {
      r_diff <- a[Comparison == paste0(cmp, '-M0')]
      interval_grid(r_diff, dx, y, co, c(1, 16, 15)[i])
      txt_grid(sprintf('%+.3f', r_diff$Estimate), 83, y, 6.2, just = 'right')
    }
  }
  r_full <- a[Comparison == 'MAll-MC3']
  txt_grid('Full split vs C3', 2, 36.5, 6)
  interval_grid(r_full, dx, 36.5, muted, 18)
  txt_grid(sprintf('%+.3f', r_full$Estimate), 83, 36.5, 6, just = 'right')
  
  axisx_grid(ax, c(.5, .8), 39, c('0.5', '0.8'))
  axisx_grid(dx, c(0, .15), 39, c('0', '0.15'))
  ln_grid(c(2, 83), c(44, 44), light)
  
  txt_grid('Adjusted associations', 2, 48, 6.5, TRUE)
  dot_grid(52, 48, ink, 16)
  txt_grid('12 mo', 54, 48, 5.8)
  dot_grid(68, 48, ink, 1)
  txt_grid('24 mo', 70, 48, 5.8)
  
  hx <- function(v) 37 + (log(v) - log(.25)) / (log(4.5) - log(.25)) * 34
  ln_grid(rep(hx(1), 2), c(51, 70), light, .6, 2)
  keys <- c('Age (per 10 years)', 'Male vs female', 'Stage III/IV vs I/II', 'C2 vs C1', 'C3 vs C1')
  for (i in 1:5) {
    y <- 52 + (i - 1) * 4.3
    txt_grid(keys[i], 2, y, 6)
    co <- if (i <= 3) ink else pal[c('C2', 'C3')[i - 3]]
    for (j in 1:2) {
      r_sub <- ah[Contrast == keys[i] & TimeMonths == c(12, 24)[j]]
      yy <- y + c(-.8, .8)[j]
      interval_grid(r_sub, hx, yy, co, c(16, 1)[j], TRUE)
    }
  }
  axisx_grid(hx, c(.25, 1, 4), 71.5, c('0.25', '1', '4'))
  txt_grid('HR', 83, 73.5, 5.7, just = 'right')
}

out_panelA_png <- file.path(out_dir, 'Figure_6_PanelA.png')
out_panelA_pdf <- file.path(out_dir, 'Figure_6_PanelA.pdf')
agg_png(out_panelA_png, width = 85, height = 76, units = 'mm', res = 600, background = 'white')
draw_panelA()
dev.off()
grDevices::cairo_pdf(out_panelA_pdf, width = 85 / 25.4, height = 76 / 25.4, family = 'Arial')
draw_panelA()
dev.off()

# ------------------------------------------------------------------------------
# Panel B: Beyond Measured Duct Increment (79 Patients)
# ------------------------------------------------------------------------------
cat("Rendering Panel B...
")
draw_panelB <- function() {
  b <- fread(file.path(data_dir, 'Figure6B_measured_duct_increment.csv'))
  
  start(85, 76); on.exit(popViewport())
  txt_grid('Beyond measured duct', 2, 3.2, 7.3, TRUE)
  txt_grid('79 patients', 83, 3.2, 6, just = 'right')
  
  dot_grid(35, 8.5, pal['C3'], 16)
  txt_grid('Primary', 37, 8.5, 5.8)
  dot_grid(55, 8.5, ink, 2)
  txt_grid('Time-varying', 57, 8.5, 5.8)
  
  txt_grid('Change in 24-month AUC', 2, 13.5, 6.5, TRUE)
  dx <- function(v) 43 + (v + .05) / .35 * 29
  ln_grid(rep(dx(0), 2), c(16, 36), light, .6, 2)
  
  for (i in 1:3) {
    y <- 19 + (i - 1) * 7
    txt_grid(c('Add C3 to clinical + duct', 'C3 vs duct (with clinical)', 'Add duct to clinical + C3')[i], 2, y, 6.0, bold = (i == 1))
    for (j in 1:2) {
      r_b <- b[Comparison == c('M3-M1', 'M2-M1', 'M3-M2')[i] & Family == c('RidgePH', 'RidgeTimeVarying')[j]]
      yy <- y + c(-1.1, 1.1)[j]
      interval_grid(r_b, dx, yy, c(pal['C3'], ink)[j], c(16, 2)[j])
      txt_grid(sprintf('%+.3f', r_b$Estimate), 83, yy, 5.8, just = 'right')
    }
  }
  axisx_grid(dx, c(0, .15, .30), 37, c('0', '0.15', '0.30'))
  ln_grid(c(2, 83), c(42, 42), light)
  
  txt_grid('Absolute 24-month AUC', 2, 46, 6.5, TRUE)
  ax <- function(v) 36 + (v - .2) / .6 * 36
  ln_grid(rep(ax(.5), 2), c(49, 69), light, .6, 2)
  
  for (i in 1:4) {
    y <- 51 + (i - 1) * 5.7
    txt_grid(c('Clinical', 'Clinical + duct', 'Clinical + C3', 'Clinical + duct + C3')[i], 2, y, 6.0)
    for (j in 1:2) {
      r_b <- b[Comparison == paste0('M', i - 1) & Family == c('RidgePH', 'RidgeTimeVarying')[j]]
      yy <- y + c(-.9, .9)[j]
      interval_grid(r_b, ax, yy, c(pal['C3'], ink)[j], c(16, 2)[j])
      if (j == 1) txt_grid(sprintf('%.3f', r_b$Estimate), 83, y, 6, just = 'right')
    }
  }
  axisx_grid(ax, c(.2, .5, .8), 71.5, c('0.2', '0.5', '0.8'))
}

out_panelB_png <- file.path(out_dir, 'Figure_6_PanelB.png')
out_panelB_pdf <- file.path(out_dir, 'Figure_6_PanelB.pdf')
agg_png(out_panelB_png, width = 85, height = 76, units = 'mm', res = 600, background = 'white')
draw_panelB()
dev.off()
grDevices::cairo_pdf(out_panelB_pdf, width = 85 / 25.4, height = 76 / 25.4, family = 'Arial')
draw_panelB()
dev.off()
# ------------------------------------------------------------------------------
# Panel C: Overall survival by observed duct morphology and C3 (image94)
# ------------------------------------------------------------------------------
cat("Rendering Panel C...\n")
draw_panelC <- function() {
  curves <- fread(file.path(data_dir, "km_curves.csv"))
  risk <- fread(file.path(data_dir, "km_risk_table.csv"))
  counts <- risk[Month == 0, .(Duct, Molecular, N = AtRisk)]
  events <- curves[, .(Events = sum(Deaths)), by = .(Duct, Molecular)]
  counts <- merge(counts, events, by = c("Duct", "Molecular"))
  
  cols <- c("Non-C3" = "#596873", "C3" = pal[["C3"]])
  
  grid.newpage()
  pushViewport(viewport(xscale = c(0, 174), yscale = c(113, 0), x = mm(0), y = mm(0),
                        width = mm(174), height = mm(113), just = c("left", "bottom")))
  
  txt_v <- function(s, x, y, size = 7.2, bold = FALSE, just = "left", color = "#263641") {
    grid.text(s, unit(x, "native"), unit(y, "native"), just = just,
              gp = gpar(fontfamily = "Arial", fontsize = size, fontface = if (bold) "bold" else "plain", col = color))
  }
  ln_v <- function(x, y, col = "#A8B1B8", lwd = 0.65, lty = 1) {
    grid.lines(unit(x, "native"), unit(y, "native"), gp = gpar(col = col, lwd = lwd, lty = lty))
  }
  pt_v <- function(x, y, col, pch = 16, sz = 1.8) {
    grid.points(unit(x, "native"), unit(y, "native"), pch = pch, size = mm(sz), gp = gpar(col = col, fill = col))
  }
  step_v <- function(d, field) {
    v <- d[[field]]
    data.table(t = c(0, rep(d$Time, each = 2)), v = c(1, as.vector(rbind(c(1, head(v, -1)), v))))
  }
  
  draw_km_sub <- function(duct, x0, x1, y0, y1, risk_y, show_y = TRUE) {
    xmax <- max(risk$Month)
    xx <- function(v) x0 + (x1 - x0) * v / xmax
    yy <- function(v) y1 - (y1 - y0) * v
    for (v in c(0, 0.5, 1)) {
      ln_v(c(x0, x1), rep(yy(v), 2), "#E7EBEE", 0.5)
      if (show_y) txt_v(sprintf("%.1f", v), x0 - 2, yy(v), 6.5, just = "right", color = muted)
    }
    for (m in c("Non-C3", "C3")) {
      d <- curves[Duct == duct & Molecular == m]
      setorder(d, Time)
      lo <- step_v(d, "Lower95")
      hi <- step_v(d, "Upper95")
      ok <- is.finite(lo$v) & is.finite(hi$v)
      grid.polygon(unit(c(xx(lo$t[ok]), rev(xx(hi$t[ok]))), "native"),
                   unit(c(yy(lo$v[ok]), rev(yy(hi$v[ok]))), "native"),
                   gp = gpar(fill = adjustcolor(cols[m], alpha.f = 0.15), col = NA))
      a_step <- step_v(d, "Survival")
      ln_v(xx(a_step$t), yy(a_step$v), cols[m], 1.25)
      c_pts <- d[Censors > 0]
      if (nrow(c_pts)) pt_v(xx(c_pts$Time), yy(c_pts$Survival), cols[m], pch = 3, sz = 1.3)
    }
    ln_v(c(x0, x1), c(y1, y1), "#65747E")
    ln_v(c(x0, x0), c(y0, y1), "#65747E")
    for (t in risk[!duplicated(Month), Month]) {
      ln_v(rep(xx(t), 2), c(y1, y1 + 1), "#65747E")
      txt_v(as.character(t), xx(t), y1 + 4, 6.5, just = "center", color = muted)
    }
    txt_v("OS (months)", mean(c(x0, x1)), y1 + 8.5, 6.8, just = "center")
    txt_v("At risk", x0 - 15, risk_y, 6.8, bold = TRUE)
    for (j in seq_along(cols)) {
      m <- names(cols)[j]
      y_pos <- risk_y + 4.5 + j * 4.2
      txt_v(m, x0 - 15, y_pos, 6.5, color = cols[m])
      a_risk <- risk[Duct == duct & Molecular == m]
      for (i in seq_len(nrow(a_risk))) txt_v(as.character(a_risk$AtRisk[i]), xx(a_risk$Month[i]), y_pos, 6.5, just = "center", color = cols[m])
    }
  }
  
  txt_v("Overall survival by observed duct morphology and C3", 4, 6, 8.5, TRUE)
  txt_v("79 patients / 45 deaths | descriptive Kaplan–Meier curves with pointwise 95% CIs", 4, 11.5, 6.8, color = muted)
  
  for (j in 1:2) {
    duct <- c("Small duct", "Large duct")[j]
    x0 <- c(22, 108)[j]; x1 <- c(78, 164)[j]
    a_cnt <- counts[Duct == duct]
    txt_v(sprintf("%s  |  n = %d", duct, sum(a_cnt$N)), x0, 19, 7.5, TRUE)
    txt_v(sprintf("C3: %d (%d deaths)   Non-C3: %d (%d deaths)",
                  a_cnt[Molecular == "C3", N], a_cnt[Molecular == "C3", Events],
                  a_cnt[Molecular == "Non-C3", N], a_cnt[Molecular == "Non-C3", Events]),
          x0, 24.5, 6.5, color = muted)
    draw_km_sub(duct, x0, x1, 31, 72, 85, TRUE)
  }
  
  grid.text("Survival", unit(6, "native"), unit(51.5, "native"), rot = 90, just = "center",
            gp = gpar(fontfamily = "Arial", fontsize = 7.0))
  txt_v("C3", 80, 106, 7.0, color = cols["C3"])
  ln_v(c(73, 78), c(106, 106), cols["C3"], 1.2)
  txt_v("Non-C3", 102, 106, 7.0, color = cols["Non-C3"])
  ln_v(c(95, 100), c(106, 106), cols["Non-C3"], 1.2)
  
  popViewport()
}

out_panelC_png <- file.path(out_dir, "Figure_6_PanelC.png")
out_panelC_pdf <- file.path(out_dir, "Figure_6_PanelC.pdf")
agg_png(out_panelC_png, width = 174, height = 113, units = "mm", res = 600, background = "white")
draw_panelC()
dev.off()
grDevices::cairo_pdf(out_panelC_pdf, width = 174 / 25.4, height = 113 / 25.4, family = "Arial")
draw_panelC()
dev.off()

# ------------------------------------------------------------------------------
# Panel D: Three-Subtype Identification (image90)
# ------------------------------------------------------------------------------
cat("Rendering Panel D...\n")
draw_panelD <- function() {
  perf <- fread(file.path(data_dir, "joint_performance_ci.csv"))
  models <- c("Anchor3", "CompactWorkflow", "JointWorkflow", "FullWorkflow")
  mlabel <- c(Anchor3 = "RNA anchors", CompactWorkflow = "Compact", JointWorkflow = "Joint", FullWorkflow = "Full workflow")
  role <- c(Anchor3 = "#85919A", CompactWorkflow = "#426C85", JointWorkflow = "#243E50", FullWorkflow = "#6F8C86")
  metrics <- c("BalancedAccuracy", "C1Recall", "C2Recall", "C3Recall")
  metricnames <- c("Balanced accuracy", "C1 recall", "C2 recall", "C3 recall")
  
  d <- perf[Model %in% models & Metric %in% metrics]
  d[, Model := factor(Model, levels = rev(models))]
  d[, Metric := factor(Metric, levels = metrics, labels = metricnames)]
  
  p <- ggplot(d, aes(100 * Estimate, Model, color = Model)) +
    geom_vline(xintercept = c(0, 50, 100), color = light, linewidth = 0.3) +
    geom_errorbar(aes(xmin = 100 * Lower95, xmax = 100 * Upper95), orientation = "y", width = 0.15, linewidth = 0.45) +
    geom_point(size = 1.5) +
    geom_text(aes(label = sprintf("%.1f", 100 * Estimate)), nudge_y = 0.28, size = 6.5 / .pt, show.legend = FALSE) +
    facet_wrap(~Metric, nrow = 1) +
    scale_y_discrete(labels = mlabel) +
    scale_x_continuous(limits = c(-5, 110), breaks = c(0, 50, 100)) +
    scale_color_manual(values = role, guide = "none") +
    labs(title = "Three-subtype identification across four workflows",
         subtitle = "Same 497 tumors · shared outer folds · conditional 95% intervals",
         x = "OOF performance (%)", y = NULL) +
    theme_classic(base_size = 7.5, base_family = "Arial") +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(size = 8.5, face = "bold", color = ink),
          plot.subtitle = element_text(size = 7.0, color = ink),
          strip.background = element_blank(), strip.text = element_text(size = 7.2, face = "bold"),
          axis.text = element_text(size = 6.8, color = ink), axis.title = element_text(size = 7.2),
          panel.spacing = unit(2.5, "mm"), plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))
  p
}

out_panelD_png <- file.path(out_dir, "Figure_6_PanelD.png")
out_panelD_pdf <- file.path(out_dir, "Figure_6_PanelD.pdf")
p_d <- draw_panelD()
agg_png(out_panelD_png, width = 174, height = 77, units = "mm", res = 600, background = "white")
print(p_d)
dev.off()
grDevices::cairo_pdf(out_panelD_pdf, width = 174 / 25.4, height = 77 / 25.4, family = "Arial")
print(p_d)
dev.off()

# ------------------------------------------------------------------------------
# Panel E: Predicted Drug Patterns Across All Three Subtypes (image92)
# ------------------------------------------------------------------------------
cat("Rendering Panel E...\n")
draw_panelE <- function() {
  ev <- fread(file.path(data_dir, "candidate_evidence_summary.csv"))
  z <- fread(file.path(data_dir, "Figure6F_drug_heterogeneity.csv"))
  ord <- ev$Drug
  short <- setNames(sub("_[0-9]+$", "", ord), ord)
  short["Sepantronium bromide_1941"] <- "Sepantronium"
  labels <- setNames(paste0(short[ord], "\n", ev$FunctionalClass), ord)
  
  z[, Drug := factor(Drug, levels = rev(ord))]
  cls <- c("C1", "C2", "C3")
  z[, CandidateGroup := factor(paste(CandidateSubtype, "nominees"), levels = paste(cls, "nominees"))]
  summ <- z[, .(N = .N, Median = median(Z), Q1 = quantile(Z, 0.25), Q3 = quantile(Z, 0.75)), by = .(Drug, CandidateGroup, FinalSubtype)]
  
  hm <- ggplot(z, aes(PatientOrder, Drug, fill = pmax(-2, pmin(2, Z)))) +
    geom_tile(width = 1, height = 1) +
    facet_grid(CandidateGroup ~ FinalSubtype, scales = "free", space = "free",
               labeller = labeller(FinalSubtype = c(C1 = "C1 · 112", C2 = "C2 · 223", C3 = "C3 · 162"))) +
    scale_y_discrete(labels = labels, expand = expansion(add = 0.5)) +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_fill_gradient2(low = "#397CA6", mid = "white", high = "#B66D7D", limits = c(-2, 2),
                         breaks = c(-2, 0, 2), labels = c("≤−2", "0", "≥2"), name = "Score (SD)",
                         guide = guide_colorbar(title.position = "top", barwidth = unit(18, "mm"), barheight = unit(1.8, "mm"))) +
    labs(x = NULL, y = NULL, title = "Individual tumors") +
    theme_classic(base_size = 7.0, base_family = "Arial") +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6.2, color = ink), strip.text.y = element_blank(),
          strip.background = element_blank(), strip.text.x = element_text(size = 6.8, face = "bold"),
          plot.title = element_text(size = 7.5, face = "bold"), panel.spacing = unit(1, "mm"),
          legend.position = "bottom", legend.title = element_text(size = 6.0), legend.text = element_text(size = 5.8),
          plot.margin = margin(1, 1, 1, 1, "mm"))
  
  q <- ggplot(summ, aes(Median, Drug, color = FinalSubtype)) +
    geom_vline(xintercept = 0, linetype = 2, color = "#C0C8CD", linewidth = 0.3) +
    geom_errorbar(aes(xmin = Q1, xmax = Q3), orientation = "y", position = position_dodge(width = 0.65), width = 0, linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.65), size = 1.1) +
    facet_grid(CandidateGroup ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    scale_x_continuous(limits = c(-2.3, 2.3), breaks = c(-2, 0, 2)) +
    scale_color_manual(values = pal, guide = "none") +
    labs(x = "Score (SD)", y = NULL, title = "Median / IQR") +
    theme_classic(base_size = 7.0, base_family = "Arial") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          strip.text.y = element_text(size = 6.2, face = "bold"), strip.background = element_blank(),
          plot.title = element_text(size = 7.5, face = "bold"), axis.text.x = element_text(size = 6.0),
          axis.title.x = element_text(size = 6.5), panel.spacing = unit(1, "mm"),
          plot.margin = margin(1, 1, 1, 1, "mm"))
  
  p <- hm + q + plot_layout(widths = c(3.4, 1)) +
    plot_annotation(title = "Predicted drug patterns across all three subtypes",
                    subtitle = "Nine nominees · 497 tumors · within-drug SD · lower = greater predicted sensitivity",
                    theme = theme(plot.title = element_text(size = 8.5, face = "bold", family = "Arial"),
                                  plot.subtitle = element_text(size = 7.0, family = "Arial")))
  p
}

out_panelE_png <- file.path(out_dir, "Figure_6_PanelE.png")
out_panelE_pdf <- file.path(out_dir, "Figure_6_PanelE.pdf")
p_e <- draw_panelE()
agg_png(out_panelE_png, width = 174, height = 112, units = "mm", res = 600, background = "white")
print(p_e)
dev.off()
grDevices::cairo_pdf(out_panelE_pdf, width = 174 / 25.4, height = 112 / 25.4, family = "Arial")
print(p_e)
dev.off()

# ------------------------------------------------------------------------------
# Panel F: Candidate Evidence Summary Matrix (image93)
# ------------------------------------------------------------------------------
cat("Rendering Panel F...\n")
draw_panelF <- function() {
  ev <- fread(file.path(data_dir, "candidate_evidence_summary.csv"))
  ord <- ev$Drug
  short <- setNames(sub("_[0-9]+$", "", ord), ord)
  short["Sepantronium bromide_1941"] <- "Sepantronium"
  
  progshort <- c("Unmapped", "JAK / STAT", "TGF-beta receptor", "B-cell receptor", "TP53 regulation", "Mitosis (broad)", "Checkpoints", "MEK / MAPK", "MAPK1 / MAPK3")
  ev[, ProgramLabel := progshort]
  
  compactexpr <- function(s) {
    s <- gsub("Up:", "↑ ", s, fixed = TRUE)
    s <- gsub("Down:", "↓ ", s, fixed = TRUE)
    s <- gsub("NA:", "Missing: ", s, fixed = TRUE)
    s <- gsub("NS:", "NS ", s, fixed = TRUE)
    s <- gsub("; ", "\n", s, fixed = TRUE)
    s <- gsub("Missing: .*", "Missing", s)
    s
  }
  ev[, RNAPlot := vapply(RNA, compactexpr, character(1))]
  ev[, ProteinPlot := vapply(Protein, compactexpr, character(1))]
  
  grid.newpage()
  pushViewport(viewport(xscale = c(0, 174), yscale = c(146, 0), x = mm(0), y = mm(0),
                        width = mm(174), height = mm(146), just = c("left", "bottom")))
  
  tx <- function(s, x, y, size = 7, bold = FALSE, col = ink, just = "left") {
    gg <- textGrob(s, unit(x, "native"), unit(y, "native"), just = c(just, "center"),
                   gp = gpar(fontfamily = "Arial", fontsize = size, fontface = if (bold) "bold" else "plain", col = col, lineheight = 1.05))
    grid.draw(gg)
  }
  
  tx("Candidate evidence: rationale, consistency and missing links", 2, 4, 8.5, TRUE)
  tx("Original nine nominees · each evidence type remains separate", 2, 8.5, 7.0, col = muted)
  
  xs <- c(2, 29, 65, 91, 119, 138, 156)
  tx("Nominee", xs[1], 16, 7.2, TRUE)
  tx("Program context", xs[2], 16, 7.2, TRUE)
  tx("RNA", xs[3], 16, 7.2, TRUE)
  tx("Protein", xs[4], 16, 7.2, TRUE)
  tx("Meta CI", xs[5], 16, 7.0, TRUE)
  tx("Top 3", xs[6], 16, 7.0, TRUE)
  tx("PRISM", xs[7], 16, 7.0, TRUE)
  
  tx("Focal subtype", xs[2], 20.5, 6.5, col = muted)
  tx("Focal subtype", xs[3], 20.5, 6.5, col = muted)
  tx("Focal subtype", xs[4], 20.5, 6.5, col = muted)
  tx("below 0", xs[5], 20.5, 6.5, col = muted)
  tx("LOSO", xs[6], 20.5, 6.5, col = muted)
  tx("C1/C2/C3", xs[7], 20.5, 6.5, col = muted)
  
  for (i in seq_len(nrow(ev))) {
    r <- ev[i]
    y <- 28 + (i - 1) * 11.9
    grid.rect(x = unit(87, "native"), y = unit(y, "native"), width = unit(172, "native"), height = unit(11.4, "native"),
              gp = gpar(fill = if (i %% 2) "#F5F7F8" else "white", col = NA))
    tx(short[r$Drug], xs[1], y - 1.5, 7.0, TRUE)
    tx(paste(r$CandidateSubtype, "nominee"), xs[1], y + 2.5, 6.5, FALSE, pal[r$CandidateSubtype])
    tx(r$ProgramLabel, xs[2], y - 2.5, 6.8)
    ps <- if (is.na(r$ProgramNES)) "Not mapped" else sprintf("%s NES %.2f\nFDR %.3g", c(Up = "↑", Down = "↓", NS = "NS")[r$ProgramState], r$ProgramNES, r$ProgramFDR)
    tx(ps, xs[2], y + 2, 6.5, FALSE, if (r$ProgramState == "Down") "#447D96" else ink)
    tx(r$RNAPlot, xs[3], y, 6.5)
    tx(r$ProteinPlot, xs[4], y, 6.5)
    tx(sprintf("%d / 2", r$PairwiseCIBelowZero), xs[5] + 5, y - 2, 7.2, TRUE, just = "center")
    tx(sprintf("I² %.0f–%.0f%%", r$MinI2, r$MaxI2), xs[5] + 5, y + 2.5, 6.0, just = "center", col = muted)
    tx(sprintf("%d / %d", r$Top3Retained, r$Omissions), xs[6] + 4, y, 7.2, TRUE, just = "center")
    tx(sprintf("%d / %d / %d", r$PRISMMeasuredC1, r$PRISMMeasuredC2, r$PRISMMeasuredC3), xs[7] + 6, y - 1.5, 7.0, TRUE, just = "center")
    tx(paste0(r$PRISMMeasured, " mapped / 3"), xs[7] + 6, y + 2.5, 6.0, just = "center", col = muted)
  }
  
  tx("↑/↓: source FDR < 0.05; NS: not significant. Missing or unmapped is not a negative assay.", 2, 135, 6.5, col = muted)
  tx("LOSO: downstream study omissions. PRISM: measured lines with projected labels; coverage is not validation.", 2, 139.5, 6.5, col = muted)
  tx("Program mapping uses the full Figure2 GSEA background; no combined evidence score.", 2, 144, 6.5, col = muted)
  
  popViewport()
}

out_panelF_png <- file.path(out_dir, "Figure_6_PanelF.png")
out_panelF_pdf <- file.path(out_dir, "Figure_6_PanelF.pdf")
agg_png(out_panelF_png, width = 174, height = 146, units = "mm", res = 600, background = "white")
draw_panelF()
dev.off()
grDevices::cairo_pdf(out_panelF_pdf, width = 174 / 25.4, height = 146 / 25.4, family = "Arial")
draw_panelF()
dev.off()

# ------------------------------------------------------------------------------
# Panel G: Ulixertinib Case Study (image91)
# ------------------------------------------------------------------------------
cat("Rendering Panel G...\n")
draw_panelG <- function() {
  z <- fread(file.path(data_dir, "Figure6F_drug_heterogeneity.csv"))
  hh <- fread(file.path(data_dir, "Figure6H_C3_hypothesis.csv"))
  hmta <- fread(file.path(data_dir, "SFig6F_drug_uncertainty.csv"))[Drug == "Ulixertinib_1908"]
  pr <- fread(file.path(data_dir, "candidate_PRISM_all_cells.csv"))[Drug == "Ulixertinib_1908"]
  
  hz <- z[as.character(Drug) == "Ulixertinib_1908"]
  hh[, Y := ifelse(Gene == "MAPK1", 2, 1) + c(C1 = 0.21, C2 = 0, C3 = -0.21)[Subtype]]
  hh[, Significant := FDR < 0.05]
  hh[, Assay := factor(Assay, levels = c("RNA", "Protein"))]
  
  x <- ggplot(hz, aes(FinalSubtype, Z, color = FinalSubtype)) +
    geom_point(position = position_jitter(width = 0.16, height = 0, seed = 20261103), size = 0.4, alpha = 0.3) +
    geom_boxplot(width = 0.45, outlier.shape = NA, fill = NA, linewidth = 0.4) +
    scale_color_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = c(C1 = "C1 · 112", C2 = "C2 · 223", C3 = "C3 · 162")) +
    labs(title = "Patient predictions", subtitle = "497 tumors · lower relative score", x = NULL, y = "Relative score (SD)") +
    theme_classic(base_size = 7.0, base_family = "Arial") +
    theme(plot.title = element_text(size = 7.5, face = "bold"), plot.subtitle = element_text(size = 6.5),
          axis.text = element_text(size = 6.5, color = ink), axis.title = element_text(size = 6.8))
  
  hmta[, Comp := factor(paste("C3 vs", Other), levels = c("C3 vs C2", "C3 vs C1"))]
  m <- ggplot(hmta, aes(Effect, Comp)) +
    geom_vline(xintercept = 0, linetype = 2, color = muted, linewidth = 0.3) +
    geom_errorbar(aes(xmin = Lower, xmax = Upper), orientation = "y", width = 0.1, color = pal["C3"], linewidth = 0.55) +
    geom_point(color = pal["C3"], size = 1.6) +
    geom_text(aes(label = sprintf("%.2f [%.2f, %.2f]  I² %.0f%%", Effect, Lower, Upper, I2_percent)), nudge_y = 0.25, size = 6.2 / .pt) +
    scale_x_continuous(limits = c(min(hmta$Lower) - 0.2, 0.15)) +
    labs(title = "Cross-study prediction effects", subtitle = "Random effects · 95% intervals · Top 3 retained 7/7",
         x = "Predicted score difference (SD)", y = NULL) +
    theme_classic(base_size = 7.0, base_family = "Arial") +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(size = 7.5, face = "bold"), plot.subtitle = element_text(size = 6.5),
          axis.text = element_text(size = 6.5, color = ink), axis.title = element_text(size = 6.8))
  
  y <- ggplot(hh, aes(Effect, Y, color = Subtype, shape = Significant)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, color = "#ACB8C0") +
    geom_point(size = 1.5) +
    facet_wrap(~Assay, nrow = 1, scales = "free_x", labeller = labeller(Assay = c(RNA = "RNA · log2 effect", Protein = "Protein · SD effect"))) +
    scale_color_manual(values = pal) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), labels = c("FALSE" = "FDR ≥ 0.05", "TRUE" = "FDR < 0.05")) +
    scale_y_continuous(breaks = c(2, 1), labels = c("MAPK1", "MAPK3"), limits = c(0.5, 2.5)) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    labs(title = "Both mapped genes", subtitle = "C3 MAPK program: NES 1.10, FDR 0.497",
         x = "Study-adjusted expression difference", y = NULL) +
    theme_classic(base_size = 7.0, base_family = "Arial") +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(size = 7.5, face = "bold"), plot.subtitle = element_text(size = 6.5),
          legend.key.width = unit(2, "mm"), legend.text = element_text(size = 6.0),
          strip.text = element_text(size = 6.5, face = "bold"), legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
          axis.text = element_text(size = 6.5, color = ink), axis.title = element_text(size = 6.8))
  
  pr[, LineLabel := factor(paste0(CellLineName, " · ", Prediction), levels = paste0(CellLineName, " · ", Prediction))]
  pp <- ggplot(pr, aes(AUC, LineLabel, color = Prediction)) +
    geom_point(size = 1.8) +
    geom_text(aes(label = sprintf("%.3f", AUC)), nudge_y = 0.24, size = 6.5 / .pt, show.legend = FALSE) +
    scale_color_manual(values = pal, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(title = "Observed PRISM responses", subtitle = "Three lines · projected subtype labels",
         x = "AUC (lower = greater sensitivity)", y = NULL) +
    theme_classic(base_size = 7.0, base_family = "Arial") +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(size = 7.5, face = "bold"), plot.subtitle = element_text(size = 6.5),
          axis.text.y = element_text(size = 6.5, color = ink), axis.text.x = element_text(size = 6.5, color = ink),
          axis.title = element_text(size = 6.8))
  
  p <- (x + m + plot_layout(widths = c(1, 1.45))) / (y + pp + plot_layout(widths = c(1.45, 1))) +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(title = "Ulixertinib: a C3 hypothesis with unresolved evidence",
                    subtitle = "MAPK3 abundance supports follow-up; pathway enrichment and three-line assays do not confirm the prediction",
                    theme = theme(plot.title = element_text(size = 8.5, face = "bold", family = "Arial"),
                                  plot.subtitle = element_text(size = 7.0, family = "Arial")))
  p
}

out_panelG_png <- file.path(out_dir, "Figure_6_PanelG.png")
out_panelG_pdf <- file.path(out_dir, "Figure_6_PanelG.pdf")
p_g <- draw_panelG()
agg_png(out_panelG_png, width = 174, height = 130, units = "mm", res = 600, background = "white")
print(p_g)
dev.off()
grDevices::cairo_pdf(out_panelG_pdf, width = 174 / 25.4, height = 130 / 25.4, family = "Arial")
print(p_g)
dev.off()

# ------------------------------------------------------------------------------
# Full Composite: Figure 6 (180 mm x 270 mm) matching PPT V2 (Slide 18)
# ------------------------------------------------------------------------------
cat("Assembling Full Figure 6 Composite...\n")
composite_png <- file.path(out_dir, "Figure_6_composite.png")
composite_pdf <- file.path(out_dir, "Figure_6_composite.pdf")

render_composite <- function() {
  grid.newpage()
  # 180 mm x 270 mm canvas
  pushViewport(viewport(x = 0, y = 0, width = mm(180), height = mm(270), just = c("left", "bottom")))
  grid.rect(gp = gpar(fill = "white", col = NA))
  
  draw_raster_panel <- function(png_path, x_mm, y_top_mm, w_mm, h_mm) {
    if (file.exists(png_path)) {
      img <- readPNG(png_path)
      y_bottom_mm <- 270 - y_top_mm - h_mm
      grid.raster(img, x = mm(x_mm + w_mm / 2), y = mm(y_bottom_mm + h_mm / 2),
                  width = mm(w_mm), height = mm(h_mm))
    }
  }
  
  draw_panel_letter <- function(letter, x_mm, y_top_mm) {
    y_from_bottom <- 270 - y_top_mm
    grid.text(letter, x = mm(x_mm), y = mm(y_from_bottom), just = c("left", "top"),
              gp = gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold", col = "#111111"))
  }
  
  # Slide 18 exact layout coordinates:
  # Row 1:
  # Panel A: [0.9, 11.3, 47.9, 42.9] | Letter at [0.0, 6.4]
  draw_raster_panel(out_panelA_png, 0.9, 11.3, 47.9, 42.9)
  draw_panel_letter("A", 0.0, 6.4)
  
  # Panel B: [53.0, 11.3, 50.3, 45.1] | Letter at [49.1, 6.2]
  draw_raster_panel(out_panelB_png, 53.0, 11.3, 50.3, 45.1)
  draw_panel_letter("B", 49.1, 6.2)
  
  # Panel C: [106.0, 10.7, 74.0, 48.0] | Letter at [104.2, 7.0]
  draw_raster_panel(out_panelC_png, 106.0, 10.7, 74.0, 48.0)
  draw_panel_letter("C", 104.2, 7.0)
  
  # Row 2:
  # Panel D: [0.3, 67.5, 80.8, 35.7] | Letter at [0.3, 56.9]
  draw_raster_panel(out_panelD_png, 0.3, 67.5, 80.8, 35.7)
  draw_panel_letter("D", 0.3, 56.9)
  
  # Panel E: [81.1, 59.0, 98.0, 63.0] | Letter at [78.2, 57.8]
  draw_raster_panel(out_panelE_png, 81.1, 59.0, 98.0, 63.0)
  draw_panel_letter("E", 78.2, 57.8)
  
  # Row 3:
  # Panel F: [2.0, 119.6, 84.4, 70.7] | Letter at [0.3, 114.4]
  draw_raster_panel(out_panelF_png, 2.0, 119.6, 84.4, 70.7)
  draw_panel_letter("F", 0.3, 114.4)
  
  # Panel G: [89.8, 122.5, 90.2, 67.4] | Letter at [85.1, 121.4]
  draw_raster_panel(out_panelG_png, 89.8, 122.5, 90.2, 67.4)
  draw_panel_letter("G", 85.1, 121.4)
  
  popViewport()
}

agg_png(composite_png, width = 180, height = 270, units = "mm", res = 300, background = "white")
render_composite()
dev.off()

grDevices::cairo_pdf(composite_pdf, width = 180 / 25.4, height = 270 / 25.4, family = "Arial")
render_composite()
dev.off()

cat("Figure 6 rendered successfully:\n")
cat("  - Individual panels:", out_dir, "\n")
cat("  - Composite PDF:", composite_pdf, "\n")
cat("  - Composite PNG:", composite_png, "\n")
