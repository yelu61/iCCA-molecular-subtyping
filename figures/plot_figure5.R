#!/usr/bin/env Rscript
# ==============================================================================
# Figure 5: Subtype Classifier Development, Performance, and Transfer
# Standalone, decoupled reproduction script.
#
# Generates Panels A-G and assembled Figure 5.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(svglite)
  library(ragg)
})

# ------------------------------------------------------------------------------
# 1. Configuration & CLI Arguments
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default_val) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) args[idx + 1] else default_val
}

data_cands <- c(get_arg("--data-dir", ""), "figures/data/figure5", "8-Visualization/data/figure5", "data/figure5", "../data/figure5")
DATA_DIR <- data_cands[dir.exists(data_cands) & data_cands != ""][1]
if (is.na(DATA_DIR)) DATA_DIR <- "8-Visualization/data/figure5"
DATA_DIR <- normalizePath(DATA_DIR)

out_cands <- c(get_arg("--out-dir", ""), if (dir.exists("8-Visualization/figures")) "8-Visualization/figures/figure5" else "", "figures/outputs/figure5", "outputs/figure5")
OUT_DIR <- out_cands[out_cands != ""][1]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)

# Color palettes
FIG5_SUBTYPE <- c(C1 = "#0072B2", C2 = "#E69F00", C3 = "#CC79A7")
FIG5_INK     <- "#293F50"
EFFECT_RAMP  <- c("#377EB8", "#FFFFFF", "#CD5258")

# Helper to read from DATA_DIR
rdata <- function(filename) {
  p <- file.path(DATA_DIR, filename)
  if (!file.exists(p)) stop(paste("Required data file not found:", p))
  fread(p)
}

# Typography and base theme
v_theme <- function() {
  theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 7.5, face = "bold", colour = FIG5_INK),
      plot.subtitle = element_text(size = 6, colour = FIG5_INK),
      plot.title.position = "plot",
      plot.caption = element_text(size = 5.8, hjust = 0, colour = "grey30", lineheight = 1.1),
      plot.caption.position = "plot",
      axis.title = element_text(size = 6, colour = FIG5_INK),
      axis.text = element_text(size = 5.8, colour = "grey20"),
      axis.line = element_line(linewidth = 0.3, colour = "grey45"),
      axis.ticks = element_line(linewidth = 0.25),
      legend.text = element_text(size = 5.8),
      legend.title = element_text(size = 6),
      legend.key.size = unit(3, "mm"),
      legend.position = "top",
      strip.background = element_rect(fill = "#EAF0F4", colour = NA),
      strip.text = element_text(size = 7, face = "bold"),
      panel.spacing = unit(3, "mm"),
      plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm")
    )
}

# Grid primitives for schematic and custom layouts
canvas <- function(w, h) {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, w), yscale = c(0, h)))
}

vt <- function(text, x, y, size = 6, bold = FALSE, just = "left", col = FIG5_INK, italic = FALSE) {
  grid.text(
    text, unit(x, "native"), unit(y, "native"), just = just,
    gp = gpar(
      fontfamily = "Arial", fontsize = size,
      fontface = if (italic) "italic" else if (bold) "bold" else "plain",
      col = col, lineheight = 1.05
    )
  )
}

vr <- function(x, y, w, h, col, edge = NA) {
  grid.rect(unit(x, "native"), unit(y, "native"), unit(w, "native"), unit(h, "native"),
            gp = gpar(fill = col, col = edge, lwd = 0.35))
}

vl <- function(x, y, col = FIG5_INK, lwd = 0.5, arrow = FALSE) {
  grid.lines(unit(x, "native"), unit(y, "native"),
             gp = gpar(col = col, lwd = lwd),
             arrow = if (arrow) grid::arrow(length = unit(1, "mm"), type = "closed") else NULL)
}

vp <- function(x, y, col = FIG5_INK, size = 0.8, pch = 16) {
  grid.points(unit(x, "native"), unit(y, "native"), pch = pch,
              size = unit(size, "mm"), gp = gpar(col = col))
}

# ------------------------------------------------------------------------------
# 2. Panel Drawing Functions
# ------------------------------------------------------------------------------

# Panel A: Workflow Schematic
draw_panel_A <- function() {
  canvas(174, 32)
  vr(14, 16, 26, 29, "#F1F5F8", "#DCE4EA")
  vt("Data", 14, 27, 7, TRUE, "centre")
  vt("497 principal tumors", 14, 21, 6, FALSE, "centre")
  vt("Expression + labels", 14, 16.5, 5.8, FALSE, "centre")
  for (j in 1:3) {
    vr(6 + (j - 1) * 8, 10.5, 7, 2.7, FIG5_SUBTYPE[j])
    vt(paste0("C", j), 6 + (j - 1) * 8, 6.5, 5.8, FALSE, "centre")
  }
  
  vr(77, 23, 87, 14, "#F1F5F8", "#DCE4EA")
  vt("Training-fold evaluation", 77, 27.5, 7, TRUE, "centre")
  vt("Study correction -> limma / RF / LASSO -> consensus -> SVM", 77, 22.5, 5.7, FALSE, "centre")
  vt("Outer 5-fold x 5 repeats; inner 4-fold", 77, 18, 5.7, FALSE, "centre")
  
  vr(77, 7.5, 87, 13, "#F1F5F8", "#DCE4EA")
  vt("All 497 tumors -> 114-gene consensus SVM", 77, 10.5, 6.8, TRUE, "centre")
  vt("Separate full-development fit", 77, 5, 5.8, FALSE, "centre")
  
  vr(150, 23, 45, 14, "#F1F5F8", "#DCE4EA")
  vt("OOF & study-held-out transfer", 150, 25.5, 6, FALSE, "centre")
  vt("Training-only gene selection", 150, 20, 5.8, FALSE, "centre")
  
  vr(150, 7.5, 45, 13, "#F1F5F8", "#DCE4EA")
  vt("Fixed-model projection", 150, 10.5, 6.3, TRUE, "centre")
  vt("506 microarray samples", 150, 5, 5.8, FALSE, "centre")
  
  for (y in c(23, 7.5)) {
    vl(c(27.5, 33), c(y, y), arrow = TRUE)
    vl(c(121, 127), c(y, y), arrow = TRUE)
  }
  popViewport()
}

# Panel B: Consensus Gene Evidence (114 genes)
draw_panel_B <- function() {
  d <- rdata("B_subtype_gene_order.csv")
  canvas(72, 76)
  vt("Consensus gene evidence", 1, 73.5, 7.5, TRUE)
  x0 <- 11
  cw <- 58 / 114
  d[, X := x0 + (DisplayOrder - 0.5) * cw + (SubtypeGroup - 1) * 0.8]
  groups <- d[, .(N = .N, L = min(X) - cw / 2, R = max(X) + cw / 2), by = SubtypeGroup]
  for (j in 1:3) {
    g <- groups[SubtypeGroup == j]
    vr((g$L + g$R) / 2, 66.5, g$R - g$L, 5, adjustcolor(FIG5_SUBTYPE[j], alpha.f = 0.15))
    vt(sprintf("C%d | %d", j, g$N), (g$L + g$R) / 2, 66.5, 6.4, TRUE, "centre")
  }
  
  vt("Method support", x0, 60, 6.4, TRUE)
  for (j in 1:3) {
    y <- 56 - (j - 1) * 3
    vt(c("limma", "RF", "LASSO")[j], x0 - 1, y, 6, FALSE, "right")
    for (i in 1:114) {
      vr(d$X[i], y, cw - 0.03, 2.4, if (d[[c("LimmaFull", "RFFull", "LASSOFull")[j]]][i]) FIG5_INK else "#E7EDF1")
    }
  }
  
  vt("Selection frequency (25 fits)", x0, 44.5, 6.4, TRUE)
  for (v in c(0, 0.5, 1)) {
    y <- 26 + 14 * v
    vl(c(x0, 71), c(y, y), "#E3E8EC")
    vt(sprintf("%d%%", 100 * v), x0 - 1, y, 5.8, FALSE, "right")
  }
  for (i in 1:114) {
    y <- 26 + 14 * d$ConsensusFraction[i]
    vl(rep(d$X[i], 2), c(26, y), "#A7BCCB", 0.4)
    vp(d$X[i], y, size = 0.45)
  }
  
  vt("Mean expression", x0, 21.5, 6.4, TRUE)
  ramp <- circlize::colorRamp2(c(-1.5, 0, 1.5), c("#377EB8", "#FFFFFF", "#CD5258"))
  for (j in 1:3) {
    y <- 17.5 - (j - 1) * 3
    vt(paste0("C", j), x0 - 1, y, 6, TRUE, "right", FIG5_SUBTYPE[j])
    for (i in 1:114) {
      vr(d$X[i], y, cw - 0.03, 2.4, ramp(d[[paste0("C", j, "MeanTrainingZ")]][i]))
    }
  }
  
  vr(2, 4, 1.3, 1.3, FIG5_INK)
  vt("Selected", 3.5, 4, 5.5)
  vr(16, 4, 1.3, 1.3, "#E7EDF1")
  vt("Not selected", 17.5, 4, 5.5)
  vt("Mean z", 39, 4, 5.5)
  vt("-1.5", 51, 4, 5.5, FALSE, "right")
  for (i in 1:30) {
    vr(52 + (i - 0.5) * 0.42, 4, 0.42, 1.8, ramp(-1.5 + (i - 0.5) * 3 / 30))
  }
  vt("1.5", 65, 4, 5.5)
  popViewport()
}

# Panel C: Subtype-Representative Genes Ranking
draw_panel_C <- function() {
  d <- rdata("C_subtype_ranking_all114.csv")
  canvas(98, 76)
  vt("Subtype-representative genes", 1, 73.5, 7.5, TRUE)
  
  for (j in 1:3) {
    cl <- names(FIG5_SUBTYPE)[j]
    z <- d[Class == cl]
    setorder(z, WithinClassRank)
    top <- z[Labelled == TRUE]
    x0 <- 8 + (j - 1) * 31
    x <- function(r) x0 + (r - 1) / 113 * 26
    y <- function(v) 30 + (v + 2.5) / 5 * 30
    
    vt(cl, x0, 68.5, 7.3, TRUE, col = FIG5_SUBTYPE[j])
    vt(sprintf("%d with delta z > 0", sum(z$HigherThanBoth)), x0, 64.8, 5.8)
    
    vl(c(x0, x0 + 26), rep(y(0), 2), "#BAC5CE", 0.4)
    vl(c(x0, x0), c(30, 60), "#7D8A92", 0.4)
    vl(c(x0, x0 + 26), c(30, 30), "#7D8A92", 0.4)
    
    vl(x(z$WithinClassRank), y(z$ExpressionMargin), "#D0D8DF", 0.65)
    vp(x(z$WithinClassRank), y(z$ExpressionMargin), "#D0D8DF", 0.55)
    pos <- z[HigherThanBoth == TRUE]
    vp(x(pos$WithinClassRank), y(pos$ExpressionMargin), FIG5_SUBTYPE[j], 0.75)
    
    for (v in c(-2, 0, 2)) {
      vl(c(x0 - 0.5, x0), rep(y(v), 2), "#7D8A92", 0.4)
      if (j == 1) vt(v, x0 - 1, y(v), 5.5, FALSE, "right")
    }
    
    ly <- seq(59, 44, length.out = 6)
    for (i in 1:6) {
      vl(c(x(top$WithinClassRank[i]), x0 + 12), c(y(top$ExpressionMargin[i]), ly[i]), "#7D91A0", 0.35)
      vt(top$Gene[i], x0 + 12.3, ly[i], 5.8, italic = TRUE)
    }
    
    for (r in c(1, 60, 114)) vt(r, x(r), 27.8, 5.5, FALSE, "centre")
    vt("Rank", x0 + 13, 25, 5.8, FALSE, "centre")
    
    t0 <- 1 + (j - 1) * 32.3
    xs <- t0 + c(0, 14, 20, 25, 30)
    heads <- c("Gene", "delta z", "Sel.%", "LRA", "RF")
    for (k in 1:5) vt(heads[k], xs[k], 21, 5.5, TRUE, if (k == 1) "left" else "centre")
    vl(c(t0, t0 + 30.7), c(19.5, 19.5), "#CCD6DD", 0.4)
    for (i in 1:6) {
      yy <- 17.8 - (i - 1) * 2.9
      if (i %% 2) vr(t0 + 15.3, yy, 30.7, 2.7, "#F1F5F7")
      vt(top$Gene[i], xs[1], yy, 5.5)
      vt(sprintf("%.2f", top$ExpressionMargin[i]), xs[2], yy, 5.5, FALSE, "centre")
      vt(round(100 * top$ConsensusFraction[i]), xs[3], yy, 5.5, FALSE, "centre")
      for (k in 1:3) {
        vr(xs[4] + (k - 2) * 1.25, yy, 0.95, 0.95, if (top[[c("LimmaFull", "RFFull", "LASSOFull")[k]]][i]) FIG5_INK else "#DCE4EA")
      }
      vt(top$RFRank[i], xs[5], yy, 5.5, FALSE, "centre")
    }
  }
  grid.text("Expression margin (delta z)", unit(2, "native"), unit(45, "native"), rot = 90,
            gp = gpar(fontfamily = "Arial", fontsize = 5.5, col = FIG5_INK))
  popViewport()
}

# Panel D: Internal Performance & Study Transfer
draw_panel_D <- function() {
  cm12 <- rdata("D_confusion.csv")
  st12 <- rdata("D_study_transfer.csv")
  canvas(100, 64)
  vt("Internal performance and study transfer", 1, 61.5, 7.5, TRUE)
  vt("OOF: n=497; BA 98.79%", 1, 56, 6.1, TRUE)
  vt("LOSO: all seven studies", 50, 56, 6.1, TRUE)
  ramp <- circlize::colorRamp2(c(0, 1), c("#F0F5F8", "#4E8DB1"))
  
  for (i in 1:3) for (j in 1:3) {
    z <- cm12[Truth == paste0("C", i) & Prediction == paste0("C", j)]
    x <- 13 + (j - 0.5) * 10
    y <- 18 + (3 - i + 0.5) * 10
    vr(x, y, 9.7, 9.7, ramp(z$RowFraction))
    vt(sprintf("%d\n%.1f%%", z$Freq, 100 * z$RowFraction), x, y, 6, FALSE, "centre",
       if (z$RowFraction > 0.5) "white" else FIG5_INK)
  }
  for (j in 1:3) {
    vt(paste0("C", j), 13 + (j - 0.5) * 10, 15.5, 5.7, FALSE, "centre")
    vt(paste0("C", j), 11, 18 + (3 - j + 0.5) * 10, 5.7, FALSE, "right")
  }
  grid.text("Reference", unit(3, "native"), unit(33, "native"), rot = 90,
            gp = gpar(fontfamily = "Arial", fontsize = 5.7, col = FIG5_INK))
  vt("Predicted", 28, 11, 5.8, FALSE, "centre")
  
  xp <- function(v) 77 + (v - 0.84) / 0.16 * 21
  for (v in c(0.85, 0.9, 0.95, 1)) {
    vl(rep(xp(v), 2), c(14, 51), "#E2E8ED")
    vt(round(100 * v), xp(v), 11.5, 5.5, FALSE, "centre")
  }
  for (i in 1:7) {
    z <- st12[i]
    y <- 49 - (i - 1) * 5.3
    vt(sprintf("%s (%d)", z$Study, z$N), 50, y, 5.6)
    vl(xp(c(z$OOF, z$LOSO)), c(y, y), "#9AAFBF", 0.7)
    vp(xp(z$OOF), y, "#8197A6", 1.1, 1)
    vp(xp(z$LOSO), y, size = 1.1)
  }
  vt("BA (%)", 98, 7.5, 5.6, FALSE, "right")
  vp(52, 4, "#8197A6", 1.1, 1)
  vt("OOF", 54, 4, 5.5)
  vp(66, 4, size = 1.1)
  vt("Held-out", 68, 4, 5.5)
  vt("Transductive centering", 1, 4, 5.5)
  popViewport()
}

# Panel E: External Frozen Programs (GSVA Heatmap)
draw_panel_E <- function() {
  eff <- rdata("program_effects.csv")
  sm  <- rdata("E_summary.csv")
  pr  <- rdata("E_composition.csv")
  
  programs <- c("HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_FATTY_ACID_METABOLISM",
                "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_NOTCH_SIGNALING",
                "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
                "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS",
                "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
  labs <- c("Bile-acid metabolism", "Fatty-acid metabolism", "Oxidative phosphorylation",
            "Notch signaling", "E2F targets", "G2M checkpoint",
            "Hypoxia", "Glycolysis", "EMT")
  
  cohs <- sort(unique(eff$Cohort))
  short <- setNames(c("Andersen", "Sia", "Jus P1", "Jus P2", "Job", "TIGER"), cohs)
  
  eff[, ProgramLabel := factor(Program, levels = rev(programs), labels = rev(labs))]
  eff[, CohortLabel := factor(short[Cohort], levels = unname(short))]
  eff[, Class := factor(Class, levels = c("C1", "C2", "C3"))]
  
  draw_summary_bar <- function() {
    canvas(174, 17)
    vt("Projected counts", 1, 4.2, 5.7, TRUE)
    vt("n; genes / 114", 1, 10, 5.7)
    for (i in 1:6) {
      z <- sm[CohortOrder == i]
      x <- 36 + (i - 1) * 24.6
      vt(short[i], x, 14.4, 6, TRUE, "centre")
      vt(sprintf("%d; %d/114", z$N, z$MatchedFeatures), x, 10, 5.5, FALSE, "centre")
      cum <- 0
      for (j in 1:3) {
        a <- pr[CohortOrder == i & Prediction == paste0("C", j)]
        n <- a$N
        wid <- 22 * a$Fraction
        vr(x - 11 + cum + wid / 2, 4.2, wid, 3, FIG5_SUBTYPE[j])
        vt(n, x - 11 + cum + wid / 2, 4.2, 5.3, FALSE, "centre", if (j == 2) FIG5_INK else "white")
        cum <- cum + wid
      }
    }
    popViewport()
  }
  
  lim <- ceiling(max(abs(eff$G), na.rm = TRUE) * 2) / 2
  heat <- ggplot(eff[Scenario == "Exclude114"], aes(CohortLabel, ProgramLabel, fill = G)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    facet_wrap(~Class, nrow = 1) +
    scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#CD5258", midpoint = 0,
                         limits = c(-lim, lim), name = "Hedges g", na.value = "#B7B7B7",
                         guide = guide_colorbar(display = "rectangles")) +
    labs(title = "Frozen programs after classifier-gene exclusion", x = NULL, y = NULL) +
    v_theme() +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(size = 5.8, angle = 40, hjust = 1),
      axis.text.y = element_text(size = 6.2),
      strip.text = element_text(size = 7, face = "bold"),
      panel.spacing = unit(3, "mm"),
      legend.position = "right", legend.key.height = unit(3, "mm"), legend.key.width = unit(2, "mm")
    )
  
  p_combined <- (wrap_elements(full = grid.grabExpr(draw_summary_bar(), width = 174 / 25.4, height = 17 / 25.4)) / heat) +
    plot_layout(heights = c(17, 67)) +
    plot_annotation(caption = "GSVA; subtype minus others; intervals and sensitivity: Supplementary Table 5", theme = v_theme())
  
  return(p_combined)
}

# Panel F: External Survival KM Curves (Jusakual & Job)
draw_panel_F <- function() {
  curves <- rdata("F_external_survival_curves.csv")
  risk   <- rdata("F_external_survival_at_risk.csv")
  tests  <- rdata("F_external_survival_tests.csv")
  
  ps <- lapply(c("2017Jusakual", "2020Job"), function(st) {
    z <- curves[Study == st]
    rr <- risk[Study == st]
    tt <- tests[Study == st]
    xmax <- max(rr$Time)
    ticks <- seq(0, xmax, by = 48)
    rr <- rr[Time %in% ticks]
    
    band <- rbindlist(lapply(split(z, z$Class), function(v) {
      setorder(v, Time)
      b <- v[rep(seq_len(nrow(v)), each = 2)]
      b[, Time := as.vector(rbind(v$Time, c(v$Time[-1], tail(v$Time, 1))))]
      b
    }))
    
    p <- ggplot(z, aes(Time, Survival, colour = Class)) +
      geom_ribbon(data = band, aes(ymin = Lower95, ymax = Upper95, fill = Class),
                  alpha = 0.08, colour = NA, show.legend = FALSE, na.rm = TRUE) +
      geom_step(linewidth = 0.4) +
      geom_point(data = z[Censored > 0], shape = 3, size = 0.8, stroke = 0.25, show.legend = FALSE) +
      scale_colour_manual(values = FIG5_SUBTYPE) +
      scale_fill_manual(values = FIG5_SUBTYPE) +
      scale_x_continuous(limits = c(0, xmax), breaks = ticks, expand = expansion(mult = c(0.035, 0.045))) +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.01, 0.02))) +
      labs(title = sprintf("%s (n=%d)", if (st == "2017Jusakual") "Jusakual" else "Job", tt$N),
           subtitle = sprintf("Holm P = %.3f", tt$P_Holm), x = "OS (months)", y = "Survival", colour = NULL) +
      v_theme() +
      theme(legend.position = if (st == "2017Jusakual") "top" else "none",
            legend.key.width = unit(2, "mm"), legend.spacing.x = unit(1, "mm"), legend.key.height = unit(1.5, "mm"),
            plot.margin = margin(1, 1, 0, 1, "mm"))
    
    rr[, Class := factor(Class, levels = c("C3", "C2", "C1"))]
    q <- ggplot(rr, aes(Time, Class, label = AtRisk, colour = Class)) +
      geom_text(size = 1.95) +
      scale_colour_manual(values = FIG5_SUBTYPE, guide = "none") +
      scale_x_continuous(limits = c(0, xmax), breaks = ticks, expand = expansion(mult = c(0.035, 0.045))) +
      labs(title = "At risk", x = NULL, y = NULL) +
      v_theme() +
      theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(margin = margin(r = 1.5, unit = "mm")),
            plot.title = element_text(size = 5.8), plot.margin = margin(0, 1, 1, 1, "mm")) +
      coord_cartesian(clip = "off")
    
    p / q + plot_layout(heights = c(3.8, 1))
  })
  
  wrap_plots(ps, ncol = 2) + plot_annotation(title = "External survival", theme = v_theme())
}

# Panel G: Feature Reduction (Candidate Weights & Cost)
draw_panel_G <- function() {
  w12   <- rdata("full_candidate_weights.csv")
  del12 <- rdata("paired_differences.csv")[Metric == "BalancedAccuracy"]
  ch12  <- rdata("choice.csv")
  
  w12[, Gene := factor(Gene, levels = rev(Gene))]
  pw <- ggplot(w12, aes(100 * RelativeWeight, Gene)) +
    geom_col(fill = FIG5_INK, width = 0.6) +
    labs(title = "Fixed 17-gene development candidate",
         subtitle = "Descriptive model weights; fixed panel not independently validated",
         x = "Relative squared weight (%)", y = NULL) +
    v_theme() +
    theme(axis.text.y = element_text(face = "italic", size = 5.8))
  
  draw_cost <- function() {
    canvas(60, 70)
    vt("Cost of feature reduction", 1, 67, 7.5, TRUE)
    vt("Adaptive compact - full workflow", 1, 61.5, 6.1, TRUE)
    vt(sprintf("Delta BA %.2f pp [%.2f, %.2f]", 100 * del12$Estimate, 100 * del12$Lower95, 100 * del12$Upper95), 1, 56.5, 6)
    
    xx <- function(v) 6 + (v + 2.6) / 2.9 * 48
    vl(rep(xx(0), 2), c(38, 51), "#8294A2", 0.5)
    grid.lines(unit(rep(xx(-1), 2), "native"), unit(c(38, 51), "native"), gp = gpar(col = "#8294A2", lwd = 0.6, lty = 2))
    vt("-1 pp target", xx(-1), 52.2, 5.7, FALSE, "centre")
    vl(xx(100 * c(del12$Lower95, del12$Upper95)), rep(44, 2), FIG5_INK, 1.1)
    vp(xx(100 * del12$Estimate), 44, size = 1.8)
    
    for (v in c(-2, -1, 0)) vt(v, xx(v), 35, 5.7, FALSE, "centre")
    vt("BA-loss target not met", 30, 30, 6.3, TRUE, "centre")
    vt("Genes selected in 25 outer fits", 1, 23.8, 6, TRUE)
    vt("Median 23 [10-51]", 1, 19.5, 5.8)
    
    freq <- ch12[, .N, by = K]
    for (i in 1:nrow(freq)) for (j in 1:freq$N[i]) vp(7 + (freq$K[i] - 10) / 41 * 45, 11 + (j - 1) * 1.5, size = 0.75)
    vl(c(7, 52), c(9, 9), "#9CABBA", 0.4)
    for (v in c(10, 30, 51)) vt(v, 7 + (v - 10) / 41 * 45, 6.5, 5.6, FALSE, "centre")
    vt("Fixed 17-gene candidate: S6K/L", 1, 2, 5.6)
    popViewport()
  }
  
  pcost <- wrap_elements(full = grid.grabExpr(draw_cost(), width = 60 / 25.4, height = 70 / 25.4))
  return(list(weights = pw, cost = pcost))
}

# ------------------------------------------------------------------------------
# 3. Execution & Rendering
# ------------------------------------------------------------------------------
if (!"Arial" %in% names(grDevices::pdfFonts())) grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])

save_panel <- function(stem, plot_func_or_obj, w, h) {
  pdf_path <- file.path(OUT_DIR, "panels", paste0(stem, ".pdf"))
  png_path <- file.path(OUT_DIR, "panels", paste0(stem, ".png"))
  
  draw <- function() {
    if (is.function(plot_func_or_obj)) plot_func_or_obj() else print(plot_func_or_obj)
  }
  
  pdf(pdf_path, width = w / 25.4, height = h / 25.4, family = "Arial")
  draw()
  dev.off()
  
  ragg::agg_png(png_path, width = w, height = h, units = "mm", res = 300, background = "white")
  draw()
  dev.off()
  
  cat(sprintf("[Rendered] %s (%.1f x %.1f mm)\n", stem, w, h))
}

cat("Rendering individual panels for Figure 5...\n")
save_panel("Panel_A_workflow", draw_panel_A, 174, 32)
save_panel("Panel_B_consensus_evidence", draw_panel_B, 72, 76)
save_panel("Panel_C_subtype_genes", draw_panel_C, 98, 76)
save_panel("Panel_D_internal_performance", draw_panel_D, 100, 64)
save_panel("Panel_E_external_frozen_programs", draw_panel_E(), 174, 90)
save_panel("Panel_F_external_survival", draw_panel_F(), 110, 70)

g_out <- draw_panel_G()
save_panel("Panel_G_weights", g_out$weights, 90, 70)
save_panel("Panel_G_cost", g_out$cost, 60, 70)

cat("Generating full assembled Figure 5 preview...\n")
composite_png <- file.path(OUT_DIR, "Figure_5_composite.png")
composite_pdf <- file.path(OUT_DIR, "Figure_5_composite.pdf")

panels_pos <- data.table(
  Letter = c("A", "B", "C", "D", "E", "F", "G", ""),
  File   = c("Panel_A_workflow.png", "Panel_B_consensus_evidence.png",
             "Panel_C_subtype_genes.png", "Panel_D_internal_performance.png",
             "Panel_E_external_frozen_programs.png", "Panel_F_external_survival.png",
             "Panel_G_weights.png", "Panel_G_cost.png"),
  X      = c(1.4, 3.9, 87.9, 4.0, 81.4, 3.4, 87.0, 140.6),
  Top    = c(5.8, 58.0, 58.9, 130.6, 127.0, 179.2, 179.2, 179.2),
  W      = c(174.0, 77.9, 88.9, 73.8, 99.3, 81.4, 45.0, 37.9),
  H      = c(48.0, 69.0, 69.0, 47.2, 51.3, 51.9, 49.9, 44.1),
  LetX   = c(1.4, 0.4, 81.7, 0.4, 75.7, 0.9, 78.5, NA),
  LetY   = c(5.8, 59.0, 58.0, 129.1, 127.9, 176.2, 176.3, NA)
)

render_fig5_composite <- function() {
  canvas(180, 270)
  grid.rect(gp = gpar(fill = "white", col = NA))
  for (i in seq_len(nrow(panels_pos))) {
    p_file <- file.path(OUT_DIR, "panels", panels_pos$File[i])
    if (file.exists(p_file)) {
      im <- png::readPNG(p_file)
      grid.raster(
        im,
        x = unit(panels_pos$X[i], "native"),
        y = unit(270 - panels_pos$Top[i], "native"),
        width = unit(panels_pos$W[i], "native"),
        height = unit(panels_pos$H[i], "native"),
        just = c("left", "top")
      )
      if (!is.na(panels_pos$LetX[i]) && nzchar(panels_pos$Letter[i])) {
        vt(panels_pos$Letter[i], panels_pos$LetX[i], 270 - panels_pos$LetY[i], 11, TRUE, col = "#111111")
      }
    }
  }
  popViewport()
}

ragg::agg_png(composite_png, width = 180, height = 270, units = "mm", res = 300, background = "white")
render_fig5_composite()
dev.off()

pdf(composite_pdf, width = 180 / 25.4, height = 270 / 25.4, family = "Arial")
render_fig5_composite()
dev.off()

cat(sprintf("[Complete] All Figure 5 panels and composite generated in %s\n", OUT_DIR))
