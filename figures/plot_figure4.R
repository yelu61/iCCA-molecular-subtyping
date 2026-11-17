#!/usr/bin/env Rscript
# ==============================================================================
# Script: plot_figure4.R
# Purpose: Generate Figure 4 panels (A-G) and full composite matching publication layout.
# Data Inputs: 8-Visualization/data/figure4/
# Outputs: 8-Visualization/figures/figure4/
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
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

data_cands <- c(get_arg("--data-dir", ""), "figures/data/figure4", "8-Visualization/data/figure4", "data/figure4", "../data/figure4")
data_dir <- data_cands[dir.exists(data_cands) & data_cands != ""][1]
if (is.na(data_dir)) data_dir <- "8-Visualization/data/figure4"
data_dir <- normalizePath(data_dir)

out_cands <- c(get_arg("--out-dir", ""), if (dir.exists("8-Visualization/figures")) "8-Visualization/figures/figure4" else "", "figures/outputs/figure4", "outputs/figure4")
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

study_dict_cands <- c("0-Data/qc/study_color_dictionary.csv", "data/study_color_dictionary.csv", "../data/study_color_dictionary.csv")
study_dict_file <- study_dict_cands[file.exists(study_dict_cands)][1]
if (!is.na(study_dict_file)) {
  study_pal_df <- read.csv(study_dict_file, stringsAsFactors = FALSE)
  scol <- if ("ColorHex" %in% names(study_pal_df)) study_pal_df$ColorHex else study_pal_df$Color
  study_colors <- setNames(scol, study_pal_df$Study)
} else {
  study_colors <- c("2019Ahn" = "#88CCEE", "2022Dong" = "#332288", "2023Deng" = "#44AA99", "2023Jeon" = "#AA4499", "2023Lee" = "#DDCC77", "2024Fan" = "#117733", "TCGA-CHOL" = "#CC6677")
}

variant_colors <- c(
  "Missense" = "#4B9B6E", "Nonsense" = "#D97933",
  "Frameshift deletion" = "#467AA8", "Frameshift insertion" = "#8C6BB1",
  "In-frame deletion" = "#D5B44B", "In-frame insertion" = "#319E9A",
  "Splice" = "#B34B5E", "Other non-silent" = "#9B8A7A", "Multi-hit" = "#333333",
  "No recorded variant" = "#E9E9E9"
)

mm <- function(x) unit(x, "mm")

txt <- function(s, x, y, size = 6.0, bold = FALSE, col = "black", just = "centre", italic = FALSE, lineheight = 1) {
  grid.text(s, mm(x), mm(y), just = just,
            gp = gpar(fontfamily = "Arial", fontsize = size,
                      fontface = if (italic) "italic" else if (bold) "bold" else "plain",
                      col = col, lineheight = lineheight))
}

rect_box <- function(x, y, w, h, fill, col = NA, lwd = 0.3) {
  grid.rect(mm(x), mm(y), mm(w), mm(h), gp = gpar(fill = fill, col = col, lwd = lwd))
}

line_seg <- function(x, y, col = "grey60", lwd = 0.5, lty = 1) {
  grid.lines(mm(x), mm(y), gp = gpar(col = col, lwd = lwd, lty = lty))
}

point_dot <- function(x, y, fill, col = fill, size = 1.4, pch = 21, lwd = 0.65) {
  grid.points(mm(x), mm(y), pch = pch, size = mm(size), gp = gpar(fill = fill, col = col, lwd = lwd))
}

# ------------------------------------------------------------------------------
# Panel A: Mutation Landscape Oncoprint (15 Genes, 308 Tumors)
# ------------------------------------------------------------------------------
cat("Rendering Panel A...\n")
draw_panelA <- function() {
  genes <- fread(file.path(data_dir, "gene_order.csv"))[order(DisplayOrder), Gene]
  d <- fread(file.path(data_dir, "mutation_display_long.csv"))
  s <- fread(file.path(data_dir, "sample_order.csv"))[order(DisplayOrder)]
  s[, Sample := as.character(Sample)]
  
  grid.newpage()
  txt("Driver mutation landscape", 3, 98, 8.5, bold = TRUE, just = "left")
  
  left <- 20; right <- 106; gap <- 1.5; usable <- right - left - 2 * gap; top <- 83; bottom <- 30
  rh <- (top - bottom) / length(genes)
  
  for (st in names(pal)) {
    ss <- s[FinalSubtype == st]
    width <- nrow(ss) / nrow(s) * usable
    cw <- width / nrow(ss)
    rect_box(left + width / 2, (top + bottom) / 2, width, top - bottom, variant_colors[["No recorded variant"]])
    
    z <- d[FinalSubtype == st & Mutated == TRUE]
    grid.rect(mm(left + (match(as.character(z$Sample), ss$Sample) - 0.5) * cw),
              mm(top - (match(z$Gene, genes) - 0.5) * rh),
              mm(cw + 0.008), mm(rh),
              gp = gpar(fill = variant_colors[z$MutationClass], col = NA))
    
    rect_box(left + width / 2, 88, width, 5, pal[[st]])
    txt(paste0(st, "  n = ", nrow(ss)), left + width / 2, 88, 7.0, bold = TRUE, col = "black")
    
    grid.rect(mm(left + (seq_len(nrow(ss)) - 0.5) * cw), mm(26),
              mm(cw + 0.008), mm(2.8),
              gp = gpar(fill = study_colors[ss$Study], col = NA))
    left <- left + width + gap
  }
  
  # Frequencies
  freq <- fread(file.path(data_dir, "mutation_frequency.csv"))
  allfreq <- freq[, .(N = sum(N), Altered = sum(Altered)), by = Gene]
  allfreq[, Fraction := Altered / N]
  centers <- c(114, 126, 138, 150)
  
  txt("Mutation frequency (%)", 132, 94, 6.8, bold = TRUE)
  for (k in seq_along(centers)) {
    txt(c("All", names(pal))[k], centers[k], 88, 6.5, bold = TRUE)
  }
  
  for (j in seq_along(genes)) {
    yy <- top - (j - 0.5) * rh
    txt(genes[j], 18, yy, 7.2, just = "right", italic = TRUE)
    if (j %% 2 == 0) rect_box(132, yy, 48, rh, "#F5F5F5")
    for (k in seq_along(centers)) {
      if (k == 1L) {
        row <- allfreq[Gene == genes[j]]; ss <- d[Gene == genes[j]]; color <- "grey50"
      } else {
        st <- names(pal)[k - 1L]; row <- freq[Gene == genes[j] & FinalSubtype == st]
        ss <- d[Gene == genes[j] & FinalSubtype == st]; color <- pal[[st]]
      }
      bw <- 12 * row$Fraction / 0.4
      if (bw > 0) rect_box(centers[k] - 6 + bw / 2, yy, bw, 3.0, adjustcolor(color, alpha.f = 0.24))
      txt(sprintf("%.1f", 100 * row$Fraction), centers[k], yy, 6.0)
    }
  }
  
  txt("Study", 18, 26, 6.5, just = "right")
  legend_order <- c("Missense", "Nonsense", "Splice", "Frameshift deletion", "Frameshift insertion",
                    "Multi-hit", "In-frame deletion", "In-frame insertion", "No recorded variant")
  types <- legend_order[legend_order %in% unique(d$MutationClass)]
  label_map <- c("Frameshift deletion" = "Frameshift del.", "Frameshift insertion" = "Frameshift ins.",
                 "In-frame deletion" = "In-frame del.", "In-frame insertion" = "In-frame ins.")
  
  txt("Mutation class", 4, 21, 6.5, bold = TRUE, just = "left")
  txt("Study", 112, 21, 6.5, bold = TRUE, just = "left")
  line_seg(c(108, 108), c(4, 22), col = "grey80", lwd = 0.4)
  
  for (i in seq_along(types)) {
    xx <- 5 + ((i - 1) %% 3) * 34; yy <- 16 - floor((i - 1) / 3) * 4.5
    rect_box(xx, yy, 2.2, 2.2, variant_colors[[types[i]]])
    lab <- if (types[i] %in% names(label_map)) label_map[[types[i]]] else types[i]
    txt(lab, xx + 2, yy, 6.0, just = "left")
  }
  
  studies <- unique(s$Study)
  for (i in seq_along(studies)) {
    yy <- 16 - (i - 1) * 4.5
    rect_box(113, yy, 2.2, 2.2, study_colors[[studies[i]]])
    txt(studies[i], 115, yy, 6.0, just = "left")
  }
}

out_panelA_png <- file.path(out_dir, "Figure_4_PanelA.png")
out_panelA_pdf <- file.path(out_dir, "Figure_4_PanelA.pdf")
agg_png(out_panelA_png, width = 158, height = 100, units = "mm", res = 600, background = "white")
draw_panelA()
dev.off()
pdf(out_panelA_pdf, width = 158 / 25.4, height = 100 / 25.4, family = "Arial")
draw_panelA()
dev.off()

# ------------------------------------------------------------------------------
# Panel B: CMH Odds Ratios for Driver Mutations
# ------------------------------------------------------------------------------
cat("Rendering Panel B...\n")
draw_panelB <- function() {
  genes <- fread(file.path(data_dir, "gene_order.csv"))[order(DisplayOrder), Gene]
  d <- fread(file.path(data_dir, "mutation_frequency.csv"))
  grid.newpage()
  
  txt("Study-stratified odds ratio (95% CI)", 62, 98, 7.5, bold = TRUE)
  top <- 83; bottom <- 20; rh <- (top - bottom) / length(genes)
  starts <- c(22, 54, 86); fw <- 28; lims <- c(0.02, 200)
  trans <- function(v, k) starts[k] + (log10(v) - log10(lims[1])) / diff(log10(lims)) * fw
  
  for (j in seq_along(genes)) {
    yy <- top - (j - 0.5) * rh
    if (j %% 2 == 0) rect_box(60, yy, 100, rh, "#F5F5F5")
    txt(genes[j], 16, yy, 7.0, just = "right", italic = TRUE)
  }
  
  for (k in seq_along(pal)) {
    st <- names(pal)[k]; z <- d[FinalSubtype == st]; n <- unique(z$N)
    txt(paste0(st, " (n = ", n, ")"), starts[k] + fw / 2, 88, 7.0, bold = TRUE, col = pal[[st]])
    line_seg(c(trans(1, k), trans(1, k)), c(bottom, top), col = "grey65", lty = 2)
    
    for (j in seq_along(genes)) {
      row <- z[Gene == genes[j]]; yy <- top - (j - 0.5) * rh
      est <- row$CMHOddsRatio; lo <- row$CMHLower95; hi <- row$CMHUpper95; q <- row$CMHFDR
      fill <- if (!is.na(q) && q < 0.05) pal[[st]] else "white"
      
      if (is.finite(est) && est > 0) {
        xx <- trans(pmax(lims[1], pmin(lims[2], est)), k)
        if (is.finite(lo) && is.finite(hi) && lo > 0 && hi > 0) {
          line_seg(trans(pmax(lims[1], pmin(lims[2], c(lo, hi))), k), c(yy, yy), col = pal[[st]], lwd = 0.8)
        }
        point_dot(xx, yy, fill = fill, col = pal[[st]], size = 1.4, pch = 21, lwd = 0.8)
      } else if (!is.na(est) && (est == 0 || is.infinite(est))) {
        xx <- if (est == 0) starts[k] else starts[k] + fw
        offsets <- if (est == 0) c(0, 1.2, 1.2) else c(0, -1.2, -1.2)
        grid.polygon(mm(xx + offsets), mm(yy + c(0, 0.7, -0.7)), gp = gpar(col = pal[[st]], fill = fill, lwd = 0.8))
      } else {
        txt("NA", starts[k] + fw / 2, yy, 6.0, col = "grey50")
      }
    }
    line_seg(c(starts[k], starts[k] + fw), c(bottom - 1, bottom - 1), col = "grey40")
    for (v in c(0.1, 1, 10, 100)) {
      xx <- trans(v, k); line_seg(c(xx, xx), c(bottom - 1, bottom - 2), col = "grey40")
      txt(as.character(v), xx, bottom - 4.5, 5.8)
    }
  }
  point_dot(30, 8, fill = "grey40", col = "grey40", size = 1.4, pch = 21)
  txt("FDR < 0.05", 33, 8, 5.8, just = "left")
  point_dot(55, 8, fill = "white", col = "grey40", size = 1.4, pch = 21)
  txt("FDR >= 0.05", 58, 8, 5.8, just = "left")
}

out_panelB_png <- file.path(out_dir, "Figure_4_PanelB.png")
out_panelB_pdf <- file.path(out_dir, "Figure_4_PanelB.pdf")
agg_png(out_panelB_png, width = 120, height = 100, units = "mm", res = 600, background = "white")
draw_panelB()
dev.off()
pdf(out_panelB_pdf, width = 120 / 25.4, height = 100 / 25.4, family = "Arial")
draw_panelB()
dev.off()

# ------------------------------------------------------------------------------
# Panel C: Proteomic Enrichment & Abundance Heatmap (248 Tumors)
# ------------------------------------------------------------------------------
cat("Rendering Panel C...\n")
draw_panelC <- function() {
  plan <- fread(file.path(data_dir, "prot_display_plan.csv"))
  main <- fread(file.path(data_dir, "prot_main_programs.csv"))
  roster <- fread(file.path(data_dir, "prot_sample_roster.csv"))
  marker <- fread(file.path(data_dir, "prot_selected_abundance_long.csv"))
  selected <- fread(file.path(data_dir, "prot_selected_proteins.csv"))
  
  grid.newpage()
  txt("Proteomic pathways & abundance", 2, 82, 8.0, bold = TRUE, just = "left")
  
  # Left: GSEA NES bubble matrix
  xs <- c(48, 58, 68); h <- 74; top <- 74; bottom <- 12
  rh <- (top - bottom) / nrow(plan)
  
  for (j in 1:3) {
    txt(names(pal)[j], xs[j], top + 3, 6.5, bold = TRUE)
  }
  
  color_ramp <- colorRamp(c("#2166AC", "#FAFAFA", "#B2182B"))
  nes_col <- function(v) rgb(color_ramp(pmax(0, pmin(1, (v / 3 + 1) / 2))), maxColorValue = 255)
  
  for (i in seq_len(nrow(plan))) {
    yy <- top - (i - 0.5) * rh
    p_row <- plan[i]
    txt(p_row$Label, 42, yy, 5.2, just = "right")
    for (j in 1:3) {
      st <- names(pal)[j]
      sub_m <- main[Pathway == p_row$Pathway & Subtype == st]
      if (nrow(sub_m) > 0) {
        rect_box(xs[j], yy, 8, rh - 0.4, nes_col(sub_m$NES))
        if (sub_m$padj < 0.05) {
          point_dot(xs[j], yy, fill = if (abs(sub_m$NES) > 2) "white" else "grey20",
                    col = if (abs(sub_m$NES) > 2) "white" else "grey20", size = 0.5, pch = 16)
        }
      }
    }
  }
  
  # Right: 30 selected proteins heatmap across 248 tumors
  genes <- unique(selected$Protein)
  heat_left <- 78; heat_right <- 128; heat_w <- heat_right - heat_left
  rh_prot <- (top - bottom) / length(genes)
  
  # Subtype track
  sub_lens <- table(roster$FinalSubtype)
  pos_x <- heat_left
  for (st in names(pal)) {
    w_st <- heat_w * sub_lens[st] / nrow(roster)
    rect_box(pos_x + w_st / 2, top + 3, w_st, 3.5, pal[st])
    txt(st, pos_x + w_st / 2, top + 3, 6.0, bold = TRUE, col = "black")
    pos_x <- pos_x + w_st
  }
  
  # Plot abundance tiles
  abund_col <- function(v) {
    zcols <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(256)
    zcols[pmax(1L, pmin(256L, 1L + floor((pmax(-2, pmin(2, v)) + 2) / 4 * 255)))]
  }
  
  cw_samp <- heat_w / nrow(roster)
  for (j in seq_along(genes)) {
    yy <- top - (j - 0.5) * rh_prot
    txt(genes[j], heat_right + 1, yy, 5.0, just = "left", italic = TRUE)
    g_dat <- marker[Protein == genes[j]][order(SampleOrder)]
    cols_v <- abund_col(g_dat$StandardizedAbundance)
    for (s_idx in seq_len(nrow(g_dat))) {
      rect_box(heat_left + (s_idx - 0.5) * cw_samp, yy, cw_samp + 0.05, rh_prot, cols_v[s_idx])
    }
  }
  
  # NES legend
  txt("NES", 10, 5, 5.5)
  for (i in 1:50) {
    rect_box(16 + (i - 0.5) * 0.3, 5, 0.32, 1.2, nes_col(-3 + (i - 1) / 49 * 6))
  }
  txt("-3", 16, 2.5, 4.8); txt("0", 23.5, 2.5, 4.8); txt("+3", 31, 2.5, 4.8)
  point_dot(38, 5, fill = "grey20", col = "grey20", size = 0.5, pch = 16)
  txt("FDR < 0.05", 40, 5, 4.8, just = "left")
}

out_panelC_png <- file.path(out_dir, "Figure_4_PanelC.png")
out_panelC_pdf <- file.path(out_dir, "Figure_4_PanelC.pdf")
agg_png(out_panelC_png, width = 145, height = 88, units = "mm", res = 600, background = "white")
draw_panelC()
dev.off()
pdf(out_panelC_pdf, width = 145 / 25.4, height = 88 / 25.4, family = "Arial")
draw_panelC()
dev.off()

# ------------------------------------------------------------------------------
# Panel D: Key Program Comparison Bridge
# ------------------------------------------------------------------------------
cat("Rendering Panel D...\n")
draw_panelD <- function() {
  main <- fread(file.path(data_dir, "D_main_programs.csv"))
  plan <- fread(file.path(data_dir, "D_fixed_program_plan.csv"))
  grid.newpage()
  
  txt("Key-program comparison", 2, 59, 7.2, bold = TRUE, just = "left")
  xs <- c(31, 44, 57, 70); stages <- c("R1", "R2", "R3", "P3")
  labels <- c("Full RNA\n497 / all", "Matched RNA\n248 / all", "Matched RNA\n248 / shared", "Protein\n248 / shared")
  for (j in 1:4) txt(labels[j], xs[j], 52.5, 5.0)
  txt("NES", 21, 52.5, 5.0, col = "grey40")
  
  for (i in 1:5) {
    y0 <- 43.5 - (i - 1) * 8.6
    z <- main[Order == i][match(stages, Stage)]
    yy <- y0 + z$NES * 0.91
    if (i %% 2 == 1) rect_box(39, y0, 76, 7.9, "#F7F8F9")
    txt(plan$Label[i], 2, y0, 5.5, bold = TRUE, col = pal[[plan$Subtype[i]]], just = "left")
    line_seg(c(27, 74), c(y0, y0), col = "grey75", lwd = 0.35, lty = 2)
    for (j in 1:3) {
      line_seg(c(xs[j], xs[j + 1]), c(yy[j], yy[j + 1]), col = pal[[plan$Subtype[i]]], lwd = 0.7, lty = if (j == 3) 2 else 1)
    }
    for (j in 1:4) {
      point_dot(xs[j], yy[j], fill = if (z$QCommon[j] < 0.05) pal[[plan$Subtype[i]]] else "white",
                col = pal[[plan$Subtype[i]]], size = 1.3)
    }
    for (v in c(-3, 3)) line_seg(c(26, 27), c(y0 + v * 0.91, y0 + v * 0.91), col = "grey55", lwd = 0.4)
    if (i == 1) {
      txt("3", 25.3, y0 + 2.73, 5.0, just = "right")
      txt("0", 25.3, y0, 5.0, just = "right")
      txt("-3", 25.3, y0 - 2.73, 5.0, just = "right")
    }
  }
  point_dot(3, 3, fill = "grey35", col = "grey35", size = 1.2)
  txt("q < .05", 5, 3, 5.2, just = "left")
  point_dot(24, 3, fill = "white", col = "grey35", size = 1.2)
  txt("q >= .05", 26, 3, 5.2, just = "left")
  txt("RNA -> protein", 59, 3, 5.2, col = "grey40")
}

out_panelD_png <- file.path(out_dir, "Figure_4_PanelD.png")
out_panelD_pdf <- file.path(out_dir, "Figure_4_PanelD.pdf")
agg_png(out_panelD_png, width = 78, height = 61, units = "mm", res = 600, background = "white")
draw_panelD()
dev.off()
pdf(out_panelD_pdf, width = 78 / 25.4, height = 61 / 25.4, family = "Arial")
draw_panelD()
dev.off()

# ------------------------------------------------------------------------------
# Panel E: Matched Pathway Concordance Quadrant
# ------------------------------------------------------------------------------
cat("Rendering Panel E...\n")
draw_panelE <- function() {
  sel <- fread(file.path(data_dir, "E_selected_pathways.csv"))
  grid.newpage()
  
  txt("Matched pathway concordance", 2, 59.2, 7.0, bold = TRUE, just = "left")
  L <- 12; R <- 56; B <- 11; T <- 55; xm <- 34; ym <- 33; s <- 5.5
  
  rect_box((L + xm) / 2, (ym + T) / 2, xm - L, T - ym, "#FCF7F2")
  rect_box((xm + R) / 2, (ym + T) / 2, R - xm, T - ym, "#F0F5EF")
  rect_box((L + xm) / 2, (B + ym) / 2, xm - L, ym - B, "#F0F2F4")
  rect_box((xm + R) / 2, (B + ym) / 2, R - xm, ym - B, "#FCF7F2")
  
  line_seg(c(L, R), c(ym, ym), col = "grey50")
  line_seg(c(xm, xm), c(B, T), col = "grey50")
  
  for (v in seq(-3, 3, 1)) {
    txt(as.character(v), xm + v * s, 8.5, 5.0)
    txt(as.character(v), L - 2.5, ym + v * s, 5.0, just = "right")
  }
  txt("RNA NES", xm, 5.0, 5.5)
  grid.text("Protein NES", mm(L - 6.5), mm(ym), rot = 90,
            gp = gpar(fontfamily = "Arial", fontsize = 5.5))
  
  # Scatter points
  for (i in seq_len(nrow(sel))) {
    r <- sel[i]
    st <- r$Subtype
    xp <- xm + r$NES_RNA * s
    yp <- ym + r$NES_Protein * s
    if (xp >= L && xp <= R && yp >= B && yp <= T) {
      point_dot(xp, yp, fill = pal[[st]], col = pal[[st]], size = 1.3)
    }
  }
  txt("Concordant (+/+)", R - 2, T - 2, 4.8, bold = TRUE, col = "#2B613A", just = "right")
  txt("Concordant (-/-)", L + 2, B + 2, 4.8, bold = TRUE, col = "#2E5575", just = "left")
}

out_panelE_png <- file.path(out_dir, "Figure_4_PanelE.png")
out_panelE_pdf <- file.path(out_dir, "Figure_4_PanelE.pdf")
agg_png(out_panelE_png, width = 65, height = 61, units = "mm", res = 600, background = "white")
draw_panelE()
dev.off()
pdf(out_panelE_pdf, width = 65 / 25.4, height = 61 / 25.4, family = "Arial")
draw_panelE()
dev.off()

# ------------------------------------------------------------------------------
# Panel F: MOFA Triple Alluvial (RNA -> RNA+Prot -> RNA+Prot+CNV)
# ------------------------------------------------------------------------------
cat("Rendering Panel F...\n")
draw_panelF <- function() {
  d <- fread(file.path(data_dir, "F_main_joined_assignments.csv"))
  paths <- fread(file.path(data_dir, "F_alluvial_paths.csv"))
  grid.newpage()
  
  txt("Multi-omics subtype concordance", 2, 70, 7.2, bold = TRUE, just = "left")
  txt("RP (RNA + Protein)", 25, 65, 5.8, bold = TRUE)
  txt("Primary (RNA)", 75, 65, 5.8, bold = TRUE)
  txt("RPC (RNA + Prot + CNV)", 125, 65, 5.8, bold = TRUE)
  
  unit_h <- 42 / 182
  # Primary bars
  prim_counts <- table(factor(d$FinalSubtype, levels = c("C1", "C2", "C3")))
  y_c <- 60
  for (st in c("C1", "C2", "C3")) {
    h_st <- prim_counts[st] * unit_h
    rect_box(75, y_c - h_st / 2, 12, h_st, pal[st])
    txt(paste0(st, " (", prim_counts[st], ")"), 75, y_c - h_st / 2, 5.5, bold = TRUE, col = "black")
    y_c <- y_c - h_st - 2.5
  }
  
  # RP and RPC bars
  rp_counts <- table(factor(d$RPCluster, levels = 1:3))
  y_c <- 60
  for (cl in 1:3) {
    h_cl <- rp_counts[cl] * unit_h
    rect_box(25, y_c - h_cl / 2, 8, h_cl, "#9FB7C8")
    txt(paste0("RP", cl, " (", rp_counts[cl], ")"), 25, y_c - h_cl / 2, 5.0, col = "black")
    y_c <- y_c - h_cl - 2.5
  }
  
  rpc_counts <- table(factor(d$RPCCluster, levels = 1:3))
  y_c <- 60
  for (cl in 1:3) {
    h_cl <- rpc_counts[cl] * unit_h
    rect_box(125, y_c - h_cl / 2, 8, h_cl, "#C4B2A0")
    txt(paste0("RPC", cl, " (", rpc_counts[cl], ")"), 125, y_c - h_cl / 2, 5.0, col = "black")
    y_c <- y_c - h_cl - 2.5
  }
  
  # Connective lines
  txt("n = 182 principal-state tumors (Dong)", 75, 8, 5.2, col = "grey40")
}

out_panelF_png <- file.path(out_dir, "Figure_4_PanelF.png")
out_panelF_pdf <- file.path(out_dir, "Figure_4_PanelF.pdf")
agg_png(out_panelF_png, width = 150, height = 75, units = "mm", res = 600, background = "white")
draw_panelF()
dev.off()
pdf(out_panelF_pdf, width = 150 / 25.4, height = 75 / 25.4, family = "Arial")
draw_panelF()
dev.off()

# ------------------------------------------------------------------------------
# Panel G: Biological Synthesis Cards & Arrows
# ------------------------------------------------------------------------------
cat("Rendering Panel G...\n")
draw_panelG <- function() {
  cells4 <- fread(file.path(data_dir, "G_all78_cells.csv"))
  summary4 <- fread(file.path(data_dir, "G_displayed_synthesis.csv"))
  grid.newpage()
  
  txt("Biological program", 2, 76, 6.4, bold = TRUE, just = "left")
  for (v in 1:2) {
    xx <- c(49, 74)[v]
    rect_box(xx, 76, 23, 7, c("#EDF3FA", "#EDF5EF")[v])
    txt(c("RNA", "Protein")[v], xx, 76, 6.7, bold = TRUE)
    for (k in 1:3) txt(names(pal)[k], xx - 7 + (k - 1) * 7, 69.5, 6.0, bold = TRUE, col = pal[k])
  }
  txt("Genomic associations", 106, 76, 6.4, bold = TRUE)
  txt("Subtype synthesis", 154.5, 76, 6.4, bold = TRUE)
  
  plan4 <- unique(cells4[, .(Pathway, Label)], by = "Pathway")
  plan4[, Row := seq_len(.N)]
  
  for (i in 1:13) {
    z <- plan4[Row == i]; y <- 65 - (i - 1) * 4.2
    if (i %% 2 == 1) rect_box(43.5, y, 85, 3.9, "#F6F7F8")
    txt(z$Label, 2, y, 5.5, just = "left")
    for (v in 1:2) for (k in 1:3) {
      d <- cells4[Pathway == z$Pathway & Subtype == names(pal)[k] & View == c("RNA", "Protein")[v]]
      fill <- if (d$FDR >= 0.05) "#EEEEEE" else if (d$NES > 0) "#F7DCCB" else "#D4E4F0"
      rect_box(c(42, 67)[v] + (k - 1) * 7, y, 6.2, 3.7, fill)
      txt(d$Glyph, c(42, 67)[v] + (k - 1) * 7, y, 7.5, bold = TRUE,
          col = if (d$FDR >= 0.05) "grey55" else if (d$NES > 0) "#A95334" else "#2C719E")
    }
  }
  
  dna_text <- c("No FDR-significant\ndriver enrichment",
                "BAP1, IDH1/2 enriched\n1p36 loss > C3",
                "KRAS, SMAD4, TP53\nenriched; 1p36 loss < C2")
  
  for (k in 1:3) {
    top <- 67 - (k - 1) * 19.2
    rect_box(106, top - 8, 35, 17.8, "#FAFAFA", col = "grey85")
    rect_box(154.5, top - 8, 53, 17.8, adjustcolor(pal[k], alpha.f = 0.05), col = adjustcolor(pal[k], alpha.f = 0.4))
    rect_box(106, top, 35, 3, adjustcolor(pal[k], alpha.f = 0.2))
    txt(names(pal)[k], 90, top, 6.0, bold = TRUE, col = pal[k], just = "left")
    txt(dna_text[k], 90, top - 8, 5.3, just = "left")
    
    sub_syn <- summary4[Subtype == names(pal)[k]]
    txt(paste(names(pal)[k], sub_syn$Text[1]), 129.5, top, 5.8, bold = TRUE, col = pal[k], just = "left")
    for (j in 2:min(4, nrow(sub_syn))) {
      txt(sub_syn$Text[j], 129.5, top - 4 - (j - 2) * 3.6, 5.0, just = "left")
    }
  }
}

out_panelG_png <- file.path(out_dir, "Figure_4_PanelG.png")
out_panelG_pdf <- file.path(out_dir, "Figure_4_PanelG.pdf")
agg_png(out_panelG_png, width = 174, height = 82, units = "mm", res = 600, background = "white")
draw_panelG()
dev.off()
pdf(out_panelG_pdf, width = 174 / 25.4, height = 82 / 25.4, family = "Arial")
draw_panelG()
dev.off()

# ------------------------------------------------------------------------------
# Full Composite: Figure 4 (180 mm x 270 mm)
# ------------------------------------------------------------------------------
cat("Assembling Full Figure 4 Composite...\n")
composite_png <- file.path(out_dir, "Figure_4_composite.png")
composite_pdf <- file.path(out_dir, "Figure_4_composite.pdf")

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
  
  # Row 1: Panel A (oncoprint) + Panel B (odds ratios)
  draw_raster_panel(out_panelA_png, 5.1, 8.5, 93.3, 57.3)
  draw_panel_letter("A", 0.2, 10.8)
  
  draw_raster_panel(out_panelB_png, 98.5, 8.7, 78.2, 57.0)
  draw_panel_letter("B", 98.7, 10.8)
  
  # Row 2: Panel C (proteomics) + Panel F (MOFA alluvial) matching Slide 12
  draw_raster_panel(out_panelC_png, 6.2, 66.5, 117.0, 75.7)
  draw_panel_letter("C", 0.2, 64.1)
  
  draw_raster_panel(out_panelF_png, 126.0, 66.5, 50.0, 75.7)
  draw_panel_letter("F", 124.0, 63.6)
  
  # Row 3: Panel D (bridge / 5 programs) + Panel E (matched pathway concordance)
  draw_raster_panel(out_panelD_png, 0.7, 144.0, 64.0, 49.8)
  draw_panel_letter("D", 0.7, 142.0)
  
  draw_raster_panel(out_panelE_png, 72.0, 144.0, 104.0, 49.8)
  draw_panel_letter("E", 68.0, 142.0)
  
  # Row 4: Panel G (synthesis)
  draw_raster_panel(out_panelG_png, 11.8, 198.0, 154.3, 66.7)
  draw_panel_letter("G", 0.2, 196.0)
  
  popViewport()
}

agg_png(composite_png, width = 180, height = 270, units = "mm", res = 300, background = "white")
render_composite()
dev.off()

pdf(composite_pdf, width = 180 / 25.4, height = 270 / 25.4, family = "Arial")
render_composite()
dev.off()

cat("Figure 4 rendered successfully:\n")
cat("  - Individual panels:", out_dir, "\n")
cat("  - Composite PDF:", composite_pdf, "\n")
cat("  - Composite PNG:", composite_png, "\n")
