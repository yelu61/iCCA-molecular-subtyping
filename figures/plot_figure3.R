#!/usr/bin/env Rscript
# ==============================================================================
# Figure 3: Native Single-Cell Expression Localizes Bulk-Associated Programs
# Standalone, decoupled reproduction script.
#
# Generates Panels A-E and assembled Figure 3.
# ==============================================================================

options(ggrastr.default.dev = "ragg")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(ggrastr)
  library(grid)
  library(data.table)
})

if (!"Arial" %in% names(grDevices::pdfFonts())) grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])

# ------------------------------------------------------------------------------
# 1. Configuration & CLI Arguments
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default_val) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) args[idx + 1] else default_val
}

data_cands <- c(get_arg("--data-dir", ""), "figures/data/figure3", "8-Visualization/data/figure3", "data/figure3", "../data/figure3")
DATA_DIR <- data_cands[dir.exists(data_cands) & data_cands != ""][1]
if (is.na(DATA_DIR)) DATA_DIR <- "8-Visualization/data/figure3"
DATA_DIR <- normalizePath(DATA_DIR)

out_cands <- c(get_arg("--out-dir", ""), if (dir.exists("8-Visualization/figures")) "8-Visualization/figures/figure3" else "", "figures/outputs/figure3", "outputs/figure3")
OUT_DIR <- out_cands[out_cands != ""][1]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)

rdata <- function(filename) {
  p <- file.path(DATA_DIR, filename)
  if (!file.exists(p)) stop(paste("Required data file not found:", p))
  read_csv(p, show_col_types = FALSE)
}

base_theme <- theme_classic(base_family = "Arial", base_size = 5.5) +
  theme(
    plot.title = element_text(size = 6.5, face = "bold", colour = "#222222"),
    plot.subtitle = element_text(size = 5, colour = "#555555"),
    plot.caption = element_text(size = 5, hjust = 0),
    axis.text = element_text(size = 5, colour = "#333333"),
    axis.title = element_text(size = 5.2),
    strip.text = element_text(size = 5.2, face = "bold"),
    strip.background = element_rect(fill = "#F0F1F2", colour = NA),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.position = "bottom",
    panel.spacing = unit(1.2, "mm"),
    plot.margin = margin(2, 1, 1, 5)
  )

theme_set(base_theme)

# ------------------------------------------------------------------------------
# 2. Panel Drawing Functions
# ------------------------------------------------------------------------------

point_rast <- function(...) geom_point_rast(..., dev = "ragg_png")

# Panel A: Opening Atlas UMAP
draw_panel_A <- function() {
  cells <- rdata("plotted_cells.csv.gz")
  paldf <- rdata("annotation_palette.csv")
  pal   <- setNames(paldf$Color, paldf$CellTypeDisplay)
  
  cells$PlotCellType <- factor(cells$PlotCellType, levels = names(pal))
  cells <- cells %>% add_count(PlotCellType, name = "GroupN") %>% arrange(desc(GroupN), CellKey)
  
  xlim <- c(-10.8, 22.1)
  ylim <- c(-13.2, 20.0)
  
  pos <- data.frame(
    PlotCellType = c("Epithelial", "Hepatocyte", "CD8 T", "Other T", "NK-like", "B cell",
                     "Plasma cell", "Mono/macrophage", "Dendritic-like", "Fibroblast",
                     "Pericyte-like", "Endothelial"),
    x = c(-6.5, -7.0, 18.1, 15.2, 11.4, 18.6, 14.8, 10.0, 17.3, 2.0, 3.7, 1.1),
    y = c(7.6, -5.7, 0.5, -5.7, 5.6, 8.2, -11.4, 18.0, 11.7, -5.7, 1.1, -10.4),
    target_x = c(-3.7, -6.4, 15.5, 13.5, 11.4, 15.6, 12, 10.2, 12.5, 4.6, 4, 4),
    target_y = c(9.5, -4.2, 1.4, -3.5, 3.8, 8.5, -9.2, 15.2, 10, -3.8, -1, -9)
  )
  
  for (i in seq_len(nrow(pos))) {
    z <- cells[cells$PlotCellType == pos$PlotCellType[i], ]
    j <- which.min((z$UMAP1 - pos$target_x[i])^2 + (z$UMAP2 - pos$target_y[i])^2)
    pos$AnchorX[i] <- z$UMAP1[j]
    pos$AnchorY[i] <- z$UMAP2[j]
  }
  
  lp <- pos
  lp$Label <- as.character(lp$PlotCellType)
  lp$Label[lp$Label == "Mono/macrophage"] <- "Monocytes/macrophages"
  
  p <- ggplot(cells, aes(UMAP1, UMAP2, colour = PlotCellType)) +
    point_rast(data = cells[!cells$PlotMixedOpenSymbol, ], shape = 16, size = 0.23, alpha = 0.88, raster.dpi = 300) +
    point_rast(data = cells[cells$PlotMixedOpenSymbol, ], shape = 1, size = 0.32, stroke = 0.10, alpha = 0.88, raster.dpi = 300) +
    scale_colour_manual(values = pal, drop = FALSE) +
    scale_x_continuous(breaks = c(-5, 5, 15)) +
    scale_y_continuous(breaks = c(-10, 0, 10)) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
    geom_segment(data = lp, aes(x = AnchorX, y = AnchorY, xend = x, yend = y), inherit.aes = FALSE, colour = "#777777", linewidth = 0.22) +
    geom_label(data = lp, aes(x = x, y = y, label = Label), inherit.aes = FALSE, family = "Arial", fontface = "bold",
               size = 6.8 / ggplot2::.pt, colour = "#333333", fill = "white", linewidth = 0.3,
               lineheight = 0.95, label.padding = unit(0.12, "lines"), label.r = unit(0.11, "lines")) +
    annotate("segment", x = -9.3, xend = -5.6, y = -11.7, yend = -11.7, linewidth = 0.35, colour = "#444444",
             arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
    annotate("segment", x = -9.3, xend = -9.3, y = -11.7, yend = -8.0, linewidth = 0.35, colour = "#444444",
             arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
    annotate("text", x = -7.2, y = -12.5, label = "UMAP 1", family = "Arial", size = 6 / ggplot2::.pt, colour = "#444444") +
    annotate("text", x = -10.1, y = -9.6, label = "UMAP 2", angle = 90, family = "Arial", size = 6 / ggplot2::.pt, colour = "#444444") +
    labs(title = "iCCA Single-Cell Atlas (34,104 plotted cells)") +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), legend.position = "none")
  
  return(p)
}

# Panel B: Grouped Program Sources Dotplot
draw_panel_B <- function() {
  ord <- c("Pattern_P1", "Pattern_P2", "WGCNA_MElightgreen", "WGCNA_MEpink",
           "Pattern_P3", "WGCNA_MEblack", "WGCNA_MEcyan", "WGCNA_MEbrown", "WGCNA_MEdarkgreen",
           "Pattern_P4", "Pattern_P5", "Pattern_P6", "WGCNA_MEmagenta", "WGCNA_MEgreenyellow",
           "WGCNA_MEred", "WGCNA_MEdarkturquoise", "WGCNA_MEtan")
  pg  <- setNames(c(rep("C1", 4), rep("C2", 5), rep("C3", 8)), ord)
  ll  <- setNames(c("P1", "P2", "lightgreen | Peroxisome", "pink | Catabolism",
                    "P3", "black | Mitochondrial", "cyan | Respiration", "brown | Morphogenesis", "darkgreen | Checkpoints",
                    "P4", "P5", "P6", "magenta | Mitotic cycle", "greenyellow | Hypoxia", "red | Glycosylation",
                    "darkturquoise | IFN", "tan | Muscle prolif."), ord)
  
  fine <- c("Epithelial", "Hepatocyte", "CD8 T", "Other T", "NK-like", "B cell",
            "Plasma cell", "Monocytes/macrophages", "Dendritic-like", "Fibroblast", "Pericyte-like", "Endothelial")
  fl   <- setNames(c("Epi.", "Hep.", "CD8 T", "Other T", "NK-like", "B", "Plasma", "Mono/Mac", "DC-like", "Fib.", "Peri.", "Endo."), fine)
  
  s <- rdata("dotplot_summary.csv") %>% filter(Level == "Fine", Method == "UCell")
  b <- expand_grid(tissue = c("Adjacent", "Tumor"), CellType = fine, ProgramID = ord) %>%
    left_join(s, by = c("tissue", "CellType", "ProgramID"))
  b$BulkGroup <- factor(pg[b$ProgramID], levels = c("C1", "C2", "C3"))
  b$ProgramID <- factor(b$ProgramID, levels = rev(ord))
  b$tissue    <- factor(b$tissue, levels = c("Adjacent", "Tumor"))
  b$CellType  <- factor(b$CellType, levels = fine)
  
  pb <- ggplot(b, aes(CellType, ProgramID)) +
    geom_tile(fill = "#FAFAFA", colour = "white", linewidth = 0.1) +
    geom_point(data = filter(b, !is.na(Score)), aes(fill = Score, size = NPatients), shape = 21, stroke = 0.1, colour = "#555555") +
    geom_text(data = filter(b, is.na(Score)), label = "-", size = 5 / ggplot2::.pt, colour = "#999999") +
    facet_grid(BulkGroup ~ tissue, scales = "free_y", space = "free_y") +
    scale_fill_gradientn(colours = c("#F3F3F3", "#B7D7DA", "#4B9CA8", "#245366"), name = "UCell") +
    scale_size_area(max_size = 1.15, breaks = c(3, 9, 12), name = "Patients") +
    scale_x_discrete(labels = fl) +
    scale_y_discrete(labels = ll) +
    labs(title = "Program sources by bulk context", x = NULL, y = NULL) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      strip.text.y = element_text(angle = 270),
      legend.key.height = unit(1.4, "mm"), legend.key.width = unit(2, "mm"),
      legend.spacing.x = unit(0.5, "mm"), legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0, "mm"),
      panel.spacing.y = unit(1.8, "mm"), plot.margin = margin(1, 1, 0, 4)
    ) +
    guides(fill = guide_colourbar(barwidth = unit(13, "mm"), barheight = unit(1.3, "mm"), title.position = "left"),
           size = guide_legend(nrow = 1))
  
  return(pb)
}

# Panel C: Three Overlap-Guided Featureplot Pairs
draw_panel_C <- function() {
  f <- rdata("featureplot_selected_cells.csv.gz")
  prs <- c("Pattern_P1", "WGCNA_MEpink", "Pattern_P3", "WGCNA_MEbrown", "Pattern_P6", "WGCNA_MEred")
  short <- c(Pattern_P1 = "P1", WGCNA_MEpink = "pink", Pattern_P3 = "P3", WGCNA_MEbrown = "brown", Pattern_P6 = "P6", WGCNA_MEred = "red")
  
  feature <- function(pr) {
    dd <- f %>% mutate(Value = .data[[pr]]) %>% arrange(Value)
    ggplot(dd, aes(UMAP1, UMAP2)) +
      point_rast(colour = "#DEDEDE", size = 0.10, raster.dpi = 300) +
      point_rast(data = filter(dd, !is.na(Value)), aes(colour = Value), size = 0.10, raster.dpi = 300) +
      scale_colour_gradientn(colours = c("#E7E7E7", "#8FC4C8", "#317F8C", "#123E50"),
                             limits = range(dd$Value, na.rm = TRUE), name = short[pr]) +
      coord_equal(xlim = c(-10.8, 22.1), ylim = c(-13.2, 20.0), expand = FALSE) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0), legend.position = "right",
            legend.key.width = unit(1.5, "mm"), legend.key.height = unit(10, "mm"),
            title.position = "top", title.hjust = 0.5)
  }
  
  heads <- c("C1 context | P1-pink\n177 shared; 73% P1 / 66% pink",
             "C2 context | P3-brown\n101 shared; 65% P3 / 21% brown",
             "C3 context | P6-red\n65 shared; 71% P6 / 21% red")
  
  pairs <- lapply(1:3, function(i) {
    (feature(prs[2 * i - 1]) / feature(prs[2 * i])) +
      plot_annotation(title = heads[i],
                      theme = theme(plot.title = element_text(size = 5.6, face = "bold", hjust = 0.5, lineheight = 1.1),
                                    plot.margin = margin(0, 0, 0, 0)))
  })
  
  pc <- wrap_plots(lapply(pairs, function(p) wrap_elements(full = p)), nrow = 1) +
    plot_annotation(title = "Overlap-guided pattern-module pairs across the three bulk contexts",
                    theme = theme(plot.title = element_text(size = 6.5, face = "bold"), plot.margin = margin(2, 1, 1, 8)))
  return(pc)
}

# Panel D: Six Patient Pairs Comparison
draw_panel_D <- function() {
  d  <- rdata("D_patient_pairs.csv")
  ss <- rdata("D_existing_statistics.csv")
  cohorts <- c("GSE151530", "GSE138709")
  shapes  <- c(GSE151530 = 16, GSE138709 = 17)
  short   <- c(Pattern_P1 = "P1", WGCNA_MEpink = "pink", Pattern_P3 = "P3", WGCNA_MEbrown = "brown", Pattern_P6 = "P6", WGCNA_MEred = "red")
  
  plots <- list()
  for (j in 1:6) {
    q <- filter(d, PanelOrder == j)
    st <- filter(ss, PanelOrder == j)
    pr <- unique(q$ProgramID)
    other <- unique(q$Other)
    ti <- unique(q$tissue)
    
    v <- bind_rows(
      transmute(q, PanelOrder, ProgramID, PatientKey, patient, cohort, tissue, Compartment = "Epithelial", X = 1, Score = EpithelialScore),
      transmute(q, PanelOrder, ProgramID, PatientKey, patient, cohort, tissue, Compartment = other, X = 2, Score = OtherScore)
    )
    v$cohort <- factor(v$cohort, levels = cohorts)
    
    qt <- vapply(cohorts, function(co) {
      r <- st[st$cohort == co, ]
      if (nrow(r) == 0) return("")
      if (is.na(r$FDR)) sprintf("%s: n=%d; q=NE", sub("GSE", "", co), r$NPatients)
      else sprintf("%s: n=%d; q=%.3f", sub("GSE", "", co), r$NPatients, r$FDR)
    }, character(1))
    qt <- paste(qt[nzchar(qt)], collapse = "\n")
    
    plots[[j]] <- ggplot(v, aes(X, Score, group = PatientKey)) +
      geom_line(colour = "#9DA8AC", linewidth = 0.25) +
      geom_point(aes(shape = cohort), colour = "#304F5C", size = 1.35, stroke = 0.25, show.legend = (j == 3)) +
      scale_shape_manual(values = shapes, limits = cohorts, drop = FALSE, name = NULL) +
      scale_x_continuous(breaks = c(1, 2), labels = c("Epi.", if (other == "Hepatocyte") "Hep." else "Mye."),
                         limits = c(0.75, 2.25), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, max(v$Score) * 1.1), breaks = scales::pretty_breaks(n = 3), expand = expansion(mult = c(0, 0.02))) +
      labs(title = short[pr], subtitle = paste0(ti, "; n=", n_distinct(v$PatientKey)), caption = qt,
           x = NULL, y = if (j == 1) "UCell score" else NULL) +
      theme(
        plot.title = element_text(size = 5.7, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 5, hjust = 0.5),
        plot.caption = element_text(size = 5, hjust = 0, lineheight = 1.05),
        legend.key.size = unit(2, "mm"), legend.margin = margin(0, 0, 0, 0), plot.margin = margin(1, 2, 1, 1)
      )
  }
  
  pd <- wrap_plots(plots, nrow = 1, guides = "collect") +
    plot_annotation(
      title = "Patient support for the same six representative programs",
      subtitle = "C1-context pair: adjacent epithelial-hepatocyte. C2/C3-context pairs: tumor epithelial-myeloid.",
      caption = "Each line joins one patient within tissue. Cohort shapes; q: original full-family BH. Epi., epithelium; Hep., hepatocyte; Mye., myeloid.",
      theme = theme(plot.title = element_text(size = 6.5, face = "bold"), plot.subtitle = element_text(size = 5),
                    plot.caption = element_text(size = 5, hjust = 0), plot.margin = margin(2, 1, 1, 8))
    ) & theme(legend.position = "bottom")
  
  return(pd)
}

# Panel E: Tumor Compartments Grid Across Patients
draw_panel_E <- function() {
  e <- rdata("E_patient_compartment_grid.csv")
  cohorts <- c("GSE151530", "GSE138709")
  short   <- c(Pattern_P1 = "P1", WGCNA_MEpink = "pink", Pattern_P3 = "P3", WGCNA_MEbrown = "brown", Pattern_P6 = "P6", WGCNA_MEred = "red")
  prs     <- c("Pattern_P1", "WGCNA_MEpink", "Pattern_P3", "WGCNA_MEbrown", "Pattern_P6", "WGCNA_MEred")
  
  pats <- e %>% distinct(PatientKey, patient, PatientOrder) %>% arrange(PatientOrder)
  e$PatientKey <- factor(e$PatientKey, levels = pats$PatientKey)
  e$ProgramID  <- factor(e$ProgramID, levels = rev(prs))
  e$CellType   <- factor(e$CellType, levels = c("Epithelial", "Myeloid", "T_NK"))
  e$cohort     <- factor(e$cohort, levels = cohorts)
  
  pe <- ggplot(e, aes(PatientKey, ProgramID)) +
    geom_tile(aes(fill = Score), colour = "white", linewidth = 0.13) +
    geom_text(data = filter(e, is.na(Score)), label = "-", size = 5 / ggplot2::.pt, colour = "#999999") +
    facet_grid(CellType ~ cohort, scales = "free_x", space = "free_x",
               labeller = labeller(CellType = c(Epithelial = "Epithelium", Myeloid = "Myeloid", T_NK = "T / NK"))) +
    scale_fill_gradientn(colours = c("#F4F4F4", "#B7D7DA", "#4B9CA8", "#245366"),
                         limits = c(0, max(e$Score, na.rm = TRUE)), name = "UCell", na.value = "#E4E4E4") +
    scale_y_discrete(labels = short) +
    scale_x_discrete(labels = setNames(pats$patient, pats$PatientKey), drop = TRUE) +
    labs(
      title = "Tumor-only distributions of the same six programs across compartments",
      subtitle = sprintf("16 patients; %d/%d estimable. Clean eligibility; equal-library patient means. Grey/dash: unsupported.",
                         sum(!is.na(e$Score)), nrow(e)),
      x = NULL, y = NULL
    ) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right", legend.key.height = unit(3, "mm"), legend.key.width = unit(2, "mm"),
      strip.text.y = element_text(angle = 0), plot.margin = margin(2, 1, 1, 8),
      panel.spacing.y = unit(1.3, "mm")
    )
  
  return(pe)
}

# ------------------------------------------------------------------------------
# 3. Execution & Rendering
# ------------------------------------------------------------------------------
save_panel <- function(stem, plot_obj, w, h) {
  pdf_path <- file.path(OUT_DIR, "panels", paste0(stem, ".pdf"))
  png_path <- file.path(OUT_DIR, "panels", paste0(stem, ".png"))
  
  pdf(pdf_path, width = w / 25.4, height = h / 25.4, family = "Arial")
  print(plot_obj)
  dev.off()
  
  ragg::agg_png(png_path, width = w, height = h, units = "mm", res = 300, background = "white")
  print(plot_obj)
  dev.off()
  
  cat(sprintf("[Rendered] %s (%.1f x %.1f mm)\n", stem, w, h))
}

cat("Rendering individual panels for Figure 3...\n")
save_panel("Panel_A_atlas", draw_panel_A(), 90, 90)
save_panel("Panel_B_grouped_sources", draw_panel_B(), 78.5, 61.4)
save_panel("Panel_C_three_pairs", draw_panel_C(), 174, 53)
save_panel("Panel_D_six_patient_pairs", draw_panel_D(), 174, 58)
save_panel("Panel_E_tumor_compartments", draw_panel_E(), 174, 63)

cat("Generating full assembled Figure 3 preview...\n")
composite_png <- file.path(OUT_DIR, "Figure_3_composite.png")
ragg::agg_png(composite_png, width = 182, height = 276, units = "mm", res = 300, background = "white")

canvas <- function(w, h) {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, w), yscale = c(0, h)))
}
vt <- function(text, x, y, size = 8, bold = TRUE, col = "#222222") {
  grid.text(text, unit(x, "native"), unit(y, "native"), just = "left",
            gp = gpar(fontfamily = "Arial", fontsize = size, fontface = if (bold) "bold" else "plain", col = col))
}

canvas(182, 276)

panels_pos <- data.table(
  Letter = c("A", "B", "C", "D", "E"),
  File   = c("Panel_A_atlas.png", "Panel_B_grouped_sources.png",
             "Panel_C_three_pairs.png", "Panel_D_six_patient_pairs.png",
             "Panel_E_tumor_compartments.png"),
  X      = c(4, 98, 4, 4, 4),
  Top    = c(4, 15, 96, 151, 211),
  W      = c(90, 78.5, 174, 174, 174),
  H      = c(90, 61.4, 53, 58, 63)
)

for (i in seq_len(nrow(panels_pos))) {
  p_file <- file.path(OUT_DIR, "panels", panels_pos$File[i])
  if (file.exists(p_file)) {
    im <- png::readPNG(p_file)
    grid.raster(
      im,
      x = unit(panels_pos$X[i], "native"),
      y = unit(276 - panels_pos$Top[i], "native"),
      width = unit(panels_pos$W[i], "native"),
      height = unit(panels_pos$H[i], "native"),
      just = c("left", "top")
    )
    vt(panels_pos$Letter[i], panels_pos$X[i] - 3, 276 - panels_pos$Top[i] + 1, 8, TRUE)
  }
}
popViewport()
dev.off()

cat(sprintf("[Complete] All Figure 3 panels and composite generated in %s\n", OUT_DIR))
