#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "1-Clustering", "publication", "tables")
figure_dir <- file.path(root, "1-Clustering", "publication", "figures")
state_levels <- c("C1", "C2", "C3", "C1_C2", "C1_C3", "C2_C3")

confusion <- fread(file.path(table_dir, "step1_locked_vs_refit_confusion_long.csv"))[Method == "SixState"]
confusion <- merge(
  CJ(Locked = state_levels, Refit = state_levels, unique = TRUE),
  confusion[, .(Locked, Refit, Samples)],
  by = c("Locked", "Refit"), all.x = TRUE
)
confusion[is.na(Samples), Samples := 0L]
confusion[, `:=`(
  Locked = factor(Locked, levels = state_levels),
  Refit = factor(Refit, levels = state_levels)
)]
p_confusion <- ggplot(confusion, aes(Refit, Locked, fill = Samples)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Samples > 0, Samples, "")), size = 3.1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C") +
  coord_equal() +
  labs(
    title = "Locked versus independently refitted six-state assignments",
    subtitle = "Exact reproduction: 50-run NMF and seed-123 KMeans changed 0 of 586 samples",
    x = "Refitted assignment", y = "Locked assignment", fill = "Samples"
  ) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

comparison <- fread(file.path(table_dir, "step1_locked_vs_refit_sample_comparison.csv"))
uncertainty <- rbindlist(list(
  comparison[, .(FinalSubtype = LockedFinalSubtype, Method = "NMF membership margin", Margin = NMFMembershipMargin)],
  comparison[, .(FinalSubtype = LockedFinalSubtype, Method = "KMeans distance margin", Margin = KMeansDistanceMargin)]
))
uncertainty[, FinalSubtype := factor(FinalSubtype, levels = state_levels)]
palette <- fread(file.path(root, "docs", "subtype_color_dictionary.csv"))
subtype_colors <- setNames(palette$ColorHex, palette$FinalSubtype)
p_uncertainty <- ggplot(uncertainty, aes(FinalSubtype, Margin, fill = FinalSubtype)) +
  geom_boxplot(width = 0.62, outlier.size = 0.55) +
  facet_wrap(~Method, scales = "free_y", nrow = 1L) +
  scale_fill_manual(values = subtype_colors, guide = "none") +
  labs(
    title = "Assignment certainty across locked molecular states",
    subtitle = "Full-data refits reproduced all 586 assignments; margins describe relative within-method certainty",
    x = NULL, y = "Assignment margin"
  ) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

for (extension in c("pdf", "png")) {
  ggsave(file.path(figure_dir, paste0("step1_locked_vs_refit_confusion.", extension)),
         p_confusion, width = 6.4, height = 5.2, dpi = 300)
  ggsave(file.path(figure_dir, paste0("step1_refit_assignment_uncertainty.", extension)),
         p_uncertainty, width = 8.2, height = 4.5, dpi = 300)
}

cat("Step 1 refit figures updated\n")
