#!/usr/bin/env Rscript
# ==============================================================================
# Script: run_all_figures.R
# Purpose: Master runner to generate all 6 primary publication figures sequentially.
# ==============================================================================

cat("==============================================================================\n")
cat(" iCCA Molecular Subtyping - Publication Figure Generation Pipeline\n")
cat("==============================================================================\n\n")

start_total <- Sys.time()

# Determine script directory
script_dir <- if (dir.exists("8-Visualization/scripts")) {
  "8-Visualization/scripts"
} else if (dir.exists("figures")) {
  "figures"
} else {
  "."
}

figures <- 1:6
results <- data.frame(
  Figure = paste0("Figure ", figures),
  Script = paste0("plot_figure", figures, ".R"),
  Status = "Pending",
  Duration_Sec = NA_real_,
  stringsAsFactors = FALSE
)

for (i in figures) {
  fig_name <- paste0("Figure ", i)
  script_name <- paste0("plot_figure", i, ".R")
  script_path <- file.path(script_dir, script_name)
  
  cat(sprintf("[ %s / 6 ] Executing %s...\n", i, script_name))
  t0 <- Sys.time()
  
  res <- tryCatch({
    source(script_path, local = new.env())
    "Success"
  }, error = function(e) {
    cat(sprintf("  [ERROR in %s]: %s\n", fig_name, e$message))
    "Failed"
  })
  
  t1 <- Sys.time()
  dur <- round(as.numeric(difftime(t1, t0, units = "secs")), 1)
  
  results$Status[i] <- res
  results$Duration_Sec[i] <- dur
  cat(sprintf("  --> %s completed in %.1f seconds (%s)\n\n", fig_name, dur, res))
}

end_total <- Sys.time()
total_dur <- round(as.numeric(difftime(end_total, start_total, units = "secs")), 1)

cat("==============================================================================\n")
cat(sprintf(" Figure Generation Finished in %.1f seconds.\n", total_dur))
cat(" Summary:\n")
print(results)
cat("==============================================================================\n")

if (any(results$Status != "Success")) {
  quit(status = 1)
} else {
  cat("\nAll 6 figures successfully reproduced.\n")
}
