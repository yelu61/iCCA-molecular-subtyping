# Environments

The preparation snapshot is redrawn with an existing R 4.5.3 arm64 system library. `environment_R_installed_packages.csv` and session records describe observed installations; they are not a solved cross-platform lock. The current figure dependencies are declared in `install_packages.R` and `environment.yml`, including ggrastr, svglite and tidyr. A fresh Conda solve or package installation was not used as validation for this snapshot.

Typical setup:

```bash
Rscript environment/install_packages.R
# Alternatively, prepare the declared Conda environment:
conda env create -f environment/environment.yml
conda activate icca-subtyping
```

Run `Rscript figures/run_all_figures.R` from the repository root. For a pinned Unix R installation, set `R_HOME` and use its `bin/exec/R` when available: some macOS framework R/Rscript wrappers redirect to the current system version. Steps5/6/8 notebook helpers support `ICCA_R_HOME` or `R_HOME`, then fall back to R on PATH.

Python 3 and Jupyter are needed for the notebook interface; the regression suite uses the Python standard library. The six-script R redraw does not require the upstream scRNA/MOFA/TIDE environments. Upstream modules still require their project inputs and modality-specific dependencies. Nothing is installed automatically during the drawing entry or tests.

The historical inventories and upstream package requirements remain useful provenance. They do not establish that all upstream analyses can execute in this public clone.
