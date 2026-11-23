# Notebook execution roles

Updated 2026-11-23 for `v2026.11.23-preparation`.

| Notebook | Role and execution scope |
|---|---|
| 0-data-qc.ipynb | Project data/QC workflow; external protected inputs required |
| 1-clustering.ipynb | Project clustering workflow; external expression objects and current publication repairs required |
| 2-transcriptomic-characterization.ipynb | Bulk characterization; project governance, corrected results and scripts required |
| 3-duct-type-analysis.ipynb | Observed duct-type analysis; project inputs required |
| 4-multiomic.ipynb | Multi-omics analyses; project inputs and scripts required |
| 5-classifier-validation.ipynb | Current classifier tables and optional refresh commands; project model and script dependencies required |
| 6-clinical-survival.ipynb | Current clinical tables and optional refresh commands; project inputs/scripts required |
| 7-single-cell-validation.ipynb | Current contract checks and project analysis outputs; external Step7 configuration, tables and scripts required |
| 8-visualization.ipynb | Public staged-data drawing entry; also supports the original research directory layout |

The current `pipelines/` Steps 0–7 expect the research directory layout and external inputs. Selected analysis owners and dependencies are retained under `_implementation/project/`; use `run.py analysis --list` and its staging interface for a portable entry. They were not executed end to end for this freeze. Step7 contract assertions match the current cohort-Harmony representation with unintegrated sensitivity, Q1 gate and patient-native CNV run. Historical cell outputs are not validation of the current code.

Step8 discovers either `figures/` or the research project's `8-Visualization/` layout, passes explicit data/output paths, and requires both a non-empty PDF and PNG before reusing output. It discovers R through `ICCA_R_HOME`, `R_HOME` or `R RHOME` and supports a pinned runtime. Python path/command checks are separate from the real six-script R redraw.

The nine public notebooks remain paired with the maintained research notebooks. Details of the retained 2026-11-12 CLI belong to [the historical notebook record](history/20261112/NOTEBOOK_STATUS.md).

The original `notebooks/Reproduce.ipynb` and the nine `_implementation/project/*ipynb` notebooks are also retained unchanged from the previous public baseline. They support the original frozen-input and analysis staging routes. `pipelines/` supplies the maintained research-facing notebooks, including the new standalone drawing interface.
