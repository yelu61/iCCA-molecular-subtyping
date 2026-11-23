# Notebook roles in the manuscript workflow

Updated 2026-11-12.

| Stage notebook | Current role | Selected execution owner |
|---|---|---|
| 0-data-qc.ipynb | Upstream source and protected interfaces | step0 |
| 1-clustering.ipynb | Upstream clustering/functions retained at original cell positions | step1 and current publication repair scripts |
| 2-transcriptomic-characterization.ipynb | Corrected bulk characterization source | step2 and the WGCNA/TIDE/current figure-specific owners |
| 3-duct-type-analysis.ipynb | Upstream duct-analysis record | duct |
| 4-multiomic.ipynb | Upstream multi-omic record | multiomic, mofa, bridge |
| 5-classifier-validation.ipynb | Current Python CLI frontend | classifier114, compact, external_programs |
| 6-clinical-survival.ipynb | Current Python CLI frontend | clinical_full_subtype, clinical_duct_stratified, clinical |
| 7-single-cell-validation.ipynb | Current Python CLI frontend | singlecell_programs, singlecell_display_sources, singlecell_cnv_cases, singlecell_score_baseline |
| 8-visualization.ipynb | Current frozen-figure frontend | run.py verify / figures, then check_results.py |

These files are under `_implementation/project`. `notebooks/Reproduce.ipynb` provides the figure workflow directly from the repository. The current analysis catalog contains 24 entries. Phase0–4 code-cell positions and source were retained. Stored historical outputs were cleared. Phase5–8 now use Python3 and the current CLI. Their original executed notebooks are retained in the research workspace archive with recorded SHA256 values.

Choose a source stage and supply its external inputs in an independent directory for upstream work. The four new frontends default to command preparation until execution is explicitly selected in their parameter cell. Stage0–4 retain their R kernels and upstream execution behavior. The selected six-figure implementation is indexed per panel in final_panel_manifest.csv and SOURCE_LINEAGE_20261110.csv.

The fixed 114-gene consensus model is primary. The fixed 17-gene development candidate retains its exploratory scope. Single-cell source and identity restrictions follow the selected run protocols. Historical notebook execution outputs are not used by the current figure reconstruction route.
