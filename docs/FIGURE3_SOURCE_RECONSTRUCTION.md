# Figure 3 source-table reconstruction

The `singlecell_display_sources` stage exposes the original selected Figure 3
table-preparation script. Its five manuscript inputs were rebuilt on 2026-11-10
in an independent directory, then rebuilt through the public `run.py` entry.
All five matched the frozen content exactly (decompressed CSV bytes for gzip).
The 12 external input hashes were unchanged afterward.

This is a transformation of existing reviewed UCell and patient tables. It does
not repeat annotation, UCell scoring, statistical fitting or biological acceptance.
The 288 patient/compartment/program slots include 199 estimates and 89 missing
values; the 42 contrast pairs and coverage masks remain unchanged.

## Running this stage

Use the recorded single-cell Python environment (Python 3.10.17, pandas 2.2.3,
NumPy 1.26.4). Supply the 12 files listed in
`config/figure3_source_inputs.json` at their exact relative paths in an independent
project directory and verify their SHA256 values. They are external upstream
interfaces, not promised downloads or all members of the frozen plotting bundle.

```sh
python run.py analysis --stage singlecell_display_sources --project-root /path/to/new_project
mkdir -p /path/to/new_project/8-Visualization/refinements/fig3_manuscript_revision_20261103_v4/main/source_data
mkdir -p /path/to/new_project/8-Visualization/refinements/fig3_manuscript_revision_20261103_v4/main/review
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run.py analysis --stage singlecell_display_sources --project-root /path/to/new_project --execute
```

The script writes its own numeric checks to `main/review/numeric_QA.json`.
`docs/FIGURE3_SOURCE_VERIFICATION_20261110.json` records the tested output hashes,
input checks and numeric checks. Reproduction from raw single-cell input still
requires the earlier stages and their eligibility restrictions.

## Selected outputs

| File in `main/source_data/` | Operation |
|---|---|
| `D_existing_statistics.csv` | Select existing UCell contrasts and assign panel order |
| `D_patient_pairs.csv` | Select the matching 42 patient pairs |
| `E_patient_compartment_grid.csv` | Summarize eligible libraries by patient; retain missing grid slots |
| `featureplot_selected_cells.csv.gz` | Apply the existing availability and support masks |
| `dotplot_summary.csv` | Copy the reviewed upstream summary without recalculation |
