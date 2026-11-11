# Reproduction scope and release preparation

Updated 2026-11-12.

## Checked execution

The current frozen-input workflow completed 37 production renderers and three final panel drivers. All 100 quantitative asset checks passed. The 756 input hashes matched before and after execution. The input-isolation regression test passed after session inputs were changed to independent file copies. Detailed current results are in `VERIFICATION_20261112.json`.

The selected implementation contains 168 files and the analysis catalog contains 24 entries. Stage5–8 notebooks use the current Python CLI. Stage0–4 preserve their upstream code and code-cell positions. Original executed notebook snapshots are retained in the research workspace.

Separate executed reconstruction checks cover the score-comparison input, five single-cell display tables and four candidate-evidence tables. Their guides and source-input contracts remain in this repository.

## Provenance scope

The source-owner review covers 271 frozen dependencies. Three historical code-version groups affect 19 records, and four original writing events remain unresolved. `PRODUCER_REVIEW_LIMITS.md` and `PROVENANCE_GAPS.csv` describe these gaps. Full upstream execution and author scientific acceptance remain outstanding.

## Data and citation

The companion frozen-input archive is `iCCA_frozen_inputs_20261104_v2.zip`, with SHA256 `0e7182c4a0b30c3cc907395096157324481de1555d78faf4168b1d66785ca762`. Its public deposit is pending. External upstream interfaces are listed separately in `upstream_external_interfaces.csv`.

The repository uses the MIT licence. Software-author metadata await completion in `CITATION.cff.template`. An archival DOI has not been assigned.

## Repository selection

This source tree contains selected manuscript code and its required dependencies, the reproduction entry points, software inventories, panel mapping, source-lineage summaries, verification results and the conceptual overview graphic. Full local review logs, manuscript drafts, slide drafts, historical result directories and large data archives are retained in the research workspace. Original public history remains accessible through Git.
