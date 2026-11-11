# iCCA molecular states

![Project overview](assets/project_overview.png)

Multi-cohort discovery and characterization of intrahepatic cholangiocarcinoma molecular states. The manuscript uses C1/C2/C3 as stable identifiers. Discovery includes 586 RNA-sequencing tumors, with 497 principal-state assignments and 89 intermediate profiles. A fixed 114-gene consensus classifier projects 506 microarray tumors from five studies. Genomic, proteomic, single-cell and clinical evidence retains its assay-specific population and patient-level interpretation.

## Reproduce the manuscript figures

The supported frozen-input workflow executes 37 production renderers and three final panel drivers, followed by 100 quantitative asset checks. Software requirements and observed package versions are in `environment/`.

```sh
python3 run.py verify --data-root /path/to/data_bundle
python3 run.py figures --data-root /path/to/data_bundle --output /path/to/new_output
python3 check_results.py --run-dir /path/to/new_output --data-root /path/to/data_bundle
```

`notebooks/Reproduce.ipynb` runs the same commands. Choose a fresh output directory outside the code and input trees. R is detected through `R RHOME`, with `--r-home` and `ICCA_R_HOME` available for explicit configuration.

The companion input archive is `iCCA_frozen_inputs_20261104_v2.zip` (474,136,575 bytes, SHA256 `0e7182c4a0b30c3cc907395096157324481de1555d78faf4168b1d66785ca762`). It unpacks to `data_bundle/`. Its public deposit is pending. The required 756 input identities are in `config/input_manifest.json`. This source repository excludes the large input archive.

## Analysis and notebook entry points

```sh
python3 run.py analysis --list
python3 run.py analysis --stage classifier114 --project-root /path/to/independent_project
```

The second command stages the selected implementation and prints the command. Supply the documented module inputs and add `--execute` for analysis. The catalog contains 24 entries. The 168 implementation files and their dependencies are identified in `IMPLEMENTATION_MANIFEST.csv` and `docs/DEPENDENCY_MAP.csv`.

Stage0–4 notebooks retain upstream code. Stage5–8 are current Python CLI frontends for classifier, clinical, single-cell and figure workflows. Notebook roles, kernels and preserved source boundaries are documented in [NOTEBOOK_STATUS.md](docs/NOTEBOOK_STATUS.md).

## Provenance and verification

- `docs/final_panel_manifest.csv` maps selected manuscript assets to drawing code and frozen sources.
- `docs/SOURCE_LINEAGE_20261110.csv` records the review of 271 frozen dependencies, including historical code-version and producing-event gaps.
- `docs/PRODUCER_REVIEW_LIMITS.md` describes the remaining provenance scope.
- `docs/VERIFICATION_20261112.json` records the synchronized notebook and figure checks.
- `docs/upstream_external_interfaces.csv` lists processed and external resources required by upstream stages.

The complete raw-input-to-manuscript pipeline was not rerun during this preparation. The verified frozen-input figure workflow and the separately checked reconstruction stages retain their stated scopes. External subtype assignments are projections. Drug predictions, PRISM comparisons and TIDE analyses support experimental hypotheses.

## Repository contents

| Location | Purpose |
|---|---|
| `run.py`, `check_results.py` | Canonical execution and validation |
| `notebooks/` | Figure-reproduction notebook |
| `config/` | Stage catalog, selected jobs and input/output contracts |
| `_implementation/project/` | Selected owners and required source dependencies |
| `_implementation/runner/`, `_implementation/final_panels/` | Figure execution adapters and final drawing drivers |
| `docs/`, `environment/` | Provenance, scope and software inventories |
| `assets/` | Conceptual project graphic |
| `tests/` | Frozen-input isolation regression check |

The graphic is an AI-generated conceptual illustration. Quantitative results and selected data figures are reproduced by the analysis code. `assets/project_overview_prompt.txt` records the generation brief.

## Citation, licence and history

Code is provided under the MIT licence in `LICENSE`. Third-party inputs retain their source-specific access and reuse conditions. `CITATION.cff.template` contains the software-citation fields awaiting real author metadata. An archival DOI has not been assigned.

The May 2026 version remains available in Git history. The research workspace retains original executed notebooks and historical result files. The current repository tree contains the manuscript implementation and reproduction materials.
