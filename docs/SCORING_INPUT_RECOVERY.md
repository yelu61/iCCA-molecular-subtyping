# Reconstructing the single-cell method-comparison input

The current `singlecell_programs` stage reads the prior rank/control-score matrix
to compare scoring methods. A later display run does not remove that dependency.
On 2026-11-10 the original scoring implementation was executed in an independent
directory. Its generated `cell_scores.npz` matched the pre-deletion artifact
byte-for-byte (11,545,751 bytes):

`5c5485d13170631eb2be493d1a420fe86ee354f03b97642a85a95c2d2470028f`

The exact artifact was restored locally. No alternate program scores or scientific
results were substituted. This validates this scoring step, not the complete
single-cell or raw-data-to-paper pipeline.

## Inputs and implementation

`config/score_baseline_inputs.json` records the three external input files, sizes
and SHA256 values: the reviewed annotated H5AD, program definitions and membership.
These inputs are not in the frozen plotting bundle. Supply them at the recorded
relative paths in an independent project directory; obtain them through the
upstream data workflow. No public input download is claimed here.

The original `score.py`, analysis contract and `r3a_annotation_evidence.py` are
retained under `_implementation/project/`. The scoring script reads the literal
`PANELS` assignment from the annotation script; it does not rerun annotation.
The input H5AD is opened read-only. An existing `scoring_complete.json` stops a
rerun. Use a fresh working directory and preserve original inputs.

## Commands

Use a Python environment matching the recorded single-cell environment. The
observed recovery used Python 3.10.17 in the original omicverse environment,
with anndata, NumPy, pandas and SciPy; installed versions are recorded in
`environment/score_recovery_python.json`.

```sh
# From the code repository, using that environment's Python:
python run.py analysis --stage singlecell_score_baseline --project-root /path/to/new_project
# Provide and verify the three external inputs at their recorded paths.
mkdir -p /path/to/new_project/7-SingleCell/runs/20261102_figure3_program_localization_v1/tables
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run.py analysis --stage singlecell_score_baseline --project-root /path/to/new_project --execute
```

Verify the generated NPZ against the hash above before using it for the historical
method comparison. A hash mismatch requires investigation of inputs, code,
environment and numeric arrays; it is not permission to replace the expected
artifact silently. The other upstream stages still require their own external
interfaces and execution validation.
