# Reconstruct Figure 6 candidate-evidence tables

The `candidate_evidence_tables` stage runs the original descriptive join script
for the selected candidate-program crosswalk, complete program evidence,
candidate evidence summary and PRISM observation grid. On 2026-11-11, both the
original-script route and the public entry reproduced all four selected CSVs
byte-for-byte. The eight historical input hashes matched. Observed Python was
3.9.6; this script uses only the Python standard library.

Supply the eight external files listed in `config/candidate_evidence_inputs.json`
at their project-relative paths in an independent directory and verify SHA256.
They are processed upstream interfaces, not bundled raw data or promised public
downloads. The frozen-figure route does not require this optional reconstruction.

```sh
python run.py analysis --stage candidate_evidence_tables --project-root /path/to/new_project
mkdir -p /path/to/new_project/8-Visualization/refinements/fig6_evidence_20261103_v21/qa/candidate_audit
python run.py analysis --stage candidate_evidence_tables --project-root /path/to/new_project --execute
```

See `candidate_source_rebuild.json` and `candidate_public_stage_verification.json`
for output/input hashes and execution scope. The operation joins existing results
and retains an explicitly retrospective mechanism crosswalk. It does not fit
models, select drugs or validate mechanisms. Missing observations, projected
cell-line identities and the PRISM discordance are preserved.
