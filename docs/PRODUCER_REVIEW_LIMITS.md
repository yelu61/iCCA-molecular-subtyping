# Source-owner review and remaining provenance limits

All 271 frozen dependencies have received a source-owner review. The index is
`SOURCE_LINEAGE_20261110.csv`. Its filename identifies the audit baseline. The full
internal review log is retained in the research workspace. Each row distinguishes code
calls, copies, historical hashes, observed notebook calls and isolated executed
reconstruction. Full original-input-to-paper execution is not certified.

The source review recovered 51 missing owners and literal code helpers, bringing
the reviewed implementation inventory to 166 files. Two fixed clinical protocols
were subsequently included, giving 168 implementation files. Inclusion records the source owner
or current counterpart; it does not assert that an unbound historical version
was the code executed originally. The public analysis catalog has 24 entries.

## Historical code-version gaps

Nineteen source records have unresolved original code versions, grouped below.
Their frozen tables remain preserved; no scientific result was substituted.

| Group | Affected records | Limitation |
|---|---:|---|
| Step1 finalization | 1 | Recorded `finalize_step1_audit.R` MD5 is c385f38983ef8908106b94ed9d61c6c7; located copies and the available Git path-history blob have 35d81e3a52bacc0ce92951c718d733f0 |
| Figure1 original v2 | 11 | Original discovery/refinement script hashes differ from current copies; same-content later typography outputs have exact code/output bindings, which do not establish the original v2 execution |
| Figure2 transcriptomic v8 | 7 | Historical suite records producer MD5 fe2738…; located script copies have 50ebadd8…; output identity does not recover the executed script version |

## Original writers not recovered

Four computed/copied artifacts have identified content or downstream uses but
no recovered original writing event:

- `governance/data/umap_coordinates_top2000.csv`.
- Figure2 v34 `tables/C_optional_heterogeneity_audit.csv`: all 11 rows and eight
  fields match independently joined upstream tables.
- Figure6 v20 `tables/joint_performance_ci.csv`: exact upstream counterpart found.
- Figure6 v20 `tables/Supplement_reference_sensitivity.csv`: exact upstream copy found.

The row-level paths, code and retention limits are in `PROVENANCE_GAPS.csv`.
Authored colors, configurations, selection contracts and display plans are
tracked as inputs/annotations; they do not require a statistical producing run.

Many other sources have reviewed direct writing calls without a complete
historical execution binding. These retain that explicit status. Completing the
review does not remove this limitation or authorize analysis-history deletion.

## Execution and portability boundaries

The checked frozen-input route, score-input recovery, five Figure3 transformations
and four Figure6 candidate-table joins retain their scoped execution checks.
The 51 added scripts were parsed and hash-checked against their original files;
their full analysis workflows were not executed again.

Two retained historical packaging sources, `package_fig6_storyline_v19.py` and
`package_fig6_translational_v17.py`, contain local figure-skill QA invocations.
They are provided byte-for-byte as source evidence and are not portable supported
public entry points. The frozen rendering route and candidate table-join stage
do not call them. Other upstream stages still require their external data,
configuration and environment contracts.
