# Public data and terminology

Updated 2026-11-23.

| Public input | Meaning |
|---|---|
| `data/subtype_master_table.csv` | 586-row assignment extract: Sample, Study, FinalSubtype, SubtypeClass. Principal states total 497; intermediate profiles total 89. Patient covariates are excluded. |
| `data/study_color_dictionary.csv` | Central study palette; plotting scripts require this file and do not replace it with a local study map. |
| `data/subtype_label_color_dictionary.csv` | Central subtype label and colour definitions. |
| `figures/data/` | 72 staged inputs for six main figure drafts. File sizes and hashes are recorded in INPUT_FILES.csv. |
| `figures/data/figure2/A_pattern_model.rds` | Pattern model and ordering data. main/all retain Sample, Study and FinalSubtype only; z, dimnames, membership, counts and pattern_order are unchanged. |

The research workspace's full governance table remains the canonical source for downstream patient analyses. The public assignment extract must not be treated as a clinical dataset. The pattern RDS and the assignment CSV were minimized in this release without changing plot inputs or their row order.

C1/C2/C3 are stable identifiers. C1-Metabolic and C3-Proliferative are bulk-derived working descriptions. C2 has no locked immune or ductal-lineage name. Observed pathology, inferred duct type and external subtype projections are separate concepts. Drug-response predictions remain hypotheses.

The current single-cell analysis set is 34,227 cells, 17 patients and 22 libraries across two cohorts; 67,898 is the broader QC candidate layer. The current patient-native CNV analysis has 42 strict candidates across three patients, with only one patient reaching 20 candidates. These are not released malignant identities. Q1 broad-compartment localization is computed; Q2 identity and Q3 estimability remain restricted.

[The historical dictionary](history/20261112/DATA_DICTIONARY.md) retains upstream interfaces and older Step7 states as dated provenance. Those paths are research-workspace paths and need not exist in a public clone.
