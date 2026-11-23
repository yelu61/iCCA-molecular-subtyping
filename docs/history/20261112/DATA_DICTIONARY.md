# Data and terminology dictionary

Last updated: 2026-10-22 (fresh Step7 source/QC fields)

## Canonical identifiers and populations

| Term | Definition | Canonical source |
|---|---|---|
| `Sample` | Harmonized analysis identifier used for cross-module joins | `governance/data/subtype_master_table.csv` |
| `Study` | Source cohort; retained as a covariate/stratum or sensitivity dimension | Step0 metadata and `0-Data/qc/study_color_dictionary.csv` |
| `FinalSubtype` | One of C1, C2, C3, C1_C2, C1_C3, C2_C3 | Governance master table |
| `SubtypeClass` | `Main` for C1/C2/C3; intermediate class otherwise | Governance master table |
| `MainSubtypeForPlot` | Governed principal-state plotting field | Governance master table |
| Discovery universe | 586 RNA-seq samples across seven studies | `0-Data/qc/step0_cohort_inclusion_manifest.csv` |
| Main-state universe | 497 discovery samples: C1=112, C2=223, C3=162 | Governance master table |
| Intermediate universe | 89 samples: 41 C1_C2, 39 C1_C3, 9 C2_C3 | Governance master table |
| Projection universe | 506 curated microarray samples | Step0 inclusion manifest |
| Curated expression universe | 1,092 discovery plus projection samples | Step0 inclusion manifest |

## Subtype naming

| Field/use | Allowed value |
|---|---|
| Stable identifiers | C1, C2, C3 |
| Short plotted descriptions | C1-Metabolic, C2, C3-Proliferative |
| Prose interpretation | C1 metabolic-adaptive; C3 stem-like/proliferative/dedifferentiated |
| Prohibited locked names | C2-Immune, C2-Immune-sensitive, C2-Ductal-lineage |

## Observed, predicted and projected variables

| Concept | Meaning | Use restriction |
|---|---|---|
| Observed duct type | Pathology-derived large-/small-duct label | Clinical association only on observed-label denominator |
| Cross-fitted duct probability | Out-of-fold prediction for eligible labeled samples | Incremental models; never substituted for pathology |
| Full-model duct projection | Prediction in unlabeled samples | Distribution/exploration only |
| Subtype OOF prediction | Nested-CV prediction within locked-label resource | Reproducibility estimate, not independent truth |
| External subtype projection | Classifier application to microarray or DepMap/CCLE | Preservation/projection, not accuracy validation |
| TIDE responder | Algorithmic prediction | Hypothesis-generating, not observed response |
| oncoPredict score | Predicted drug-response score; lower means greater predicted sensitivity | Hypothesis-generating |

## Key analytical interfaces

| Interface | Shape/scope | Role |
|---|---|---|
| `0-Data/qc/0-Data-rnaseq.Rdata` | Genes × samples | Protected Step0 discovery expression input |
| `1-Clustering/data/combined_expr_after_correction.rds` | Genes × 586 samples | Locked corrected discovery expression |
| `1-Clustering/data/sample_info_with_clusters.xlsx` | 586 sample assignments | Step1 locked clustering interface |
| `governance/data/subtype_master_table.csv` | 586 governed rows | Single source of subtype truth downstream |
| `3-DuctType/data/pred.csv` | Duct observed/predicted interface | Fields must remain semantically separate |
| `4-Multiomic/publication/models/step4_mofa_model.hdf5` | Saved MOFA model | Supportive 2022Dong multi-omic model |
| `5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_consensus_v1/models/candidate_models.rds` | 114-gene consensus linear SVM | Primary subtype classifier model (historical 123-gene in `iCCA_Subtype_SVM_Classifier.rds` is superseded) |
| `5-ClassifierValidation/runs/bulk-rnaseq__20261031_step5_minpanel_v1/models/compact_candidate.rds` | 17-gene compact SVM | Exploratory panel candidate |
| `7-SingleCell/config/current_run.json` | Active G0/G1 IDs and stage status | No implicit latest-run or archive fallback |
| `7-SingleCell/runs/20261022_g0_source_rebuild_r3/tables/library_manifest.csv` | Dataset × source library | Source counts, Patient/Specimen/Tissue, explicit hashes; timepoint and treatment Unknown pending review |
| `7-SingleCell/runs/20261023_g1_omicverse_scdblfinder_v1/tables/library_qc_summary.csv` | Dataset × library | Provisional QC views and doublet flags; not accepted identities |
| `docs/Step7_Clean_Restart_20261022.md` | Current Step7 contract | Old source-dependent conclusions retired |
| `8-Visualization/tables/figure_manifest.csv` | Panel contract | Step8 source-to-target map |

## Analysis-specific denominators

| Analysis | Denominator |
|---|---|
| Main bulk transcriptomics/TIDE/oncoPredict | 497 |
| Observed duct classifier | 115: 29 large, 86 small |
| Observed duct by main subtype | 93 |
| DNA principal-state analysis | 308 of 359 assay-profiled tumors |
| Mass-spectrometry concordance | 248 principal states of 292 samples |
| TCGA-CHOL RPPA | 23 principal states of 25 matched tumors; C3 n=3 |
| MOFA | 182 principal states of 204 triple-matched tumors |
| Survival/incremental prognosis | 331 patients, 153 deaths |
| Compact C3 panel | 341 patients, C3 n=106 |
| Step7 source counts | 36 libraries / 77,204 source-called cells; 20 source patients any tissue / 19 with tumor tissue; patient-inference scope pending |
| Step7 QC v1 | 74,174 metric-pass / 67,898 after primary doublet exclusion; pending biological acceptance |
| Step7 doublet diagnostics | scDblFinder primary; 6,276 excluded; Scrublet diagnostic; unresolved status preserved |
| Step7 malignancy / program / composition inference | NOT_RUN; no accepted candidate denominator or association N |
| PRISM iCCA join | 3 cell lines: C1=0, C2=1, C3=2 |

## Current single-cell fields

| Field | Meaning and boundary |
|---|---|
| Dataset / LibraryID / GSM / PatientKey / SpecimenKey / Tissue / Region | Source identities; treatment/timepoint remain Unknown where unverified; tissue origin is not cell malignancy |
| nUMIs / detected_genes / mito_perc / cell_complexity | Metrics from untransformed counts; mito_perc is a fraction, not percentage |
| QC_mt10 / QC_mt15 / QC_mt20 | UMI>500, genes>250 and mitochondrial threshold masks; 10/20 comparisons do not include rerun doublet classification |
| scDblFinderScore / scDblFinderCall / scDblFinderAvailable / DoubletStatus | Primary algorithm evidence; low-support/error states do not mean singlet |
| ScrubletScore / ScrubletCall / ScrubletAvailable / ScrubletStatus | Separately stored diagnostic method; scores not directly comparable across algorithms |
| ExcludedPrimaryDoublet / QC_v1_retained | Explicit exclusion and retained-cell mask for technical QC v1 |
| DoubletMethodDisagreement / DoubletReviewRequired | Uncertainty to inspect with independent lineage evidence |
| AtlasRole / StandaloneCompositionRole | Library-level roles; >=100 is an operational support screen, not a power guarantee |

No malignant, subtype program, or integrated embedding labels are inherited. Full native gene counts are retained separately from representation gene masks.

## Statistical vocabulary

- `Confirmed`: passed a named release gate and directly supports the stated claim.
- `Supported`: source-backed, but subset-specific, sensitivity-level or not independent.
- `Exploratory`: hypothesis-generating without clinical/causal validation.
- `Qualified`: technically auditable evidence with a stated scientific limitation that blocks a stronger claim.
- `Negative`: the prespecified analysis was estimable but did not support the tested claim.
- `Not estimable`: the required denominator or contrast is absent; no effect is imputed.
- `Superseded`: retained for provenance but replaced by a canonical analysis.
- `Blocked`: required work cannot proceed under the current evidence or environment.

## Color governance

Subtype colors come from `docs/subtype_color_dictionary.csv`. Study colors come from
`0-Data/qc/study_color_dictionary.csv`. Do not redefine these mappings in notebooks
or plotting scripts.
