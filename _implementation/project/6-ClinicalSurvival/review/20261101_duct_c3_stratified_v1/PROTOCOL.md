# Duct-by-C3 stratified prognosis: prespecified exploratory extension

Recorded 2026-11-01 before fitting this extension. Status: exploratory review;
not a release of the clinical pipeline or an independent validation.

## Question and fixed population

Does C3 mark adverse overall survival among patients with observed small-duct
morphology? Use exactly the 79 complete cases from the executed measured-duct
comparison (45 deaths). The small-duct subset contains 61 patients (15 C3,
46 non-C3); the large-duct subset contains 18 (11 C3, 7 non-C3). No predicted
duct labels, population expansion, new cutoff, or subtype relabeling is allowed.
The patient is the independent unit. These are post hoc subgroups selected
because of the present manuscript question, not confirmatory strata.

## Estimands and analysis

1. Describe all four observed duct-by-C3 groups using Kaplan–Meier OS, pointwise
   log-log 95% CIs, censor marks, and number-at-risk tables. Show the complete
   observed time range and report 12/24-month support. Stratify panels by duct
   with C3/non-C3 colors. No pairwise log-rank screening or selection.
2. In the 61 small-duct patients, estimate C3 versus non-C3 OS association by
   ordinary Cox regression adjusted for age per 10 years, sex, and stage
   III/IV versus I/II, with study-specific baseline hazards. Save all
   coefficients, Wald 95% CIs and exploratory P values, and Schoenfeld
   proportional-hazards diagnostics.
3. Regardless of diagnostic P values, fit a prespecified sensitivity model
   with log1p(t/12) interactions for age and C3, retaining the same covariates
   and study strata. Report C3 HR at 12 and 24 months using the full covariance
   matrix and pointwise Wald 95% CIs. Report convergence, information-matrix
   checks and warning diagnostics. If the model is not estimable, retain that
   failure instead of selecting another model for a favorable estimate.
4. The existing 79-patient joint-model prediction comparison is reused as
   contextual evidence and is not refitted. This extension evaluates
   association and descriptive survival, not incremental prediction accuracy.

## Interpretation and graphical contract

Null and discordant results remain visible. Small-duct contrasts cannot
establish an interaction with duct; the seven large-duct/non-C3 cases do not
support a definitive absence-of-effect conclusion. No unadjusted event-rate
comparison, therapy-predictive interpretation, causal claim, multiplicity-
controlled confirmation, or clinical-utility claim is allowed. All intervals
are pointwise and exploratory, with no multiple-testing significance claim.

Quantitative grid; R-only rendering; independent editable-text PDFs at
174 mm width, without panel letters; R-generated 300-dpi PNG review previews.
Canonical C3 magenta and a neutral non-C3 color; study colors loaded centrally
but not used where study is not the visual grouping. Source data, scripts,
model objects, exact patient/risk counts, legends, hashes, and verification
accompany the PDFs. Preserve all source and historical outputs.
