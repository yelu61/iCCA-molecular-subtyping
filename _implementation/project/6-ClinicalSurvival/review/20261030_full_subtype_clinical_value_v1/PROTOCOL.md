# Full eligible-cohort clinical value: protocol before new performance results

2026-10-30. Required for the user's redraw of “molecular subtyping adds clinical information → relationship to existing duct pathology.” The previous331-patient three-class and79-patient binary analyses cannot be mixed as if they used the same model/population. This run fits three models fairly on the full existing eligible cohort.

## Question and population

- Use all331 principal-state patients with positive OS, valid death status and complete age/sex/stage; no requirement for measured duct. Expected153 deaths, four studies, with all C1/C2/C3 retained. Join by unique Sample and reproduce existing Step6 clinical population exactly.
- M0: age, sex, stage I/II vs III/IV; MC3: M0+C3/non-C3; MAll: M0+C2+C3 indicator variables (C1 reference). MAll includes all three states; MAll−MC3 asks whether separating C1 and C2 adds information beyond the C3 split.
- Measured/predicted duct does not enter this fit. The79-patient comparison remains a separate pathology sensitivity; it cannot be pooled with this run.

## Methods fixed before execution

- Reuse the numerical definitions and independently verified predictor of the existing79-patient Step6 extension, structurally adapted to the331-patient population and the full three-state model. Preserve that script and all existing owner outputs.
- Primary model family: Study-stratified Cox with fixed ridge theta1 and training-fold standardization. Prespecified sensitivity: same ridge model with time-varying age and included C2/C3 effects interacting with log(1+t/12), using the existing counting-process predictor. C2 and C3 interactions are both present for MAll; no post-hoc choice of the better model family.
- Same5 repeats ×5folds for all three models, stratified by Study and death status, seeds20261031–20261035. Store folds and average each patient's five OOF probabilities. No refitting of discovery labels, feature selection or expression preprocessing. These are internal predictions conditional on the existing discovered states.
- Primary endpoint:24-month cumulative/dynamic AUC. Prespecified primary contrasts MC3−M0 and MAll−M0; secondary MAll−MC3 quantifies the added C1/C2 distinction. Also retain12-month AUC, C-index from12-month risk, within-Study discrimination, Brier12/24 and point calibration.
- IPCW censoring distribution estimated by Study; same definitions as the preceding run, fixed numerical floor1e−4, with weight/support diagnostics.
-1000 paired within-Study patient bootstraps of fixed aggregatedOOF predictions, including re-estimation of censoring weights. Pointwise95% intervals;97.5% Bonferroni intervals apply only to the two primaryAUC24 contrasts. This is conditional evaluation uncertainty, not pipeline-refit, discovery or independent external-validation uncertainty.
- Native fixed- and time-varying survival prediction parity, independent AUC/Brier calculations, patient/fold alignment, source hashes and explicit numerical warning review are required. No convergence/metric failures are silently removed; missingness/estimability and relevant source-level limitations stay visible.

## Figure contract

- Question1: how much risk-discrimination information is gained from theC3 split and the complete three-state classification? Full-model performance/paired gain are shown; C1/C2 do not receive artificial separate “ΔC-index per subtype.”
- Question2: how do the molecular states overlap with duct phenotypes? Use all497 main samples and the explicit hybrid prediction field (93OOF,404projection). Show individual scores and class fractions under the existing0.5 threshold. The93 measured labels remain a distinct reference, with no inflated pathology denominator or independent-validation claim.
- Archetype: quantitative comparison with the prognostic panel as the anchor, then a distribution/association panel. Existing R/grid/Cairo/ragg workflow only. Each panel is a standalone PDF, editable Arial>=7pt, no panel letters or whole-figure assembly; review previews are PNG. Widths are chosen for standalone review, not shrunk automatically into half-width slots.
- All497 inference rows retained; no point downsampling. Source subgroup denominators, classifier provenance, fields and derived display transformations are exported.
- All work remains review-only in new directories. Do not overwrite publication tables, originalFig6 PDFs, the79-patient analysis, source data orFig1–5. No registry selection/submission lock is advanced.

## Numerical execution note

The first counting-process attempt encountered near-equal floating-point month values that became effectively zero intervals under survival's default near-tie reconciliation. Before constructing training intervals, the runner now applies survival::aeqSurv to that training fold's outcome times, matching the ordinary-Cox default tie convention. Source times and evaluation outcomes remain unchanged; no patient is removed, and no new comparative performance results were inspected to choose this repair. The failed log is retained separately.
