# Preparation snapshot: scope and verification

Updated 2026-11-23. Tag: `v2026.11.23-preparation`.

## Retained analysis and visualization code

All 177 files under the previous public `_implementation/` tree are retained byte-for-byte, including the 168 selected project files indexed in `IMPLEMENTATION_MANIFEST.csv`. The 24 analysis entries, 37 selected main/supplementary drawing jobs, three final drivers, seven configuration files, `run.py`, `check_results.py`, reproduction notebook and input-isolation test are present. A regression check verifies the selected source hashes and every analysis entry. The new standalone six-script workflow is additional code.

The preserved full selected-panel route was verified at baseline e4af73a: 37 production renderers, three final drivers, 100 quantitative assets and 756 unchanged inputs. Its [dated execution record](history/20261112/VERIFICATION_20261112.json) remains scoped to that run. We did not repeat the full selected-panel route or upstream fitting for this freeze.

## Checks executed for this snapshot

- Seven Python regression tests pass, covering analysis-source retention, analysis entry presence, input isolation and the portable plotting notebook entry.
- The six standalone R scripts execute in an independent temporary copy using R 4.5.3. They produce 92 output files, including six one-page composite PDF/PNG pairs. All 75 staged input/dictionary hashes are unchanged after drawing.
- Notebook syntax, runtime selection and command/path checks pass. A pinned R 4.5 runtime is selected correctly despite the macOS framework wrapper redirect. Direct plotting dependencies include ggrastr, svglite and tidyr.
- PDF text checks confirm the corrected S6J cross-reference and retained direction arrows in individual panels. Cairo PDF export avoids the observed Unicode substitutions. The previous overview and its prompt are unchanged.
- The governed study palette passes the project validator. File identities and scoped machine-readable results are in `INPUT_FILES.csv` and `VERIFICATION_20261123.json`.

The run reports 24 legend-theme warnings (`title.position`/`title.hjust`); these concern styling and remain recorded. A fresh environment installation and cross-platform execution were not tested.

## Draft layout and scientific boundaries

Visual inspection of the composite PDF contact sheet identifies layout work still required, including overlapping colour-bar/text labels in the single-cell program projection panel and crowded multi-part labels in the discovery/multi-omics drafts. The composite PDFs embed raster panels; individual PDFs are also generated. These drafts are not approved journal-ready artwork.

The local selected manuscript figures and the new standalone drafts still require final version binding, caption alignment and author acceptance. The snapshot does not assert that the current manuscript is complete. Q2 malignant identity, Q3 estimability, external projection interpretation, fixed 17-gene exploratory scope and therapeutic-hypothesis limits remain unchanged.

The provenance review retains 19 historical code-version gaps and four original writer-event gaps. The companion frozen-input archive awaits public deposit; real software-author metadata and a persistent archival identifier remain pending. Upstream data/model dependencies are described in `upstream_external_interfaces.csv`.
