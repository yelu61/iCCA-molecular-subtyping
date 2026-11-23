#!/usr/bin/env python3
"""Compute single-transform PCA/UMAP/Harmony without replacing native expression."""
from pathlib import Path
import json
import os
import warnings

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

warnings.filterwarnings('ignore', category=FutureWarning)
ROOT = Path(os.environ.get('ICCA_PROJECT_ROOT', str(Path.cwd())))
RUN = ROOT / '7-SingleCell/runs/20261102_two_cohort_repair_v1'
OUT = RUN / 'r2_integration'
FINAL = OUT / 'integrated_representation.h5ad'
if FINAL.exists():
    raise SystemExit('R2 integrated object already exists; refusing to overwrite')
OUT.mkdir(exist_ok=True)


def main():
    a = ad.read_h5ad(RUN / 'r1_source_qc/r1_native_counts.h5ad')
    counts = a.X.tocsr() if sparse.issparse(a.X) else sparse.csr_matrix(a.X)
    assert np.all(counts.data >= 0) and np.all(counts.data == np.floor(counts.data))
    a.layers['counts'] = counts.copy()
    # The only normalization/log transform in the repair workflow.
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.uns['x_semantics'] = 'normalize_total(target_sum=1e4) followed once by log1p; native availability in varm'
    common = a.var.available_all_libraries.to_numpy()
    exclude = a.var_names.str.match(r'^(MT-|RPL|RPS|HB[ABDEGMQZ][0-9]?)')
    eligible = common & ~exclude
    work = a[:, eligible].copy()
    sc.pp.highly_variable_genes(work, n_top_genes=min(3000, work.n_vars), flavor='seurat', batch_key='cohort')
    work = work[:, work.var.highly_variable].copy()
    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=min(50, work.n_vars - 1), svd_solver='arpack', random_state=123)
    a.obsm['X_pca'] = work.obsm['X_pca'].copy()

    sc.pp.neighbors(work, n_neighbors=15, n_pcs=30, use_rep='X_pca')
    a.obsp['unintegrated_connectivities'] = work.obsp['connectivities'].copy()
    a.obsp['unintegrated_distances'] = work.obsp['distances'].copy()
    sc.tl.umap(work, random_state=123)
    a.obsm['X_umap_unintegrated'] = work.obsm['X_umap'].copy()
    sc.tl.leiden(work, resolution=0.5, random_state=123, key_added='leiden_unintegrated')
    a.obs['leiden_unintegrated'] = work.obs.leiden_unintegrated.to_numpy()

    sc.external.pp.harmony_integrate(work, key='cohort', max_iter_harmony=30)
    a.obsm['X_pca_harmony'] = work.obsm['X_pca_harmony'].copy()
    sc.pp.neighbors(work, n_neighbors=15, n_pcs=30, use_rep='X_pca_harmony')
    a.obsp['harmony_connectivities'] = work.obsp['connectivities'].copy()
    a.obsp['harmony_distances'] = work.obsp['distances'].copy()
    a.obsp['connectivities'] = work.obsp['connectivities'].copy()
    a.obsp['distances'] = work.obsp['distances'].copy()
    a.uns['neighbors'] = work.uns['neighbors']
    sc.tl.umap(work, random_state=123)
    a.obsm['X_umap_harmony'] = work.obsm['X_umap'].copy()
    a.obsm['X_umap'] = a.obsm['X_umap_harmony'].copy()
    sc.tl.leiden(work, resolution=0.5, random_state=123, key_added='leiden_harmony')
    a.obs['leiden_harmony'] = work.obs.leiden_harmony.to_numpy()
    a.uns['integration_contract'] = 'Harmony cohort only; native X and counts unchanged; patient biology not an integration batch'
    a.uns['hvg_exclusion'] = 'MT/RPL/RPS/HB prefixes excluded from integration HVGs only'
    a.write(FINAL, compression='gzip')

    rows = []
    for representation in ['unintegrated', 'harmony']:
        col = 'leiden_' + representation
        tab = a.obs.groupby([col, 'cohort'], observed=True).size().rename('n').reset_index()
        tab['representation'] = representation
        rows.append(tab)
    pd.concat(rows).to_csv(OUT / 'cluster_cohort_composition.csv', index=False)
    a.obs.to_csv(OUT / 'integrated_obs.csv.gz')
    summary = {
        'cells': a.n_obs, 'genes': a.n_vars,
        'common_native_genes_all_libraries': int(common.sum()),
        'integration_eligible_genes': int(eligible.sum()), 'hvg_used': work.n_vars,
        'harmony_batch': 'cohort', 'n_pcs_neighbors': 30,
        'unintegrated_clusters': int(a.obs.leiden_unintegrated.nunique()),
        'harmony_clusters': int(a.obs.leiden_harmony.nunique()),
        'counts_layer_preserved': True, 'x_transform': 'normalize_total_then_log1p_once',
        'status': 'REPRESENTATION_COMPUTED_ANNOTATION_GATE_PENDING',
    }
    (OUT / 'r2_summary.json').write_text(json.dumps(summary, indent=2) + '\n')

    # Diagnostic overview; these are review plots, not release-gated manuscript figures.
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    for j, rep in enumerate(['unintegrated', 'harmony']):
        a.obsm['X_umap'] = a.obsm['X_umap_' + rep]
        for i, key in enumerate(['cohort', 'patient', 'leiden_' + rep]):
            sc.pl.umap(a, color=key, ax=axes[j, i], show=False, frameon=False,
                       title=f'{rep}: {key}', legend_loc='right margin' if key != 'patient' else None)
    plt.tight_layout()
    fig.savefig(OUT / 'integration_overview.png', dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
