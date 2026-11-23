#!/usr/bin/env python3
"""Generate native-expression marker, reference and per-patient annotation evidence."""
from pathlib import Path
import json
import warnings

import anndata as ad
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

warnings.filterwarnings('ignore', category=FutureWarning)
ROOT = Path('/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/Ongoing/202504ZL_iCCA')
RUN = ROOT / '7-SingleCell/runs/20261102_two_cohort_repair_v1'
OUT = RUN / 'r3_annotation'
FINAL = OUT / 'r3a_complete.json'
if FINAL.exists():
    raise SystemExit('R3a is already complete; refusing to overwrite')
OUT.mkdir(exist_ok=True)

PANELS = {
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'KRT7', 'CDH1', 'MUC1'],
    'Hepatocyte': ['ALB', 'APOA1', 'APOA2', 'APOB', 'TTR', 'HNF4A', 'SERPINA1', 'CYP3A4'],
    'T_NK': ['PTPRC', 'CD3D', 'CD3E', 'CD3G', 'TRAC', 'CD8A', 'CD4', 'NKG7', 'GNLY', 'KLRD1'],
    'B_Plasma': ['CD79A', 'CD79B', 'MS4A1', 'BANK1', 'MZB1', 'JCHAIN', 'DERL3'],
    'Myeloid': ['LST1', 'LYZ', 'FCER1G', 'TYROBP', 'CD14', 'FCGR3A', 'CD68', 'AIF1', 'C1QA', 'C1QB'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'COL6A1', 'PDGFRA', 'ACTA2'],
    'Endothelial': ['VWF', 'PECAM1', 'ENG', 'CDH5', 'KDR', 'FLT1', 'CLDN5', 'RAMP2'],
}
MODELS = {
    'celltypist_liver': '/Users/luye/.celltypist/data/models/Healthy_Human_Liver.pkl',
    'celltypist_immune': '/Users/luye/.celltypist/data/models/Immune_All_High.pkl',
}


def main():
    a = ad.read_h5ad(RUN / 'r2_integration/integrated_representation.h5ad')
    availability = a.varm['native_availability']
    common = a.var.available_all_libraries.to_numpy()
    counts = a.layers['counts'].tocsr()
    obs = a.obs.copy()
    marker_rows = []
    for panel, genes in PANELS.items():
        for gene in genes:
            marker_rows.append({'panel': panel, 'gene': gene, 'in_union': gene in a.var_names,
                                'available_library_count': int(availability.loc[gene].sum()) if gene in a.var_names else 0,
                                'available_all_libraries': bool(availability.loc[gene].all()) if gene in a.var_names else False})
        use = [g for g in genes if g in a.var_names and availability.loc[g].all()]
        if not use:
            raise RuntimeError(f'No universally available markers for {panel}')
        X = a[:, use].X
        C = counts[:, a.var_names.get_indexer(use)]
        obs[f'marker_mean_{panel}'] = np.asarray(X.mean(axis=1)).ravel()
        obs[f'marker_detection_{panel}'] = np.asarray((C > 0).mean(axis=1)).ravel()
        print(panel, 'universally available markers:', use, flush=True)
    pd.DataFrame(marker_rows).to_csv(OUT / 'marker_native_availability.csv', index=False)

    # Reference predictions are evidence columns, never automatic final labels.
    model_audit = []
    for cohort in ['GSE151530', 'GSE138709']:
        mask = (a.obs.cohort == cohort).to_numpy()
        libs = a.obs.loc[mask, 'sample'].unique()
        cohort_genes = availability[libs].all(axis=1).to_numpy()
        view = a[mask, cohort_genes].copy()
        for key, model_path in MODELS.items():
            pred = celltypist.annotate(view, model=model_path, majority_voting=False)
            labels = pred.predicted_labels
            obs.loc[mask, key] = labels.predicted_labels.to_numpy()
            obs.loc[mask, key + '_confidence'] = pred.probability_matrix.max(axis=1).to_numpy()
            model = celltypist.models.Model.load(model_path)
            missing = sorted(set(model.features) - set(view.var_names))
            model_audit.append({'cohort': cohort, 'model': key, 'model_features': len(model.features),
                                'available_features': len(model.features) - len(missing),
                                'missing_features': ';'.join(missing)})
            print(cohort, key, 'missing model features', len(missing), flush=True)
            obs.loc[mask, [key, key + '_confidence']].to_csv(
                OUT / f'{cohort}_{key}_predictions.csv.gz')
    pd.DataFrame(model_audit).to_csv(OUT / 'celltypist_feature_audit.csv', index=False)
    obs.to_csv(OUT / 'annotation_evidence_obs.csv.gz')

    cluster = 'leiden_harmony'
    obs['cluster'] = a.obs[cluster].astype(str).to_numpy()
    score_cols = [f'marker_mean_{p}' for p in PANELS] + [f'marker_detection_{p}' for p in PANELS]
    obs.groupby(['cluster', 'cohort'], observed=True)[score_cols].mean().to_csv(
        OUT / 'cluster_cohort_marker_evidence.csv')
    obs.groupby(['cluster', 'cohort', 'patient'], observed=True)[score_cols].mean().to_csv(
        OUT / 'cluster_patient_marker_evidence.csv')
    obs.groupby(['cluster', 'cohort', 'patient'], observed=True).size().rename('n_cells').to_csv(
        OUT / 'cluster_patient_counts.csv')
    obs.groupby(['cluster', 'cohort', 'celltypist_liver'], observed=True).size().rename('n_cells').to_csv(
        OUT / 'cluster_celltypist_liver_counts.csv')
    obs.groupby(['cluster', 'cohort', 'celltypist_immune'], observed=True).size().rename('n_cells').to_csv(
        OUT / 'cluster_celltypist_immune_counts.csv')
    obs.groupby(['cluster', 'cohort', 'author_celltype'], observed=True).size().rename('n_cells').to_csv(
        OUT / 'cluster_author_label_counts.csv')

    # Native-expression cluster markers within each cohort.
    marker_tables = []
    for cohort in ['GSE151530', 'GSE138709']:
        mask = (a.obs.cohort == cohort).to_numpy()
        view = a[mask, common].copy()
        view.obs['cluster_review'] = pd.Categorical(view.obs[cluster].astype(str))
        sc.tl.rank_genes_groups(view, 'cluster_review', method='wilcoxon', n_genes=50, use_raw=False)
        result = view.uns['rank_genes_groups']
        for group in result['names'].dtype.names:
            frame = pd.DataFrame({key: result[key][group] for key in ['names', 'scores', 'logfoldchanges', 'pvals_adj']})
            frame['cohort'] = cohort
            frame['cluster'] = group
            marker_tables.append(frame)
    pd.concat(marker_tables).to_csv(OUT / 'cluster_native_markers.csv', index=False)
    summary = {
        'cells': len(obs), 'clusters': int(obs.cluster.nunique()),
        'panels': list(PANELS), 'reference_models': list(MODELS),
        'annotation_status': 'EVIDENCE_GENERATED_MANUAL_ADJUDICATION_PENDING',
        'bulk_programs_used_for_annotation': False,
    }
    FINAL.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
