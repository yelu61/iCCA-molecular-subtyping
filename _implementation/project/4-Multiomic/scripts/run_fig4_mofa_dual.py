"""Additive two-view MOFA ablation on the locked triple-view sample/feature inputs."""
from __future__ import annotations

import contextlib
import hashlib
import importlib.metadata
import itertools
import json
import os
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from mofapy2.run.entry_point import entry_point
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '4-Multiomic/refinements/fig4_mofa_dual_20261030'
SEED = 20260715
VIEWS = ['RNA', 'Protein']


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def bootstrap(a, b, seed):
    rng = np.random.default_rng(seed)
    result = []
    for _ in range(2000):
        idx = rng.integers(0, len(a), len(a))
        result.append(adjusted_rand_score(a[idx], b[idx]))
    return np.quantile(result, [.025, .975])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    paths = [f'4-Multiomic/mofa_inputs/step4_mofa_input_{v}.csv' for v in VIEWS] + [
        '4-Multiomic/mofa_inputs/step4_mofa_sample_metadata.csv',
        '4-Multiomic/publication/models/step4_mofa_model.hdf5',
        '4-Multiomic/publication/models/step4_mofa_run_manifest.json',
        '4-Multiomic/publication/tables/step4_multiomic_cluster_assignments.csv',
        '4-Multiomic/publication/tables/step4_multiomic_cluster_comparison.csv',
        'governance/data/subtype_master_table.csv',
        '4-Multiomic/scripts/run_fig4_mofa_dual.py']
    inputs = pd.DataFrame({'Path': paths, 'MD5': [md5(ROOT / p) for p in paths]})
    if (OUT / 'input_hashes.csv').exists() and inputs.equals(pd.read_csv(OUT / 'input_hashes.csv')):
        if (OUT / 'run_manifest.json').exists() and json.loads((OUT / 'run_manifest.json').read_text()).get('converged'):
            print('Reusing verified additive dual-view MOFA outputs.', flush=True)
            return
    meta = pd.read_csv(ROOT / paths[2], dtype={'Sample': str})
    assert len(meta) == 204 and meta.Sample.is_unique
    triple = pd.read_csv(ROOT / paths[5], dtype={'Sample': str})
    assert list(triple.Sample) == list(meta.Sample)
    assert triple.FinalSubtype.equals(meta.FinalSubtype)
    assert triple.SubtypeClass.equals(meta.SubtypeClass)
    master = pd.read_csv(ROOT / paths[7], dtype={'Sample': str}).set_index('Sample')
    assert meta.FinalSubtype.tolist() == master.loc[meta.Sample, 'FinalSubtype'].tolist()
    assert meta.SubtypeClass.tolist() == master.loc[meta.Sample, 'SubtypeClass'].tolist()
    matrices = {}
    for v, p in zip(VIEWS, paths[:2]):
        frame = pd.read_csv(ROOT / p, dtype={'Sample': str}).set_index('Sample')
        assert list(frame.index) == list(meta.Sample) and frame.shape == (204, 1000)
        assert frame.columns.is_unique and np.isfinite(frame.to_numpy()).all()
        matrices[v] = frame
    print('Fitting RNA+Protein MOFA on unchanged 204 cases, 1000 features/view; seed 20260715.', flush=True)
    model = entry_point()
    with open(OUT / 'training.log', 'w') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        model.set_data_options(scale_views=False, scale_groups=False, center_groups=True, use_float32=False)
        model.set_data_matrix(data=[[matrices[v].to_numpy(dtype=float)] for v in VIEWS],
                              likelihoods=['gaussian'] * 2, views_names=VIEWS, groups_names=['2022Dong'],
                              samples_names=[list(meta.Sample)], features_names=[list(matrices[v].columns) for v in VIEWS])
        model.set_model_options(factors=10, spikeslab_factors=False, spikeslab_weights=True,
                                ard_factors=True, ard_weights=True)
        model.set_train_options(iter=1000, convergence_mode='medium', startELBO=1, freqELBO=1,
                                seed=SEED, verbose=False, quiet=False, dropR2=None)
        model.build()
        model.run()
        model.save(str(OUT / 'mofa_dual_model.hdf5'), save_data=False, save_parameters=False)
    converged = 'Converged!' in (OUT / 'training.log').read_text()
    if not converged:
        raise RuntimeError('Dual-view MOFA did not reach its predeclared convergence criterion; see training.log')
    ex = model.model.getExpectations(True)
    factors = np.asarray(ex['Z'], dtype=float)
    factor_names = [f'Factor{i+1}' for i in range(factors.shape[1])]
    scores = pd.DataFrame(factors, columns=factor_names)
    scores.insert(0, 'Sample', list(meta.Sample))
    scores = scores.merge(meta, on='Sample', validate='one_to_one')
    scores.to_csv(OUT / 'factor_scores.csv', index=False)
    r2 = np.asarray(model.model.calculate_variance_explained(total=False)[0], dtype=float) * 100
    total = np.asarray(model.model.calculate_variance_explained(total=True)[0], dtype=float) * 100
    variance = pd.DataFrame([dict(View=v, Factor=f, VarianceExplainedPercent=r2[i, k])
                            for i, v in enumerate(VIEWS) for k, f in enumerate(factor_names)])
    variance.to_csv(OUT / 'variance_explained_per_factor.csv', index=False)
    pd.DataFrame({'View': VIEWS, 'TotalVarianceExplainedPercent': total}).to_csv(OUT / 'total_variance_explained.csv', index=False)
    weights = pd.DataFrame([dict(View=v, Feature=feature, Factor=f, Weight=float(ex['W'][i][j, k]))
                            for i, v in enumerate(VIEWS) for j, feature in enumerate(matrices[v].columns)
                            for k, f in enumerate(factor_names)])
    weights.to_csv(OUT / 'feature_weights.csv', index=False)
    factor_use = [f for f in factor_names if variance.loc[variance.Factor.eq(f), 'VarianceExplainedPercent'].sum() >= .1]
    if len(factor_use) < 2:
        factor_use = factor_names[:2]
    scaled = StandardScaler().fit_transform(scores[factor_use])
    labels = KMeans(n_clusters=3, n_init=100, random_state=SEED).fit_predict(scaled) + 1
    joined = meta.copy()
    joined['RPCluster'] = labels
    joined['RPCCluster'] = triple.MOFAKMeansCluster.to_numpy()
    joined.to_csv(OUT / 'full_fit_assignments.csv', index=False)
    matched = joined.loc[joined.SubtypeClass.eq('Main')].copy()
    assert len(matched) == 182 and matched.FinalSubtype.value_counts().to_dict() == {'C2': 92, 'C3': 48, 'C1': 42}
    matched.to_csv(OUT / 'main_joined_assignments.csv', index=False)
    paths_rows = []
    for rp, subtype, rpc in itertools.product(range(1, 4), ['C1', 'C2', 'C3'], range(1, 4)):
        samples = matched.loc[matched.RPCluster.eq(rp) & matched.FinalSubtype.eq(subtype) & matched.RPCCluster.eq(rpc), 'Sample']
        paths_rows.append(dict(PathID=f'RP{rp}|{subtype}|RPC{rpc}', RPCluster=rp, FinalSubtype=subtype,
                               RPCCluster=rpc, N=len(samples), Samples=';'.join(sorted(samples))))
    paths_table = pd.DataFrame(paths_rows)
    assert paths_table.N.sum() == len(matched)
    paths_table.to_csv(OUT / 'alluvial_paths.csv', index=False)
    comparisons = []
    for k, (name, a, b) in enumerate([
        ('Primary RNA vs RNA+Protein MOFA factor k-means', matched.FinalSubtype.to_numpy(), matched.RPCluster.to_numpy()),
        ('Primary RNA vs RNA+Protein+CNV MOFA factor k-means', matched.FinalSubtype.to_numpy(), matched.RPCCluster.to_numpy()),
        ('RNA+Protein vs RNA+Protein+CNV MOFA factor k-means', matched.RPCluster.to_numpy(), matched.RPCCluster.to_numpy())]):
        lo, hi = bootstrap(a, b, SEED + 303 + k)
        comparisons.append(dict(Comparison=name, N=len(matched), ARI=adjusted_rand_score(a, b),
                                BootstrapLower95=lo, BootstrapUpper95=hi, BootstrapReplicates=2000))
    pd.DataFrame(comparisons).to_csv(OUT / 'agreement.csv', index=False)
    stats = model.model.getTrainingStats()
    stats = {k: np.asarray(v).reshape(-1) for k, v in stats.items() if np.asarray(v).ndim == 1}
    maxlen = max(map(len, stats.values()))
    pd.DataFrame({k: np.pad(v, (0, maxlen-len(v)), constant_values=np.nan) for k, v in stats.items()}).to_csv(OUT / 'training_statistics.csv', index=False)
    # The existing triple-view model, its assignments and all protected inputs must stay byte-identical.
    assert inputs.MD5.tolist() == [md5(ROOT / p) for p in paths]
    manifest = dict(seed=SEED, views=VIEWS, fit_samples=204, main_samples=182,
                    features_per_view={v: 1000 for v in VIEWS}, initial_factors=10,
                    retained_factors=factors.shape[1], factor_use_for_kmeans=factor_use,
                    kmeans_clusters=3, kmeans_n_init=100, iterations_requested=1000,
                    iterations_completed=maxlen, convergence_mode='medium', converged=converged,
                    python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                    scipy=scipy.__version__, mofapy2=importlib.metadata.version('mofapy2'),
                    thread_limits={k: os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},
                    note='New two-view ablation on exact locked triple-view cases/features. Frozen triple model reused, not refitted. One deterministic seed; fixed-k correspondence, not optimal-k or multi-seed validation.')
    (OUT / 'run_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    inputs.to_csv(OUT / 'input_hashes.csv', index=False)
    artifacts = [p for p in OUT.iterdir() if p.is_file() and p.name != 'output_hashes.csv']
    pd.DataFrame({'File': [p.name for p in artifacts], 'MD5': [md5(p) for p in artifacts]}).to_csv(OUT / 'output_hashes.csv', index=False)
    print(pd.DataFrame(comparisons).to_string(index=False), flush=True)
    print(f'Converged in {maxlen} iterations. All frozen inputs and the triple model are unchanged.', flush=True)


if __name__ == '__main__':
    main()
