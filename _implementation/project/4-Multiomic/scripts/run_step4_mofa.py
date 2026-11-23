#!/usr/bin/env python3
"""Train a deterministic MOFA model for the Step4 matched 2022Dong cohort."""

from __future__ import annotations

import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import kruskal
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler

from mofapy2.run.entry_point import entry_point


ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "4-Multiomic" / "mofa_inputs"
TABLE_DIR = ROOT / "4-Multiomic" / "publication" / "tables"
MODEL_DIR = ROOT / "4-Multiomic" / "publication" / "models"
SEED = 20260715
VIEWS = ["RNA", "Protein", "CNV"]


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan)
    ok = np.isfinite(values)
    p = values[ok]
    if p.size == 0:
        return out
    order = np.argsort(p)
    ranked = p[order] * p.size / np.arange(1, p.size + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    out[ok] = adjusted
    return out


def bootstrap_ari(labels_a: np.ndarray, labels_b: np.ndarray, seed: int, replicates: int = 2000) -> tuple[float, float]:
    labels_a = np.asarray(labels_a)
    labels_b = np.asarray(labels_b)
    rng = np.random.default_rng(seed)
    estimates = np.empty(replicates)
    for index in range(replicates):
        sampled = rng.integers(0, labels_a.size, labels_a.size)
        estimates[index] = adjusted_rand_score(labels_a[sampled], labels_b[sampled])
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def load_inputs() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    meta = pd.read_csv(INPUT_DIR / "step4_mofa_sample_metadata.csv", dtype={"Sample": str})
    if meta["Sample"].duplicated().any():
        raise ValueError("Duplicated MOFA sample identifiers")
    matrices: dict[str, pd.DataFrame] = {}
    for view in VIEWS:
        frame = pd.read_csv(INPUT_DIR / f"step4_mofa_input_{view}.csv", dtype={"Sample": str})
        if frame["Sample"].duplicated().any():
            raise ValueError(f"Duplicated samples in {view}")
        frame = frame.set_index("Sample")
        if list(frame.index) != list(meta["Sample"]):
            raise ValueError(f"Sample order mismatch in {view}")
        if frame.columns.duplicated().any():
            raise ValueError(f"Duplicated features in {view}")
        values = frame.to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"Non-finite values in {view}")
        matrices[view] = frame
    all_features = [feature for view in VIEWS for feature in matrices[view].columns]
    if len(all_features) != len(set(all_features)):
        raise ValueError("MOFA feature identifiers must be globally unique across views")
    return meta, matrices


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    meta, matrices = load_inputs()

    data = [[matrices[view].to_numpy(dtype=float)] for view in VIEWS]
    feature_names = [list(matrices[view].columns) for view in VIEWS]
    samples = list(meta["Sample"])

    model = entry_point()
    model.set_data_options(scale_views=False, scale_groups=False, center_groups=True, use_float32=False)
    model.set_data_matrix(
        data=data,
        likelihoods=["gaussian"] * len(VIEWS),
        views_names=VIEWS,
        groups_names=["2022Dong"],
        samples_names=[samples],
        features_names=feature_names,
    )
    model.set_model_options(
        factors=10,
        spikeslab_factors=False,
        spikeslab_weights=True,
        ard_factors=True,
        ard_weights=True,
    )
    model.set_train_options(
        iter=1000,
        convergence_mode="medium",
        startELBO=1,
        freqELBO=1,
        seed=SEED,
        verbose=False,
        quiet=False,
        dropR2=None,
    )
    model.build()
    model.run()

    model_file = MODEL_DIR / "step4_mofa_model.hdf5"
    model.save(str(model_file), save_data=False, save_parameters=False)

    expectations = model.model.getExpectations(True)
    factors = np.asarray(expectations["Z"], dtype=float)
    weights = expectations["W"]
    factor_names = [f"Factor{i + 1}" for i in range(factors.shape[1])]

    factor_frame = pd.DataFrame(factors, columns=factor_names)
    factor_frame.insert(0, "Sample", samples)
    factor_frame = factor_frame.merge(meta, on="Sample", validate="one_to_one")
    factor_frame.to_csv(TABLE_DIR / "step4_mofa_factor_scores.csv", index=False)

    weight_rows = []
    for view_index, view in enumerate(VIEWS):
        view_weights = np.asarray(weights[view_index], dtype=float)
        for feature_index, feature in enumerate(feature_names[view_index]):
            original = feature.split("::", 1)[1] if "::" in feature else feature
            for factor_index, factor in enumerate(factor_names):
                weight_rows.append(
                    {
                        "View": view,
                        "Feature": feature,
                        "OriginalFeature": original,
                        "Factor": factor,
                        "Weight": view_weights[feature_index, factor_index],
                    }
                )
    weight_frame = pd.DataFrame(weight_rows)
    weight_frame["AbsoluteWeight"] = weight_frame["Weight"].abs()
    weight_frame["WithinViewFactorRank"] = (
        weight_frame.groupby(["View", "Factor"])["AbsoluteWeight"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    weight_frame.to_csv(TABLE_DIR / "step4_mofa_feature_weights.csv", index=False)
    weight_frame.loc[weight_frame["WithinViewFactorRank"] <= 20].to_csv(
        TABLE_DIR / "step4_mofa_top_feature_weights.csv", index=False
    )

    r2 = np.asarray(model.model.calculate_variance_explained(total=False)[0], dtype=float) * 100
    total_r2 = np.asarray(model.model.calculate_variance_explained(total=True)[0], dtype=float) * 100
    variance_rows = []
    for view_index, view in enumerate(VIEWS):
        for factor_index, factor in enumerate(factor_names):
            variance_rows.append(
                {
                    "View": view,
                    "Factor": factor,
                    "VarianceExplainedPercent": r2[view_index, factor_index],
                }
            )
    variance_frame = pd.DataFrame(variance_rows)
    variance_frame.to_csv(TABLE_DIR / "step4_mofa_variance_explained_per_factor.csv", index=False)
    pd.DataFrame(
        {"View": VIEWS, "TotalVarianceExplainedPercent": total_r2}
    ).to_csv(TABLE_DIR / "step4_mofa_total_variance_explained.csv", index=False)

    association_rows = []
    main_mask = factor_frame["SubtypeClass"].eq("Main")
    for factor in factor_names:
        groups = [
            factor_frame.loc[main_mask & factor_frame["FinalSubtype"].eq(subtype), factor].to_numpy()
            for subtype in ["C1", "C2", "C3"]
        ]
        statistic, pvalue = kruskal(*groups)
        association_rows.append(
            {
                "Factor": factor,
                "MainStateN": int(main_mask.sum()),
                "KruskalWallisStatistic": statistic,
                "PValue": pvalue,
            }
        )
    associations = pd.DataFrame(association_rows)
    associations["FDR"] = bh_adjust(associations["PValue"].to_numpy())
    associations = associations.sort_values(["FDR", "PValue", "Factor"])
    associations.to_csv(TABLE_DIR / "step4_mofa_factor_subtype_associations.csv", index=False)

    factor_means = (
        factor_frame.groupby(["FinalSubtype", "SubtypeClass"], observed=True)[factor_names]
        .agg(["mean", "median", "count"])
    )
    factor_means.columns = [f"{factor}_{stat}" for factor, stat in factor_means.columns]
    factor_means.reset_index().to_csv(TABLE_DIR / "step4_mofa_factor_summary_by_state.csv", index=False)

    factor_use = [
        factor
        for factor in factor_names
        if variance_frame.loc[variance_frame["Factor"].eq(factor), "VarianceExplainedPercent"].sum() >= 0.1
    ]
    if len(factor_use) < 2:
        factor_use = factor_names[:2]
    factor_matrix = StandardScaler().fit_transform(factor_frame[factor_use])
    mofa_cluster = KMeans(n_clusters=3, n_init=100, random_state=SEED).fit_predict(factor_matrix) + 1

    concatenated = np.concatenate([matrices[view].to_numpy(dtype=float) for view in VIEWS], axis=1)
    pca = PCA(n_components=min(20, concatenated.shape[0] - 1), random_state=SEED)
    integrated_pcs = pca.fit_transform(concatenated)
    pca_cluster = KMeans(n_clusters=3, n_init=100, random_state=SEED).fit_predict(integrated_pcs) + 1

    cluster_frame = factor_frame[["Sample", "FinalSubtype", "SubtypeClass"]].copy()
    cluster_frame["MOFAKMeansCluster"] = mofa_cluster
    cluster_frame["IntegratedPCAKMeansCluster"] = pca_cluster
    cluster_frame.to_csv(TABLE_DIR / "step4_multiomic_cluster_assignments.csv", index=False)

    locked_labels = cluster_frame.loc[main_mask, "FinalSubtype"]
    locked_vs_mofa_ci = bootstrap_ari(
        locked_labels.to_numpy(), mofa_cluster[main_mask.to_numpy()], SEED + 101
    )
    locked_vs_pca_ci = bootstrap_ari(
        locked_labels.to_numpy(), pca_cluster[main_mask.to_numpy()], SEED + 102
    )
    mofa_vs_pca_ci = bootstrap_ari(mofa_cluster, pca_cluster, SEED + 103)
    cluster_comparison = pd.DataFrame(
        [
            {
                "Comparison": "Locked main subtype vs MOFA-factor k-means",
                "Population": "Main states",
                "N": int(main_mask.sum()),
                "ARI": adjusted_rand_score(locked_labels, mofa_cluster[main_mask.to_numpy()]),
                "BootstrapLower95": locked_vs_mofa_ci[0],
                "BootstrapUpper95": locked_vs_mofa_ci[1],
                "BootstrapReplicates": 2000,
            },
            {
                "Comparison": "Locked main subtype vs integrated-PCA k-means",
                "Population": "Main states",
                "N": int(main_mask.sum()),
                "ARI": adjusted_rand_score(locked_labels, pca_cluster[main_mask.to_numpy()]),
                "BootstrapLower95": locked_vs_pca_ci[0],
                "BootstrapUpper95": locked_vs_pca_ci[1],
                "BootstrapReplicates": 2000,
            },
            {
                "Comparison": "MOFA-factor vs integrated-PCA k-means",
                "Population": "All matched states",
                "N": len(cluster_frame),
                "ARI": adjusted_rand_score(mofa_cluster, pca_cluster),
                "BootstrapLower95": mofa_vs_pca_ci[0],
                "BootstrapUpper95": mofa_vs_pca_ci[1],
                "BootstrapReplicates": 2000,
            },
        ]
    )
    cluster_comparison.to_csv(TABLE_DIR / "step4_multiomic_cluster_comparison.csv", index=False)

    training_stats = model.model.getTrainingStats()
    one_dimensional_stats = {
        key: np.asarray(value).reshape(-1)
        for key, value in training_stats.items()
        if np.asarray(value).ndim <= 1
    }
    max_len = max(len(value) for value in one_dimensional_stats.values())
    stats_frame = pd.DataFrame({"IterationIndex": np.arange(1, max_len + 1)})
    for key, values in one_dimensional_stats.items():
        padded = np.full(max_len, np.nan)
        padded[: len(values)] = values
        stats_frame[key] = padded
    stats_frame.to_csv(TABLE_DIR / "step4_mofa_training_statistics.csv", index=False)

    run_manifest = {
        "seed": SEED,
        "samples": len(samples),
        "main_state_samples": int(main_mask.sum()),
        "views": VIEWS,
        "features_per_view": {view: int(matrices[view].shape[1]) for view in VIEWS},
        "initial_factors": 10,
        "retained_factors": int(factors.shape[1]),
        "factor_use_for_kmeans": factor_use,
        "iterations_requested": 1000,
        "convergence_mode": "medium",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "model_file": str(model_file.relative_to(ROOT)),
    }
    (MODEL_DIR / "step4_mofa_run_manifest.json").write_text(json.dumps(run_manifest, indent=2) + "\n")

    print(cluster_comparison.to_string(index=False))
    print(associations.head().to_string(index=False))


if __name__ == "__main__":
    main()
