#!/usr/bin/env python3
"""Reproducible, study-aware TIDE refresh for Step2.

TIDE is used as a hypothesis-generating immune-evasion score, not as a
validated clinical response endpoint in iCCA. Expression is centered within
study before TIDEpy is called, and downstream subtype tests adjust for Study.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "2-TranscriptomicCharacterization"
DATA = MODULE / "data"
TABLES = MODULE / "publication" / "tables"
FIGURES = MODULE / "publication" / "figures"
TIDEPY = Path(os.environ.get("TIDEPY_BIN") or shutil.which("tidepy") or "tidepy")
METRICS = ["TIDE", "Dysfunction", "Exclusion", "IFNG", "CD274", "CD8", "MDSC", "CAF", "TAM M2", "CTL"]
PALETTE = {"C1": "#0072B2", "C2": "#E69F00", "C3": "#CC79A7"}


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tide(path: Path) -> pd.DataFrame:
    for sep in ("\t", ","):
        out = pd.read_csv(path, sep=sep)
        if out.shape[1] > 3:
            return out
    raise ValueError(f"Could not parse TIDE output: {path}")


def adjusted_tests(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    omnibus, pairwise = [], []
    for metric in [m for m in METRICS if m in df.columns]:
        work = df[[metric, "FinalSubtype", "Study"]].dropna().copy()
        if work[metric].nunique() < 2:
            omnibus.append({"Metric": metric, "N": len(work), "Studies": work.Study.nunique(),
                            "F": np.nan, "P": np.nan, "Note": "not tested: zero variance"})
            continue
        fit = smf.ols(f'Q("{metric}") ~ C(FinalSubtype) + C(Study)', data=work).fit(cov_type="HC3")
        reduced = smf.ols(f'Q("{metric}") ~ C(Study)', data=work).fit()
        full_plain = smf.ols(f'Q("{metric}") ~ C(FinalSubtype) + C(Study)', data=work).fit()
        fstat, fp, _ = full_plain.compare_f_test(reduced)
        omnibus.append({"Metric": metric, "N": len(work), "Studies": work.Study.nunique(),
                        "F": fstat, "P": fp, "Note": "Study-adjusted nested-model F test"})
        for a, b in (("C2", "C1"), ("C3", "C1"), ("C3", "C2")):
            key_a = f"C(FinalSubtype)[T.{a}]"
            key_b = f"C(FinalSubtype)[T.{b}]"
            names = list(fit.params.index)
            contrast = np.zeros(len(names))
            if a != "C1":
                contrast[names.index(key_a)] += 1
            if b != "C1":
                contrast[names.index(key_b)] -= 1
            test = fit.t_test(contrast)
            pairwise.append({
                "Metric": metric, "Contrast": f"{a}-{b}", "N": len(work),
                "Estimate": test.effect.item(), "SE_HC3": test.sd.item(),
                "CI_low": float(test.conf_int()[0, 0]), "CI_high": float(test.conf_int()[0, 1]),
                "P": float(test.pvalue),
            })
    om = pd.DataFrame(omnibus)
    pw = pd.DataFrame(pairwise)
    om["FDR"] = false_discovery_rate(om["P"])
    pw["FDR_within_metric"] = pw.groupby("Metric")["P"].transform(false_discovery_rate)
    pw["FDR_global"] = false_discovery_rate(pw["P"])
    return om, pw


def false_discovery_rate(values: pd.Series) -> pd.Series:
    x = pd.Series(values, dtype=float)
    valid = np.isfinite(x.to_numpy())
    answer = np.full(len(x), np.nan)
    if not valid.any():
        return pd.Series(answer, index=x.index)
    positions = np.where(valid)[0]
    order = np.argsort(x.to_numpy()[valid])
    ranked = x.to_numpy()[valid][order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q); out[order] = np.minimum(q, 1)
    answer[positions] = out
    return pd.Series(answer, index=x.index)


def responder_model(df: pd.DataFrame) -> pd.DataFrame:
    work = df[["Responder", "FinalSubtype", "Study"]].dropna().copy()
    work["Responder01"] = work["Responder"].astype(str).str.lower().isin(["true", "1", "yes"]).astype(int)
    fit = smf.logit("Responder01 ~ C(FinalSubtype) + C(Study)", data=work).fit(disp=False, cov_type="HC3")
    rows = []
    for subtype in ("C2", "C3"):
        term = f"C(FinalSubtype)[T.{subtype}]"
        beta, se, p = fit.params[term], fit.bse[term], fit.pvalues[term]
        rows.append({"Contrast": f"{subtype}-C1", "N": len(work), "OddsRatio": np.exp(beta),
                     "CI_low": np.exp(beta - 1.96 * se), "CI_high": np.exp(beta + 1.96 * se), "P": p})
    out = pd.DataFrame(rows)
    out["FDR"] = false_discovery_rate(out["P"])
    return out


def random_effects_meta(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    per_study, pooled = [], []
    for metric in ("TIDE", "Dysfunction", "Exclusion"):
        if metric not in df:
            continue
        for a, b in (("C2", "C1"), ("C3", "C1"), ("C3", "C2")):
            effects = []
            for study, d in df.groupby("Study"):
                xa = d.loc[d.FinalSubtype.eq(a), metric].dropna().to_numpy()
                xb = d.loc[d.FinalSubtype.eq(b), metric].dropna().to_numpy()
                if min(len(xa), len(xb)) < 3:
                    continue
                est = xa.mean() - xb.mean()
                se = np.sqrt(xa.var(ddof=1) / len(xa) + xb.var(ddof=1) / len(xb))
                if not np.isfinite(se) or se <= 0:
                    continue
                effects.append((study, est, se, len(xa), len(xb)))
                per_study.append({"Metric": metric, "Contrast": f"{a}-{b}", "Study": study,
                                  "Estimate": est, "SE": se, "N_A": len(xa), "N_B": len(xb)})
            if len(effects) < 2:
                continue
            yi = np.array([x[1] for x in effects]); vi = np.array([x[2] ** 2 for x in effects])
            wi = 1 / vi; fixed = np.sum(wi * yi) / np.sum(wi)
            q = np.sum(wi * (yi - fixed) ** 2); df_q = len(yi) - 1
            c = np.sum(wi) - np.sum(wi ** 2) / np.sum(wi)
            tau2 = max(0.0, (q - df_q) / c) if c > 0 else 0.0
            wr = 1 / (vi + tau2); est = np.sum(wr * yi) / np.sum(wr); se = np.sqrt(1 / np.sum(wr))
            z = est / se
            pooled.append({"Metric": metric, "Contrast": f"{a}-{b}", "Studies": len(yi),
                           "RandomEffect": est, "SE": se, "CI_low": est - 1.96 * se,
                           "CI_high": est + 1.96 * se, "P": 2 * norm.sf(abs(z)), "Tau2": tau2,
                           "I2_percent": max(0.0, (q - df_q) / q * 100) if q > 0 else 0.0})
    study = pd.DataFrame(per_study)
    meta = pd.DataFrame(pooled)
    if len(meta):
        meta["FDR"] = false_discovery_rate(meta["P"])
    return study, meta


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    expression_path = DATA / "tide_expression_main_subtypes.tsv"
    master_path = ROOT / "governance" / "data" / "subtype_master_table.csv"
    if not TIDEPY.exists():
        raise FileNotFoundError(f"TIDEpy executable not found: {TIDEPY}")
    expr = pd.read_csv(expression_path, sep="\t", index_col=0)
    master = pd.read_csv(master_path, dtype={"Sample": str})
    master = master.loc[master.SubtypeClass.eq("Main"), ["Sample", "FinalSubtype", "Study"]].copy()
    expr.columns = expr.columns.astype(str)
    common = master.Sample[master.Sample.isin(expr.columns)].tolist()
    if len(common) != len(master) or len(common) != 497:
        raise ValueError(f"Expected all 497 main samples; found {len(common)}")
    meta = master.set_index("Sample").loc[common]
    expr = expr.loc[:, common].apply(pd.to_numeric, errors="coerce")
    if expr.isna().any().any():
        raise ValueError("Expression contains non-numeric or missing values")
    centered = expr.copy()
    for _, samples in meta.groupby("Study").groups.items():
        centered.loc[:, list(samples)] = expr.loc[:, list(samples)].sub(expr.loc[:, list(samples)].mean(axis=1), axis=0)
    centered_path = DATA / "tide_expression_main_subtypes_study_centered.tsv.gz"
    centered.to_csv(centered_path, sep="\t", compression="gzip", float_format="%.6g")
    raw_path = TABLES / "step2_tidepy_raw_results.tsv"
    cmd = [str(TIDEPY), str(centered_path), "-o", str(raw_path), "-c", "Other", "--ignore_norm", "--vthres", "0"]
    run = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if run.returncode:
        raise RuntimeError(f"TIDEpy failed ({run.returncode}):\n{run.stdout}\n{run.stderr}")
    tide = read_tide(raw_path)
    first = tide.columns[0]
    if first not in ("Sample", "RowName"):
        tide = tide.rename(columns={first: "Sample"})
    elif first == "RowName":
        tide = tide.rename(columns={first: "Sample"})
    tide["Sample"] = tide["Sample"].astype(str)
    if tide.Sample.duplicated().any() or set(tide.Sample) != set(common):
        raise ValueError("TIDE output sample IDs do not exactly match the 497 main samples")
    result = meta.reset_index().merge(tide, on="Sample", how="left", validate="one_to_one")
    result.to_csv(TABLES / "step2_tide_sample_results.csv", index=False)
    summary = result.groupby("FinalSubtype").agg(N=("Sample", "size"), TIDE_median=("TIDE", "median"),
        TIDE_IQR=("TIDE", lambda x: x.quantile(.75)-x.quantile(.25)),
        Responder_rate=("Responder", lambda x: x.astype(str).str.lower().isin(["true", "1", "yes"]).mean()),
        Dysfunction_median=("Dysfunction", "median"), Exclusion_median=("Exclusion", "median")).reset_index()
    summary.to_csv(TABLES / "step2_tide_subtype_summary.csv", index=False)
    omnibus, pairwise = adjusted_tests(result)
    omnibus.to_csv(TABLES / "step2_tide_adjusted_omnibus.csv", index=False)
    pairwise.to_csv(TABLES / "step2_tide_adjusted_pairwise.csv", index=False)
    responder_model(result).to_csv(TABLES / "step2_tide_responder_adjusted.csv", index=False)
    by_study, meta_results = random_effects_meta(result)
    by_study.to_csv(TABLES / "step2_tide_study_specific_effects.csv", index=False)
    meta_results.to_csv(TABLES / "step2_tide_meta_analysis.csv", index=False)
    legacy_path = DATA / "tide_results_main_subtypes_annotated.csv"
    legacy_rows = []
    if legacy_path.exists():
        legacy = pd.read_csv(legacy_path, dtype={"Sample": str})
        joined = result.merge(legacy, on="Sample", suffixes=("_refresh", "_legacy"))
        for metric in ("TIDE", "Dysfunction", "Exclusion"):
            legacy_rows.append({"Metric": metric, "N": joined.shape[0],
                                "Pearson_r": joined[f"{metric}_refresh"].corr(joined[f"{metric}_legacy"]),
                                "Median_abs_difference": (joined[f"{metric}_refresh"]-joined[f"{metric}_legacy"]).abs().median()})
    pd.DataFrame(legacy_rows).to_csv(TABLES / "step2_tide_legacy_concordance.csv", index=False)

    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    for ax, metric in zip(axes, ("TIDE", "Dysfunction", "Exclusion")):
        sns.boxplot(data=result, x="FinalSubtype", y=metric, hue="FinalSubtype", order=["C1","C2","C3"],
                    palette=PALETTE, legend=False, showfliers=False, ax=ax)
        sns.stripplot(data=result, x="FinalSubtype", y=metric, order=["C1","C2","C3"], color="black", alpha=.2, size=2, ax=ax)
        ax.set(xlabel="", title=metric)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGURES / f"step2_tide_study_centered_scores.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    try:
        tide_version = importlib.metadata.version("TIDEpy")
    except importlib.metadata.PackageNotFoundError:
        tide_version = "unknown"
    manifest = {
        "analysis": "Step2 study-aware TIDE refresh", "tidepy": str(TIDEPY),
        "tidepy_version": tide_version, "python": platform.python_version(),
        "command": cmd, "input_expression": str(expression_path.relative_to(ROOT)),
        "input_expression_md5": md5(expression_path), "master_table_md5": md5(master_path),
        "n_samples": len(result), "n_studies": result.Study.nunique(), "n_genes": expr.shape[0],
        "normalization": "gene-wise mean centering within Study before TIDEpy; --ignore_norm",
        "cancer": "Other", "interpretation": "hypothesis-generating immune-evasion score; not clinical response validation",
    }
    (TABLES / "step2_tide_run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    gates = pd.DataFrame([
        ["main_sample_count_497", len(result) == 497, len(result)],
        ["exact_unique_sample_ids", result.Sample.nunique() == 497, result.Sample.nunique()],
        ["all_three_main_subtypes", set(result.FinalSubtype) == {"C1","C2","C3"}, ",".join(sorted(result.FinalSubtype.unique()))],
        ["study_adjusted_tests_written", len(pairwise) > 0, len(pairwise)],
        ["study_effect_meta_written", len(meta_results) > 0, len(meta_results)],
        ["tide_is_hypothesis_generating", True, "explicit in manifest and notebook"],
    ], columns=["Check", "Pass", "Detail"])
    gates.to_csv(TABLES / "step2_tide_release_gate_status.csv", index=False)
    if not gates.Pass.all():
        raise RuntimeError("One or more Step2 TIDE release checks failed")
    print(json.dumps({"status": "PASS", "samples": len(result), "studies": result.Study.nunique(), "output": str(TABLES)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
