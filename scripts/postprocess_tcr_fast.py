#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def find_first(run_dir: Path, pattern: str) -> Optional[Path]:
    hits = sorted(run_dir.glob(pattern))
    return hits[0] if hits else None


def extract_best_pair(contig_df: pd.DataFrame) -> pd.DataFrame:
    df = contig_df.copy()
    for col in ["productive", "high_confidence", "full_length"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.lower().eq("true")

    if "productive" not in df.columns or "high_confidence" not in df.columns:
        return pd.DataFrame(
            columns=[
                "barcode",
                "cdr3_TRA",
                "cdr3_TRB",
                "v_TRA",
                "j_TRA",
                "v_TRB",
                "d_TRB",
                "j_TRB",
                "umi_TRA",
                "umi_TRB",
                "reads_TRA",
                "reads_TRB",
            ]
        )

    df = df[df["productive"] & df["high_confidence"]]
    df = df[df["chain"].isin(["TRA", "TRB"])]

    if df.empty:
        return pd.DataFrame(
            columns=[
                "barcode",
                "cdr3_TRA",
                "cdr3_TRB",
                "v_TRA",
                "j_TRA",
                "v_TRB",
                "d_TRB",
                "j_TRB",
                "umi_TRA",
                "umi_TRB",
                "reads_TRA",
                "reads_TRB",
            ]
        )

    # Keep top chain per barcode by UMI support for speed and determinism.
    df["umis"] = pd.to_numeric(df.get("umis", 0), errors="coerce").fillna(0)
    df["reads"] = pd.to_numeric(df.get("reads", 0), errors="coerce").fillna(0)

    tra = (
        df[df["chain"] == "TRA"]
        .sort_values(["barcode", "umis", "reads"], ascending=[True, False, False])
        .drop_duplicates("barcode")
        .rename(
            columns={
                "cdr3": "cdr3_TRA",
                "v_gene": "v_TRA",
                "j_gene": "j_TRA",
                "umis": "umi_TRA",
                "reads": "reads_TRA",
            }
        )
    )
    trb = (
        df[df["chain"] == "TRB"]
        .sort_values(["barcode", "umis", "reads"], ascending=[True, False, False])
        .drop_duplicates("barcode")
        .rename(
            columns={
                "cdr3": "cdr3_TRB",
                "v_gene": "v_TRB",
                "d_gene": "d_TRB",
                "j_gene": "j_TRB",
                "umis": "umi_TRB",
                "reads": "reads_TRB",
            }
        )
    )

    keep_tra = ["barcode", "cdr3_TRA", "v_TRA", "j_TRA", "umi_TRA", "reads_TRA"]
    keep_trb = ["barcode", "cdr3_TRB", "v_TRB", "d_TRB", "j_TRB", "umi_TRB", "reads_TRB"]

    paired = pd.merge(tra[keep_tra], trb[keep_trb], on="barcode", how="inner")
    return paired


def fast_scanpy_qc_norm_harmony(adatas: List[ad.AnnData]) -> ad.AnnData:
    adata = ad.concat(adatas, axis=0, join="outer", merge="same")

    mito = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["mt"] = mito
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Speed-focused QC thresholds.
    adata = adata[(adata.obs["n_genes_by_counts"] >= 200) & (adata.obs["pct_counts_mt"] <= 20)].copy()
    sc.pp.filter_genes(adata, min_cells=3)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000, batch_key="sample")
    adata = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.pca(adata, n_comps=50, svd_solver="randomized")

    try:
        import scanpy.external as sce

        sce.pp.harmony_integrate(adata, key="sample", basis="X_pca")
    except Exception:
        # Keep pipeline robust if harmonypy is unavailable.
        adata.obsm["X_pca_harmony"] = adata.obsm["X_pca"].copy()

    return adata


def main() -> None:
    p = argparse.ArgumentParser(description="Postprocess cellranger outputs for fast TCR-centric analysis.")
    p.add_argument("--samples-tsv", required=True)
    p.add_argument("--cellranger-root", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    samples = pd.read_csv(args.samples_tsv, sep="\t")
    required = {"sample_id"}
    missing = required - set(samples.columns)
    if missing:
        raise ValueError(f"Missing required columns in samples TSV: {sorted(missing)}")

    root = Path(args.cellranger_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adatas: List[ad.AnnData] = []
    paired_rows: List[pd.DataFrame] = []
    summary: Dict[str, Dict[str, int]] = {}

    for _, row in samples.iterrows():
        sample_id = str(row["sample_id"])
        run_dir = root / sample_id

        gex_h5 = find_first(run_dir, "outs/**/sample_filtered_feature_bc_matrix.h5")
        contigs = find_first(run_dir, "outs/**/filtered_contig_annotations.csv")

        summary[sample_id] = {"has_gex": int(gex_h5 is not None), "has_contigs": int(contigs is not None)}

        if gex_h5 is not None:
            a = sc.read_10x_h5(gex_h5.as_posix(), gex_only=True)
            a.var_names_make_unique()
            a.obs["sample"] = sample_id
            a.obs_names = [f"{sample_id}:{bc}" for bc in a.obs_names]
            adatas.append(a)

        if contigs is not None:
            c = pd.read_csv(contigs)
            paired = extract_best_pair(c)
            if not paired.empty:
                paired["sample"] = sample_id
                paired["cell_id"] = paired["sample"] + ":" + paired["barcode"]
                paired["pair_id"] = paired["cdr3_TRA"].astype(str) + "|" + paired["cdr3_TRB"].astype(str)
                paired_rows.append(paired)
                summary[sample_id]["paired_cells"] = int(paired.shape[0])
            else:
                summary[sample_id]["paired_cells"] = 0
        else:
            summary[sample_id]["paired_cells"] = 0

    if not paired_rows:
        raise RuntimeError("No paired TRA/TRB cells found in filtered_contig_annotations.csv files.")

    paired_df = pd.concat(paired_rows, ignore_index=True)
    paired_df.to_csv(outdir / "paired_alpha_beta_cells.csv", index=False)

    per_sample = (
        paired_df.groupby(["sample", "pair_id"], as_index=False)
        .size()
        .rename(columns={"size": "n_cells"})
    )
    per_sample["relative_abundance"] = per_sample["n_cells"] / per_sample.groupby("sample")["n_cells"].transform("sum")

    global_ab = (
        paired_df.groupby(["pair_id"], as_index=False).size().rename(columns={"size": "n_cells"})
    )
    global_ab["relative_abundance"] = global_ab["n_cells"] / global_ab["n_cells"].sum()
    global_ab["sample"] = "ALL"

    abundance = pd.concat([per_sample, global_ab], ignore_index=True)
    abundance.to_csv(outdir / "paired_alpha_beta_clonotype_abundance.csv", index=False)

    if adatas:
        adata = fast_scanpy_qc_norm_harmony(adatas)
        tcr_cell_ids = set(paired_df["cell_id"].tolist())
        adata.obs["has_paired_tcr_ab"] = adata.obs_names.isin(tcr_cell_ids)
        adata.write_h5ad(outdir / "adata_qc_norm_harmony.h5ad", compression="gzip")
        summary["_global"] = {
            "n_cells_after_qc": int(adata.n_obs),
            "n_genes_after_hvg": int(adata.n_vars),
            "n_paired_tcr_cells": int(adata.obs["has_paired_tcr_ab"].sum()),
        }

    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
