#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


def _pick_count_col(df: pd.DataFrame) -> str:
    for c in ["duplicate_count", "consensus_count", "clone_count", "count"]:
        if c in df.columns:
            return c
    return ""


def _chain_from_row(row: pd.Series) -> str:
    for key in ["locus", "chain", "chain_type"]:
        if key in row and pd.notna(row[key]):
            v = str(row[key]).upper()
            if "TRA" in v:
                return "TRA"
            if "TRB" in v:
                return "TRB"
    return ""


def _cdr3_from_row(row: pd.Series) -> str:
    for key in ["junction_aa", "cdr3_aa", "cdr3"]:
        if key in row and pd.notna(row[key]):
            return str(row[key])
    return ""


def _best_pairs_from_airr(airr_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    df = airr_df.copy()
    if "productive" in df.columns:
        df = df[df["productive"].astype(str).str.lower().isin(["true", "t", "1"])]

    if df.empty:
        return pd.DataFrame()

    df["chain"] = df.apply(_chain_from_row, axis=1)
    df = df[df["chain"].isin(["TRA", "TRB"])]
    if df.empty:
        return pd.DataFrame()

    count_col = _pick_count_col(df)
    if not count_col:
        df["count"] = 1
        count_col = "count"
    df[count_col] = pd.to_numeric(df[count_col], errors="coerce").fillna(1)

    # Use available cell id column; fallback makes this bulk-like and cannot do true pairing.
    cell_col = ""
    for c in ["cell_id", "cell", "barcode", "sequence_id"]:
        if c in df.columns:
            cell_col = c
            break

    if not cell_col:
        return pd.DataFrame()

    tra = (
        df[df["chain"] == "TRA"]
        .sort_values([cell_col, count_col], ascending=[True, False])
        .drop_duplicates(cell_col)
    )
    trb = (
        df[df["chain"] == "TRB"]
        .sort_values([cell_col, count_col], ascending=[True, False])
        .drop_duplicates(cell_col)
    )

    keep_common = ["v_call", "d_call", "j_call"]
    keep_tra = [c for c in [cell_col, "junction_aa", "cdr3_aa", "cdr3", count_col] + keep_common if c in tra.columns]
    keep_trb = [c for c in [cell_col, "junction_aa", "cdr3_aa", "cdr3", count_col] + keep_common if c in trb.columns]

    tra = tra[keep_tra].copy()
    trb = trb[keep_trb].copy()

    rename_map_tra = {count_col: "count_TRA", "v_call": "v_TRA", "j_call": "j_TRA", "d_call": "d_TRA"}
    rename_map_trb = {count_col: "count_TRB", "v_call": "v_TRB", "j_call": "j_TRB", "d_call": "d_TRB"}
    if "junction_aa" in tra.columns:
        rename_map_tra["junction_aa"] = "cdr3_TRA"
    elif "cdr3_aa" in tra.columns:
        rename_map_tra["cdr3_aa"] = "cdr3_TRA"
    elif "cdr3" in tra.columns:
        rename_map_tra["cdr3"] = "cdr3_TRA"
    if "junction_aa" in trb.columns:
        rename_map_trb["junction_aa"] = "cdr3_TRB"
    elif "cdr3_aa" in trb.columns:
        rename_map_trb["cdr3_aa"] = "cdr3_TRB"
    elif "cdr3" in trb.columns:
        rename_map_trb["cdr3"] = "cdr3_TRB"

    tra = tra.rename(columns=rename_map_tra)
    trb = trb.rename(columns=rename_map_trb)
    paired = pd.merge(tra, trb, on=cell_col, how="inner")

    if paired.empty:
        return pd.DataFrame()

    paired = paired.rename(columns={cell_col: "cell_id"})
    paired["sample"] = sample_id
    paired["pair_id"] = paired["cdr3_TRA"].astype(str) + "|" + paired["cdr3_TRB"].astype(str)
    return paired


def _bulk_chain_abundance(airr_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    df = airr_df.copy()
    if "productive" in df.columns:
        df = df[df["productive"].astype(str).str.lower().isin(["true", "t", "1"])]
    if df.empty:
        return pd.DataFrame(columns=["sample", "chain", "cdr3", "n", "relative_abundance"])

    df["chain"] = df.apply(_chain_from_row, axis=1)
    df = df[df["chain"].isin(["TRA", "TRB"])]
    if df.empty:
        return pd.DataFrame(columns=["sample", "chain", "cdr3", "n", "relative_abundance"])

    df["cdr3"] = df.apply(_cdr3_from_row, axis=1)
    df = df[df["cdr3"] != ""]

    count_col = _pick_count_col(df)
    if not count_col:
        df["n"] = 1
    else:
        df["n"] = pd.to_numeric(df[count_col], errors="coerce").fillna(1)

    out = df.groupby(["chain", "cdr3"], as_index=False)["n"].sum()
    out["sample"] = sample_id
    out["relative_abundance"] = out["n"] / out["n"].sum()
    return out[["sample", "chain", "cdr3", "n", "relative_abundance"]]


def main() -> None:
    p = argparse.ArgumentParser(description="Combine triage TCR outputs into paired-chain abundance tables.")
    p.add_argument("--samples-tsv", required=True)
    p.add_argument("--result-root", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    samples = pd.read_csv(args.samples_tsv, sep="\t")
    if "sample_id" not in samples.columns:
        raise ValueError("samples TSV must contain 'sample_id'")

    root = Path(args.result_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    paired_rows: List[pd.DataFrame] = []
    bulk_rows: List[pd.DataFrame] = []
    summary: Dict[str, Dict[str, int]] = {}

    for _, row in samples.iterrows():
        sample = str(row["sample_id"])
        sample_dir = root / sample
        airr_files = sorted(sample_dir.glob("*.airr.tsv")) + sorted(sample_dir.glob("*_airr.tsv"))
        summary[sample] = {"airr_files": len(airr_files)}
        if not airr_files:
            continue

        # If multiple exist, use the largest by rows.
        best_df: Optional[pd.DataFrame] = None
        best_n = -1
        for af in airr_files:
            try:
                d = pd.read_csv(af, sep="\t")
                if d.shape[0] > best_n:
                    best_df = d
                    best_n = d.shape[0]
            except Exception:
                continue

        if best_df is None or best_df.empty:
            continue

        paired = _best_pairs_from_airr(best_df, sample)
        if not paired.empty:
            paired_rows.append(paired)
            summary[sample]["paired_cells"] = int(paired.shape[0])
        else:
            summary[sample]["paired_cells"] = 0

        bulk = _bulk_chain_abundance(best_df, sample)
        if not bulk.empty:
            bulk_rows.append(bulk)
            summary[sample]["chain_clonotypes"] = int(bulk.shape[0])
        else:
            summary[sample]["chain_clonotypes"] = 0

    if paired_rows:
        paired_df = pd.concat(paired_rows, ignore_index=True)
        paired_df.to_csv(outdir / "paired_alpha_beta_cells.csv", index=False)
        ab = (
            paired_df.groupby(["sample", "pair_id"], as_index=False)
            .size()
            .rename(columns={"size": "n_cells"})
        )
        ab["relative_abundance"] = ab["n_cells"] / ab.groupby("sample")["n_cells"].transform("sum")

        global_ab = paired_df.groupby(["pair_id"], as_index=False).size().rename(columns={"size": "n_cells"})
        global_ab["relative_abundance"] = global_ab["n_cells"] / global_ab["n_cells"].sum()
        global_ab["sample"] = "ALL"
        pd.concat([ab, global_ab], ignore_index=True).to_csv(
            outdir / "paired_alpha_beta_clonotype_abundance.csv", index=False
        )

    if bulk_rows:
        bulk_df = pd.concat(bulk_rows, ignore_index=True)
        bulk_df.to_csv(outdir / "bulk_chain_abundance.csv", index=False)

    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
