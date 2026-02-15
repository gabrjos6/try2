# Fast scRNA + TCR Pipeline (100-sample throughput)

This pipeline is optimized for **speed-first** processing of 10x 5' Gene Expression + VDJ FASTQs.

It performs:
- FASTQ preprocessing/alignment/quantification via `cellranger multi`
- Fast QC filtering
- Fast normalization (`normalize_total` + `log1p`)
- Batch correction with Harmony on PCA space
- TCR filtering for confident productive paired `TRA`/`TRB`
- Output of paired alpha/beta chains with relative abundance

## Why this design

For 10x data, `cellranger multi` is the most robust path to keep GEX and VDJ in sync and quickly produce high-confidence paired chains. Postprocessing is intentionally lightweight to stay within a short wall clock budget.

## Requirements

- `cellranger` (current recommended release)
- Python 3.10+
- Python packages:
  - `scanpy`
  - `anndata`
  - `pandas`
  - `numpy`
  - `harmonypy`

## Input format

Create a TSV with one sample per line:

`samples.tsv`
```tsv
sample_id	gex_fastqs	gex_fastq_id	vdj_fastqs	vdj_fastq_id	chemistry
S1	/data/fastqs	S1_GEX	/data/fastqs	S1_VDJ	auto
S2	/data/fastqs	S2_GEX	/data/fastqs	S2_VDJ	auto
```

Columns:
- `sample_id`: unique run/sample ID
- `gex_fastqs`: folder containing GEX FASTQs
- `gex_fastq_id`: FASTQ sample prefix for GEX
- `vdj_fastqs`: folder containing VDJ FASTQs
- `vdj_fastq_id`: FASTQ sample prefix for VDJ
- `chemistry`: optional (`auto` recommended)

## Run

```bash
bash scripts/run_fast_tcr_pipeline.sh \
  --samples /abs/path/samples.tsv \
  --transcriptome /abs/path/refdata-gex \
  --vdj-reference /abs/path/refdata-cellranger-vdj \
  --outdir /abs/path/results_fast_tcr \
  --total-cores 32 \
  --total-mem-gb 384 \
  --jobs 8
```

Notes:
- `jobs * cores_per_job` is capped at `total-cores`.
- For fastest total wall time with 32 cores, start with `jobs=8` (4 cores/job) and tune by sample size.

## Outputs

- `results_fast_tcr/cellranger_runs/<sample_id>/...` raw pipeline outputs per sample
- `results_fast_tcr/postprocess/paired_alpha_beta_cells.csv`
  - One row per cell with best productive high-confidence `TRA` + `TRB`
- `results_fast_tcr/postprocess/paired_alpha_beta_clonotype_abundance.csv`
  - Paired alpha/beta clonotypes with count and relative abundance per sample + global
- `results_fast_tcr/postprocess/adata_qc_norm_harmony.h5ad`
  - QC-filtered, normalized object with Harmony embedding in `obsm["X_pca_harmony"]`
- `results_fast_tcr/postprocess/run_summary.json`
  - Basic counts and processing summary

## Performance-first choices (intentional)

- Skips expensive optional tasks (full marker annotation, trajectory, differential testing)
- Uses lightweight QC thresholds and Harmony (fast and practical for many batches)
- Avoids slow imputation/integration methods that are not essential for downstream ML feature generation

## Caveats

- Assumes 10x Chromium 5' GEX + VDJ data model.
- If your FASTQs are not 10x-compatible, use a separate upstream aligner/caller path (e.g., STARsolo + MiXCR/TRUST4) and adapt postprocessing accordingly.
