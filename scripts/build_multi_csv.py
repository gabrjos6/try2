#!/usr/bin/env python3
import argparse
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description="Build cellranger multi config CSV for a single sample.")
    p.add_argument("--sample-id", required=True)
    p.add_argument("--transcriptome", required=True)
    p.add_argument("--vdj-reference", required=True)
    p.add_argument("--gex-fastqs", required=True)
    p.add_argument("--gex-fastq-id", required=True)
    p.add_argument("--vdj-fastqs", required=True)
    p.add_argument("--vdj-fastq-id", required=True)
    p.add_argument("--chemistry", default="auto")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Minimal multi config for 5' GEX + VDJ-T.
    lines = [
        "[gene-expression]",
        f"reference,{Path(args.transcriptome).resolve().as_posix()}",
        f"chemistry,{args.chemistry}",
        "",
        "[vdj]",
        f"reference,{Path(args.vdj_reference).resolve().as_posix()}",
        "",
        "[libraries]",
        "fastq_id,fastqs,feature_types",
        f"{args.gex_fastq_id},{Path(args.gex_fastqs).resolve().as_posix()},Gene Expression",
        f"{args.vdj_fastq_id},{Path(args.vdj_fastqs).resolve().as_posix()},VDJ-T",
    ]

    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
