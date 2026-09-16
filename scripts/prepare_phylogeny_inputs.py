#!/usr/bin/env python3
"""Stage QC-selected assemblies for a core-genome phylogeny.

The cohort manifest remains the source of truth for FASTA paths. This helper
selects samples by ``cohort.qc.tsv`` status, creates stable ``sample_id.fasta``
files (symlinks by default), and writes an audit table for downstream tools.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil

import pandas as pd


def resolve_manifest_path(manifest: Path, value: object) -> Path:
    path = Path(str(value))
    if not path.is_absolute():
        path = (manifest.parent / path).resolve()
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage MOMENTO QC-selected FASTAs for phylogeny")
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path, help="cohort.qc.tsv")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--qc-status", default="PASS", help="comma-separated statuses; default PASS")
    parser.add_argument("--mode", choices=["symlink", "copy"], default="symlink")
    parser.add_argument("--force", action="store_true", help="replace existing staged FASTA links/files")
    args = parser.parse_args()

    manifest = args.manifest.resolve()
    qc_path = args.qc.resolve()
    outdir = args.output_dir.resolve()
    genomes = outdir / "genomes"
    genomes.mkdir(parents=True, exist_ok=True)

    m = pd.read_csv(manifest, sep="\t")
    q = pd.read_csv(qc_path, sep="\t")
    required_manifest = {"sample_id", "fasta"}
    required_qc = {"sample_id", "qc_status"}
    if not required_manifest.issubset(m.columns):
        raise SystemExit("Manifest must contain sample_id and fasta columns")
    if not required_qc.issubset(q.columns):
        raise SystemExit("QC table must contain sample_id and qc_status columns")
    if m["sample_id"].duplicated().any():
        raise SystemExit("Manifest contains duplicate sample_id values")
    if q["sample_id"].duplicated().any():
        raise SystemExit("QC table contains duplicate sample_id values")

    keep_statuses = {x.strip() for x in args.qc_status.split(",") if x.strip()}
    selected = q.loc[q["qc_status"].astype(str).isin(keep_statuses), ["sample_id", "qc_status"]].copy()
    selected = selected.merge(m[["sample_id", "fasta"]], on="sample_id", how="left", validate="one_to_one")
    if selected["fasta"].isna().any():
        missing = selected.loc[selected["fasta"].isna(), "sample_id"].astype(str).tolist()
        raise SystemExit("Selected samples missing from manifest: " + ", ".join(missing))

    rows = []
    for row in selected.itertuples(index=False):
        sample = str(row.sample_id)
        source = resolve_manifest_path(manifest, row.fasta)
        if not source.exists():
            raise SystemExit(f"Missing FASTA for {sample}: {source}")
        staged = genomes / f"{sample}.fasta"
        if staged.exists() or staged.is_symlink():
            if not args.force:
                raise SystemExit(f"Staged FASTA already exists (use --force): {staged}")
            staged.unlink()
        if args.mode == "symlink":
            staged.symlink_to(source)
        else:
            shutil.copy2(source, staged)
        rows.append(
            {
                "sample_id": sample,
                "qc_status": row.qc_status,
                "source_fasta": str(source),
                "staged_fasta": str(staged),
                "mode": args.mode,
            }
        )

    audit = pd.DataFrame(rows).sort_values("sample_id")
    audit.to_csv(outdir / "phylogeny_samples.tsv", sep="\t", index=False)
    (outdir / "sample_ids.txt").write_text("\n".join(audit["sample_id"].astype(str)) + "\n")
    print(
        f"PHYLOGENY INPUTS: {len(audit)} samples with qc_status in "
        f"{','.join(sorted(keep_statuses))} -> {genomes}"
    )
    print(f"Audit: {outdir / 'phylogeny_samples.tsv'}")


if __name__ == "__main__":
    main()
