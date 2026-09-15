import argparse
from pathlib import Path
import pandas as pd

from .schema import REQUIRED_METHYLATION_COLUMNS, missing_columns
from .qc import motif_prevalence, flag_samples_by_core_motif
from .matrix import build_matrix
from .cluster import pca
from .cohort import assess_import_status, build_cohort
from .io.pacbio import import_pacbio_gff
from .methylotype import run_methylotype_analysis


def read_table(path):
    return pd.read_csv(path, sep="\t")


def cmd_validate(a):
    df = read_table(a.input)
    missing = missing_columns(df.columns, REQUIRED_METHYLATION_COLUMNS)
    if missing:
        raise SystemExit("Missing required columns: " + ", ".join(missing))
    print(f"PASS: {len(df):,} rows; {df.sample_id.nunique()} samples; {df.motif.nunique()} motifs")


def cmd_import_pacbio(a):
    try:
        summary, sites = import_pacbio_gff(
            a.gff,
            a.fasta,
            a.sample,
            platform=a.platform,
            min_score=a.min_score,
            min_identification_qv=a.min_identification_qv,
            motif_attribute=a.motif_attribute,
            modification_attribute=a.modification_attribute,
            filter_cognate_positions=not a.keep_all_motif_positions,
        )
    except (OSError, ValueError) as exc:
        raise SystemExit(f"FAIL: {a.sample}: {exc}") from exc

    status = assess_import_status(summary, sites, a.low_data_threshold)
    summary = summary.copy()
    summary["qc_status"] = status

    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(a.output, sep="\t", index=False)
    if a.sites_output:
        Path(a.sites_output).parent.mkdir(parents=True, exist_ok=True)
        sites.to_csv(a.sites_output, sep="\t", index=False)

    included = int(sites["included_in_summary"].sum()) if len(sites) else 0
    prefix = "PASS" if status == "PASS" else status
    print(
        f"{prefix}: {a.sample}; {len(summary)} motif/modification rows; "
        f"{included:,}/{len(sites):,} site-by-motif assignments included -> {a.output}"
    )
    if status == "LOW_DATA":
        print(
            f"WARNING: only {len(sites):,} site-by-motif assignments were parsed "
            f"(< --low-data-threshold {a.low_data_threshold:,}); review this sample "
            "before cohort analysis."
        )


def cmd_build_cohort(a):
    try:
        combined, audit = build_cohort(
            a.manifest,
            a.output_dir,
            platform=a.platform,
            min_score=a.min_score,
            min_identification_qv=a.min_identification_qv,
            motif_attribute=a.motif_attribute,
            modification_attribute=a.modification_attribute,
            filter_cognate_positions=not a.keep_all_motif_positions,
            low_data_threshold=a.low_data_threshold,
            min_context_coverage=a.min_context_coverage,
            min_context_concordance=a.min_context_concordance,
            min_gff_span_fraction=a.min_gff_span_fraction,
            min_gff_bin_coverage=a.min_gff_bin_coverage,
            min_gff_records_for_span=a.min_gff_records_for_span,
        )
    except (OSError, ValueError) as exc:
        raise SystemExit(f"FAIL: cohort: {exc}") from exc

    counts = audit["qc_status"].value_counts().to_dict() if len(audit) else {}
    status_text = ", ".join(f"{key}={value}" for key, value in sorted(counts.items()))
    imported_statuses = ["PASS", "LOW_DATA", "NO_SUMMARY", "PARTIAL_GFF"]
    imported = int((audit["qc_status"].isin(imported_statuses)).sum()) if len(audit) else 0
    print(
        f"COHORT: {imported}/{len(audit)} samples imported; {len(combined):,} motif/modification rows"
        + (f"; {status_text}" if status_text else "")
    )
    print(f"Summary: {Path(a.output_dir) / 'cohort.momento.tsv'}")
    print(f"Audit:   {Path(a.output_dir) / 'cohort.qc.tsv'}")


def cmd_matrix(a):
    df = read_table(a.input)
    exclude = [x for x in a.exclude.split(",") if x] if a.exclude else None
    mat = build_matrix(df, a.value, a.min_called_sites, exclude)
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    mat.to_csv(a.output, sep="\t")
    print(a.output)


def cmd_qc(a):
    df = read_table(a.input)
    motif_prevalence(df, a.min_called_sites).to_csv(a.prevalence_output, sep="\t", index=False)
    if a.core_motif:
        flag_samples_by_core_motif(df, a.core_motif).to_csv(a.core_output, sep="\t", index=False)


def cmd_cluster(a):
    mat = pd.read_csv(a.matrix, sep="\t", index_col=0)
    pca(mat).to_csv(a.output, sep="\t")


def _csv_values(text):
    return tuple(x.strip() for x in text.split(",") if x.strip())


def cmd_methylotypes(a):
    df = read_table(a.input)
    try:
        result = run_methylotype_analysis(
            df,
            a.output_dir,
            qc_statuses=_csv_values(a.qc_statuses),
            modifications=_csv_values(a.modifications),
            exclude_families=_csv_values(a.exclude_families),
            min_called_sites=a.min_called_sites,
        )
    except (OSError, ValueError) as exc:
        raise SystemExit(f"FAIL: methylotypes: {exc}") from exc

    family_long = result["family_long"]
    presence = result["presence"]
    nested = result["nested"]
    variance = result["variance"]
    raw_features = set()
    for row in family_long.itertuples():
        for motif in str(row.raw_motifs).split(";"):
            if motif:
                raw_features.add((row.modification, motif))
    raw_count = len(raw_features)
    family_count = family_long[["motif_family", "modification"]].drop_duplicates().shape[0]
    print(
        f"METHYLOTYPES: {presence.shape[0]} samples; {raw_count} raw motif/modification features; "
        f"{family_count} reverse-complement families; {presence.shape[1]} accessory features; "
        f"{len(nested)} nested-family relationships"
    )
    if len(variance):
        axes = ", ".join(
            f"{r.axis}={r.variance_fraction:.3f}" for r in variance.itertuples()
        )
        print(f"PCoA positive-eigenvalue variance fractions: {axes}")
    print(f"Outputs: {Path(a.output_dir)}")


def add_pacbio_import_options(parser):
    parser.add_argument("--platform", default="pacbio")
    parser.add_argument("--min-score", type=float, help="optional minimum GFF column-6 score")
    parser.add_argument(
        "--min-identification-qv",
        type=float,
        help=(
            "optional minimum PacBio identificationQv; no universal default is "
            "assumed because historical pipelines differ"
        ),
    )
    parser.add_argument(
        "--keep-all-motif-positions",
        action="store_true",
        help=(
            "disable context-based cognate modified-position filtering; useful "
            "for diagnostics and legacy comparisons"
        ),
    )
    parser.add_argument(
        "--low-data-threshold",
        type=int,
        default=100,
        help=(
            "warn and set qc_status=LOW_DATA when fewer than this many "
            "site-by-motif assignments are parsed (default: 100; warning only)"
        ),
    )
    parser.add_argument("--motif-attribute", help="force a non-standard GFF motif attribute key")
    parser.add_argument("--modification-attribute", help="force a non-standard modification attribute key")


def build_parser():
    p = argparse.ArgumentParser(prog="momento", description="MOMENTO: comparative bacterial methylomics")
    sub = p.add_subparsers(dest="command", required=True)

    v = sub.add_parser("validate", help="validate a canonical MOMENTO long table")
    v.add_argument("--input", required=True)
    v.set_defaults(func=cmd_validate)

    ip = sub.add_parser("import-pacbio", help="import PacBio methylation GFF + matching assembly FASTA")
    ip.add_argument("--gff", required=True, help="PacBio methylation/motif GFF")
    ip.add_argument("--fasta", required=True, help="matching assembly FASTA")
    ip.add_argument("--sample", required=True, help="stable MOMENTO sample ID")
    ip.add_argument("--output", required=True, help="motif-level MOMENTO TSV")
    ip.add_argument("--sites-output", help="optional site-level TSV with QC/filter flags")
    add_pacbio_import_options(ip)
    ip.set_defaults(func=cmd_import_pacbio)

    bc = sub.add_parser(
        "build-cohort",
        help="import a manifest-defined PacBio cohort with FASTA/GFF preflight QC",
    )
    bc.add_argument(
        "--manifest",
        required=True,
        help="TSV with required columns sample_id, fasta, gff; relative paths resolve from the manifest",
    )
    bc.add_argument("--output-dir", required=True, help="cohort output directory")
    bc.add_argument(
        "--min-context-coverage",
        type=float,
        default=0.95,
        help="minimum fraction of GFF contexts comparable to the FASTA (default: 0.95)",
    )
    bc.add_argument(
        "--min-context-concordance",
        type=float,
        default=0.99,
        help="minimum exact-match fraction among compared contexts (default: 0.99)",
    )
    bc.add_argument(
        "--min-gff-span-fraction",
        type=float,
        default=0.90,
        help=(
            "minimum genome-weighted min-to-max coordinate span for sufficiently "
            "large GFFs before qc_status=PARTIAL_GFF (default: 0.90)"
        ),
    )
    bc.add_argument(
        "--min-gff-bin-coverage",
        type=float,
        default=0.80,
        help=(
            "minimum genome fraction in 20-bin windows containing GFF features "
            "before qc_status=PARTIAL_GFF (default: 0.80)"
        ),
    )
    bc.add_argument(
        "--min-gff-records-for-span",
        type=int,
        default=1000,
        help=(
            "minimum GFF feature records required before coordinate-span QC is "
            "enforced (default: 1000)"
        ),
    )
    add_pacbio_import_options(bc)
    bc.set_defaults(func=cmd_build_cohort)

    m = sub.add_parser("matrix")
    m.add_argument("--input", required=True)
    m.add_argument("--output", required=True)
    m.add_argument("--value", choices=["presence", "called_sites", "methylated_fraction"], default="presence")
    m.add_argument("--min-called-sites", type=int, default=1)
    m.add_argument("--exclude", default="")
    m.set_defaults(func=cmd_matrix)

    q = sub.add_parser("qc")
    q.add_argument("--input", required=True)
    q.add_argument("--prevalence-output", required=True)
    q.add_argument("--min-called-sites", type=int, default=1)
    q.add_argument("--core-motif")
    q.add_argument("--core-output", default="core_motif_qc.tsv")
    q.set_defaults(func=cmd_qc)

    c = sub.add_parser("cluster")
    c.add_argument("--matrix", required=True)
    c.add_argument("--output", required=True)
    c.set_defaults(func=cmd_cluster)

    mt = sub.add_parser(
        "methylotypes",
        help="collapse reverse-complement motif families and run binary accessory methylotype analysis",
    )
    mt.add_argument("--input", required=True, help="cohort.momento.tsv")
    mt.add_argument("--output-dir", required=True)
    mt.add_argument(
        "--qc-statuses",
        default="PASS",
        help="comma-separated sample QC labels to include (default: PASS)",
    )
    mt.add_argument(
        "--modifications",
        default="m6A,m4C",
        help="comma-separated modification classes to include (default: m6A,m4C)",
    )
    mt.add_argument(
        "--exclude-families",
        default="RAATTY",
        help="comma-separated motif families to exclude from accessory matrix (default: RAATTY)",
    )
    mt.add_argument("--min-called-sites", type=int, default=1)
    mt.set_defaults(func=cmd_methylotypes)
    return p


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
