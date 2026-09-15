import argparse
from pathlib import Path
import pandas as pd

from .schema import REQUIRED_METHYLATION_COLUMNS, missing_columns
from .qc import motif_prevalence, flag_samples_by_core_motif
from .matrix import build_matrix
from .cluster import pca
from .io.pacbio import import_pacbio_gff


def read_table(path):
    return pd.read_csv(path, sep="\t")


def assess_import_status(summary, sites, low_data_threshold=100):
    """Return a conservative import-level QC label.

    ``LOW_DATA`` is a warning, not a failure: very small PacBio GFFs may be
    technically valid but should not enter cohort analyses unnoticed. The
    threshold is intentionally configurable because expected call counts depend
    on organism and dataset.
    """
    n_sites = len(sites)
    if n_sites < low_data_threshold:
        return "LOW_DATA"
    if len(summary) == 0:
        return "NO_SUMMARY"
    return "PASS"


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
    ip.add_argument("--platform", default="pacbio")
    ip.add_argument("--min-score", type=float, help="optional minimum GFF column-6 score")
    ip.add_argument(
        "--min-identification-qv",
        type=float,
        help=(
            "optional minimum PacBio identificationQv; no universal default is "
            "assumed because historical pipelines differ"
        ),
    )
    ip.add_argument(
        "--keep-all-motif-positions",
        action="store_true",
        help=(
            "disable context-based cognate modified-position filtering; useful "
            "for diagnostics and legacy comparisons"
        ),
    )
    ip.add_argument(
        "--low-data-threshold",
        type=int,
        default=100,
        help=(
            "warn and set qc_status=LOW_DATA when fewer than this many "
            "site-by-motif assignments are parsed (default: 100; warning only)"
        ),
    )
    ip.add_argument("--motif-attribute", help="force a non-standard GFF motif attribute key")
    ip.add_argument("--modification-attribute", help="force a non-standard modification attribute key")
    ip.set_defaults(func=cmd_import_pacbio)

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
    return p


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
