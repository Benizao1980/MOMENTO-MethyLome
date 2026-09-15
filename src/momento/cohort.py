"""Manifest-driven cohort import and FASTA/GFF preflight QC for MOMENTO."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from .io.pacbio import import_pacbio_gff, iter_gff, read_fasta, reverse_complement_iupac

MANIFEST_REQUIRED_COLUMNS = ("sample_id", "fasta", "gff")
AUDIT_COLUMNS = [
    "sample_id", "fasta", "gff", "preflight_status", "context_usable",
    "context_compared", "context_exact", "context_coverage",
    "context_concordance", "qc_status", "motif_rows", "site_assignments",
    "included_assignments", "message",
]


def assess_import_status(summary: pd.DataFrame, sites: pd.DataFrame, low_data_threshold: int = 100) -> str:
    """Return a conservative import-level QC label."""
    if len(sites) < low_data_threshold:
        return "LOW_DATA"
    if len(summary) == 0:
        return "NO_SUMMARY"
    return "PASS"


def load_manifest(path: str | Path) -> pd.DataFrame:
    """Load a TSV cohort manifest and resolve relative input paths.

    Required columns are ``sample_id``, ``fasta`` and ``gff``. Relative paths
    are interpreted relative to the manifest file, making manifests portable
    within a project directory.
    """
    manifest_path = Path(path)
    df = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    missing = [c for c in MANIFEST_REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError("Manifest missing required columns: " + ", ".join(missing))
    if df.empty:
        raise ValueError("Manifest contains no samples")
    if (df["sample_id"].str.strip() == "").any():
        raise ValueError("Manifest contains an empty sample_id")
    duplicated = df.loc[df["sample_id"].duplicated(keep=False), "sample_id"].tolist()
    if duplicated:
        raise ValueError("Manifest contains duplicate sample_id values: " + ", ".join(sorted(set(duplicated))))

    base = manifest_path.parent
    out = df.copy()
    for column in ("fasta", "gff"):
        resolved = []
        for value in out[column]:
            p = Path(value)
            if not p.is_absolute():
                p = base / p
            resolved.append(str(p))
        out[column] = resolved
    return out


def context_preflight(
    gff_path: str | Path,
    fasta_path: str | Path,
    *,
    min_coverage: float = 0.95,
    min_concordance: float = 0.99,
) -> dict:
    """Check whether PacBio context sequences agree with the nominated FASTA.

    ``context_coverage`` is the fraction of context-bearing GFF records that can
    be compared to the FASTA. ``context_concordance`` is the fraction of those
    compared records that match exactly. FASTA parsing is intentionally strict,
    so duplicate record IDs are a hard preflight error rather than being
    silently overwritten.
    """
    if not 0 <= min_coverage <= 1:
        raise ValueError("min_coverage must be between 0 and 1")
    if not 0 <= min_concordance <= 1:
        raise ValueError("min_concordance must be between 0 and 1")

    seqs = read_fasta(fasta_path)
    usable = compared = exact = 0

    for record in iter_gff(gff_path):
        context = (record["attributes"].get("context") or "").upper().strip()
        if not context or len(context) % 2 == 0:
            continue
        usable += 1

        seq = seqs.get(record["seqid"])
        if seq is None:
            continue

        flank = len(context) // 2
        left = record["start"] - 1 - flank
        right = record["start"] + flank
        if left < 0 or right > len(seq):
            continue

        observed = seq[left:right]
        if record["strand"] == "-":
            observed = reverse_complement_iupac(observed)

        compared += 1
        if observed == context:
            exact += 1

    coverage: Optional[float] = compared / usable if usable else None
    concordance: Optional[float] = exact / compared if compared else None

    if usable == 0:
        status = "NO_CONTEXT"
        message = "No odd-length context attributes available; sequence concordance not assessed"
    elif coverage is None or coverage < min_coverage:
        status = "FAIL"
        message = f"Context coverage {coverage or 0:.4f} is below minimum {min_coverage:.4f}"
    elif concordance is None or concordance < min_concordance:
        status = "FAIL"
        message = f"Context concordance {concordance or 0:.4f} is below minimum {min_concordance:.4f}"
    else:
        status = "PASS"
        message = ""

    return {
        "preflight_status": status,
        "context_usable": usable,
        "context_compared": compared,
        "context_exact": exact,
        "context_coverage": coverage,
        "context_concordance": concordance,
        "message": message,
    }


def build_cohort(
    manifest_path: str | Path,
    output_dir: str | Path,
    *,
    platform: str = "pacbio",
    min_score: Optional[float] = None,
    min_identification_qv: Optional[float] = None,
    motif_attribute: Optional[str] = None,
    modification_attribute: Optional[str] = None,
    filter_cognate_positions: bool = True,
    low_data_threshold: int = 100,
    min_context_coverage: float = 0.95,
    min_context_concordance: float = 0.99,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Import all manifest samples, continuing past sample-level failures.

    Per-sample outputs are written under ``samples/<sample_id>/``. Cohort-level
    outputs are ``cohort.momento.tsv`` and ``cohort.qc.tsv``. A failed sample is
    recorded in the audit table but does not abort the remaining cohort.
    """
    manifest = load_manifest(manifest_path)
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    summaries: list[pd.DataFrame] = []
    audit_rows: list[dict] = []

    for row in manifest.to_dict(orient="records"):
        sample = row["sample_id"]
        fasta = row["fasta"]
        gff = row["gff"]
        audit = {
            "sample_id": sample,
            "fasta": fasta,
            "gff": gff,
            "preflight_status": "NOT_RUN",
            "context_usable": 0,
            "context_compared": 0,
            "context_exact": 0,
            "context_coverage": None,
            "context_concordance": None,
            "qc_status": "FAIL",
            "motif_rows": 0,
            "site_assignments": 0,
            "included_assignments": 0,
            "message": "",
        }

        try:
            preflight = context_preflight(
                gff,
                fasta,
                min_coverage=min_context_coverage,
                min_concordance=min_context_concordance,
            )
            audit.update(preflight)
        except (OSError, ValueError) as exc:
            audit["preflight_status"] = "FAIL"
            audit["message"] = str(exc)
            audit_rows.append(audit)
            continue

        if audit["preflight_status"] == "FAIL":
            audit_rows.append(audit)
            continue

        try:
            summary, sites = import_pacbio_gff(
                gff,
                fasta,
                sample,
                platform=platform,
                min_score=min_score,
                min_identification_qv=min_identification_qv,
                motif_attribute=motif_attribute,
                modification_attribute=modification_attribute,
                filter_cognate_positions=filter_cognate_positions,
            )
        except (OSError, ValueError) as exc:
            audit["message"] = str(exc)
            audit_rows.append(audit)
            continue

        status = assess_import_status(summary, sites, low_data_threshold)
        summary = summary.copy()
        summary["qc_status"] = status

        sample_dir = outdir / "samples" / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        summary.to_csv(sample_dir / f"{sample}.momento.tsv", sep="\t", index=False)
        sites.to_csv(sample_dir / f"{sample}.sites.tsv", sep="\t", index=False)

        audit["qc_status"] = status
        audit["motif_rows"] = len(summary)
        audit["site_assignments"] = len(sites)
        audit["included_assignments"] = int(sites["included_in_summary"].sum()) if len(sites) else 0
        if status == "LOW_DATA" and not audit["message"]:
            audit["message"] = (
                f"Only {len(sites)} site-by-motif assignments parsed "
                f"(< low-data threshold {low_data_threshold})"
            )

        summaries.append(summary)
        audit_rows.append(audit)

    if summaries:
        combined = pd.concat(summaries, ignore_index=True)
        combined = combined.sort_values(["sample_id", "motif", "modification"]).reset_index(drop=True)
    else:
        combined = pd.DataFrame()

    audit_df = pd.DataFrame(audit_rows, columns=AUDIT_COLUMNS)
    combined.to_csv(outdir / "cohort.momento.tsv", sep="\t", index=False)
    audit_df.to_csv(outdir / "cohort.qc.tsv", sep="\t", index=False)
    return combined, audit_df
