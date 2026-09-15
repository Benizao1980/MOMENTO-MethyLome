"""Manifest-driven cohort import and FASTA/GFF preflight QC for MOMENTO."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from .io.pacbio import (
    CANONICAL_COLUMNS,
    import_pacbio_gff,
    iter_gff,
    read_fasta,
)

MANIFEST_REQUIRED_COLUMNS = ("sample_id", "fasta", "gff")
AUDIT_COLUMNS = [
    "sample_id", "fasta", "gff", "manifest_note", "preflight_status",
    "context_usable", "context_compared", "context_exact", "context_coverage",
    "context_concordance", "gff_feature_records", "gff_matched_records",
    "gff_span_fraction", "gff_bin_coverage_fraction", "qc_status",
    "motif_rows", "site_assignments", "included_assignments", "message",
]
_SEQUENCE_COMPLEMENT = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMSWVHDBNtgcayrkmswvhdbn",
)


def _reverse_complement_sequence(sequence: str) -> str:
    """Reverse-complement assembly sequence while preserving unknown symbols.

    Legacy assemblies occasionally contain placeholder characters such as
    ``?``. Those characters should contribute a context mismatch if they differ
    from the PacBio GFF, rather than crashing preflight. Known IUPAC bases are
    complemented; unsupported symbols are retained after reversal.
    """
    return sequence.translate(_SEQUENCE_COMPLEMENT)[::-1].upper()


def _coordinate_coverage_metrics(
    seqs: dict[str, str],
    positions_by_seqid: dict[str, list[int]],
    *,
    bins: int = 20,
) -> tuple[float, float]:
    """Return genome-weighted GFF span and occupied-bin coverage fractions.

    ``gff_span_fraction`` measures the summed min-to-max coordinate span of GFF
    feature records on each FASTA record, weighted by FASTA length.
    ``gff_bin_coverage_fraction`` divides each FASTA record into equal genomic
    bins and asks what fraction of genome length lies in bins containing at
    least one GFF feature record. The latter catches large internal gaps that a
    simple min/max span can miss.
    """
    if bins < 1:
        raise ValueError("bins must be at least 1")

    genome_bases = sum(len(seq) for seq in seqs.values())
    if genome_bases == 0:
        return 0.0, 0.0

    span_bases = 0
    occupied_bases = 0

    for seqid, seq in seqs.items():
        length = len(seq)
        positions = positions_by_seqid.get(seqid, [])
        if not positions or length == 0:
            continue

        valid_positions = [pos for pos in positions if 1 <= pos <= length]
        if not valid_positions:
            continue

        span_bases += max(valid_positions) - min(valid_positions) + 1

        occupied_bins = {
            min(int((pos - 1) / length * bins), bins - 1)
            for pos in valid_positions
        }
        for bin_index in occupied_bins:
            left = (bin_index * length) // bins
            right = ((bin_index + 1) * length) // bins
            occupied_bases += right - left

    return span_bases / genome_bases, occupied_bases / genome_bases


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
    within a project directory. Extra columns such as ``note`` are preserved.
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
    min_gff_span_fraction: float = 0.90,
    min_gff_bin_coverage: float = 0.80,
    min_gff_records_for_span: int = 1000,
    gff_span_bins: int = 20,
) -> dict:
    """Check PacBio GFF/FASTA identity and whether the GFF spans the genome.

    ``context_coverage`` is the fraction of context-bearing GFF records that can
    be compared to the FASTA. ``context_concordance`` is the fraction of those
    compared records that match exactly. ``gff_span_fraction`` and
    ``gff_bin_coverage_fraction`` use all GFF feature records to identify
    partial/truncated modification files that can otherwise make unassessed
    genomic regions look biologically unmethylated.

    Coordinate-span QC is applied only when at least
    ``min_gff_records_for_span`` feature records are present, so genuinely sparse
    legacy GFFs remain governed by LOW_DATA rather than being misclassified as
    truncated. FASTA parsing is intentionally strict, so duplicate record IDs
    are a hard preflight error rather than being silently overwritten.
    """
    if not 0 <= min_coverage <= 1:
        raise ValueError("min_coverage must be between 0 and 1")
    if not 0 <= min_concordance <= 1:
        raise ValueError("min_concordance must be between 0 and 1")
    if not 0 <= min_gff_span_fraction <= 1:
        raise ValueError("min_gff_span_fraction must be between 0 and 1")
    if not 0 <= min_gff_bin_coverage <= 1:
        raise ValueError("min_gff_bin_coverage must be between 0 and 1")
    if min_gff_records_for_span < 1:
        raise ValueError("min_gff_records_for_span must be at least 1")
    if gff_span_bins < 1:
        raise ValueError("gff_span_bins must be at least 1")

    seqs = read_fasta(fasta_path)
    usable = compared = exact = 0
    gff_feature_records = 0
    gff_matched_records = 0
    positions_by_seqid: dict[str, list[int]] = {seqid: [] for seqid in seqs}

    for record in iter_gff(gff_path):
        gff_feature_records += 1
        seq = seqs.get(record["seqid"])
        if seq is not None:
            gff_matched_records += 1
            positions_by_seqid[record["seqid"]].append(record["start"])

        context = (record["attributes"].get("context") or "").upper().strip()
        if not context or len(context) % 2 == 0:
            continue
        usable += 1

        if seq is None:
            continue

        flank = len(context) // 2
        left = record["start"] - 1 - flank
        right = record["start"] + flank
        if left < 0 or right > len(seq):
            continue

        observed = seq[left:right]
        if record["strand"] == "-":
            observed = _reverse_complement_sequence(observed)

        compared += 1
        if observed == context:
            exact += 1

    coverage: Optional[float] = compared / usable if usable else None
    concordance: Optional[float] = exact / compared if compared else None
    gff_span_fraction, gff_bin_coverage_fraction = _coordinate_coverage_metrics(
        seqs,
        positions_by_seqid,
        bins=gff_span_bins,
    )

    span_assessed = gff_feature_records >= min_gff_records_for_span

    if gff_feature_records > 0 and gff_matched_records == 0:
        status = "FAIL"
        message = "No GFF feature seqids matched any FASTA record ID"
    elif usable and (coverage is None or coverage < min_coverage):
        status = "FAIL"
        message = f"Context coverage {coverage or 0:.4f} is below minimum {min_coverage:.4f}"
    elif usable and (concordance is None or concordance < min_concordance):
        status = "FAIL"
        message = f"Context concordance {concordance or 0:.4f} is below minimum {min_concordance:.4f}"
    elif span_assessed and (
        gff_span_fraction < min_gff_span_fraction
        or gff_bin_coverage_fraction < min_gff_bin_coverage
    ):
        status = "PARTIAL_GFF"
        reasons = []
        if gff_span_fraction < min_gff_span_fraction:
            reasons.append(
                f"coordinate span {gff_span_fraction:.4f} < {min_gff_span_fraction:.4f}"
            )
        if gff_bin_coverage_fraction < min_gff_bin_coverage:
            reasons.append(
                "occupied-bin genome fraction "
                f"{gff_bin_coverage_fraction:.4f} < {min_gff_bin_coverage:.4f}"
            )
        message = "Partial/truncated GFF suspected: " + "; ".join(reasons)
    elif usable == 0:
        status = "NO_CONTEXT"
        message = "No odd-length context attributes available; sequence concordance not assessed"
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
        "gff_feature_records": gff_feature_records,
        "gff_matched_records": gff_matched_records,
        "gff_span_fraction": gff_span_fraction,
        "gff_bin_coverage_fraction": gff_bin_coverage_fraction,
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
    min_gff_span_fraction: float = 0.90,
    min_gff_bin_coverage: float = 0.80,
    min_gff_records_for_span: int = 1000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Import all manifest samples, continuing past sample-level failures.

    Per-sample outputs are written under ``samples/<sample_id>/``. Cohort-level
    outputs are ``cohort.momento.tsv`` and ``cohort.qc.tsv``. A failed sample is
    recorded in the audit table but does not abort the remaining cohort.
    ``PARTIAL_GFF`` samples are imported for auditing but retain that QC label so
    they are excluded automatically by analyses restricted to PASS samples.
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
            "manifest_note": row.get("note", ""),
            "preflight_status": "NOT_RUN",
            "context_usable": 0,
            "context_compared": 0,
            "context_exact": 0,
            "context_coverage": None,
            "context_concordance": None,
            "gff_feature_records": 0,
            "gff_matched_records": 0,
            "gff_span_fraction": None,
            "gff_bin_coverage_fraction": None,
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
                min_gff_span_fraction=min_gff_span_fraction,
                min_gff_bin_coverage=min_gff_bin_coverage,
                min_gff_records_for_span=min_gff_records_for_span,
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

        import_status = assess_import_status(summary, sites, low_data_threshold)
        status = (
            "PARTIAL_GFF"
            if audit["preflight_status"] == "PARTIAL_GFF"
            else import_status
        )
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
        if import_status == "LOW_DATA" and not audit["message"]:
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
        combined = pd.DataFrame(columns=CANONICAL_COLUMNS)

    audit_df = pd.DataFrame(audit_rows, columns=AUDIT_COLUMNS)
    combined.to_csv(outdir / "cohort.momento.tsv", sep="\t", index=False)
    audit_df.to_csv(outdir / "cohort.qc.tsv", sep="\t", index=False)
    return combined, audit_df
