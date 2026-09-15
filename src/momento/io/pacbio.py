"""PacBio methylation-site GFF importer for MOMENTO.

This module does not call methylation from raw signal. It standardises
platform-specific PacBio/SMRT modification outputs into MOMENTO's canonical
isolate-by-motif representation and uses the matching assembly to count genomic
motif opportunities.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import math
import re
import statistics
from typing import Dict, Iterable, Optional, Sequence

import pandas as pd

IUPAC_REGEX = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
    "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
    "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}
IUPAC_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R",
    "S": "S", "W": "W", "K": "M", "M": "K", "B": "V", "D": "H",
    "H": "D", "V": "B", "N": "N",
}
MOTIF_KEY_PRIORITY = (
    "motif", "motifstring", "motif_string", "motifseq", "motif_seq",
    "recognitionmotif", "recognition_motif", "recognitionsequence",
    "recognition_sequence",
)
MODIFICATION_KEY_PRIORITY = (
    "modification", "modificationtype", "modification_type", "modifiedbase",
    "modified_base", "modbase", "mod_base",
)
CANONICAL_COLUMNS = [
    "sample_id", "platform", "motif", "modification", "called_sites",
    "genomic_sites", "methylated_fraction", "median_score", "source_file",
    "qc_status",
]
SITE_COLUMNS = [
    "sample_id", "seqid", "start", "end", "strand", "motif",
    "modification", "score", "identification_qv", "context",
    "motif_position", "cognate_position", "cognate_position_source",
    "passes_score", "passes_identification_qv", "passes_motif_position",
    "included_in_summary", "motif_attribute", "source", "feature_type",
    "attributes_raw",
]
NUMERIC_TOKEN_RE = re.compile(r"[+-]?\d+(?:\.\d+)?")


def read_fasta(path: str | Path) -> Dict[str, str]:
    """Read a FASTA file using only the standard library."""
    seqs: Dict[str, list[str]] = {}
    current: Optional[str] = None
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in seqs:
                    raise ValueError(f"Duplicate FASTA record ID: {current}")
                seqs[current] = []
            else:
                if current is None:
                    raise ValueError("FASTA sequence encountered before first header")
                seqs[current].append(line.upper())
    if not seqs:
        raise ValueError(f"No FASTA records found in {path}")
    return {name: "".join(parts) for name, parts in seqs.items()}


def parse_attributes(text: str) -> Dict[str, str]:
    """Parse common GFF3/GFF key=value and key value attribute styles."""
    attrs: Dict[str, str] = {}
    for field in text.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(None, 1)
        else:
            key, value = field, ""
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def _normalise_key(key: str) -> str:
    return re.sub(r"[^a-z0-9_]", "", key.lower())


def choose_attribute(
    attrs: Dict[str, str], forced: Optional[str], priority: Sequence[str]
) -> Optional[str]:
    """Choose an explicit motif/modification attribute without guessing context."""
    if forced:
        for key in attrs:
            if key == forced or key.lower() == forced.lower():
                return key
        return None
    keymap = {_normalise_key(k): k for k in attrs}
    for candidate in priority:
        norm = _normalise_key(candidate)
        if norm in keymap:
            return keymap[norm]
    if priority is MOTIF_KEY_PRIORITY:
        hits = [k for k in attrs if "motif" in k.lower()]
    else:
        hits = [
            k for k in attrs
            if "modification" in k.lower() or "modifiedbase" in k.lower()
        ]
    return hits[0] if len(hits) == 1 else None


def _strip_pacbio_motif_position(text: str) -> str:
    """Strip an old SMRT/MotifMaker numeric modified-base suffix.

    Some PacBio motif GFF generations encode one motif as e.g. ``RAATTY,2``.
    The numeric suffix records modified-base-position metadata and is not part
    of the recognition sequence.
    """
    fields = [field.strip() for field in text.split(",")]
    if len(fields) == 1:
        return text
    if fields[0] and all(
        NUMERIC_TOKEN_RE.fullmatch(field) for field in fields[1:] if field
    ):
        return fields[0]
    raise ValueError(f"Unsupported comma-delimited motif value: {text!r}")


def clean_motif_component(text: str) -> str:
    """Convert common PacBio/modified-base motif notation to plain IUPAC."""
    motif = _strip_pacbio_motif_position(text.strip()).upper()
    for old, new in {
        "[6MA]": "A", "[M6A]": "A", "[4MC]": "C", "[M4C]": "C",
        "[5MC]": "C", "[M5C]": "C",
    }.items():
        motif = motif.replace(old, new)
    motif = re.sub(r"[\s^*\-]", "", motif)
    if not motif:
        raise ValueError("Empty motif")
    invalid = sorted({base for base in motif if base not in IUPAC_REGEX})
    if invalid:
        raise ValueError(f"Unsupported motif characters: {invalid}")
    return motif


def canonicalise_motif(raw: str) -> str:
    """Normalise one motif while preserving slash-separated partner motifs."""
    parts = [part for part in raw.split("/") if part.strip()]
    if not parts:
        raise ValueError(f"Empty motif value: {raw!r}")
    return "/".join(clean_motif_component(part) for part in parts)


def canonicalise_motif_values(raw: str) -> list[str]:
    """Return one or more canonical motifs encoded in one PacBio GFF value.

    Older PacBio/MotifMaker GFFs use commas in two different ways:

    * ``RAATTY,2`` -- one motif plus numeric modified-base-position metadata;
    * ``GCCAA,RAATTY`` -- one modification call assigned to overlapping motifs.

    Partner motifs containing ``/`` retain their paired representation.
    """
    raw = raw.strip()
    if not raw:
        raise ValueError("Empty motif value")

    if "/" in raw:
        return [canonicalise_motif(raw)]

    fields = [field.strip() for field in raw.split(",") if field.strip()]
    if not fields:
        raise ValueError(f"Empty motif value: {raw!r}")

    if len(fields) == 1 or all(
        NUMERIC_TOKEN_RE.fullmatch(field) for field in fields[1:]
    ):
        return [canonicalise_motif(raw)]

    motifs: list[str] = []
    for field in fields:
        if NUMERIC_TOKEN_RE.fullmatch(field):
            continue
        motif = canonicalise_motif(field)
        if motif not in motifs:
            motifs.append(motif)
    if not motifs:
        raise ValueError(f"No motif sequence found in value: {raw!r}")
    return motifs


def reverse_complement_iupac(motif: str) -> str:
    return "".join(IUPAC_COMPLEMENT[base] for base in reversed(motif))


def distinct_orientation_patterns(canonical_motif: str) -> list[str]:
    """Return unique recognition strings plus their reverse complements."""
    patterns: set[str] = set()
    for component in canonical_motif.split("/"):
        patterns.add(component)
        patterns.add(reverse_complement_iupac(component))
    return sorted(patterns)


def _iupac_pattern(motif: str) -> re.Pattern[str]:
    body = "".join(IUPAC_REGEX[base] for base in motif)
    return re.compile(f"(?=({body}))")


def count_genomic_opportunities(
    seqs: Dict[str, str], canonical_motif: str
) -> int:
    """Count physical motif occurrences across all assembly sequences.

    Reverse-complement-equivalent patterns are de-duplicated. For a
    self-reverse-complementary motif this therefore counts each physical motif
    once. Use :func:`count_genomic_target_sites` when a strand-specific
    modified-base position is known.
    """
    patterns = [_iupac_pattern(x) for x in distinct_orientation_patterns(canonical_motif)]
    return sum(
        sum(1 for _ in pattern.finditer(seq))
        for seq in seqs.values()
        for pattern in patterns
    )


def count_genomic_target_sites(
    seqs: Dict[str, str],
    canonical_motif: str,
    cognate_position: Optional[int],
) -> int:
    """Count strand-specific methylatable target sites.

    For a single self-reverse-complementary motif, a physical motif occurrence
    contains one target on each DNA strand once a cognate modified position is
    defined. Therefore the strand-specific opportunity count is twice the
    physical motif count. Non-palindromic motifs are already represented by
    their two distinct orientations in :func:`count_genomic_opportunities`.
    """
    physical = count_genomic_opportunities(seqs, canonical_motif)
    components = canonical_motif.split("/")
    if (
        cognate_position is not None
        and len(components) == 1
        and reverse_complement_iupac(components[0]) == components[0]
    ):
        return physical * 2
    return physical


def motif_positions_in_context(context: str, canonical_motif: str) -> list[int]:
    """Return 1-based positions of the context-centre base within a motif.

    PacBio/MotifMaker ``context`` strings place the called base at the centre.
    The function tests each slash-separated motif component and returns every
    distinct motif-relative position whose match covers the centre base.
    """
    context = (context or "").upper().strip()
    if not context:
        return []
    centre = len(context) // 2
    positions: set[int] = set()
    for component in canonical_motif.split("/"):
        pattern = _iupac_pattern(component)
        for match in pattern.finditer(context):
            start = match.start()
            if start <= centre < start + len(component):
                positions.add(centre - start + 1)
    return sorted(positions)


def _finite_float(value: object) -> Optional[float]:
    if value in (None, "", "."):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def infer_cognate_position(
    rows: Sequence[dict],
    *,
    min_median_qv_delta: float = 20.0,
) -> tuple[Optional[int], str]:
    """Infer the cognate modified-base position for one motif/modification.

    The inference is deliberately conservative:

    * one observed motif position -> use it;
    * multiple positions -> use the position with the highest median
      ``identificationQv`` only when it exceeds the runner-up by at least
      ``min_median_qv_delta``;
    * otherwise leave the position unresolved and do not position-filter.

    This prevents a motif annotation alone from being treated as proof that
    every called base inside that motif is the biologically cognate target.
    """
    by_position: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        position = row.get("motif_position")
        if isinstance(position, int):
            by_position[position].append(row)

    if not by_position:
        return None, "no_context_position"
    if len(by_position) == 1:
        return next(iter(by_position)), "single_observed_position"

    qv_medians: list[tuple[float, int]] = []
    for position, pos_rows in by_position.items():
        qvs = [
            float(row["identification_qv"])
            for row in pos_rows
            if row.get("identification_qv") is not None
            and math.isfinite(float(row["identification_qv"]))
        ]
        if qvs:
            qv_medians.append((statistics.median(qvs), position))

    qv_medians.sort(reverse=True)
    if len(qv_medians) >= 2:
        best_qv, best_position = qv_medians[0]
        second_qv = qv_medians[1][0]
        if best_qv - second_qv >= min_median_qv_delta:
            return best_position, "highest_median_identification_qv"
    elif len(qv_medians) == 1:
        return qv_medians[0][1], "only_position_with_identification_qv"

    return None, "ambiguous_multiple_positions"


def iter_gff(path: str | Path) -> Iterable[dict]:
    """Yield parsed non-comment GFF records."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                raise ValueError(
                    f"{path}:{line_no}: expected 9 GFF columns, found {len(cols)}"
                )
            yield {
                "line_no": line_no,
                "seqid": cols[0],
                "source": cols[1],
                "feature_type": cols[2],
                "start": int(cols[3]),
                "end": int(cols[4]),
                "score": None if cols[5] == "." else float(cols[5]),
                "strand": cols[6],
                "phase": cols[7],
                "attributes_raw": cols[8],
                "attributes": parse_attributes(cols[8]),
            }


def import_pacbio_gff(
    gff_path: str | Path,
    fasta_path: str | Path,
    sample_id: str,
    *,
    platform: str = "pacbio",
    min_score: Optional[float] = None,
    min_identification_qv: Optional[float] = None,
    motif_attribute: Optional[str] = None,
    modification_attribute: Optional[str] = None,
    filter_cognate_positions: bool = True,
    min_cognate_qv_delta: float = 20.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Import one PacBio methylation GFF and matching assembly.

    The importer first preserves all motif assignments in a site table. When a
    PacBio ``context`` attribute is available it determines the called base's
    motif-relative position, then infers the cognate modified position for each
    motif/modification combination. Position filtering is applied only when the
    inference is unambiguous.

    ``min_identification_qv`` is optional and deliberately has no universal
    default. Different historical PacBio pipelines can use different confidence
    conventions. Calls excluded from the motif summary remain in the site table
    with explicit pass/fail flags.

    Returns ``(motif_summary, site_table)``.
    """
    seqs = read_fasta(fasta_path)
    sites: list[dict] = []
    observed_keys: Counter[str] = Counter()
    examples: list[str] = []
    feature_count = 0

    for rec in iter_gff(gff_path):
        feature_count += 1
        attrs = rec["attributes"]
        observed_keys.update(attrs.keys())
        if len(examples) < 3:
            examples.append(str(rec["attributes_raw"]))

        motif_key = choose_attribute(attrs, motif_attribute, MOTIF_KEY_PRIORITY)
        if motif_key is None or not attrs.get(motif_key):
            continue

        motifs = canonicalise_motif_values(attrs[motif_key])
        mod_key = choose_attribute(
            attrs, modification_attribute, MODIFICATION_KEY_PRIORITY
        )
        modification = (
            attrs[mod_key]
            if mod_key and attrs.get(mod_key)
            else str(rec["feature_type"])
        )
        score = rec["score"]
        identification_qv = _finite_float(attrs.get("identificationQv"))
        context = attrs.get("context", "")

        for motif in motifs:
            positions = motif_positions_in_context(context, motif)
            motif_position = positions[0] if len(positions) == 1 else None
            sites.append({
                "sample_id": sample_id,
                "seqid": rec["seqid"],
                "start": rec["start"],
                "end": rec["end"],
                "strand": rec["strand"],
                "motif": motif,
                "modification": modification,
                "score": score,
                "identification_qv": identification_qv,
                "context": context,
                "motif_position": motif_position,
                "cognate_position": None,
                "cognate_position_source": "",
                "passes_score": True,
                "passes_identification_qv": True,
                "passes_motif_position": True,
                "included_in_summary": True,
                "motif_attribute": motif_key,
                "source": rec["source"],
                "feature_type": rec["feature_type"],
                "attributes_raw": rec["attributes_raw"],
            })

    if feature_count == 0:
        raise ValueError(f"No feature records found in {gff_path}")
    if not sites:
        keys = ", ".join(key for key, _ in observed_keys.most_common()) or "<none>"
        preview = " | ".join(examples) or "<none>"
        raise ValueError(
            "Could not identify an explicit motif attribute in PacBio GFF. "
            f"Observed attribute keys: {keys}. Example attributes: {preview}. "
            "If the motif uses another field, pass motif_attribute explicitly."
        )

    rows_by_key: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in sites:
        rows_by_key[(row["motif"], row["modification"])].append(row)

    cognate_by_key: dict[tuple[str, str], Optional[int]] = {}
    cognate_source_by_key: dict[tuple[str, str], str] = {}
    for key, key_rows in rows_by_key.items():
        cognate, source = infer_cognate_position(
            key_rows, min_median_qv_delta=min_cognate_qv_delta
        )
        cognate_by_key[key] = cognate
        cognate_source_by_key[key] = source

    for row in sites:
        key = (row["motif"], row["modification"])
        cognate = cognate_by_key[key]
        row["cognate_position"] = cognate
        row["cognate_position_source"] = cognate_source_by_key[key]

        score = row["score"]
        row["passes_score"] = (
            min_score is None
            or score is None
            or (math.isfinite(float(score)) and float(score) >= min_score)
        )

        identification_qv = row["identification_qv"]
        row["passes_identification_qv"] = (
            min_identification_qv is None
            or (
                identification_qv is not None
                and math.isfinite(float(identification_qv))
                and float(identification_qv) >= min_identification_qv
            )
        )

        motif_position = row["motif_position"]
        row["passes_motif_position"] = (
            not filter_cognate_positions
            or cognate is None
            or motif_position is None
            or motif_position == cognate
        )

        row["included_in_summary"] = bool(
            row["passes_score"]
            and row["passes_identification_qv"]
            and row["passes_motif_position"]
        )

    grouped: Counter[tuple[str, str]] = Counter()
    scores: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in sites:
        if not row["included_in_summary"]:
            continue
        key = (row["motif"], row["modification"])
        grouped[key] += 1
        score = row["score"]
        if score is not None and math.isfinite(float(score)):
            scores[key].append(float(score))

    opportunity_cache: dict[tuple[str, Optional[int]], int] = {}
    summary_rows: list[dict] = []
    for (motif, modification), called_sites in sorted(grouped.items()):
        cognate = cognate_by_key[(motif, modification)]
        opportunity_key = (motif, cognate)
        if opportunity_key not in opportunity_cache:
            opportunity_cache[opportunity_key] = count_genomic_target_sites(
                seqs, motif, cognate
            )
        genomic_sites = opportunity_cache[opportunity_key]
        fraction = called_sites / genomic_sites if genomic_sites else math.nan
        median_score = (
            statistics.median(scores[(motif, modification)])
            if scores[(motif, modification)]
            else math.nan
        )
        summary_rows.append({
            "sample_id": sample_id,
            "platform": platform,
            "motif": motif,
            "modification": modification,
            "called_sites": int(called_sites),
            "genomic_sites": int(genomic_sites),
            "methylated_fraction": fraction,
            "median_score": median_score,
            "source_file": Path(gff_path).name,
            "qc_status": "UNASSESSED",
        })

    return (
        pd.DataFrame(summary_rows, columns=CANONICAL_COLUMNS),
        pd.DataFrame(sites, columns=SITE_COLUMNS),
    )


def standardise_summary_table(
    path: str | Path,
    sample_id_col: str = "sample_id",
    motif_col: str = "motif",
    count_col: str = "called_sites",
) -> pd.DataFrame:
    """Compatibility adapter for already summarised PacBio tables."""
    df = pd.read_csv(path, sep=None, engine="python")
    return pd.DataFrame({
        "sample_id": df[sample_id_col].astype(str),
        "platform": "pacbio",
        "motif": df[motif_col].astype(str),
        "modification": df.get("modification", "unknown"),
        "called_sites": pd.to_numeric(df[count_col], errors="coerce"),
        "genomic_sites": pd.to_numeric(df.get("genomic_sites"), errors="coerce"),
        "methylated_fraction": pd.to_numeric(
            df.get("methylated_fraction"), errors="coerce"
        ),
        "median_score": pd.to_numeric(df.get("median_score"), errors="coerce"),
        "source_file": str(path),
        "qc_status": "UNASSESSED",
    })
