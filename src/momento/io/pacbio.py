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
    "modification", "score", "motif_attribute", "source", "feature_type",
    "attributes_raw",
]


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


def clean_motif_component(text: str) -> str:
    """Convert common modified-base notation to a plain IUPAC motif."""
    motif = text.strip().upper()
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
    """Normalise a motif while preserving slash-separated paired motifs."""
    parts = [part for part in raw.split("/") if part.strip()]
    if not parts:
        raise ValueError(f"Empty motif value: {raw!r}")
    return "/".join(clean_motif_component(part) for part in parts)


def reverse_complement_iupac(motif: str) -> str:
    return "".join(IUPAC_COMPLEMENT[base] for base in reversed(motif))


def distinct_orientation_patterns(canonical_motif: str) -> list[str]:
    """Return unique explicit recognition strings plus reverse complements."""
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
    """Count overlapping motif opportunities across all assembly sequences."""
    patterns = [_iupac_pattern(x) for x in distinct_orientation_patterns(canonical_motif)]
    return sum(sum(1 for _ in pattern.finditer(seq)) for seq in seqs.values() for pattern in patterns)


def iter_gff(path: str | Path) -> Iterable[dict]:
    """Yield parsed non-comment GFF records."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                raise ValueError(f"{path}:{line_no}: expected 9 GFF columns, found {len(cols)}")
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
    motif_attribute: Optional[str] = None,
    modification_attribute: Optional[str] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Import one PacBio methylation GFF and matching assembly.

    Returns ``(motif_summary, site_table)``. The site table preserves the raw GFF
    attribute field for reproducibility and later position-level analyses.
    """
    seqs = read_fasta(fasta_path)
    grouped: Counter[tuple[str, str]] = Counter()
    scores: dict[tuple[str, str], list[float]] = defaultdict(list)
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
        score = rec["score"]
        if min_score is not None and score is not None and score < min_score:
            continue
        motif_key = choose_attribute(attrs, motif_attribute, MOTIF_KEY_PRIORITY)
        if motif_key is None or not attrs.get(motif_key):
            continue
        motif = canonicalise_motif(attrs[motif_key])
        mod_key = choose_attribute(attrs, modification_attribute, MODIFICATION_KEY_PRIORITY)
        modification = attrs[mod_key] if mod_key and attrs.get(mod_key) else str(rec["feature_type"])
        key = (motif, modification)
        grouped[key] += 1
        if score is not None and math.isfinite(score):
            scores[key].append(float(score))
        sites.append({
            "sample_id": sample_id,
            "seqid": rec["seqid"], "start": rec["start"], "end": rec["end"],
            "strand": rec["strand"], "motif": motif, "modification": modification,
            "score": score, "motif_attribute": motif_key, "source": rec["source"],
            "feature_type": rec["feature_type"], "attributes_raw": rec["attributes_raw"],
        })

    if feature_count == 0:
        raise ValueError(f"No feature records found in {gff_path}")
    if not grouped:
        keys = ", ".join(key for key, _ in observed_keys.most_common()) or "<none>"
        preview = " | ".join(examples) or "<none>"
        raise ValueError(
            "Could not identify an explicit motif attribute in PacBio GFF. "
            f"Observed attribute keys: {keys}. Example attributes: {preview}. "
            "If the motif uses another field, pass motif_attribute explicitly."
        )

    opportunity_cache: dict[str, int] = {}
    summary_rows: list[dict] = []
    for (motif, modification), called_sites in sorted(grouped.items()):
        genomic_sites = opportunity_cache.setdefault(
            motif, count_genomic_opportunities(seqs, motif)
        )
        fraction = called_sites / genomic_sites if genomic_sites else math.nan
        median_score = statistics.median(scores[(motif, modification)]) if scores[(motif, modification)] else math.nan
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
        "methylated_fraction": pd.to_numeric(df.get("methylated_fraction"), errors="coerce"),
        "median_score": pd.to_numeric(df.get("median_score"), errors="coerce"),
        "source_file": str(path),
        "qc_status": "UNASSESSED",
    })
