from pathlib import Path

import pandas as pd
import pytest

from momento.methylotype import (
    build_family_presence,
    canonical_motif_family,
    collapse_motif_families,
    find_nested_motif_pairs,
    nested_alignment,
    run_methylotype_analysis,
)


def _summary(rows):
    columns = [
        "sample_id", "platform", "motif", "modification", "called_sites",
        "genomic_sites", "methylated_fraction", "median_score", "source_file",
        "qc_status",
    ]
    return pd.DataFrame(rows, columns=columns)


def test_reverse_complement_motifs_collapse_to_one_family():
    assert canonical_motif_family("GACNNNNNNNTAG") == "CTANNNNNNNGTC"
    assert canonical_motif_family("CTANNNNNNNGTC") == "CTANNNNNNNGTC"

    df = _summary([
        ["S1", "pacbio", "GACNNNNNNNTAG", "m6A", 10, 10, 1.0, 100, "x.gff", "PASS"],
        ["S1", "pacbio", "CTANNNNNNNGTC", "m6A", 11, 11, 1.0, 100, "x.gff", "PASS"],
    ])
    collapsed = collapse_motif_families(df)
    assert len(collapsed) == 1
    assert collapsed.iloc[0]["motif_family"] == "CTANNNNNNNGTC"
    assert collapsed.iloc[0]["n_raw_motifs"] == 2


def test_iupac_nesting_detects_aaaattY_within_raatty():
    alignment = nested_alignment("RAATTY", "AAAATTY")
    assert alignment is not None
    assert alignment["child_start_1based"] in {1, 2}

    # The broader core motif is not itself nested in the more specific motif.
    assert nested_alignment("AAAATTY", "RAATTY") is None


def test_nested_pairs_report_cooccurrence_and_core_exclusion():
    df = _summary([
        ["S1", "pacbio", "RAATTY", "m6A", 100, 100, 1.0, 100, "x", "PASS"],
        ["S2", "pacbio", "RAATTY", "m6A", 100, 100, 1.0, 100, "x", "PASS"],
        ["S1", "pacbio", "AAAATTY", "m6A", 20, 40, 0.5, 100, "x", "PASS"],
    ])
    fam = collapse_motif_families(df)
    nested = find_nested_motif_pairs(fam)
    hit = nested.loc[
        (nested["parent_family"] == "RAATTY")
        & (nested["child_family"] == "AAAATTY")
    ].iloc[0]
    assert hit["child_isolates"] == 1
    assert hit["cooccurring_isolates"] == 1
    assert hit["child_with_parent_fraction"] == 1.0
    assert bool(hit["parent_excluded_from_accessory"])


def test_presence_matrix_excludes_core_but_keeps_all_samples():
    df = _summary([
        ["S1", "pacbio", "RAATTY", "m6A", 100, 100, 1.0, 100, "x", "PASS"],
        ["S2", "pacbio", "RAATTY", "m6A", 100, 100, 1.0, 100, "x", "PASS"],
        ["S2", "pacbio", "GATC", "m6A", 10, 10, 1.0, 100, "x", "PASS"],
    ])
    fam = collapse_motif_families(df)
    matrix = build_family_presence(fam)
    assert list(matrix.index) == ["S1", "S2"]
    assert "m6A|RAATTY" not in matrix.columns
    assert matrix.loc["S1", "m6A|GATC"] == 0
    assert matrix.loc["S2", "m6A|GATC"] == 1


def test_methylotype_analysis_uses_pass_samples_and_writes_outputs(tmp_path: Path):
    df = _summary([
        ["S1", "pacbio", "RAATTY", "m6A", 100, 100, 1.0, 100, "x", "PASS"],
        ["S1", "pacbio", "GATC", "m6A", 10, 10, 1.0, 100, "x", "PASS"],
        ["S2", "pacbio", "RAATTY", "m6A", 100, 100, 1.0, 100, "x", "PASS"],
        ["S2", "pacbio", "GTAC", "m6A", 10, 10, 1.0, 100, "x", "PASS"],
        ["S3", "pacbio", "RAATTY", "m6A", 100, 100, 1.0, 100, "x", "PASS"],
        ["S3", "pacbio", "GATC", "m6A", 10, 10, 1.0, 100, "x", "PASS"],
        ["BAD", "pacbio", "GCCAA", "m6A", 10, 10, 1.0, 100, "x", "PARTIAL_GFF"],
    ])
    result = run_methylotype_analysis(df, tmp_path)
    assert set(result["presence"].index) == {"S1", "S2", "S3"}
    assert "BAD" not in result["presence"].index
    assert result["presence"].shape == (3, 2)
    assert result["distance"].shape == (3, 3)
    assert result["pcoa"].shape[0] == 3

    expected = {
        "motif_families_long.tsv",
        "motif_family_prevalence.tsv",
        "nested_motif_pairs.tsv",
        "accessory_presence_absence.tsv",
        "jaccard_distance.tsv",
        "pcoa.tsv",
        "pcoa_variance.tsv",
        "cluster_order.tsv",
        "linkage.tsv",
    }
    assert expected.issubset({p.name for p in tmp_path.iterdir()})


def test_analysis_rejects_all_core_only_cohort(tmp_path: Path):
    df = _summary([
        ["S1", "pacbio", "RAATTY", "m6A", 100, 100, 1.0, 100, "x", "PASS"],
    ])
    with pytest.raises(ValueError, match="No accessory motif-family features"):
        run_methylotype_analysis(df, tmp_path)
