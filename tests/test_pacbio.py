import math

import pytest

from momento.io.pacbio import (
    canonicalise_motif,
    canonicalise_motif_values,
    count_genomic_opportunities,
    count_genomic_target_sites,
    import_pacbio_gff,
    motif_positions_in_context,
    parse_attributes,
    reverse_complement_iupac,
)


def _context_with_call(motif: str, position: int) -> str:
    """Place a concrete motif in a 41-bp context with position at the centre."""
    centre = 20
    start = centre - (position - 1)
    context = list("C" * 41)
    context[start:start + len(motif)] = list(motif)
    return "".join(context)


def test_parse_attributes_gff3():
    attrs = parse_attributes("motif=RAATTY;coverage=55;modification=m6A")
    assert attrs["motif"] == "RAATTY"
    assert attrs["coverage"] == "55"
    assert attrs["modification"] == "m6A"


def test_iupac_reverse_complement_and_target_counting():
    seqs = {"chr": "GAATTCAAATTTGATC", "plasmid": "GAATTC"}
    assert canonicalise_motif("RAATTY") == "RAATTY"
    assert reverse_complement_iupac("RAATTY") == "RAATTY"
    assert count_genomic_opportunities(seqs, "RAATTY") == 3
    assert count_genomic_target_sites(seqs, "RAATTY", 3) == 6
    assert count_genomic_opportunities(seqs, "GATC") == 1


def test_motif_position_from_pacbio_context():
    position3 = _context_with_call("GAATTC", 3)
    position2 = _context_with_call("GAATTC", 2)
    assert motif_positions_in_context(position3, "RAATTY") == [3]
    assert motif_positions_in_context(position2, "RAATTY") == [2]
    assert motif_positions_in_context("", "RAATTY") == []


def test_old_pacbio_numeric_motif_position_suffix():
    assert canonicalise_motif("RAATTY,2") == "RAATTY"
    assert canonicalise_motif("CATG,1") == "CATG"
    assert canonicalise_motif(
        "AGTNNNNNNRTTG,0/CAAYNNNNNNACT,12"
    ) == "AGTNNNNNNRTTG/CAAYNNNNNNACT"
    with pytest.raises(ValueError, match="Unsupported comma-delimited motif value"):
        canonicalise_motif("RAATTY,unexpected")


def test_overlapping_comma_separated_motif_values():
    assert canonicalise_motif_values("GCCAA,RAATTY") == ["GCCAA", "RAATTY"]
    assert canonicalise_motif_values("GCCAA,RAATTY,2") == ["GCCAA", "RAATTY"]
    assert canonicalise_motif_values("RAATTY,2") == ["RAATTY"]


def test_import_infers_cognate_position_and_keeps_excluded_sites(tmp_path):
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(">chr\nGAATTCGAATTC\n")

    pos2 = _context_with_call("GAATTC", 2)
    pos3 = _context_with_call("GAATTC", 3)
    gff = tmp_path / "toy.gff"
    gff.write_text(
        "##gff-version 3\n"
        f"chr\tkinModCall\tm6A\t1\t1\t90\t+\t.\tmotif=RAATTY;context={pos2};identificationQv=10\n"
        f"chr\tkinModCall\tm6A\t2\t2\t90\t+\t.\tmotif=RAATTY;context={pos2};identificationQv=15\n"
        f"chr\tkinModCall\tm6A\t3\t3\t200\t+\t.\tmotif=RAATTY;context={pos3};identificationQv=70\n"
        f"chr\tkinModCall\tm6A\t4\t4\t200\t+\t.\tmotif=RAATTY;context={pos3};identificationQv=100\n"
        f"chr\tkinModCall\tm6A\t5\t5\t200\t+\t.\tmotif=RAATTY;context={pos3};identificationQv=120\n"
    )

    summary, sites = import_pacbio_gff(
        gff,
        fasta,
        "TOY",
        min_identification_qv=80,
    )

    assert list(summary.columns) == [
        "sample_id", "platform", "motif", "modification", "called_sites",
        "genomic_sites", "methylated_fraction", "median_score", "source_file",
        "qc_status",
    ]
    raatty = summary.loc[summary.motif == "RAATTY"].iloc[0]
    assert raatty.called_sites == 2
    assert raatty.genomic_sites == 4
    assert math.isclose(raatty.methylated_fraction, 0.5)

    # All five site assignments are retained for audit/QC even though only two
    # pass both cognate-position and identification-QV filters.
    assert len(sites) == 5
    assert set(sites.cognate_position.dropna().astype(int)) == {3}
    assert sites.included_in_summary.sum() == 2
    assert (~sites.passes_motif_position).sum() == 2
    assert (~sites.passes_identification_qv).sum() == 3
    assert set(sites.cognate_position_source) == {
        "highest_median_identification_qv"
    }


def test_import_without_context_falls_back_to_legacy_counting(tmp_path):
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(">chr\nGAATTCAAATTT\n")
    gff = tmp_path / "toy.gff"
    gff.write_text(
        "chr\tkinModCall\tm6A\t2\t2\t42\t+\t.\tmotif=RAATTY,2;coverage=55\n"
        "chr\tkinModCall\tm6A\t8\t8\t44\t+\t.\tmotif=RAATTY,2;coverage=60\n"
    )
    summary, sites = import_pacbio_gff(gff, fasta, "TOY")
    raatty = summary.loc[summary.motif == "RAATTY"].iloc[0]
    assert raatty.called_sites == 2
    assert raatty.genomic_sites == 2
    assert sites.cognate_position.isna().all()
    assert sites.included_in_summary.all()


def test_import_fails_loudly_without_motif_attribute(tmp_path):
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(">chr\nGAATTC\n")
    gff = tmp_path / "toy.gff"
    gff.write_text(
        "chr\tkinModCall\tm6A\t2\t2\t42\t+\t.\t"
        "coverage=55;context=NNGAATTCNN\n"
    )
    with pytest.raises(
        ValueError, match="Could not identify an explicit motif attribute"
    ):
        import_pacbio_gff(gff, fasta, "TOY")
