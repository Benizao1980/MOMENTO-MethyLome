import math

import pytest

from momento.io.pacbio import (
    canonicalise_motif,
    canonicalise_motif_values,
    count_genomic_opportunities,
    import_pacbio_gff,
    parse_attributes,
    reverse_complement_iupac,
)


def test_parse_attributes_gff3():
    attrs = parse_attributes("motif=RAATTY;coverage=55;modification=m6A")
    assert attrs["motif"] == "RAATTY"
    assert attrs["coverage"] == "55"
    assert attrs["modification"] == "m6A"


def test_iupac_and_reverse_complement_counting():
    seqs = {"chr": "GAATTCAAATTTGATC", "plasmid": "GAATTC"}
    assert canonicalise_motif("RAATTY") == "RAATTY"
    assert reverse_complement_iupac("RAATTY") == "RAATTY"
    assert count_genomic_opportunities(seqs, "RAATTY") == 3
    assert count_genomic_opportunities(seqs, "GATC") == 1


def test_old_pacbio_numeric_motif_position_suffix():
    assert canonicalise_motif("RAATTY,2") == "RAATTY"
    assert canonicalise_motif("CATG,1") == "CATG"
    assert canonicalise_motif("AGTNNNNNNRTTG,0/CAAYNNNNNNACT,12") == "AGTNNNNNNRTTG/CAAYNNNNNNACT"
    with pytest.raises(ValueError, match="Unsupported comma-delimited motif value"):
        canonicalise_motif("RAATTY,unexpected")


def test_overlapping_comma_separated_motif_values():
    assert canonicalise_motif_values("GCCAA,RAATTY") == ["GCCAA", "RAATTY"]
    assert canonicalise_motif_values("GCCAA,RAATTY,2") == ["GCCAA", "RAATTY"]
    assert canonicalise_motif_values("RAATTY,2") == ["RAATTY"]


def test_import_pacbio_gff(tmp_path):
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(">chr\nGAATTCAAATTTGATCGATCGCCAA\n")
    gff = tmp_path / "toy.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chr\tkinModCall\tm6A\t2\t2\t42\t+\t.\tmotif=RAATTY,2;coverage=55\n"
        "chr\tkinModCall\tm6A\t8\t8\t44\t+\t.\tmotif=RAATTY,2;coverage=60\n"
        "chr\tkinModCall\tm6A\t14\t14\t50\t+\t.\tmotif=GATC,1;coverage=62\n"
        "chr\tkinModCall\tm6A\t22\t22\t48\t+\t.\tmotif=GCCAA,RAATTY;coverage=58\n"
    )
    summary, sites = import_pacbio_gff(gff, fasta, "TOY")
    assert list(summary.columns) == [
        "sample_id", "platform", "motif", "modification", "called_sites",
        "genomic_sites", "methylated_fraction", "median_score", "source_file",
        "qc_status",
    ]
    raatty = summary.loc[summary.motif == "RAATTY"].iloc[0]
    assert raatty.called_sites == 3
    assert raatty.genomic_sites == 2
    assert math.isclose(raatty.methylated_fraction, 1.5)
    gccaa = summary.loc[summary.motif == "GCCAA"].iloc[0]
    assert gccaa.called_sites == 1
    assert gccaa.genomic_sites == 1
    # Four GFF records, but the overlapping GCCAA/RAATTY record is represented
    # twice in the site-by-motif table.
    assert len(sites) == 5
    assert sites.iloc[0].attributes_raw.startswith("motif=RAATTY,2")


def test_import_fails_loudly_without_motif_attribute(tmp_path):
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(">chr\nGAATTC\n")
    gff = tmp_path / "toy.gff"
    gff.write_text("chr\tkinModCall\tm6A\t2\t2\t42\t+\t.\tcoverage=55;context=NNGAATTCNN\n")
    with pytest.raises(ValueError, match="Could not identify an explicit motif attribute"):
        import_pacbio_gff(gff, fasta, "TOY")
