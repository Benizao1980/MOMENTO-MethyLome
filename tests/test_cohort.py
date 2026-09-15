from pathlib import Path

import pandas as pd
import pytest

from momento.cohort import build_cohort, context_preflight, load_manifest


def _write_context_gff(path: Path, seqid: str, context: str, position: int = 21):
    path.write_text(
        "##gff-version 3\n"
        f"{seqid}\tkinModCall\tm6A\t{position}\t{position}\t100\t+\t.\t"
        f"motif=RAATTY;context={context};identificationQv=100\n"
    )


def test_context_preflight_exact_match(tmp_path):
    sequence = "A" * 20 + "G" + "C" * 20
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(f">chr\n{sequence}\n")
    gff = tmp_path / "calls.gff"
    _write_context_gff(gff, "chr", sequence)

    result = context_preflight(gff, fasta)
    assert result["preflight_status"] == "PASS"
    assert result["context_usable"] == 1
    assert result["context_compared"] == 1
    assert result["context_exact"] == 1
    assert result["context_coverage"] == 1.0
    assert result["context_concordance"] == 1.0


def test_context_preflight_rejects_wrong_assembly(tmp_path):
    context = "A" * 20 + "G" + "C" * 20
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">chr\n" + "T" * 41 + "\n")
    gff = tmp_path / "calls.gff"
    _write_context_gff(gff, "chr", context)

    result = context_preflight(gff, fasta)
    assert result["preflight_status"] == "FAIL"
    assert result["context_coverage"] == 1.0
    assert result["context_concordance"] == 0.0


def test_context_preflight_preserves_unknown_assembly_symbols(tmp_path):
    sequence = "A" * 20 + "?" + "C" * 20
    expected_reverse_context = "G" * 20 + "?" + "T" * 20
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(f">chr\n{sequence}\n")
    gff = tmp_path / "calls.gff"
    gff.write_text(
        "##gff-version 3\n"
        f"chr\tkinModCall\tm6A\t21\t21\t100\t-\t.\t"
        f"motif=RAATTY;context={expected_reverse_context};identificationQv=100\n"
    )

    result = context_preflight(gff, fasta)
    assert result["preflight_status"] == "PASS"
    assert result["context_compared"] == 1
    assert result["context_exact"] == 1


def test_context_preflight_rejects_duplicate_fasta_ids(tmp_path):
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">chr\nAAAA\n>chr\nTTTT\n")
    gff = tmp_path / "calls.gff"
    gff.write_text("##gff-version 3\n")

    with pytest.raises(ValueError, match="Duplicate FASTA record ID"):
        context_preflight(gff, fasta)


def test_load_manifest_rejects_duplicate_sample_ids(tmp_path):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample_id\tfasta\tgff\n"
        "S1\ta.fasta\ta.gff\n"
        "S1\tb.fasta\tb.gff\n"
    )
    with pytest.raises(ValueError, match="duplicate sample_id"):
        load_manifest(manifest)


def test_build_cohort_continues_after_sample_failure(tmp_path):
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(">chr\nGAATTCAAATTT\n")

    good_gff = tmp_path / "good.gff"
    good_gff.write_text(
        "chr\tkinModCall\tm6A\t2\t2\t42\t+\t.\tmotif=RAATTY,2;coverage=55\n"
        "chr\tkinModCall\tm6A\t8\t8\t44\t+\t.\tmotif=RAATTY,2;coverage=60\n"
    )
    bad_gff = tmp_path / "bad.gff"
    bad_gff.write_text("##gff-version 3\n")

    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample_id\tfasta\tgff\n"
        "GOOD\ttoy.fasta\tgood.gff\n"
        "BAD\ttoy.fasta\tbad.gff\n"
    )

    outdir = tmp_path / "results"
    combined, audit = build_cohort(manifest, outdir, low_data_threshold=1)

    assert set(combined["sample_id"]) == {"GOOD"}
    status = dict(zip(audit["sample_id"], audit["qc_status"]))
    assert status["GOOD"] == "PASS"
    assert status["BAD"] == "FAIL"
    assert (outdir / "cohort.momento.tsv").exists()
    assert (outdir / "cohort.qc.tsv").exists()
    assert (outdir / "samples" / "GOOD" / "GOOD.momento.tsv").exists()

    audit_on_disk = pd.read_csv(outdir / "cohort.qc.tsv", sep="\t")
    assert set(audit_on_disk["sample_id"]) == {"GOOD", "BAD"}
