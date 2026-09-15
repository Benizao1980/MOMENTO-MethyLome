from types import SimpleNamespace

import pandas as pd
import pytest

import momento.cli as cli


def test_assess_import_status_flags_tiny_imports():
    summary = pd.DataFrame([{"motif": "RAATTY"}])
    tiny_sites = pd.DataFrame({"motif": ["RAATTY"] * 6})
    normal_sites = pd.DataFrame({"motif": ["RAATTY"] * 100})

    assert cli.assess_import_status(summary, tiny_sites, 100) == "LOW_DATA"
    assert cli.assess_import_status(summary, normal_sites, 100) == "PASS"


def test_import_pacbio_cli_reports_clean_failure(monkeypatch):
    def fail_import(*args, **kwargs):
        raise ValueError("No feature records found in broken.gff")

    monkeypatch.setattr(cli, "import_pacbio_gff", fail_import)
    args = SimpleNamespace(
        gff="broken.gff",
        fasta="genome.fasta",
        sample="BROKEN",
        platform="pacbio",
        min_score=None,
        min_identification_qv=None,
        motif_attribute=None,
        modification_attribute=None,
        keep_all_motif_positions=False,
    )

    with pytest.raises(SystemExit, match="FAIL: BROKEN: No feature records found"):
        cli.cmd_import_pacbio(args)
