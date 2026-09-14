import pytest
from momento.motifs import *
def test_motif(): assert canonicalise_motif(" raatty ")=="RAATTY"
def test_bad():
    with pytest.raises(ValueError): canonicalise_motif("ABCXYZ")
