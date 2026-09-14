import re
IUPAC=set("ACGTRYSWKMBDHVN")

def canonicalise_motif(motif):
    motif=motif.strip().upper()
    for part in motif.split("/"):
        clean=re.sub(r"[^A-Z]","",part)
        if not clean or any(b not in IUPAC for b in clean):
            raise ValueError(f"Invalid motif: {motif!r}")
    return motif

def is_core_motif(prevalence, threshold=0.90):
    return prevalence>=threshold
