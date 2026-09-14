from pathlib import Path
import pandas as pd
from momento.schema import REQUIRED_METHYLATION_COLUMNS,missing_columns
R=Path(__file__).resolve().parents[1]; df=pd.read_csv(R/"data/toy_methylation_long.tsv",sep="\t")
m=missing_columns(df.columns,REQUIRED_METHYLATION_COLUMNS)
if m: raise SystemExit(m)
print(f"PASS: {df.sample_id.nunique()} samples, {df.motif.nunique()} motifs")
