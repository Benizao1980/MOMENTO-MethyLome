from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[1]; m=pd.read_csv(R/"metadata/toy_samples.tsv",sep="\t")
print(m[["sample_id","host_group","lineage"]]); print("Use grouped/nested CV; compare genome-only, methylome-only, and combined models.")
