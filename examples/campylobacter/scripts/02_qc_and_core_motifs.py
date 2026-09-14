from pathlib import Path
import pandas as pd
from momento.qc import motif_prevalence,flag_samples_by_core_motif
R=Path(__file__).resolve().parents[1]; O=R/"results"; O.mkdir(exist_ok=True); df=pd.read_csv(R/"data/toy_methylation_long.tsv",sep="\t")
motif_prevalence(df,100).to_csv(O/"motif_prevalence.tsv",sep="\t",index=False)
flag_samples_by_core_motif(df,"RAATTY").to_csv(O/"raatty_qc.tsv",sep="\t",index=False)
