from pathlib import Path
import pandas as pd
from momento.matrix import build_matrix
from momento.cluster import pca,jaccard_distance
R=Path(__file__).resolve().parents[1]; O=R/"results"; O.mkdir(exist_ok=True); df=pd.read_csv(R/"data/toy_methylation_long.tsv",sep="\t")
m=build_matrix(df,"presence",100,["RAATTY"]); m.to_csv(O/"motif_presence_no_core.tsv",sep="\t"); pca(m).to_csv(O/"pca_presence_no_core.tsv",sep="\t"); jaccard_distance(m).to_csv(O/"jaccard_distance.tsv",sep="\t")
