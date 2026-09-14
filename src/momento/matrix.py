import numpy as np
import pandas as pd

def build_matrix(df,value="presence",min_called_sites=1,exclude_motifs=None):
    work=df.copy()
    if exclude_motifs: work=work.loc[~work["motif"].isin(exclude_motifs)]
    if value=="presence": work["_value"]=(pd.to_numeric(work["called_sites"],errors="coerce").fillna(0)>=min_called_sites).astype(int)
    elif value=="called_sites": work["_value"]=pd.to_numeric(work["called_sites"],errors="coerce").fillna(0)
    elif value=="methylated_fraction": work["_value"]=pd.to_numeric(work["methylated_fraction"],errors="coerce")
    else: raise ValueError("value must be presence, called_sites, or methylated_fraction")
    mat=work.pivot_table(index="sample_id",columns="motif",values="_value",aggfunc="max")
    if value in {"presence","called_sites"}: mat=mat.fillna(0)
    return mat.sort_index().sort_index(axis=1)

def log1p_counts(matrix): return np.log1p(matrix)
