import numpy as np
import pandas as pd

def add_fraction(df):
    out=df.copy()
    called=pd.to_numeric(out["called_sites"],errors="coerce")
    genomic=pd.to_numeric(out["genomic_sites"],errors="coerce")
    out["methylated_fraction"]=np.where(genomic>0,called/genomic,np.nan)
    return out

def motif_prevalence(df,min_called_sites=1):
    work=df.copy(); work["called_sites"]=pd.to_numeric(work["called_sites"],errors="coerce").fillna(0)
    total=work["sample_id"].nunique()
    present=(work.loc[work["called_sites"]>=min_called_sites].groupby("motif")["sample_id"].nunique().rename("n_samples").reset_index())
    present["total_samples"]=total
    present["prevalence"]=present["n_samples"]/total if total else np.nan
    return present.sort_values(["prevalence","n_samples"],ascending=False)

def flag_samples_by_core_motif(df,core_motif,min_fraction=0.90,min_called_sites_if_no_denominator=1000):
    work=df.loc[df["motif"]==core_motif].copy()
    work["called_sites"]=pd.to_numeric(work["called_sites"],errors="coerce")
    work["methylated_fraction"]=pd.to_numeric(work["methylated_fraction"],errors="coerce")
    def status(row):
        frac=row["methylated_fraction"]
        if pd.notna(frac): return "PASS" if frac>=min_fraction else "WARN"
        count=row["called_sites"]
        return "PASS" if pd.notna(count) and count>=min_called_sites_if_no_denominator else "WARN"
    work["core_motif_qc"]=work.apply(status,axis=1)
    return work[["sample_id","motif","called_sites","genomic_sites","methylated_fraction","core_motif_qc"]]
