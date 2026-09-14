"""PacBio importer scaffold. Real SMRT Analysis exports vary by version; parser validation against representative GFFs is required."""
import pandas as pd

def standardise_summary_table(path,sample_id_col="sample_id",motif_col="motif",count_col="called_sites"):
    df=pd.read_csv(path,sep=None,engine="python")
    return pd.DataFrame({
      "sample_id":df[sample_id_col].astype(str),"platform":"pacbio","motif":df[motif_col].astype(str),"modification":df.get("modification","unknown"),
      "called_sites":pd.to_numeric(df[count_col],errors="coerce"),"genomic_sites":pd.to_numeric(df.get("genomic_sites"),errors="coerce"),
      "methylated_fraction":pd.to_numeric(df.get("methylated_fraction"),errors="coerce"),"median_score":pd.to_numeric(df.get("median_score"),errors="coerce"),
      "source_file":str(path),"qc_status":"UNASSESSED"})
