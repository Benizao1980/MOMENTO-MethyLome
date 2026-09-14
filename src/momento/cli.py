import argparse
from pathlib import Path
import pandas as pd
from .schema import REQUIRED_METHYLATION_COLUMNS,missing_columns
from .qc import motif_prevalence,flag_samples_by_core_motif
from .matrix import build_matrix
from .cluster import pca

def read_table(path): return pd.read_csv(path,sep="\t")
def cmd_validate(a):
    df=read_table(a.input); m=missing_columns(df.columns,REQUIRED_METHYLATION_COLUMNS)
    if m: raise SystemExit("Missing required columns: "+", ".join(m))
    print(f"PASS: {len(df):,} rows; {df.sample_id.nunique()} samples; {df.motif.nunique()} motifs")
def cmd_matrix(a):
    df=read_table(a.input); ex=[x for x in a.exclude.split(",") if x] if a.exclude else None
    mat=build_matrix(df,a.value,a.min_called_sites,ex); Path(a.output).parent.mkdir(parents=True,exist_ok=True); mat.to_csv(a.output,sep="\t"); print(a.output)
def cmd_qc(a):
    df=read_table(a.input); motif_prevalence(df,a.min_called_sites).to_csv(a.prevalence_output,sep="\t",index=False)
    if a.core_motif: flag_samples_by_core_motif(df,a.core_motif).to_csv(a.core_output,sep="\t",index=False)
def cmd_cluster(a):
    mat=pd.read_csv(a.matrix,sep="\t",index_col=0); pca(mat).to_csv(a.output,sep="\t")
def build_parser():
    p=argparse.ArgumentParser(prog="momento",description="MOMENTO: comparative bacterial methylomics"); s=p.add_subparsers(dest="command",required=True)
    v=s.add_parser("validate"); v.add_argument("--input",required=True); v.set_defaults(func=cmd_validate)
    m=s.add_parser("matrix"); m.add_argument("--input",required=True); m.add_argument("--output",required=True); m.add_argument("--value",choices=["presence","called_sites","methylated_fraction"],default="presence"); m.add_argument("--min-called-sites",type=int,default=1); m.add_argument("--exclude",default=""); m.set_defaults(func=cmd_matrix)
    q=s.add_parser("qc"); q.add_argument("--input",required=True); q.add_argument("--prevalence-output",required=True); q.add_argument("--min-called-sites",type=int,default=1); q.add_argument("--core-motif"); q.add_argument("--core-output",default="core_motif_qc.tsv"); q.set_defaults(func=cmd_qc)
    c=s.add_parser("cluster"); c.add_argument("--matrix",required=True); c.add_argument("--output",required=True); c.set_defaults(func=cmd_cluster); return p
def main():
    a=build_parser().parse_args(); a.func(a)
if __name__=="__main__": main()
