import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances

def pca(matrix,n_components=2):
    n_components=min(n_components,matrix.shape[0],matrix.shape[1])
    model=PCA(n_components=n_components); coords=model.fit_transform(matrix.values)
    out=pd.DataFrame(coords,index=matrix.index,columns=[f"PC{i+1}" for i in range(n_components)])
    out.index.name="sample_id"; out.attrs["explained_variance_ratio"]=model.explained_variance_ratio_.tolist(); return out

def jaccard_distance(matrix):
    d=pairwise_distances(matrix.astype(bool).values,metric="jaccard")
    return pd.DataFrame(d,index=matrix.index,columns=matrix.index)
