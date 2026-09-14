import pandas as pd
from momento.matrix import build_matrix
def test_presence():
 d=pd.DataFrame({"sample_id":["A","A","B","B"],"motif":["X","Y","X","Y"],"called_sites":[10,0,2,8],"methylated_fraction":[.9,0,.2,.8]}); m=build_matrix(d,"presence",5); assert m.loc["A","X"]==1 and m.loc["B","X"]==0
