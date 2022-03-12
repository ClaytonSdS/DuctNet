from importlib import reload
import pandas as pd
import numpy as np

from DuctNet import DuctNet as dn
reload(dn)

file_path = "C:\\Users\\clayt\\AppData\\Local\\Programs\\Python\\Python39\\Lib\\site-packages\\DuctNet\data"

net = dn.Net(f'{file_path}\\node.inp',
            f'{file_path}\\link.inp',
             f'{file_path}\\connect.inp')



def fun(x):
    x.Set_Dynamic()
    return x.DP

df = pd.DataFrame({'DP': [fun(lk[i]) for i in range(len(lk))]})
df['cum'] = df.cumsum()
DP_total = df.tail(1)['cum'].values
p = [DP_total]+[DP_total-df.loc[i,'cum'] for i in df.index]

b = np.zeros(len(lk)-1)
              
