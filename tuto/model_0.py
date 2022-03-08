from importlib import reload
import pandas as pd
import numpy as np

from DuctNet import DuctNet as dn
reload(dn)

lk = [dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.5)]
lk += [dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.1)]
lk += [dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.3)]
lk += [dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.7)]

no = [dn.Node(5000)]
no += [dn.Node(3400)]
no += [dn.Node(3150)]
no += [dn.Node(2170)]
no += [dn.Node(10)]

net = dn.Net(no, lk)
net.net[0][1] = 0
net.net[1][2] = 1
net.net[2][3] = 2
net.net[3][4] = 3
net.net



def fun(x):
    x.Set_Dynamic()
    return x.DP

df = pd.DataFrame({'DP': [fun(lk[i]) for i in range(len(lk))]})
df['cum'] = df.cumsum()
DP_total = df.tail(1)['cum'].values
p = [DP_total]+[DP_total-df.loc[i,'cum'] for i in df.index]

b = np.zeros(len(lk)-1)
              
