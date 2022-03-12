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



lk = {0: dn.Link(5, 0.002, 0.002, 0.002, 1000, 0.5),
      1:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.1),
      2:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.3),
      3:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.7),
      5: dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.7)}

no = {0: dn.Node(5000),
      1: dn.Node(3400),
      2: dn.Node(3150),
      3: dn.Node(2170),
      4: dn.Node(10)}

lk = {0: dn.Link(5, 0.002, 0.002, 0.002, 1000, 0.5),
      1:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.1),
      2:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.3),
      3:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.7),
      4:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.7),
      5: dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.7),
      6:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.7),
      7:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.7),
      8:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.7),
      9:dn.Link(5, 0.002, 0.002,  0.002, 1000, 0.7),
      }

no = {0: dn.Node(5000),
      1: dn.Node(3400),
      2: dn.Node(3150),
      3: dn.Node(2170),
      4: dn.Node(10),
      5: dn.Node(10),
      6: dn.Node(10),
      7: dn.Node(10),
      8: dn.Node(10),
      9: dn.Node(10),
      10: dn.Node(10)}

# caso 1 - #######################################################
net = dn.Net(connectivity= np.array([[0,1,2,3,5],
                                   [0,1,2,3,5],
                                   [1,2,3,4,2]]), ligacoes=lk, nos=no)
# caso 2 - #######################################################
net = dn.Net(connectivity= np.array([[0,1,2,3,4,5,6,7,8,9],
                                   [0,1,2,3,5,6,7,8,9,10],
                                   [1,2,3,4,1,2,3,1,2,3]]), ligacoes=lk, nos=no)
#caso 3 serie - ##################################################
net = dn.Net(connectivity= np.array([[0,1,2,3],
                                   [0,1,2,3],
                                   [1,2,3,4]]), ligacoes=lk, nos=no)



def fun(x):
    x.Set_Dynamic()
    return x.DP

df = pd.DataFrame({'DP': [fun(lk[i]) for i in range(len(lk))]})
df['cum'] = df.cumsum()
DP_total = df.tail(1)['cum'].values
p = [DP_total]+[DP_total-df.loc[i,'cum'] for i in df.index]

b = np.zeros(len(lk)-1)
              
