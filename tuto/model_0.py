from importlib import reload
import pandas as pd
import numpy as np

from DuctNet import DuctNet as dn
reload(dn)

file_path = "C:\\Users\\clayt\\AppData\\Local\\Programs\\Python\\Python39\\Lib\\site-packages\\DuctNet\data"

net = dn.Net(f'{file_path}\\node.inp',
            f'{file_path}\\link.inp',
             f'{file_path}\\connect.inp')


net.Start_Iteration(300)
