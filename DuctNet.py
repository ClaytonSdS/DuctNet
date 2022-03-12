import pandas as pd
import numpy as np

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline

from scipy.interpolate import interp1d as interp1d
from scipy.interpolate import interp2d as interp2d
from scipy.optimize import minimize_scalar as minimize_scalar
from scipy.optimize import minimize as minimize
from scipy.optimize import bisect as bisect


# from UnitConv import UnitConv as uc
class Link():
    def __init__(self, m, A_r, A_i, A_o, rho, zeta):
        self.m = m
        self.A_r = A_r
        self.A_i = A_i
        self.A_o = A_o
        self.rho = rho
        self.zeta = zeta
        self.Set_a()


    def Set_Dynamic(self):
        self.U = self.m / (self.A_r * self.rho)
        self.p_d = self.rho * self.U ** 2 / 2
        self.DP = self.zeta * self.p_d

    def Set_a(self):
        self.a = 1 / ((self.m / (2 * self.rho)) * (self.zeta / self.A_r ** 2 + 1 / self.A_o ** 2 - 1 / self.A_i ** 2))


class Node():
    def __init__(self, p):
        self.p = p


class Net():
    def __init__(self, fnode, flink, fconnect):
        self.node = pd.read_csv(fnode, sep='\t', index_col='idx')
        self.dflink = pd.read_csv(flink, sep='\t', index_col='idx')
        self.connect = pd.read_csv(fconnect, sep='\t', index_col='lk')
        self.DictLink()

        #self.Set_Connection()

        self.mass_balance = {}
        self.matrix_a = {}
        self.matrix_b = {}
        self.result = {}

    def DictLink(self):
        self.link = {ID: Link(self.dflink.loc[ID, 'm'],
                              self.dflink.loc[ID, 'A_r'],
                              self.dflink.loc[ID, 'A_i'],
                              self.dflink.loc[ID, 'A_o'],
                              self.dflink.loc[ID, 'rho'],
                              self.dflink.loc[ID, 'zeta'])
                     for ID in self.dflink.index}

    def Set_Connection(self):
        pass
        #self.node_variable = list(set(self.from_).intersection(self.to_))
        #self.node_boundary = list(set(self.from_).symmetric_difference(set(self.to_)))

        #self.nodes = list(self.node_variable + self.node_boundary)
        #self.nodes.sort()

        # checar
        #self.nodes_links = {self.nodes[i]: [] for i in range(0, len(self.nodes))}
        #self.links_nodes = {f"L_{self.link[i]}": [self.from_[i], self.to_[i]] for i in range(len(self.link))}
        #self.nodes_near = {}

    # m_line = {0: a0(p0) - a0(p1); 1: a1(p1) - a1(p2); 2: a2(p2 + p5) - a2(p3); 3: a3(p3) - a3(p4); 5:a5(p5) - a5(p2)
    # {2: (a1 + a5) - a3}
    # preciso saber o que entra e o que sai do nó
    # m_entra - m_sai = 0
    # if to -> chegando
    # if from -> saindo

