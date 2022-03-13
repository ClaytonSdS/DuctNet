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
        self.Set_Boundary_Variable_Nodes()

        self.matrix_a = {}
        self.matrix_b = {}
        self.result = {}
        self.Set_Equations_by_Node()

        # COLOCAR OS NOS ASSOCIADOS AS LIGACOES
        self.nodes_in_links = {self.connect.index[node]:  list(self.connect.loc[self.connect.index[node]]) for node in range(len(self.connect.index))}

        # COLOCAR AS LIGACOES ASSOCIADAS AOS NÓS VARIAVEIS
        self.links_in_nodes = {self.node_variable[node]: self.Set_from_to(self.node_variable[node]) for node  in range(len(self.node_variable))}

        # ATIVAR FUNÇÃO PARA ENCONTRAR OS VIZINHOS DOS NÓS VARIAVEIS
        #self.Set_Neighbors()
        self.Set_mass_balance()

    def DictLink(self):
        self.link = {ID: Link(self.dflink.loc[ID, 'm'],
                              self.dflink.loc[ID, 'A_r'],
                              self.dflink.loc[ID, 'A_i'],
                              self.dflink.loc[ID, 'A_o'],
                              self.dflink.loc[ID, 'rho'],
                              self.dflink.loc[ID, 'zeta'])
                     for ID in self.dflink.index}

    def Set_from_to(self, no):
        return {'in': list(self.connect.loc[self.connect['to']==no].index), 'out': list(self.connect.loc[self.connect['from']==no].index)}

    def Set_Boundary_Variable_Nodes(self):
        self.node_variable = list(set(self.connect['from'].values).intersection(self.connect['to'].values))
        self.node_boundary = list(set(self.connect['from'].values).symmetric_difference(set(self.connect['to'].values)))

    def find_different_node(self, lista, no):
        for x in range(len(lista)):
            if lista[x] != no:
                return lista[x]

    def Set_mass_balance(self):
        nodes = self.node_variable + self.node_boundary
        nodes.sort()

        self.mass_balance = pd.DataFrame({"node": nodes}).set_index("node")
        #in_ = sum([net.equations[net.links_in_nodes[2]["in"][i]].values for i in range(len(net.links_in_nodes[2]["in"]))])
        #out_ = sum([net.equations[net.links_in_nodes[2]["out"][i]].values for i in range(len(net.links_in_nodes[2]["out"]))])
        for node in range(len(self.node_variable)):
            self.mass_balance[self.node_variable[node]] = sum([self.equations[self.links_in_nodes[self.node_variable[node]]["in"][i]].values for i in range(len(self.links_in_nodes[self.node_variable[node]]["in"]))]) - sum([self.equations[self.links_in_nodes[self.node_variable[node]]["out"][i]].values for i in range(len(self.links_in_nodes[self.node_variable[node]]["out"]))])

    def Set_Neighbors(self):
        self.nodes_neighbors = {}
        for no in range(len(self.node_variable)):
            self.nodes_neighbors[self.node_variable[no]] = []
            for link in range(len(self.links_in_nodes[self.node_variable[no]])):
                self.nodes_neighbors[self.node_variable[no]].append(self.find_different_node(self.nodes_in_links[self.links_in_nodes[self.node_variable[no]][link]], self.node_variable[no]))

    def Set_Equations_by_Node(self):
        nodes = self.node_variable+self.node_boundary
        nodes.sort()
        links = {self.connect.index[i]:[] for i in range(len(self.connect.index))}
        self.equations = pd.DataFrame({"node": nodes}).set_index("node")
        for link in range(len(self.connect.index)):
            self.equations[self.connect.index[link]] = np.zeros(len(nodes))

        for link in range(len(self.connect.index.values)):
            self.equations[self.connect.index.values[link]][self.connect.loc[self.connect.index.values[link]].values[0]] = self.link[self.connect.index.values[link]].a
            self.equations[self.connect.index.values[link]][self.connect.loc[self.connect.index.values[link]].values[1]] = - self.link[self.connect.index.values[link]].a

        #net.connect.loc["d1"]
    # m_line = {0: a0(p0) - a0(p1); 1: a1(p1) - a1(p2); 2: a2(p2 + p5) - a2(p3); 3: a3(p3) - a3(p4); 5:a5(p5) - a5(p2)
    # {2: (a1 + a5) - a3}
    # preciso saber o que entra e o que sai do nó
    # m_entra - m_sai = 0
    # if to -> chegando
    # if from -> saindo

