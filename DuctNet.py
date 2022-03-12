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
    def __init__(self, connectivity, ligacoes, nos):
        self.ligacoes = ligacoes
        self.nos = nos

        self.link = connectivity[0]
        self.from_ = connectivity[1]
        self.to_ = connectivity[2]
        self.Set_Connection()

        self.mass_balance = {}
        self.matrix_a = {}
        self.matrix_b = {}
        self.result = {}

        # checar
        for node in range(len(self.nodes)):
            self.find_links_nodes(self.nodes[node])

        self.Set_Nodes_Neighbors()

        for y in range(len(self.node_variable)):
            self.Set_Mass_Balance_Equations(self.node_variable[y])

        for z in range(len(self.node_variable)):
            self.matrix_a[self.node_variable[z]] = self.Set_Matrix_A_Coefs(self.node_variable[z])
            self.matrix_b[self.node_variable[z]] = self.Set_Matrix_B_Coefs(self.node_variable[z])

        self.A = np.array(list(self.matrix_a.values()))
        self.B = np.array(list(self.matrix_b.values()))
        self.C = np.linalg.solve(self.A, self.B)

        for j in range(len(self.node_variable)):
            self.result[self.node_variable[j]] = self.C[j]

    def Set_Connection(self):
        self.connection = pd.DataFrame({"link": self.link, "from": self.from_, "to": self.to_})


        self.node_variable = list(set(self.from_).intersection(self.to_))
        self.node_boundary = list(set(self.from_).symmetric_difference(set(self.to_)))

        self.nodes = list(self.node_variable + self.node_boundary)
        self.nodes.sort()

        # checar
        self.nodes_links = {self.nodes[i]: [] for i in range(0, len(self.nodes))}
        self.links_nodes = {f"L_{self.link[i]}": [self.from_[i], self.to_[i]] for i in range(len(self.link))}
        self.nodes_near = {}

    def find_links_nodes(self, node_to_search):
        for x in range(len(self.to_)):
            if self.from_[x] == node_to_search or self.to_[x] == node_to_search:
                self.nodes_links[node_to_search].append(f"L_{self.link[x]}")

    def Set_Mass_Balance_Equations(self, node):
        self.mass_balance[node] = {node: [- self.ligacoes[int(self.nodes_links[node][i].split('_')[-1])].a for i in
                                          range(len(self.nodes_links[node]))]}
        self.mass_balance[node][node] = sum(self.mass_balance[node][node])
        for x in range(len(self.nodes_near[node])):
            self.mass_balance[node][self.nodes_near[node][x]] = self.ligacoes[
                int(self.nodes_links[node][x].split('_')[-1])].a

    def Set_Nodes_Neighbors(self):
        self.neighbors_by_link = {}
        for x in range(len(self.node_variable)):
            self.neighbors_by_link[self.node_variable[x]] = []
        for x in range(len(self.node_variable)):
            for j in range(len(self.nodes_links[self.node_variable[x]])):
                self.nodes_links[self.node_variable[x]][j]
                self.neighbors_by_link[self.node_variable[x]].append(self.nodes_links[self.node_variable[x]][j])
        for x in range(len(self.node_variable)):
            for i in range(len(self.neighbors_by_link[self.node_variable[x]])):
                self.neighbors_by_link[self.node_variable[x]][i] = self.links_nodes[
                    self.neighbors_by_link[self.node_variable[x]][i]]
        for x in range(len(self.node_variable)):
            self.nodes_near[self.node_variable[x]] = []
            for i in range(len(self.neighbors_by_link[self.node_variable[x]])):
                for z in range(len(self.neighbors_by_link[self.node_variable[x]][i])):
                    if self.neighbors_by_link[self.node_variable[x]][i][z] != self.node_variable[x]:
                        self.nodes_near[self.node_variable[x]].append(
                            self.neighbors_by_link[self.node_variable[x]][i][z])

    def Set_Matrix_A_Coefs(self, no):
        encontrados = []
        for x in range(len(self.node_variable)):
            if no in self.mass_balance[self.node_variable[x]].keys():
                encontrados.append(self.mass_balance[self.node_variable[x]][no])
            else:
                encontrados.append(0)
        return encontrados

    def Set_Matrix_B_Coefs(self, no):
        encontrados = []
        nos_boundary = self.node_boundary
        nos = list(self.mass_balance[no].keys())
        for x in range(len(nos)):
            if nos[x] in nos_boundary:
                encontrados.append(self.mass_balance[no][nos[x]] * self.nos[int(nos[x])].p)
        if len(encontrados) == 0:
            encontrados = [0]
        return encontrados

    # m_line = {0: a0(p0) - a0(p1); 1: a1(p1) - a1(p2); 2: a2(p2 + p5) - a2(p3); 3: a3(p3) - a3(p4); 5:a5(p5) - a5(p2)
    # {2: (a1 + a5) - a3}
    # preciso saber o que entra e o que sai do nó
    # m_entra - m_sai = 0
    # if to -> chegando
    # if from -> saindo

