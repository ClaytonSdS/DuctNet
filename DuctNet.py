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
        self.node = pd.read_csv(fnode, sep='|', index_col='idx')
        self.dflink = pd.read_csv(flink, sep='|', index_col='idx')
        self.connect = pd.read_csv(fconnect, sep='|', index_col='lk')
        self.DictLink()

        # NÓS DE CONTORNO E NÓS DE BALANÇO DE MASSA (VARIAVEIS)
        self.node_variable = list(set(self.connect['from'].values).intersection(self.connect['to'].values))
        self.node_boundary = list(set(self.connect['from'].values).symmetric_difference(set(self.connect['to'].values)))


        self.Set_Mass_Equations()

        # COLOCAR OS NOS ASSOCIADOS AS LIGACOES
        self.nodes_in_links = {self.connect.index[node]:  list(self.connect.loc[self.connect.index[node]]) for node in range(len(self.connect.index))}

        # COLOCAR AS LIGACOES ASSOCIADAS AOS NÓS VARIAVEIS
        self.links_in_nodes = {self.node_variable[node]: self.Set_In_Out(self.node_variable[node]) for node  in range(len(self.node_variable))}

        # ATIVAR FUNÇÃO PARA ENCONTRAR OS VIZINHOS DOS NÓS VARIAVEIS
        self.Set_mass_balance()

        # CRIAR AS MATRIX DO SOLVER
        self.matrix_a = self.mass_balance[:].loc[self.node_variable].values
        self.Add_Pressure_in_Nodes_Boundary()   # ADICIONAR PRESSAO AOS TERMOS DE CONTORNO
        self.matrix_b = self.mass_balance[:].loc[self.node_boundary].cumsum().values[-1] * -1
        self.t = np.linalg.solve(self.matrix_a,self.matrix_b)
        self.result = {self.node_variable[i]:np.linalg.solve(self.matrix_a,self.matrix_b)[i] for i in range(len(self.node_variable))}

    def DictLink(self):
        self.link = {ID: Link(self.dflink.loc[ID, 'm'],
                              self.dflink.loc[ID, 'A_r'],
                              self.dflink.loc[ID, 'A_i'],
                              self.dflink.loc[ID, 'A_o'],
                              self.dflink.loc[ID, 'rho'],
                              self.dflink.loc[ID, 'zeta'])
                     for ID in self.dflink.index}

    # SEPARAR FLUXOS CHEGANDO E SAINDO DE NÓS (IN, OUT)
    def Set_In_Out(self, no):
        # se exisir negativo em in -> colocar em out (em modulo)
        # se existir negativo em out -> colocar em in (em modulo)

        return {'in': list(self.connect.loc[self.connect['to']==no].index), 'out': list(self.connect.loc[self.connect['from']==no].index)}

    # FUNÇÃO PARA CALCULO DOS BALANÇO DE MASSAS NOS NÓS VARIAVEIS
    def Set_mass_balance(self):
        nodes = self.node_variable + self.node_boundary
        nodes.sort()
        self.mass_balance = pd.DataFrame({"node": nodes}).set_index("node")
        #in_ = sum([net.equations[net.links_in_nodes[2]["in"][i]].values for i in range(len(net.links_in_nodes[2]["in"]))])
        #out_ = sum([net.equations[net.links_in_nodes[2]["out"][i]].values for i in range(len(net.links_in_nodes[2]["out"]))])
        for node in range(len(self.node_variable)):
            mass_entering = [self.equations[self.links_in_nodes[self.node_variable[node]]["in"][i]].values for i in range(len(self.links_in_nodes[self.node_variable[node]]["in"]))]
            mass_leaving =
            self.mass_balance[self.node_variable[node]] = sum() - sum()

    # ADICIONAR OS TERMOS DE PRESSÃO AOS TERMOS CONSTANTES
    def Add_Pressure_in_Nodes_Boundary(self):
        for n_contorno in range(len(self.node_boundary)):
            self.mass_balance.loc[self.mass_balance.index[self.node_boundary[n_contorno]]] = self.mass_balance.loc[self.mass_balance.index[self.node_boundary[n_contorno]]].values * self.node.loc[self.node.index[self.node_boundary[n_contorno]]].values

    # FUNÇÃO PARA ASSOCIAR OS TERMOS a_n COM AS PRESSOES p_saida E p_from
    def Set_Mass_Equations(self):
        nodes = self.node_variable+self.node_boundary
        nodes.sort()
        links = {self.connect.index[i]:[] for i in range(len(self.connect.index))}
        self.equations = pd.DataFrame({"node": nodes}).set_index("node")

        # CRIAR COLUNAS COM OS NOMES DAS LIGAÇÕES E ATRIBUIR UM ARRAY DE ZEROS
        for link_name in range(len(self.connect.index)):
            self.equations[self.connect.index[link_name]] = np.zeros(len(nodes))

        # COLOCAR a_n E -a_n NAS LIGAÇÕES DE ACORDO COM AS PRESSOES ENTRADA E SAIDA (SENTIDO DO ESCOAMENTO)
        for link in range(len(self.connect.index.values)):
            # se exisir negativo em in -> colocar em out (em modulo)
            # se existir negativo em out -> colocar em in (em modulo)

            self.equations[self.connect.index.values[link]][self.connect.loc[self.connect.index.values[link]].values[0]] = self.link[self.connect.index.values[link]].a
            self.equations[self.connect.index.values[link]][self.connect.loc[self.connect.index.values[link]].values[1]] = - self.link[self.connect.index.values[link]].a