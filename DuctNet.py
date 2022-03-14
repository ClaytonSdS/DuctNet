import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import re

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

    class Tubo(Link):
        def __init__(self, m, A_r, A_i, A_o, rho, zeta, l):
            Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
            self.mu =  1.12e-3
            self.l = l
            self.D = np.sqrt((4 * self.A_r) / np.pi)
            self.U = self.m / (self.A_r * self.rho)
            self.nu = self.mu / self.rho
            self.df = pd.read_excel(
                r'tubo_chartB.xlsx',
                engine='openpyxl').set_index('Reynolds')
            self.Set_Zeta()

        def lambda_(self, Re):
            interp = interp1d(self.df.index, self.df['lambda'].values)
            if Re <= 2000:
                return 64 / Re
            if (Re > 2000 and Re <= 4000):
                return float(interp(Re))
            if (Re > 4000 and Re < 100000):
                return 0.3164 / (Re ** 0.25)
            if Re >= 100000:
                return 1 / ((1.8 * np.log10(Re) - 1.64) ** 2)

        def Set_Zeta(self):
            Re = round(self.U * self.D / self.nu)
            lambdaValue = float(self.lambda_(Re))
            self.zeta = lambdaValue * (self.l / self.D)

    class Difusor(Link):
        def __init__(self, m, A_r, A_i, A_o, rho, zeta, l):
            Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
            self.mu =  1.12e-3
            self.l = l
            self.D = np.sqrt((4 * self.A_r) / np.pi)
            self.U = self.m / (self.A_r * self.rho)
            self.nu = self.mu / self.rho
            df_z = {2: pd.read_excel(
                r'C:\Users\clayt\OneDrive\Área de Trabalho\IC\Diagramas\Tabelas Excel\diagrama 5.2\difusor_zetad_Nair2.xlsx',
                skiprows=1, engine='openpyxl').set_index('ang'),  #
                    4: pd.read_excel(
                        r'C:\Users\clayt\OneDrive\Área de Trabalho\IC\Diagramas\Tabelas Excel\diagrama 5.2\difusor_zetad_Nair4.xlsx',
                        skiprows=1, engine='openpyxl').set_index('ang'),  #
                    6: pd.read_excel(
                        r'C:\Users\clayt\OneDrive\Área de Trabalho\IC\Diagramas\Tabelas Excel\diagrama 5.2\difusor_zetad_Nair6.xlsx',
                        skiprows=1, engine='openpyxl').set_index('ang'),  #
                    10: pd.read_excel(
                        r'C:\Users\clayt\OneDrive\Área de Trabalho\IC\Diagramas\Tabelas Excel\diagrama 5.2\difusor_zetad_Nair10.xlsx',
                        skiprows=1, engine='openpyxl').set_index('ang'),  #
                    16: pd.read_excel(
                        r'C:\Users\clayt\OneDrive\Área de Trabalho\IC\Diagramas\Tabelas Excel\diagrama 5.2\difusor_zetad_Nair16.xlsx',
                        skiprows=1, engine='openpyxl').set_index('ang')}

            interpz = {i: interp2d(df_z[i].index, df_z[i].columns, np.array(df_z[i]).T) for i in [2, 4, 6, 10, 16]}

            self.Set_Zeta()



        def Set_Zeta(self):
            self.zeta = lambdaValue * (self.l / self.D)

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
        #self.t = np.linalg.solve(self.matrix_a,self.matrix_b)
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
        return {'in': list(self.connect.loc[self.connect['to']==no].index), 'out': list(self.connect.loc[self.connect['from']==no].index)}

    # FUNÇÃO PARA CALCULO DOS BALANÇO DE MASSAS NOS NÓS VARIAVEIS
    def Set_mass_balance(self):
        nodes = self.node_variable + self.node_boundary
        nodes.sort()
        self.mass_balance = pd.DataFrame({"node": nodes}).set_index("node")

        for node_index in range(len(self.node_variable)):
            _node_ = self.node_variable[node_index]
            _links_in_ = self.links_in_nodes[_node_]["in"]
            _links_out_ = self.links_in_nodes[_node_]["out"]

            mass_entering = [self.equations[_links_in_[i]].values for i in range(len(_links_in_))]
            mass_leaving = [self.equations[_links_out_[i]].values for i in range(len(_links_out_))]

            self.mass_balance[_node_] = sum(mass_entering) - sum(mass_leaving)

    # ADICIONAR OS TERMOS DE PRESSÃO AOS TERMOS CONSTANTES
    def Add_Pressure_in_Nodes_Boundary(self):
        for node_index in range(len(self.node_boundary)):
            _node_ = self.node_boundary[node_index]
            _pressure_node_ = self.node.loc[_node_].values
            _a_node_ = self.mass_balance.loc[_node_].values
            self.mass_balance.loc[_node_] = _a_node_ * _pressure_node_

    # FUNÇÃO PARA ASSOCIAR OS TERMOS a_n COM AS PRESSOES p_saida E p_from
    def Set_Mass_Equations(self):
        nodes = self.node_variable+self.node_boundary
        nodes.sort()
        self.equations = pd.DataFrame({"node": nodes}).set_index("node")

        # CRIAR COLUNAS COM OS NOMES DAS LIGAÇÕES E ATRIBUIR UM ARRAY DE ZEROS
        for link_name in range(len(self.connect.index)):
            self.equations[self.connect.index[link_name]] = np.zeros(len(nodes))

        # COLOCAR a_n E -a_n NAS LIGAÇÕES DE ACORDO COM AS PRESSOES ENTRADA E SAIDA (SENTIDO DO ESCOAMENTO)
        for link_index in range(len(self.connect.index.values)):
            _ligacao_ = self.connect.index.values[link_index]
            _todas_ligacoes_ = self.connect.index.values
            _find_ = self.connect.loc
            self.equations[_todas_ligacoes_[link_index]][_find_[_todas_ligacoes_[link_index]].values[0]] = self.link[_ligacao_].a
            self.equations[_todas_ligacoes_[link_index]][_find_[_todas_ligacoes_[link_index]].values[1]] = - self.link[_ligacao_].a

    def Verify_Flux_Direction(self, ligacao):
        pass


    def Identify_Device_Type(self):
        for index in range(len(self.connect.index.values)):
            device = self.connect.index.values[index]
            type_device = re.split("\d", device)[0]

            # DEVICE == TUBO
            if type_device == 'd':
                length = float(input(f"Digite o Valor do Comprimento do Tubo [m] {device}"))
                self.link[self.connect.index.values[index]] = Link.Tubo(self.link[device].m,
                                                                        self.link[device].A_r,
                                                                        self.link[device].A_i,
                                                                        self.link[device].A_o,
                                                                        self.link[device].rho,
                                                                        self.link[device].zeta,
                                                                        length)

    def Set_M_Line_Equations(self):
        nodes = self.node_variable + self.node_boundary
        nodes.sort()
        self.m_line = self.equations[:]
        for node_index in range(len(self.node.index)):

            _node_ = self.node.index[node_index]
            _pressure_node_ = self.node.loc[_node_].values
            _a_node_ = self.m_line.loc[_node_].values
            self.m_line.loc[_node_] = _a_node_ * _pressure_node_

    def Start_Iteration(self):
        self.alpha_p = 0.5
        self.alpha_m = 0.5
        self.Identify_Device_Type()

        for x in range(400):
            self.Set_M_Line_Equations()

            # ATUALIZAR OS VALORES DAS PRESSOES
            self.shape = self.node.loc[self.node_variable].values.shape[0]
            pressure_line = np.array(list(self.result.values()))
            pressure_dot = self.node.loc[self.node_variable].values.reshape(-1)
            delta_p = pressure_line - pressure_dot
            self.node.loc[self.node_variable] = (pressure_dot + self.alpha_p * delta_p).reshape((self.shape, 1))

            # ATUALIZAR OS VALORES DAS VAZOES MASSICAS
            ligacoes = self.connect.index.values
            self.mass_line = {ligacoes[index]: self.m_line.cumsum().values[-1][index] for index in range(len(ligacoes))}
            self.mass_dot = {ligacoes[index]: self.link[ligacoes[index]].m for index in range(len(ligacoes))}
            self.delta_m = {ligacoes[index]: self.mass_line[ligacoes[index]] - self.mass_dot[ligacoes[index]] for index
                            in range(len(ligacoes))}
            for link_index in range(len(self.connect.index.values)):
                _ligacao_ = self.connect.index.values[link_index]
                self.link[_ligacao_].m = self.link[_ligacao_].m + self.alpha_m * self.delta_m[_ligacao_]
                self.link[_ligacao_].Set_a()
            self.Set_Mass_Equations()
            self.Set_mass_balance()

            # CRIAR AS MATRIX DO SOLVER
            self.matrix_a = self.mass_balance[:].loc[self.node_variable].values
            self.Add_Pressure_in_Nodes_Boundary()  # ADICIONAR PRESSAO AOS TERMOS DE CONTORNO
            self.matrix_b = self.mass_balance[:].loc[self.node_boundary].cumsum().values[-1] * -1
            self.result = {self.node_variable[i]: np.linalg.solve(self.matrix_a, self.matrix_b)[i] for i in
                           range(len(self.node_variable))}


