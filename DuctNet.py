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



class Node():
    def __init__(self, p):
        self.p = p


class Net():
    def __init__(self, fnode, flink, fconnect):
        self.node = pd.read_csv(fnode, sep='|', index_col='idx')
        self.dflink = pd.read_csv(flink, sep='|', index_col='idx')
        self.connect = pd.read_csv(fconnect, sep='|', index_col='lk')
        self.DictLink()
        print('aaa')
        # NÓS DE CONTORNO E NÓS DE BALANÇO DE MASSA (VARIAVEIS)
        self.node_variable = list(set(self.connect['from'].values).intersection(self.connect['to'].values))
        self.node_boundary = list(set(self.connect['from'].values).symmetric_difference(set(self.connect['to'].values)))

        # ADICIONAR M_EXTRA SE LEN > 1:
        for x in range(len(self.node.index.values)):
            no = self.node.index.values[x]
            if len(self.connect.loc[self.connect['to']==no].index.values) > 1:
                self.dflink.loc[self.connect.loc[self.connect['to'] == no].index.values, "m_extra"] = self.dflink.loc[self.connect.loc[self.connect['to'] == no].index.values, "m"].sum()

        # COLOCAR OS NOS ASSOCIADOS AS LIGACOES
        #self.nodes_in_links = {self.connect.index[node]:  list(self.connect.loc[self.connect.index[node]]) for node in range(len(self.connect.index))}

        # COLOCAR AS LIGACOES ASSOCIADAS AOS NÓS VARIAVEIS
        self.links_in_nodes = {self.node_variable[node]: self.Set_In_Out(self.node_variable[node]) for node  in range(len(self.node_variable))}

        # GERAR AS EQUAÇÕES PARA CADA LIGAÇÃO
        self.Set_Mass_Equations()

        # GERAR OS BALANÇOS DE MASSA DOS NÓS VARIAVEIS
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

    # VALOR ABSOLUTO PARA VAZÃO MASSICA / INVERTE FROM e TO VALORES
    def Algorithm_0(self):
        for node_index in range(len(self.node_variable)):
            _node_ = self.node_variable[node_index]

            # dict with in and outs by nodes
            for in_index in range(len(self.links_in_nodes[_node_]["in"])):
                link_ = self.links_in_nodes[_node_]["in"][in_index]
                if self.link[link_].m < 0:
                    self.link[link_].m = abs(self.link[link_].m)
                    self.connect.loc[link_] = self.connect.loc[link_].values[::-1]
            for out_index in range(len(self.links_in_nodes[_node_]["out"])):
                link_ = self.links_in_nodes[_node_]["out"][out_index]
                if self.link[link_].m < 0:
                    self.link[link_].m = abs(self.link[link_].m)
                    self.connect.loc[link_] = self.connect.loc[link_].values[::-1]

    # TROCAR AS DIREÇÕES IN E OUT
    def Algorithm_1(self):
        for node_index in range(len(self.node_variable)):
            _node_ = self.node_variable[node_index]
            # dict with in and outs by nodes
            for in_index in range(len(self.links_in_nodes[_node_]["in"])):
                try:
                    link_ = self.links_in_nodes[_node_]["in"][in_index]
                    if self.link[link_].m < 0:
                        self.links_in_nodes[_node_]["out"].append(link_)
                        self.links_in_nodes[_node_]["in"].remove(link_)
                except IndexError: pass

            for out_index in range(len(self.links_in_nodes[_node_]["out"])):
                try:
                    link_ = self.links_in_nodes[_node_]["out"][out_index]
                    if self.link[link_].m < 0:
                        self.links_in_nodes[_node_]["in"].append(link_)
                        self.links_in_nodes[_node_]["out"].remove(link_)
                except IndexError:
                    pass

    # TROCAR AS DIREÇÕES IN E OUT / INVERTE FROM e TO VALORES
    def Algorithm_1_2(self):
        for node_index in range(len(self.node_variable)):
            _node_ = self.node_variable[node_index]
            # dict with in and outs by nodes
            for in_index in range(len(self.links_in_nodes[_node_]["in"])):
                try:
                    link_ = self.links_in_nodes[_node_]["in"][in_index]
                    if self.link[link_].m < 0:
                        self.links_in_nodes[_node_]["out"].append(link_)
                        self.links_in_nodes[_node_]["in"].remove(link_)
                        self.connect.loc[link_] = self.connect.loc[link_].values[::-1]

                except IndexError:
                    pass

            for out_index in range(len(self.links_in_nodes[_node_]["out"])):
                try:
                    link_ = self.links_in_nodes[_node_]["out"][out_index]
                    if self.link[link_].m < 0:
                        self.links_in_nodes[_node_]["in"].append(link_)
                        self.links_in_nodes[_node_]["out"].remove(link_)
                        self.connect.loc[link_] = self.connect.loc[link_].values[::-1]
                except IndexError:
                    pass


    def Identify_Device_Type(self):
        for index in range(len(self.connect.index.values)):
            device = self.connect.index.values[index]
            type_device = re.split("\d", device)[0]

            # DEVICE == TUBO
            if type_device == 'd':
                length = float(input(f"[TUBO] Digite o Valor do Comprimento do Tubo {device} [m]"))
                self.link[self.connect.index.values[index]] = Tubo(self.link[device].m,
                                                                        self.link[device].A_r,
                                                                        self.link[device].A_i,
                                                                        self.link[device].A_o,
                                                                        self.link[device].rho,
                                                                        self.link[device].zeta,
                                                                        length)
            # DEVICE == CONEXÃO EM T
            if type_device == 'ct':
                juncao = "Yes"
                referencia = "Yes"
                particao = int(input(f"[CONEXAO T - {device}] POSSUI PARTICAO?\n(1) Yes\n(2) No"))
                if particao == 1: particao = "Yes"
                if particao == 2: particao = "No"

                self.link[self.connect.index.values[index]] = Conexao_T(self.link[device].m,
                                                                        self.link[device].A_r,
                                                                        self.link[device].A_i,
                                                                        self.link[device].A_o,
                                                                        self.link[device].rho,
                                                                        self.link[device].zeta,
                                                                        juncao,referencia, particao)
            # DEVICE == DIFUSOR
            if type_device == 'df':
                angulo = float(input(f"[DIFUSOR] Digite o Valor do Ângulo do Difusor {device} [º]"))
                length = float(input(f"[DIFUSOR] Digite o Valor do Comprimento Anterior A Entrada do Difusor {device} [m]"))
                self.link[self.connect.index.values[index]] = Difusor(self.link[device].m,
                                                                        self.link[device].A_r,
                                                                        self.link[device].A_i,
                                                                        self.link[device].A_o,
                                                                        self.link[device].rho,
                                                                        self.link[device].zeta,
                                                                        angulo,
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

    def Start_Iteration(self, iterations):
        self.alpha_p = 0.5
        self.alpha_m = 0.5
        self.Identify_Device_Type()
        # criar os index
        self.plot_pressure = []
        self.plot_mass = []

        for x in range(iterations):
            # GERAR DATAFRAME COM OS TERMOS (a_n * p_n) - (a_n-1 * p_n-1)
            self.Set_M_Line_Equations()

            # ATUALIZAR OS VALORES DAS PRESSOES
            self.shape = self.node.loc[self.node_variable].values.shape[0]
            pressure_line = np.array(list(self.result.values()))
            pressure_dot = self.node.loc[self.node_variable].values.reshape(-1)
            delta_p = pressure_line - pressure_dot
            self.node.loc[self.node_variable] = (pressure_dot + self.alpha_p * delta_p).reshape((self.shape, 1))

            # ARMAZENAR DADOS DAS ITERACOES DAS PRESSOES E DAS VAZOES
            self.plot_pressure.append(list(self.node.loc[:].values.reshape(-1)))
            self.plot_mass.append([self.link[self.connect.index.values[index]].m for index in range(len(self.connect.index.values))])

            # ATUALIZAR OS VALORES DAS VAZOES MASSICAS
            ligacoes = self.connect.index.values
            self.mass_line = {ligacoes[index]: self.m_line.cumsum().values[-1][index] for index in range(len(ligacoes))}
            self.mass_dot = {ligacoes[index]: self.link[ligacoes[index]].m for index in range(len(ligacoes))}
            self.delta_m = {ligacoes[index]: self.mass_line[ligacoes[index]] - self.mass_dot[ligacoes[index]] for index in range(len(ligacoes))}

            for link_index in range(len(self.connect.index.values)):
                _ligacao_ = self.connect.index.values[link_index]
                self.link[_ligacao_].m = self.link[_ligacao_].m + self.alpha_m * self.delta_m[_ligacao_]
                self.link[_ligacao_].Set_Zeta()
                self.link[_ligacao_].Set_a()

            self.Set_Mass_Equations()
            self.Set_mass_balance()

            # CRIAR AS MATRIX DO SOLVER
            self.matrix_a = self.mass_balance[:].loc[self.node_variable].values
            self.Add_Pressure_in_Nodes_Boundary()  # ADICIONAR PRESSAO AOS TERMOS DE CONTORNO
            self.matrix_b = self.mass_balance[:].loc[self.node_boundary].cumsum().values[-1] * -1
            self.result = {self.node_variable[i]: np.linalg.solve(self.matrix_a, self.matrix_b)[i] for i in
                           range(len(self.node_variable))}

        # DADOS PARA PLOTAGEM
        self.pressure_df = pd.DataFrame(self.plot_pressure)
        self.mass_df = pd.DataFrame(self.plot_mass)
        self.pressure_df.columns = [f"p_{self.node.index.values[i]}" for i in range(len(self.node.index.values))]
        self.mass_df.columns = self.connect.index.values



class Tubo(Link):
    def __init__(self, m, A_r, A_i, A_o, rho, zeta, l):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu = 1.12e-3
        self.l = l
        self.D = np.sqrt((4 * self.A_r) / np.pi)
        self.U = self.m / (self.A_r * self.rho)
        self.nu = self.mu / self.rho
        self.df = pd.read_excel(
            r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\tubo_chartB.xlsx',
            engine='openpyxl').set_index('Reynolds')

        self.Set_Zeta()
        self.Set_a()

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
    def __init__(self, m, A_r, A_i, A_o, rho, zeta, angle, l, UniformVelocityProfile=True):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu = 1.12e-3
        self.UniformVelocityProfile = UniformVelocityProfile
        self.angle = angle
        self.l = l
        self.Nair = self.A_o / self.A_i
        self.D = np.sqrt((4 * self.A_r) / np.pi)
        self.U = self.m / (self.A_r * self.rho)
        self.nu = self.mu / self.rho
        self.Re = abs(round(self.U * self.D / self.nu))
        self.df_z = {2: pd.read_excel(
            r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_zetad_Nair2.xlsx',
            skiprows=1, engine='openpyxl').set_index('ang'),  #
                     4: pd.read_excel(
                         r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_zetad_Nair4.xlsx',
                         skiprows=1, engine='openpyxl').set_index('ang'),  #
                     6: pd.read_excel(
                         r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_zetad_Nair6.xlsx',
                         skiprows=1, engine='openpyxl').set_index('ang'),  #
                     10: pd.read_excel(
                         r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_zetad_Nair10.xlsx',
                         skiprows=1, engine='openpyxl').set_index('ang'),  #
                     16: pd.read_excel(
                         r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_zetad_Nair16.xlsx',
                         skiprows=1, engine='openpyxl').set_index('ang')}
        self.interpz = {i: interp2d(self.df_z[i].index, self.df_z[i].columns, np.array(self.df_z[i]).T) for i in
                        [2, 4, 6, 10, 16]}
        self.df_kd1 = {50000: pd.read_excel(
            r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n2_Re50000.xlsx',
            skiprows=2, engine='openpyxl').set_index('ang'),  #
                       100000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n2_Re100000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       300000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n2_Re300000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       400000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n2_Re300000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       2000000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n2_Re2000000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       5000000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n2_Re2000000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       6000000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n2_Re6000000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       8000000000000000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n2_Re6000000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang')}

        self.df_kd2 = {50000: pd.read_excel(
            r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n4_Re50000.xlsx',
            skiprows=2, engine='openpyxl').set_index('ang'),  #
                       100000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n4_Re100000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       300000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n4_Re300000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       400000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n4_Re300000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       2000000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n4_Re2000000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       5000000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n4_Re2000000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       6000000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n4_Re6000000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang'),  #
                       8000000000000000: pd.read_excel(
                           r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\difusor_Kd_n4_Re6000000.xlsx',
                           skiprows=2, engine='openpyxl').set_index('ang')}
        self.interpKd1 = {i: interp2d(self.df_kd1[i].index, self.df_kd1[i].columns, np.array(self.df_kd1[i]).T) for i in
                          [50000, 100000, 300000, 400000, 2000000, 5000000, 6000000, 8000000000000000]}
        self.interpKd2 = {i: interp2d(self.df_kd2[i].index, self.df_kd2[i].columns, np.array(self.df_kd2[i]).T) for i in
                          [50000, 100000, 300000, 400000, 2000000, 5000000, 6000000, 8000000000000000]}
        self.Set_Zeta()

    def Zetad(self, alpha, nair, Re):
        nair_values = [2, 4, 6, 10, 16]
        zetad_lista = []
        for i in range(len(nair_values)):
            zetad_lista = zetad_lista + [float(self.interpz[nair_values[i]](alpha, Re))]
        zetad_interp = interp1d(nair_values, zetad_lista)
        return float(zetad_interp(nair))

    def Kd(self, alpha, nair, Re, l_0, D_0):
        Re_values = [50000, 100000, 300000, 400000, 2000000, 5000000, 6000000, 8000000000000000]
        nair_values = [2, 4, 16]
        values_n2 = []
        values_n4 = []
        Kd_values = []
        for i in range(len(Re_values)):
            values_n2 = values_n2 + [float(self.interpKd1[Re_values[i]](alpha, l_0 / D_0))]
            values_n4 = values_n4 + [float(self.interpKd2[Re_values[i]](alpha, l_0 / D_0))]
        for x in range(len(nair_values)):
            Kd_values = Kd_values + [float(interpReynolds[nair_values[x]](Re))]
        interpNair = interp1d(nair_values, Kd_values)

        return float(interpNair(nair))

    def Set_Zeta(self):
        self.zeta_d = self.Zetad(alpha=self.angle, nair=self.Nair, Re=self.Re)
        if self.UniformVelocityProfile == True:
            self.zeta = self.zeta_d
        else:
            self.k_d = self.Kd(alpha=self.angle, nair=self.Nair, Re=self.Re, l_0=self.l, D_0=self.D)
            self.zeta = self.zeta_d * self.k_d

class Conexao_T(Link):
    def __init__(self, m, A_r, A_i, A_o, rho, zeta, Merging, Zeta_cs, Partition):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu = 1.12e-3
        self.Merging = Merging
        self.Zeta_cs = Zeta_cs
        self.Partition = Partition
        self.D = np.sqrt((4 * self.A_r) / np.pi)
        self.U = self.m / (self.A_r * self.rho)
        self.nu = self.mu / self.rho

        self.df_graphA_1 = pd.read_excel(
            r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\conexao_t_zetaLinha1_cs.xlsx',
            skiprows=2, engine='openpyxl').set_index('Qs/Qc')
        self.interp_graphA_1 = interp2d(self.df_graphA_1.index, self.df_graphA_1.columns, np.array(self.df_graphA_1).T)

        self.df_graphA_2 = pd.read_excel(
            r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\conexao_t_zeta1_cs.xlsx',
            skiprows=2, engine='openpyxl').set_index('Qs/Qc')
        self.interp_graphA_2 = interp1d(self.df_graphA_2.index, self.df_graphA_2[1])

        self.df_K1 = pd.read_excel(
            r'C:\Users\clayt\AppData\Local\Programs\Python\Python39\Lib\site-packages\DuctNet\loss_in_devices\conexao_t_df_k1.xlsx',
            skiprows=1, engine='openpyxl').set_index('Qs/Qc')
        self.interp_k1 = interp1d(self.df_K1.index, self.df_K1[90])

        self.A_s = 1
        self.A_c = 1
        if self.Merging == True:
            self.A_s = self.A_i
            self.A_c = self.A_o
        if self.Merging == False:
            self.A_s = self.A_o
            self.A_c = self.A_i

        self.U_s = self.m / (self.A_s * self.rho)
        self.U_c = self.m / (self.A_c * self.rho)
        self.Qs = self.A_s * self.U_s
        self.Qc = self.A_c * self.U_c

        self.A = self.A_parameter(self.A_s, self.A_c, self.Qs, self.Qc)
        self.k_1 = self.k1_parameter(self.A_s, self.A_c, self.Qs, self.Qc)
        self.Merging_func2 = self.A * (1 + (self.A_c / self.A_s) ** 2 + 3 * (self.A_c / self.A_s) ** 2 * (
                    (self.Qs / self.Qc) ** 2 - self.Qs / self.Qc))
        self.Dividing_func1 = 1 + self.k_1 * (self.U_s / self.U_c) ** 2

        self.Set_Zeta()

    def k1_parameter(self, Fs, Fc, Qs, Qc):
        if Fs / Fc <= 1:
            return float(self.interp_k1(Qs / Qc))
        if Fs / Fc > 1:
            if Qs / Qc <= 0.4:
                return 0.9
            if Qs / Qc > 0.4:
                return 0

    def A_parameter(self, Fs, Fc, Qs, Qc):
        if Fs / Fc <= 0.35 and Qs / Qc <= 1:
            return 1.0
        if Fs / Fc > 0.35 and Qs / Qc <= 0.4:
            return 0.9 * (1 - Qs / Qc)
        if Fs / Fc > 0.35 and Qs / Qc > 0.4:
            return 0.55

    def ConvertZeta_cs(self, Fs, Fc, Qs, Qc, Zeta_cs):
        return Zeta_cs / ((Qs * Fc) / (Qc * Fs)) ** 2

    def Verificador(self, Fs, Fc):
        if Fs / Fc == 1:
            Merging_func1 = float(self.interp_graphA_2(self.Qs / self.Qc))
            return Merging_func1
        else:
            print('Nao existe dados de Juncao com Particaoo para Fs/Fc = {}'.format(Fs / Fc))

    def Set_Zeta(self):
        if self.Merging:
            WithPartition_Zeta_cs = lambda Zeta_cs: self.Verificador(self.A_s, self.A_c) if str(
                self.Zeta_cs) == 'Yes' else self.ConvertZeta_cs(self.A_s, self.A_c, self.Qs, self.Qc,
                                                                Zeta_cs=self.Verificador(self.A_s, self.A_c))
            WithoutPartition_Zeta_cs = lambda Zeta_cs: self.Merging_func2 if str(
                self.Zeta_cs) == 'Yes' else self.ConvertZeta_cs(self.A_s, self.A_c, self.Qs, self.Qc,
                                                                Zeta_cs=self.Merging_func2)
            VerifierPartition = lambda Partition: WithPartition_Zeta_cs(self.Zeta_cs) if str(
                self.Partition) == 'Yes' else WithoutPartition_Zeta_cs(self.Zeta_cs)
            self.zeta = VerifierPartition(self.Partition)

        if not self.Merging:
            WithoutPartition_Zeta_cs = lambda Zeta_cs: self.Dividing_func1 if str(Zeta_cs) == 'Yes' else ConvertZeta_cs(
                self.A_s, self.A_c, self.Qs, self.Qc, Zeta_cs=self.Dividing_func1)
            WithPartition_Zeta_cs = lambda Zeta_cs: print('Nao existe dados de Divisao com Particao.') if str(
                Zeta_cs) == 'Yes' else print('Nao existe dados de Divisao com Particao.')
            VerifierPartition = lambda Partition: WithoutPartition_Zeta_cs(self.Zeta_cs) if str(
                self.Partition) == 'No' else WithPartition_Zeta_cs(self.Zeta_cs)
            self.zeta = VerifierPartition(self.Partition)