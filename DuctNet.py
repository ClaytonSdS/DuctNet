import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import re
import os
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline

from scipy.interpolate import interp1d as interp1d
from scipy.interpolate import interp2d as interp2d
from scipy.optimize import minimize_scalar as minimize_scalar
from scipy.optimize import minimize as minimize
from scipy.optimize import bisect as bisect
import pathlib
import os
import random





win_path = r'{}'.format(str(os.path.dirname(__file__)))
loss_path = pathlib.PurePath(str(win_path),'loss_in_devices')
path_example = str(pathlib.PurePath(str(win_path), 'loss_in_devices', 'conexao_t_zetaLinha1_cs.xlsx'))


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
        try:
            if self.Merging:
                self.U = self.m_extra / (self.A_r * self.rho)
            if not self.Merging:
                self.U = self.m_extra / (self.A_r*self.rho)
        except AttributeError:
            self.U = self.m / (self.A_r * self.rho)

        self.p_d = self.rho * self.U ** 2 / 2
        self.DP = self.zeta * self.p_d

    def Set_a(self):
        try:
            # CONEXÃO T NÃO SIMÉTRICA

            # CONEXÃO T SIMÉTRICA
            if self.Merging:
                Ar_dot = self.A_o * (self.m / self.m_extra)
                Ao_dot = self.A_o * (self.m / self.m_extra)
                self.a = 1 / ((self.m/ (2 * self.rho)) * (self.zeta/Ar_dot**2 + 1/Ao_dot**2 - 1/ self.A_i**2))

            if not self.Merging:
                Ar_dot = self.A_i * (self.m / (self.m_extra))
                Ai_dot = self.A_i * (self.m / (self.m_extra))
                self.a = 1 / ((self.m / (2 * self.rho)) * (self.zeta/Ar_dot**2 + 1/self.A_o**2 - 1/Ai_dot**2))

        except AttributeError:
            self.a = 1 / ((self.m / (2 * self.rho)) * (self.zeta / self.A_r ** 2 + 1/self.A_o ** 2 - 1 / self.A_i ** 2))

class Node():
    def __init__(self, p):
        self.p = p

class Net():
    def __init__(self, fnode, flink, fconnect, pipes, pipes_c, diffusers, tees, reductions, ball_valves, cotovelos_abruptos, placas_orificios):
        self.fnode = fnode
        self.flink = flink
        self.fconnect = fconnect
        self.ALL_pipes = pipes
        self.ALL_diffusers = diffusers
        self.ALL_tees = tees
        self.ALL_pipes_corrugated = pipes_c
        self.ALL_reductions = reductions
        self.ALL_ballvalves = ball_valves
        self.ALL_cotovelosabruptos = cotovelos_abruptos
        self.ALL_placasorificios = placas_orificios


        self.node = pd.read_csv(fnode, sep='|', index_col='idx')
        self.dflink = pd.read_csv(flink, sep='|', index_col='idx')
        self.connect = pd.read_csv(fconnect, sep='|', index_col='lk')

        # NÓS VARIAVEIS (BALANÇO DE MASSA)
        #self.node_variable = list(set(self.connect['from'].values).intersection(self.connect['to'].values))
        self.node_variable = self.node.loc[self.node["Condicao_Contorno"]!=True].index

        # CONDIÇÕES DE CONTORNO DE VAZÃO MASSICA
        self.mass_boundary = self.dflink.loc[self.dflink["Condicao_Contorno"] == True].index
        self.MassBoundary_Analyser()

        # CONDIÇÕES DE CONTORNO DE PRESSÃO
        self.node_boundary = self.node.loc[self.node["Condicao_Contorno"]==True].index


        # ATRIBUIR VALORES GUESS PARA M E P
        self.DictLink()

        # IDENTIFICAR DISPOSITIVOS E ADICIONAR SEUS RESPECTIVOS PARAMETROS
        self.Identify_Device_Type()

        # CRIAÇÃO DE LISTAS VAZIAS PARA COLOCAR LIGAÇÕES DE JUNÇÕES E DIVISÕES
        self.merging_links = []
        self.division_links = []

        # ADICIONAR M_EXTRA SE LEN > 1:
        # CASO DE JUNÇÃO E DIVISÃO
        for x in range(len(self.node.index.values)):
            no = self.node.index.values[x]
            # JUNCAO
            if len(self.connect.loc[self.connect['to']==no].index.values) > 1:
                self.addList_Merging(list(self.connect.loc[self.connect['to']==no].index.values))
            # DIVISAO
            if len(self.connect.loc[self.connect['from']==no].index.values) > 1:
                self.addList_Division(list(self.connect.loc[self.connect['from']==no].index.values))

        # ATUALIZAR INFORMAÇÕES DE MERGING
        for juncao in range(len(self.merging_links)):
            self.link[self.merging_links[juncao]].Merging = True


        # CASO DE DIVISÃO
        for divisao in range(len(self.division_links)):
            _divisao_ = self.division_links[divisao]
            self.link[self.division_links[divisao]].Merging = False
            self.link[_divisao_].Qc = self.link[_divisao_].m*2/(self.link[_divisao_].rho)
            self.link[_divisao_].Qs = self.link[_divisao_].m/ (self.link[_divisao_].rho)
            self.link[_divisao_].Set_Parameters()
            self.link[_divisao_].Set_Zeta()

        # CALCULAR VAZOES EXTRAS
        self.Refresh_M_Extra()


        # ATIVAR PARAMETROS PARA JUNÇÃO
        for juncao in range(len(self.merging_links)):
            _ligacao_ = self.merging_links[juncao]
            self.link[_ligacao_].m_extra = self.dflink.loc[_ligacao_, "m_extra"]
            self.link[_ligacao_].Set_Parameters()

        # ATIVAR PARAMETROS PARA DIVISAO
        for divisao in range(len(self.division_links)):
            _ligacao_ = self.division_links[divisao]
            self.link[_ligacao_].m_extra = self.dflink.loc[_ligacao_, "m_extra"]
            self.link[_ligacao_].Set_Parameters()

        # ATIVAR A FUNÇÃO PARA OS CALCULOS DOS VALORES DE ZETA E a
        for link_index in range(len(self.connect.index.values)):
            _ligacao_ = self.connect.index.values[link_index]
            self.link[_ligacao_].Set_Zeta()
            self.link[_ligacao_].Set_a()

        # COLOCAR AS LIGACOES ASSOCIADAS AOS NÓS VARIAVEIS
        self.links_in_nodes = {self.node_variable[node]: self.Set_In_Out(self.node_variable[node]) for node  in range(len(self.node_variable))}

        self.Refresh_Px()

        # GERAR AS EQUAÇÕES PARA CADA LIGAÇÃO
        self.Set_Mass_Equations()

        # GERAR OS BALANÇOS DE MASSA DOS NÓS VARIAVEIS
        self.Set_Mass_Balance()

        # CRIAR AS MATRIX DO SOLVER
        self.matrix_a = self.mass_balance[:].loc[self.node_variable].values
        self.Add_Pressure_in_Nodes_Boundary()   # ADICIONAR PRESSAO AOS TERMOS DE CONTORNO
        self.matrix_b = self.mass_balance[:].loc[self.node_boundary].cumsum().values[-1] * -1

        self.result = {self.node_variable[i]:np.linalg.solve(self.matrix_a,self.matrix_b)[i] for i in range(len(self.node_variable))}

    def DictLink(self):
        self.link = {ID: Link(self.dflink.loc[ID, 'm'],
                              self.dflink.loc[ID, 'A_r'],
                              self.dflink.loc[ID, 'A_i'],
                              self.dflink.loc[ID, 'A_o'],
                              self.dflink.loc[ID, 'rho'],
                              self.dflink.loc[ID, 'zeta'])
                     for ID in self.dflink.index}

    def addList_Merging(self, lista):
        for x in range(len(lista)):
            if not lista[x] in self.merging_links:
                self.merging_links.append(lista[x])

    def addList_Division(self, lista):
        for x in range(len(lista)):
            if not lista[x] in self.division_links:
                self.division_links.append(lista[x])

    # SEPARAR FLUXOS CHEGANDO E SAINDO DE NÓS (IN, OUT)
    def Set_In_Out(self, no):
        return {'in': list(self.connect.loc[self.connect['to']==no].index), 'out': list(self.connect.loc[self.connect['from']==no].index)}

    # FUNÇÃO PARA CALCULO DOS BALANÇO DE MASSAS NOS NÓS VARIAVEIS
    def Set_Mass_Balance(self):
        nodes = self.node.index
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
            _pressure_node_ = self.node.loc[_node_, "p"]
            _a_node_ = self.mass_balance.loc[_node_].values
            self.mass_balance.loc[_node_] = _a_node_ * _pressure_node_

    # FUNÇÃO PARA ASSOCIAR OS TERMOS a_n COM AS PRESSOES p_saida E p_from
    def Set_Mass_Equations(self):
        nodes = self.node.index
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

    def Identify_Device_Type(self):
        # pd.read_csv(fnode, sep='|', index_col='idx')

        self.pipes =  pd.read_csv(self.ALL_pipes, sep='|', index_col='idx')
        self.pipes_corrugated = pd.read_csv(self.ALL_pipes_corrugated, sep='|', index_col='idx')
        self.connecs_t = pd.read_csv(self.ALL_tees, sep='|', index_col='idx')
        self.dfs = pd.read_csv(self.ALL_diffusers, sep='|', index_col='idx')
        self.reductions = pd.read_csv(self.ALL_reductions, sep='|', index_col='idx')
        self.ball_valves = pd.read_csv(self.ALL_ballvalves, sep='|', index_col='idx')
        self.cotovelos_abrupto = pd.read_csv(self.ALL_cotovelosabruptos, sep='|', index_col='idx')
        self.placas_orificio = pd.read_csv(self.ALL_placasorificios, sep='|', index_col='idx')

        # antigos
        #self.pipes = pd.read_csv(pathlib.PurePath(str(win_path), 'data', 'd_parameters.inp'), sep='|', index_col='idx')
        #self.connecs_t = pd.read_csv(pathlib.PurePath(str(win_path), 'data', 'ct_parameters.inp'), sep='|', index_col='idx')
        #self.dfs = pd.read_csv(pathlib.PurePath(str(win_path), 'data', 'd_parameters.inp'), sep='|', index_col='idx')

        for index in range(len(self.connect.index.values)):
            device = self.connect.index.values[index]
            type_device = re.split("\d", device)[0]

            # (self, m, A_r, A_i, A_o, rho, zeta, l, l_c, h)
            # DEVICE == TUBO CORRUGADO
            if type_device == 'dc':
                self.link[self.connect.index.values[index]] = Tubo_Corrugado(self.link[device].m,
                                                                   self.link[device].A_r,
                                                                   self.link[device].A_i,
                                                                   self.link[device].A_o,
                                                                   self.link[device].rho,
                                                                   self.link[device].zeta,
                                                                   self.pipes_corrugated.loc[self.connect.index.values[index]].values[0],
                                                                   self.pipes_corrugated.loc[self.connect.index.values[index]].values[2],
                                                                   self.pipes_corrugated.loc[self.connect.index.values[index]].values[1])

            # DEVICE == TUBO
            if type_device == 'd':
                self.link[self.connect.index.values[index]] = Tubo(self.link[device].m,
                                                                        self.link[device].A_r,
                                                                        self.link[device].A_i,
                                                                        self.link[device].A_o,
                                                                        self.link[device].rho,
                                                                        self.link[device].zeta,
                                                                        self.pipes.loc[self.connect.index.values[index]].values[0])

            # DEVICE == COTOVELO ABRUPTO
            if type_device == 'ca':
                self.link[self.connect.index.values[index]] = Cotovelo_Abrupto(self.link[device].m,
                                                                   self.link[device].A_r,
                                                                   self.link[device].A_i,
                                                                   self.link[device].A_o,
                                                                   self.link[device].rho,
                                                                   self.link[device].zeta,
                                                                   self.cotovelos_abrupto.loc[
                                                                       self.connect.index.values[index]].values[0])

            # DEVICE == VALVULA ESFERA
            if type_device == 've':
                self.link[self.connect.index.values[index]] = Valvula_Esfera(self.link[device].m,
                                                                   self.link[device].A_r,
                                                                   self.link[device].A_i,
                                                                   self.link[device].A_o,
                                                                   self.link[device].rho,
                                                                   self.link[device].zeta,
                                                                   self.ball_valves.loc[
                                                                       self.connect.index.values[index]].values[0],
                                                                   self.ball_valves.loc[
                                                                       self.connect.index.values[index]].values[1],
                                                                   self.ball_valves.loc[
                                                                       self.connect.index.values[index]].values[2]
                                                                   )

            # DEVICE == REDUCAO ABRUPTA
            if type_device == 'rd':
                self.link[self.connect.index.values[index]] = Reducao(self.link[device].m,
                                                                   self.link[device].A_r,
                                                                   self.link[device].A_i,
                                                                   self.link[device].A_o,
                                                                   self.link[device].rho,
                                                                   self.link[device].zeta,
                                                                   self.reductions.loc[
                                                                       self.connect.index.values[index]].values[0])

            # DEVICE == CONEXÃO EM T
            if type_device == 'ct':
                juncao = "Yes"

                simetrico = self.connecs_t.loc[device, "Simetric"]

                if simetrico != False and simetrico != True:
                    simetrico = True

                straight = self.connecs_t.loc[device, "Straight"]


                self.link[self.connect.index.values[index]] = Conexao_T(self.link[device].m,
                                                                        self.link[device].m*2,
                                                                        self.link[device].A_r,
                                                                        self.link[device].A_i,
                                                                        self.link[device].A_o,
                                                                        self.link[device].rho,
                                                                        self.link[device].zeta,
                                                                        juncao, self.connecs_t.loc[self.connect.index.values[index]].values[0], simetrico, straight)

            a = 2
            # DEVICE == PLACA DE ORIFICIO
            if type_device == 'po':
                self.link[self.connect.index.values[index]] = PlacaOrificio(self.link[device].m,
                                                                        self.link[device].A_r,
                                                                        self.link[device].A_i,
                                                                        self.link[device].A_o,
                                                                        self.link[device].rho,
                                                                        self.link[device].zeta,
                                                                        self.placas_orificio.loc[self.connect.index.values[index]].values[0])

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
        nodes = self.node.index
        self.m_line = self.equations[:]
        for node_index in range(len(self.node.index)):
            _node_ = self.node.index[node_index]
            _pressure_node_ = self.node.loc[_node_, "p"]
            _a_node_ = self.m_line.loc[_node_].values
            self.m_line.loc[_node_] = _a_node_ * _pressure_node_

    def Refresh_M_Extra(self):
        for x in range(len(self.node.index.values)):
            no = self.node.index.values[x]
            # JUNÇÃO
            if len(self.connect.loc[self.connect['to'] == no].index.values) > 1:
                self.dflink.loc[self.connect.loc[self.connect['to'] == no].index.values, "m_extra"] = self.dflink.loc[self.connect.loc[self.connect['to'] == no].index.values, "m"].sum()

            # DIVISÃO
            if len(self.connect.loc[self.connect['from'] == no].index.values) > 1:
                self.dflink.loc[self.connect.loc[self.connect['from'] == no].index.values, "m_extra"] = (self.dflink.loc[self.connect.loc[self.connect['from'] == no].index.values, "m"].sum())

    def Refresh_Px(self):
        for x in range(len(self.mass_boundary)):
            ligacao = self.mass_boundary[x]
            node_right = self.connect.loc[ligacao]["to"]
            node_left = self.connect.loc[ligacao]["from"]

            # NÓ A ESQUERDA É FIXADO -  ENTRADA
            if self.connect.loc[ligacao]["from"] in self.node_boundary:
                K = self.link[ligacao].m  / self.link[ligacao].a
                p1 = self.node.loc[node_right].values[0]
                self.px = K + p1
                self.node.loc[node_left, "p"] = self.px

            # NÓ A DIREITA É FIXADO - SAIDA
            if self.connect.loc[ligacao]["to"] in self.node_boundary:
                K = self.link[ligacao].m / self.link[ligacao].a
                p1 = self.node.loc[node_left].values[0]
                self.px =  p1 - K
                self.node.loc[node_right, "p"] = self.px

    # FIXAR COMO NÓ DE CONTORNO PARA PRESSÃO A JUSANTE/ A MONTANTE PARA CONDIÇÕES DE CONTORNO DE VAZÃO MÁSSICA
    def MassBoundary_Analyser(self):
        Possible_Nodes_Boundary = list(set(self.connect['from'].values).symmetric_difference(set(self.connect['to'].values)))
        for x in range(len(Possible_Nodes_Boundary)):
            if len(self.connect.loc[self.mass_boundary].loc[self.connect['from'] == Possible_Nodes_Boundary[x]]) > 0 or len(self.connect.loc[self.mass_boundary].loc[self.connect['to'] == Possible_Nodes_Boundary[x]]) > 0:
                Real_Node_Boundary = Possible_Nodes_Boundary[x]
                self.node.loc[Real_Node_Boundary, "Condicao_Contorno"] = True

    def Start_Iteration(self, alpha=0.2, iterations=1000, tol=1e-16):
        self.alpha_p = alpha
        self.alpha_m = alpha

        self.plot_pressure = []
        self.plot_mass = []
        err = 0


        for x in range(iterations):
            self.Refresh_Px()

            # GERAR DATAFRAME COM OS TERMOS (a_n * p_n) - (a_n-1 * p_n-1)
            self.Set_M_Line_Equations()

            # ATUALIZAR OS VALORES DAS PRESSOES
            self.shape = self.node.loc[self.node_variable, "p"].shape[0]
            pressure_line = np.array(list(self.result.values()))
            pressure_dot = self.node.loc[self.node_variable, "p"].to_numpy().reshape(-1)
            delta_p = pressure_line - pressure_dot

            err = np.sqrt(sum([delta**2 for delta in delta_p])/len(delta_p))

            self.node.loc[self.node_variable, "p"] = (pressure_dot + self.alpha_p * delta_p).reshape((self.shape, 1))

            # ATUALIZAR OS VALORES DAS VAZOES MASSICAS
            ligacoes = self.connect.index.values
            self.mass_line = {ligacoes[index]: self.m_line.cumsum().values[-1][index] for index in range(len(ligacoes))}
            self.mass_dot = {ligacoes[index]: self.link[ligacoes[index]].m for index in range(len(ligacoes))}
            self.delta_m = {ligacoes[index]: self.mass_line[ligacoes[index]] - self.mass_dot[ligacoes[index]] for index in range(len(ligacoes))}

            # ATUALIZAR VAZÕES MASSICAS NAS LIGAÇÕES E NO DATAFRAME DFLINK
            for link_index in range(len(self.connect.index.values)):
                _ligacao_ = self.connect.index.values[link_index]
                self.link[_ligacao_].m = self.link[_ligacao_].m + self.alpha_m * self.delta_m[_ligacao_]
                self.dflink.loc[_ligacao_, "m"] = self.link[_ligacao_].m

            # CALCULAR NOVAS VAZOES EXTRAS E ALTERAR NO DATAFRAME
            self.Refresh_M_Extra()

            # ATIVAR PARAMETROS PARA JUNÇÃO
            for juncao in range(len(self.merging_links)):
                _ligacao_ = self.merging_links[juncao]
                self.link[_ligacao_].m_extra = self.dflink.loc[_ligacao_, "m_extra"]
                self.link[_ligacao_].m_line = self.mass_line[_ligacao_]
                self.link[_ligacao_].m_dot = self.mass_dot[_ligacao_]
                self.link[_ligacao_].Set_Parameters()

            # ATIVAR PARAMETROS PARA DIVISAO
            for divisao in range(len(self.division_links)):
                _ligacao_ = self.division_links[divisao]
                self.link[_ligacao_].m_extra = self.dflink.loc[_ligacao_, "m_extra"]
                self.link[_ligacao_].Set_Parameters()

            # ATIVAR A FUNÇÃO PARA OS CALCULOS DOS VALORES DE ZETA E a
            for link_index in range(len(self.connect.index.values)):
                _ligacao_ = self.connect.index.values[link_index]
                self.link[_ligacao_].Set_Zeta()
                self.link[_ligacao_].Set_a()

            self.Set_Mass_Equations()
            self.Set_Mass_Balance()

            # CRIAR AS MATRIX DO SOLVER
            self.matrix_a = self.mass_balance[:].loc[self.node_variable].values
            self.Add_Pressure_in_Nodes_Boundary()  # ADICIONAR PRESSAO AOS TERMOS DE CONTORNO PARA FICAR a_n * p_n
            self.matrix_b = self.mass_balance[:].loc[self.node_boundary].cumsum().values[-1] * -1
            self.result = {self.node_variable[i]: np.linalg.solve(self.matrix_a, self.matrix_b)[i] for i in range(len(self.node_variable))}

            # ARMAZENAR DADOS DAS ITERACOES DAS PRESSOES E DAS VAZOES
            self.plot_pressure.append(list(self.node.loc[:,"p"].values.reshape(-1)))
            self.plot_mass.append([self.link[self.connect.index.values[index]].m for index in range(len(self.connect.index.values))])

            #if iterations <= 100:
               # plt.plot(self.plot_mass)
               # plt.show()
               # plt.pause(0.03)

            if err < tol:
                print(f"[FINISHED BY TOLERANCE]-ITERATIONS [{x}/{iterations}]")
                break

            if x == iterations:
                print(f"[FINISHED BY ITERARATIONS]-ITERATIONS [{x}/{iterations}]")




        # DADOS PARA PLOTAGEM
        #self.pressure_df = pd.DataFrame(self.plot_pressure)
        #self.mass_df = pd.DataFrame(self.plot_mass)
        #self.pressure_df.columns = [f"p_{self.node.index.values[i]}" for i in range(len(self.node.index.values))]
        #self.mass_df.columns = self.connect.index.values

class Tubo(Link):
    def __init__(self, m, A_r, A_i, A_o, rho, zeta, l):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu =  1.849e-5
        self.l = l
        self.D = np.sqrt((4 * self.A_r) / np.pi)
        self.U = self.m / (self.A_r * self.rho)
        self.nu = self.mu / self.rho

        self.df = pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'tubo_chartB.xlsx'),engine='openpyxl').set_index('Reynolds')

        self.Set_Zeta()
        self.Set_a()

    def lambda_(self, Re):
        interp = interp1d(self.df.index, self.df['lambda'].values)
        if Re <= 2000:
            return 64 / Re
        if (Re > 2000) and (Re <= 4000):
            return float(interp(Re))
        if (Re > 4000) and (Re < 100000):
            return 0.3164 / (Re ** 0.25)
        if Re >= 100000:
            return 1 / ((1.8 * np.log10(Re) - 1.64) ** 2)

    def Set_Zeta(self):
        self.U = self.m / (self.A_r * self.rho)
        self.Re = round(self.U * self.D / self.nu)
        lambdaValue = float(self.lambda_(self.Re))
        self.zeta = lambdaValue * (self.l / self.D)

class PlacaOrificio(Link):
    def __init__(self, m, A_r, A_i, A_o, rho, zeta, C):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu =  1.849e-5
        self.C = C

        self.D_i = np.sqrt((4 * self.A_i) / np.pi)
        self.D_o = np.sqrt((4 * self.A_o) / np.pi)
        self.beta = self.D_o / self.D_i

        self.U = self.m / (self.A_i * self.rho)
        self.nu = self.mu / self.rho

        self.Set_Zeta()
        self.Set_a()


    def Set_Zeta(self):
        self.U = self.m / (self.A_i * self.rho)
        self.Re = round(self.U * self.D_i / self.nu)
        self.zeta = (np.sqrt(1- (self.beta**4) * (1 - self.C**2))/(self.C * self.beta ** 2) - 1) ** 2

class Cotovelo_Abrupto(Link):
    def __init__(self, m, A_r, A_i, A_o, rho, zeta, delta):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu =  1.849e-5
        self.D = np.sqrt((4 * self.A_r) / np.pi)
        self.U = self.m / (self.A_r * self.rho)
        self.nu = self.mu / self.rho
        self.delta = delta
        self.delta_rad = self.delta * np.pi/180

        self.Set_Zeta()
        self.Set_a()

    def Set_Zeta(self):
        self.A = 0.95 + 33.5/self.delta
        self.zeta_loc = 0.95*(np.sin(self.delta_rad/2)**2) + 2.05*(np.sin(self.delta_rad/2)**4)
        self.zeta = self.A * self.zeta_loc

class Valvula_Esfera(Link):
    def __init__(self, m, A_r, A_i, A_o, rho, zeta, D_1, D_2, theta):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu =  1.849e-5
        self.D_1 = D_1
        self.D_2 = D_2
        self.theta = theta
        self.beta = self.D_1/self.D_2
        self.U = self.m / (self.A_r * self.rho)
        self.nu = self.mu / self.rho

        self.df = pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'ball_valve_friction_factor.xlsx'),engine='openpyxl').set_index('D')
        self.interp = interp1d(self.df.index, self.df['zeta_fr'].values)

        self.Set_Zeta()
        self.Set_a()

    def Set_Zeta(self):
        self.f_t = self.interp(self.D_2)
        if self.beta == 1 and self.theta == 0:
            self.K = 3 * self.f_t
            #self.zeta = self.K
            self.zeta = self.K/(self.rho*9.81)

class Tubo_Corrugado(Link):
    def __init__(self, m, A_r, A_i, A_o, rho, zeta, l, l_c, h):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu =1.849e-5
        self.l = l
        self.l_c = l_c
        self.h = h
        self.D = np.sqrt((4 * self.A_r) / np.pi)
        self.U = self.m / (self.A_r * self.rho)
        self.nu = self.mu / self.rho


        self.df = pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'corrugated_pipe.xlsx'), engine='openpyxl').set_index('Re')
        #print(self.df)

        self.Set_Zeta()
        self.Set_a()

    def lambda_(self, Re):
        interp = interp2d(self.df.index, self.df.columns, np.array(self.df).T)
        return float(interp(Re, self.h/self.l_c))


    def Set_Zeta(self):
        self.U = self.m / (self.A_r * self.rho)
        self.Re = round(self.U * self.D / self.nu)
        lambdaValue = float(self.lambda_(self.Re))
        self.zeta = lambdaValue * (self.l / self.D)

class Reducao(Link):
    def __init__(self, m, A_r, A_i, A_o, rho, zeta, l):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu = 1.849e-5
        self.l_0 = l
        self.D_0 = np.sqrt((4 * self.A_o) / np.pi)  # Area ref - F0 (saida)
        self.D_1 = np.sqrt((4 * self.A_i) / np.pi)  # Area ref - F1 (entrada)

        self.F_0 = self.A_r
        self.F_1 = self.A_i

        self.U = self.m / (self.A_r * self.rho)
        self.nu = self.mu / self.rho

        self.df = pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'tubo_chartB.xlsx'),engine='openpyxl').set_index('Reynolds')

        self.Set_Zeta()
        self.Set_a()

    def lambda_(self, Re):
        interp = interp1d(self.df.index, self.df['lambda'].values)
        if Re <= 2000:
            return 64 / Re
        if (Re > 2000) and (Re <= 4000):
            return float(interp(Re))
        if (Re > 4000) and (Re < 100000):
            return 0.3164 / (Re ** 0.25)
        if Re >= 100000:
            return 1 / ((1.8 * np.log10(Re) - 1.64) ** 2)

    def Set_Zeta(self):
        self.U = self.m / (self.A_r * self.rho)
        self.Re = round(self.U * self.D_0 / self.nu)
        lambdaValue = float(self.lambda_(self.Re))
        self.zeta_fr = lambdaValue * (self.l_0 / self.D_0)

        self.zeta = 0.5*(1 - self.F_0/self.F_1)**(3/4) + self.zeta_fr

class Difusor(Link):
    def __init__(self, m, A_r, A_i, A_o, rho, zeta, angle, l, UniformVelocityProfile=True):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu = 1.849e-5
        self.UniformVelocityProfile = UniformVelocityProfile
        self.angle = angle
        self.l = l
        self.Nair = self.A_o / self.A_i
        self.D = np.sqrt((4 * self.A_r) / np.pi)
        self.U = self.m / (self.A_r * self.rho)
        self.nu = self.mu / self.rho
        self.Re = abs(round(self.U * self.D / self.nu))

        self.df_z = {2: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_zetad_Nair2.xlsx'),skiprows=1, engine='openpyxl').set_index('ang'),  #
                     4: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_zetad_Nair4.xlsx'),skiprows=1, engine='openpyxl').set_index('ang'),  #
                     6: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_zetad_Nair6.xlsx'),skiprows=1, engine='openpyxl').set_index('ang'),  #
                     10: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_zetad_Nair10.xlsx'),skiprows=1, engine='openpyxl').set_index('ang'),  #
                     16: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_zetad_Nair16.xlsx'),skiprows=1, engine='openpyxl').set_index('ang')}
        self.interpz = {i: interp2d(self.df_z[i].index, self.df_z[i].columns, np.array(self.df_z[i]).T) for i in
                        [2, 4, 6, 10, 16]}

        self.df_kd1 = {50000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n2_Re50000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       100000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n2_Re100000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       300000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n2_Re300000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       400000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n2_Re300000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       2000000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n2_Re2000000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       5000000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n2_Re2000000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       6000000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n2_Re6000000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       8000000000000000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n2_Re6000000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang')}

        self.df_kd2 = {50000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n4_Re50000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       100000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n4_Re100000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       300000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n4_Re300000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       400000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n4_Re300000.xlsx'), skiprows=2, engine='openpyxl').set_index('ang'),  #
                       2000000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n4_Re2000000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       5000000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n4_Re2000000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang'),  #
                       6000000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n4_Re6000000.xlsx'), skiprows=2, engine='openpyxl').set_index('ang'),  #
                       8000000000000000: pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'difusor_Kd_n4_Re6000000.xlsx'),skiprows=2, engine='openpyxl').set_index('ang')}

        self.interpKd1 = {i: interp2d(self.df_kd1[i].index, self.df_kd1[i].columns, np.array(self.df_kd1[i]).T) for i in
                          [50000, 100000, 300000, 400000, 2000000, 5000000, 6000000, 8000000000000000]}
        self.interpKd2 = {i: interp2d(self.df_kd2[i].index, self.df_kd2[i].columns, np.array(self.df_kd2[i]).T) for i in
                          [50000, 100000, 300000, 400000, 2000000, 5000000, 6000000, 8000000000000000]}
        self.Set_Zeta()
        self.Set_a()

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
    def __init__(self, m, m_extra, A_r, A_i, A_o, rho, zeta, Merging, Partition, Simetric, Straight):
        Link.__init__(self, m, A_r, A_i, A_o, rho, zeta)
        self.mu = 1.849e-5
        self.Merging = Merging

        self.Simetric = Simetric
        self.Straight = Straight

        self.m_extra = m_extra
        self.m_line = m
        self.m_dot = m

        # CONDIÇÕES PARA PARTIÇÃO
        self.Partition = Partition
        if Partition == True:
            self.partition = "Yes"

        if Partition == False:
            self.partition = "No"

        # PARAMETROS CONEXÃO T - SIMÉTRICA
        self.df_graphA_1 = pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'conexao_t_zetaLinha1_cs.xlsx'),skiprows=2, engine='openpyxl').set_index('Qs/Qc')
        self.interp_graphA_1 = interp2d(self.df_graphA_1.index, self.df_graphA_1.columns, np.array(self.df_graphA_1).T)

        self.df_graphA_2 = pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'conexao_t_zeta1_cs.xlsx'), skiprows=2, engine='openpyxl').set_index('Qs/Qc')
        self.interp_graphA_2 = interp1d(self.df_graphA_2.index, self.df_graphA_2[1])

        self.df_K1 = pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'conexao_t_df_k1.xlsx'),skiprows=1, engine='openpyxl').set_index('Qs/Qc')
        self.interp_k1 = interp1d(self.df_K1.index, self.df_K1[90])

        # PARAMETROS CONEXÃO T - NÃO SIMÉTRICA
        self.ZetaLinha_cst = pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'conexao_t_DivisaoNaoSimetrica-TrechoReto.xlsx'), engine='openpyxl').set_index('Qs/Qc')
        self.interp_ZetaLinha_cst = interp1d(self.ZetaLinha_cst.index, self.ZetaLinha_cst["Fs/Fc"])

        self.ZetaLinha_cs = pd.read_excel(pathlib.PurePath(str(win_path), 'loss_in_devices', 'conexao_t_DivisaoNaoSimetrica-TrechoLateral.xlsx'), engine='openpyxl').set_index('Qs/Qc')
        self.interp_ZetaLinha_cs = interp2d(self.ZetaLinha_cs.index, self.ZetaLinha_cs.columns, np.array(self.ZetaLinha_cs).T)

        self.A_s = self.A_i
        self.A_c = self.A_r

        if not self.Merging:
            self.U_c = 2*self.m/(self.rho)
            self.U_s = self.m/self.rho

        self.Set_Parameters()
        self.Set_Zeta()
        self.Set_a()

    def k1_parameter(self, Fs, Fc, Qs, Qc):
        Qs = abs(Qs)
        Qc = abs(Qc)
        if Fs / Fc <= 1:
            return float(self.interp_k1(Qs/Qc))
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

    # FUNÇÃO DE INTERPOLAÇÃO PARA - JUNÇÃO COM PARTIÇÃO
    def Verificador(self, Fs, Fc):
        if Fs / Fc == 1:
           # Merging_func1 = float(self.interp_graphA_2(self.Qs / self.Qc))
            Merging_func1 = 1
            return Merging_func1

        else:
            print('Nao existe dados de Juncao com Particao para Fs/Fc = {}'.format(Fs / Fc))

    def Set_Parameters(self):
        self.Qs = self.m / self.rho
        self.Qc = self.m_extra / self.rho

        self.U_s = self.Qs / self.A_s
        self.U_c = self.Qc / self.A_c

        self.R_Q = self.Qs / self.Qc
        self.R_A = self.A_c / self.A_s


        #############################################################################
        # CALCULO DOS PARAMETROS PARA CONEXÃO EM T - NÃO SIMÉTRICA
        #############################################################################
        # JUNÇÃO NÃO SIMÉTRICA - CALCULO DE ZETA NO TRECHO RETO
        if self.Merging and self.Simetric == False and self.Straight == True:
            self.Merging_func4 = 1.55*self.R_Q - self.R_Q**2

        # JUNÇÃO NÃO SIMÉTRICA - CALCULO DE ZETA NO TRECHO LATERAL
        if self.Merging and self.Simetric == False and self.Straight == False:
            self.A = self.A_parameter(self.A_s, self.A_c, self.Qs, self.Qc)
            self.Merging_func3 = self.A * (1 + (self.R_Q * self.R_A) ** 2 - 2 * (1 - self.R_A) ** 2)
            #self.Merging_func3 = (1 + (self.R_Q * self.R_A) ** 2 - 2 * (1 - self.R_A) ** 2)

        # DIVISÃO NÃO SIMÉTRICA - CALCULO DE ZETA NO TRECHO RETO
        if not self.Merging and self.Simetric == False  and self.Straight == True:
            Zeta_cst = self.interp_ZetaLinha_cst(self.Qs/self.Qc)
            self.Dividing_func2 = Zeta_cst / (((1 - self.Qs / self.Qc) ** 2) * ((self.A_c / self.A_s) ** 2))


        # DIVISÃO NÃO SIMÉTRICA - CALCULO DE ZETA NO TRECHO LATERAL
        if not self.Merging and self.Simetric == False and self.Straight == False:
            # PROBLEMA ESTÁ AQUI

            Zeta_cs = self.interp_ZetaLinha_cs(self.Qs/self.Qc, self.A_s/self.A_c)
            self.Dividing_func3 = Zeta_cs / (((self.Qs * self.A_c) / (self.Qc * self.A_s)) ** 2)

        #############################################################################
        # CALCULO DOS PARAMETROS PARA CONEXÃO EM T - SIMÉTRICA
        #############################################################################
        if self.Merging and self.Simetric == True:
            #self.A = self.A_parameter(self.A_s, self.A_c, self.Qs, self.Qc)
            #self.Merging_func2 = self.A * (1 + (self.R_A)**2 + 3*(self.R_A)**2 * ((self.R_Q)** 2 - self.R_Q)) # resultado não converge considerando o parametro A
            self.Merging_func2 = (1 + (self.R_A)**2 + 3 * (self.R_A)**2 * ((self.R_Q)**2 - self.R_Q))

        if not self.Merging and self.Simetric == True:
            self.k_1 = self.k1_parameter(self.A_s, self.A_c, self.Qs, self.Qc)  # usado para DIVISÃO sem partição
            self.Dividing_func1 = 1 + self.k_1 * (self.U_s / self.U_c)**2

    def Set_Zeta(self):

        # CONEXÃO T - NÃO SIMÉTRICA
        if self.Merging and self.Simetric == False and self.Straight == True:
            self.zeta = self.Merging_func4

        if self.Merging and self.Simetric == False and self.Straight == False:
            self.zeta = self.Merging_func3

        if not self.Merging and self.Simetric == False and self.Straight == True:
            self.zeta = self.Dividing_func2

        if not self.Merging and self.Simetric == False and self.Straight == False:
            self.zeta = self.Dividing_func3


        # CONEXÃO T - SIMÉTRICA
        if self.Merging and self.Simetric == True:
            WithPartition = self.Verificador(self.A_s, self.A_c)
            WithoutPartition = self.Merging_func2
            VerifierPartition = lambda Partition: WithPartition if str(self.Partition) == 'Yes' else WithoutPartition
            self.zeta = VerifierPartition(self.Partition)

        if not self.Merging and self.Simetric == True:
            self.zeta = self.Dividing_func1

