
from DuctNet import DuctNet as dn
from importlib import reload


# PATH - Caminho da onde estão seus arquivos .inp
path = "C:\\Users\\clayt\\OneDrive\\Área de Trabalho\\Projeto_Rede"

# REDE_DUTOS - Variavel para criar sua rede de dutos conforme os parametros dos arquivos .inp
net = dn.Net(f'{path}\\node.inp',
            f'{path}\\link.inp',
            f'{path}\\connect.inp',
            f'{path}\\d_parameters.inp',
            f'{path}\\dc_parameters.inp',
            f'{path}\\df_parameters.inp',
            f'{path}\\ct_parameters.inp',
             f'{path}\\rd_parameters.inp',
             f'{path}\\ve_parameters.inp',
             f'{path}\\ca_parameters.inp',
             f'{path}\\po_parameters.inp')

# FUNÇÃO - Ativar iterações com os respectivos parâmetros.
Rede_Dutos.Start_Iteration(iterations=1000, tol=1, alpha=0.4)









