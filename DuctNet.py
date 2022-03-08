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

#from UnitConv import UnitConv as uc
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
        self.U = self.m/(self.A_r * self.rho)
        self.p_d = self.rho*self.U**2/2
        self.DP = self.zeta * self.p_d

    def Set_a(self):
        self.a = 1/((self.m/(2*self.rho)) * (self.zeta/self.A_r**2 + 1/self.A_o**2 - 1/self.A_i**2))
   
        
class Node():
    def __init__(self, p):
        self.p = p


class Net():
    def __init__(self, no, lk):
        self.no = no
        self.lk = lk
        self.Set_Net()
        
    def Set_Net(self):
        self.net = np.empty((len(self.no), len(self.no)))
        self.net[:] = np.nan

    def Set_Mass_Balance(self):
        print('oi')


        
    


        
