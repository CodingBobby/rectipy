import numpy as np
import pandas as pd


def A(name):
   if name is 'water':
      return 23.1965
   elif name is 'ethanol':
      return 23.8048
   elif name is 'ethyl acetate':
      return 21.0445
   elif name is 'benzene':
      return 20.7937
   elif name is 'p-xylene':
      return 20.9892
   elif name is 'toluene':
      return 20.9065


def B(name):
   if name is 'water':
      return 3816.44
   elif name is 'ethanol':
      return 3803.98
   elif name is 'ethyl acetate':
      return 2790.5
   elif name is 'benzene':
      return 2788.51
   elif name is 'p-xylene':
      return 3346.65
   elif name is 'toluene':
      return 3096.52


def C(name):
   if name is 'water':
      return -46.13
   elif name is 'ethanol':
      return -41.68
   elif name is 'ethyl acetate':
      return -57.15
   elif name is 'benzene':
      return -52.36
   elif name is 'p-xylene':
      return -57.84
   elif name is 'toluene':
      return -53.67


def P(name, T):
   p = np.exp(A(name) - B(name)/(C(name) + T + 237.15))
   return p


def T(name, p):
   return B(name)/(A(name) - np.log(p)) - C(name) - 237.15


class vL(object):
   A_12 = 0
   A_21 = 0


def vanLaar(name_1, name_2):
   c = vL()

   data = pd.read_csv('modules/data/vanLaar.csv', sep=',', na_values='.')
   
   for i, j in enumerate(data['comp1'].array):
      if j == name_1:
         if data['comp2'].array[i] == name_2:
            c.A_12 = float(data['A12'].array[i])
            c.A_21 = float(data['A21'].array[i])

   return c
