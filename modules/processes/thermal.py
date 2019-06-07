import modules.data.antoine as Antoine
import numpy as np
import matplotlib.pyplot as plt


class Data(object):
   x = []
   y = []


def lin_m( x1, y1, x2, y2 ):
   return (y2 - y1) / (x2 - x1)


def lin_b( x1, y1, x2, y2 ):
   m = lin_m(x1, y1, x2, y2)
   return y1 - (m * x1)


def lin( x, x1, y1, x2, y2 ):
   m = lin_m(x1, y1, x2, y2)
   b = lin_b(x1, y1, x2, y2)
   return m * x + b


def interpolate(xdata, ydata, x):
   xa = 0
   xb = 0
   ya = 0
   yb = 0

   found = False

   for i in range(0, len(xdata)):
      if found is False:
         if x > xdata[i]:
            yb = ydata[i]
            xb = xdata[i]
         else:
            ya = ydata[i]
            xa = xdata[i]
            found = True

   m = (ya-yb)/(xa-xb)
   b = yb

   return m*(x-xb)+b


class Component(object):
   def __init__(self, name):
      self.name = name
      self.A = Antoine.A(name)
      self.B = Antoine.B(name)
      self.C = Antoine.C(name)


   def pressure(self, temperature): 
      return Antoine.P(self.name, temperature)


   def temperature(self, pressure):
      return Antoine.T(self.name, pressure)


class Mixture(object):
   def __init__(self, component_1, component_2):
      self.comp_1 = component_1
      self.comp_2 = component_2

      # van Laar mixture coefficients
      vL = Antoine.vanLaar(self.comp_1.name, self.comp_2.name)
      self.A_12 = vL.A_12
      self.A_21 = vL.A_21

      # barometric pressure in Pa
      self.p_u = 101325

      # datapoints of vapour-liquid-equilibrium
      self.vle_data = Data()
      self.vle_data.x = np.linspace(0, 1, 200)
      for x in self.vle_data.x:
         y = self.VLE(x)
         self.vle_data.y.append(y)
   

   def VLE(self, x):
      gamma1 = self.gamma_one(x)
      gamma2 = self.gamma_two(x)
      T = self.temperature(x)
      p1 = x * self.comp_1.pressure(T) * gamma1
      p2 = (1-x) * self.comp_2.pressure(T) * gamma2
      p_tot = p1 + p2
      return p1/p_tot


   def VLE_B(self, y):
      return interpolate(self.vle_data.y, self.vle_data.x, y)


   def temperature(self, x):
      # maximum allowed error
      err = 10**(-12)

      Ta = Antoine.T(self.comp_1.name, self.p_u)
      Tb = Antoine.T(self.comp_2.name, self.p_u)

      T1 = Ta - (Ta-Tb)/3
      T2 = Ta - 2*(Ta-Tb)/3

      p1 = self.pressure(x, T1)
      p2 = self.pressure(x, T2)

      pnew = 0
      Tnew = 0

      while np.abs(pnew - self.p_u) > (err * self.p_u):
         Tnew = T1 + (T2 - T1)*(self.p_u - p1)/(p2 - p1)
         pnew = self.pressure(x, Tnew)

         if np.abs(p1 - self.p_u) < np.abs(p2 - self.p_u):
            p2 = pnew
            T2 = Tnew
         else:
            p1 = pnew
            T1 = Tnew

      return Tnew


   def pressure(self, x, T):
      g1 = self.gamma_one(x)
      g2 = self.gamma_two(x)
      a = x * g1 * Antoine.P(self.comp_1.name, T)
      b = (1 - x) * g2 * Antoine.P(self.comp_2.name, T)
      return a + b
   
   
   def gamma_one(self, x):
      num = self.A_12 * (1 - x)**2 * self.A_21**2
      return np.exp(num/self.denom(x))


   def gamma_two(self, x):
      num = self.A_21 * x**2 * self.A_12**2
      return np.exp(num/self.denom(x))


   def denom(self, x):
      d = x * self.A_12 + (1 - x) * self.A_21
      return d**2
   

   def plot(self, **kwargs):
      data = self.vle_data

      fig, ax = plt.subplots(1, 1, figsize = (7, 7))
      ax.set_aspect('equal', 'box')
      ax.axis([0, 1, 0, 1])
      ax.set_xlabel(r'$X_1$')
      ax.set_ylabel(r'$X_2$')
      ax.grid(True)

      ax.plot([0, 1], [0, 1], 'black', lw = 0.7)
      ax.plot(data.x, data.y, 'black', lw = 0.7, label = 'VLE')

      if kwargs.get('noplot', None) is True:
         return fig, ax, plt

      else:
         ax.legend()
         plt.show()
         return 1


class Rectification(object):
   # xA: mol concentration of comp 1 in Rest
   # xD: mol concentration of comp 1 in Product
   def __init__(self, mixture, xA, xF, xD):
      self.mix = mixture
      self.xA  = xA
      self.xF  = xF
      self.xD  = xD
   

   # xF: backfeed
   def min_operation(self):
      x_B = self.mix.VLE(self.xF)
      intercept = self.xD/x_B - 1

      x = [0, self.xD]
      y = [intercept, self.xD]

      vmin = self.xD/intercept - 1

      return x, y, vmin

   
   def opt_operation(self, vmin):
      vopt = vmin * 1.5
      intercept = self.xD/(vopt+1)

      x = [0, self.xD]
      y = [intercept, self.xD]

      return x, y, vopt


   def sub_operation(self, vopt):
      intercept = self.xD/(vopt+1)
      slope = (self.xD-intercept)/self.xD

      height = slope * self.xF + intercept

      x = [self.xA, self.xF]
      y = [self.xA, height]

      return x, y


   def opt_steps(self, xF, inter):
      currX = self.xD
      currY = self.xD
      x = []
      y = []

      # Rectification steps
      while currX > self.xA:
         x.append(currX)
         y.append(currY)
         
         y.append(currY)
         currX = self.mix.VLE_B(currY)
         
         if currX > xF:
            currY = lin(currX, xF, inter, self.xD, self.xD)
         elif currX > self.xA:
            currY = lin(currX, self.xA, self.xA, xF, inter)
         else:
            currY = currX
         
         x.append(currX)
      
      x.append(currX)
      y.append(currY)

      return x, y

   
   def plot(self):
      fig, ax, plt = self.mix.plot(noplot = True)

      mx, my, vmin = self.min_operation()
      ox, oy, vopt = self.opt_operation(vmin)
      sx, sy = self.sub_operation(vopt)

      rx, ry = self.opt_steps(sx[1], sy[1])

      # minimum operation line
      ax.plot(mx, my, 'blue', linestyle = '--', dashes = (6, 8), lw = 0.7)
      # optimum operation line part 1
      ax.plot([ox[0], sx[1]], [oy[0], sy[1]], 'blue', linestyle = '--', dashes = (6, 8), lw = 0.7)
      # optimum operation line part 2
      ax.plot([sx[1], ox[1]], [sy[1], oy[1]], 'blue', lw = 0.7, label = 'Operation Line')
      # optimum sub operation line
      ax.plot(sx, sy, 'green', lw = 0.7, label = 'Output Line')

      # rectification steps
      ax.plot(rx, ry, 'red', lw = 0.7, label = 'Steps')

      ax.legend()
      plt.show()
      return 1
