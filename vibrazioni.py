# Author: MohsenGhalbi 
# Date: 07/01/2020
# Description: Sistema vibrante

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import integrate

# parametri
fn = 1; 
wn = 2*np.pi*fn;
xi = 1.5; # fattore di smorzamento
m = 1;
F0 = 0; # forzante
W = 0.2*wn; # pulsazione della forzante
tf = 10; # tempo di integrazione
dt = 0.001; # passo di integrazione
t = np.arange(0,tf,dt)
y1 = [0.1]; # posizione
y2 = [-10]; # velocità
y1p = [];
y2p = [];

def ode45_step(f, x, t, dt, *args):
    """
    One step of 4th Order Runge-Kutta method
    """
    k = dt
    k1 = k * f(t, x, *args)
    k2 = k * f(t + 0.5*k, x + 0.5*k1, *args)
    k3 = k * f(t + 0.5*k, x + 0.5*k2, *args)
    k4 = k * f(t + dt, x + k3, *args)
    return x + 1/6. * (k1 + 2*k2 + 2*k3 + k4)

def ode45(f, t, x0, *args):
    """
    4th Order Runge-Kutta method
    """
    n = len(t)
    x = np.zeros((n, len(x0)))
    x[0] = x0
    for i in range(n-1):
        dt = t[i+1] - t[i] 
        x[i+1] = ode45_step(f, x[i], t[i], dt, *args)
    return x

def sist_1gdl(t, y, xi, wn, F0, W, m):	
	dy0 = y[1];
	dy1=-2*xi*wn*y[1]-pow(wn,2)*y[0]+F0/m*np.sin(W*t);
	dy = np.array([dy0, dy1]).transpose()
	return dy;

y = ode45(sist_1gdl, t, [y1[0], y2[0]], xi, wn, F0, W, m)

for i in  np.arange(len(t)):
    y1p.append(y2[i]);
    y2p.append(-2*xi*wn*y2[i]-pow(wn,2)*y1[i]+F0/m*np.sin(W*t[i]));
    y1.append(y1[i]+y1p[i]*dt);
    y2.append(y2[i]+y2p[i]*dt);


plt.plot(t,y1[:-1],t,y[:,0]);
plt.show();
