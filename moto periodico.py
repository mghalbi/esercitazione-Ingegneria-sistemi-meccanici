# Author: Monse Ghalbi
# Date: 30/12/2019
# Description: moto periodico

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import integrate 

pa=-0.5e5 #pressione di aspirazione [Pa] 
pm=4.8e5 #pressione di mandata [Pa]
c=280e-3 #corsa stantuffo [m]
r=c/2  #lunghezza manovella;
D=210e-3 #diametro stantuffo
A=np.pi*pow(D,2)/4
Jm=0.1 #momento d'inerzia motore [kgm^2]
m=54 #massa del piede di biella [kg]
n=195 #191.69 %(202)velocità di rotazione media albero manovella [rpm]
w=n*2*np.pi/60 #velocità di rotazione media albero manovella [rad/s]
i=0.03 #irregolarità periodica
tau=1/7.5
eta=0.85
M0=308 # [Nm]
K=-0.1225 # [Nm/rpm]

def motore(phip):
	nm = phip/(2*np.pi)*60/tau;
	Mm = M0 + K*nm
	Mmrid = Mm*eta/tau
	return Mmrid

# Coppia motrice
Fa = np.abs(pa)*A
Fm = -pm*A
Lr = Fa*c+np.abs(Fm)*c
Mmrid = Lr/(2*np.pi)

# Coppia resistente
dphi = np.pi/100
vphia = np.arange(0,np.pi,dphi)
vphim = np.arange(np.pi, 2*np.pi, dphi)
vphi =  np.concatenate((vphia, vphim), axis=None)
vFa = Fa*np.ones(len(vphia))
vFm = Fm*np.ones(len(vphim))
vFr = np.concatenate((vFa, vFm), axis=None)
taum = -r*np.sin(vphi)
vMr = vFr*np.array(taum).transpose()

# coppia d'inerzia
vMi = -m*pow(r,2)*pow(w,2)*np.array(np.sin(vphi))*np.array(np.cos(vphi)).transpose()
# coppia totale
vMmrid = Mmrid*np.ones(len(vphi))
vMtot = vMmrid + vMr +vMi
plt.plot(vphi,vMtot)
plt.plot(vphi,vMmrid)
plt.plot(vphi,vMi)
plt.plot(vphi,vMr)
plt.show()

vE=integrate.cumtrapz(vMtot,vphi, initial=0)

plt.plot(vphi,vE)
plt.plot(vphi,vMtot)
plt.show()

DEmax=np.max(vE)-np.min(vE)
Jtrid=DEmax/(i*pow(w,2))
Jv=Jtrid*pow(tau,2)/eta-Jm
print(['Jv [kg*m^2]=', Jv])
Jvrid=Jv/pow(tau,2)*eta;
Jmrid=Jm/pow(tau,2)*eta;
dt=0.01; tf=30;
t=np.arange(0,tf,dt);
phi = [0] 
phip = [0]

for j in np.arange(0,len(t)):
	Mmrid = motore(phip[j])
	phig = 2*np.pi*(phi[j]/(2*np.pi)-np.fix(phi[j]/(2*np.pi)))
	Mr = np.interp(phig,vphi,vMr)	
	print((Jvrid+Jmrid+m*pow(r,2)*pow(np.sin(phi[j]),2)))
	phipp=(Mmrid+Mr-m*pow(r,2)*pow(phip[j],2)*np.sin(phi[j])*np.cos(phi[j]))/(Jvrid+Jmrid+m*pow(r,2)*pow(np.sin(phi[j]),2))
	phip.append(phip[j]+phipp*dt)
	phi.append(phi[j]+phip[j]*dt)


plt.plot(t,np.array(phip[:-1])/(2*np.pi/60))
plt.show()


