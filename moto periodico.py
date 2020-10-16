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
vhpigrad = vphi*18*np.pi
vFa = Fa*np.ones(len(vphia))
vFm = Fm*np.ones(len(vphim))
vFr = np.concatenate((vFa, vFm), axis=None)
taum = -r*np.sin(vphi)
vMr = vFr*np.array(taum).transpose()

#Andamento del momento resistente ridotto
plt.plot(vhpigrad,vMr, label='$M_r^* [Nm]$')

plt.xlabel(r'$\varphi [°]$')
plt.ylabel('$M_r^*$ [Nm]')
plt.grid()
plt.show()

# coppia d'inerzia
vMi = -m*pow(r,2)*pow(w,2)*np.array(np.sin(vphi))*np.array(np.cos(vphi)).transpose()
# coppia totale
vMmrid = Mmrid*np.ones(len(vphi))
vMtot = vMmrid + vMr +vMi
plt.plot(vhpigrad,vMmrid, label='$M_m^* [Nm]$')
plt.plot(vhpigrad,vMi, label='$M_i^* [Nm]$')
plt.plot(vhpigrad,vMr, label='$M_r^* [Nm]$')
plt.legend(loc="lower left")
plt.xlabel(r'$\varphi [°]$')
plt.ylabel('Nm')
plt.grid()
plt.show()

vE=integrate.cumtrapz(vMtot,vphi, initial=0)

plt.plot(vhpigrad,vE,label='$E_r$')
plt.plot(vhpigrad,vMtot,label='$M_m^*+M_r^*+M_i^* [Nm]$')
plt.legend(loc="lower left")
plt.xlabel(r'$\varphi [°]$')
plt.title('Andamento di $M^∗$ e di $E_r$')
plt.grid()
fig = plt.gcf()
fig.canvas.set_window_title('Andamento di $M^∗$ e di $E_r$')
plt.show()

DEmax=np.max(vE)-np.min(vE)
Jtrid=DEmax/(i*pow(w,2))
Jv=Jtrid*pow(tau,2)/eta-Jm
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
	phipp=(Mmrid+Mr-m*pow(r,2)*pow(phip[j],2)*np.sin(phi[j])*np.cos(phi[j]))/(Jvrid+Jmrid+m*pow(r,2)*pow(np.sin(phi[j]),2))
	phip.append(phip[j]+phipp*dt)
	phi.append(phi[j]+phip[j]*dt)


plt.plot(t,np.array(phip[:-1])/(2*np.pi/60))
plt.grid()
plt.ylabel('n [rpm]')
plt.xlabel('tempo [s]')
fig = plt.gcf()
plt.title('Andamento della velocità dell\'albero di manovella')
fig.canvas.set_window_title("Andamento della velocità dell'albero di manovella")
plt.show()


