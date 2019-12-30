# Author: MohsenGhalbi 
# Date: 08/10/2019
# Description: DInamica di un veicolo ferrovario

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# inizializzazione delle variabili
m = 8000
S = 7.5
Cr = 0.7
D =  0.6
Jm = 0.1
tau =5.0/24.0
eta =  0.97
p = 30.0/1000.0
Cs = 900
n0 = 1000
g = 9.81
alpha = np.arctan(p)
rho = 1.25

# Determinazione della funzione caratteristica del motore
k = -Cs/n0
n = np.arange(0,n0,10)
# Curva caratteristica del motore
Cm = Cs + k*n

# Coppia motore ridotta all'utilizzatore
Cmrid = Cm*eta/(tau*D/2)
wm = n*2*np.pi/60 # conversione della velocita in rad/s
wmrid = wm*tau*D/2 # [m/s]
wmridkmh = wmrid * 3.6
Fp = m*g*np.sin(alpha) # forza peso

# curva caratteristica del carico
vkmh = np.arange(0,40) # velocita [km/h]
v = vkmh/3.6 # [m/s]
vrid = v/(tau*D/2)/(2*np.pi/60)

Fa = list(map(lambda x:0.5*S*Cr*rho*pow(x,2),v))

Fr = Fa + Fp

Frrid = Fr*(tau*D/2)/eta

# rappresentazione del grafico
mpl.style.use('seaborn')

fig, (ax1, ax2) = plt.subplots(2);

ax1.plot(wmridkmh, Cmrid, 'C1', label='Forza motrice')
ax1.plot(vkmh, Fr,  'C0', label='Forza resistente')
ax2.plot(n, Cm, 'C1', label='Forza motrice')
ax2.plot(vrid, Frrid,  'C0', label='Forza resistente')
#ax1.xlabel('velocità [km/h]')
#ax1.ylabel('Forza [N]')
#ax1.title('Curva caratteristiche del carico')
#ax1.legend()
plt.show()
