# Author: Mohsen Ghalbi 
# Date: 31/12/2019
# Description: Qudrilatero articolato

import numpy as np
import matplotlib.pyplot as plt
import cmath
from numpy.linalg import inv
import matplotlib as mpl

r=100e-3
b=315e-3
c=500e-3
d=425e-3
delta=7/4*np.pi
m2=2
J2=0.016
m3=3.5
J3=0.066
n=150; # [rpm]
w=n*2*np.pi/60  # [rad/s]

dalpha = 2*np.pi/100
alphaf = 2*np.pi
valpha = np.arange(0,alphaf, dalpha)
toll = 0.001
beta0 = 15*np.pi/180
gamma0 = 100*np.pi/180
Mm = []

def resto(alpha, x):
	beta = x[0]
	gamma = x[1]
	f = r*cmath.exp(complex(0,1)*alpha) + b*cmath.exp(complex(0,1)*beta)-c*cmath.exp(complex(0,1)*gamma)- d*cmath.exp(complex(0,1)*delta)
	return np.array([f.real, f.imag]).transpose()

def jacob(x):
	beta = x[0]
	gamma = x[1]
	dfbeta = complex(0,1)*b*cmath.exp(complex(0,1)*beta)
	dfgamma = -complex(0,1)*b*cmath.exp(complex(0,1)*gamma)
	return np.array([[dfbeta.real, dfgamma], [dfbeta.imag, dfgamma.imag]])

x = [beta0, gamma0]
matx = []
matx1 = []
for i in np.arange(0, len(valpha)):
	# anailisi di posizione
	Dx = [10*toll, 10*toll]
	R = [10*toll, 10*toll]
	
	while (abs(Dx[0]) > toll or abs(Dx[1]) > toll or abs(R[0]) > toll or abs(R[1]) > toll):
		R = resto(valpha[i],x)		
		J = jacob(x)
		Dx = np.dot(-inv(J),R)
		x = x + Dx
	beta = x[0]
	gamma = x[1]
		
	# analisi della velocità
	cbetap = -complex(0,1)*b*cmath.exp(complex(0,1)*beta)
	cgammap = complex(0,1)*c*cmath.exp(complex(0,1)*gamma)
	Tn = complex(0,1)*r*w*cmath.exp(complex(0,1)*valpha[i])
	xp = np.dot(inv(np.array([[cbetap.real, cgammap.real], [cbetap.imag, cgammap.imag]])),np.array([Tn.real, Tn.imag]).transpose())
	betap = xp[0]
	gammap = xp[1]
	matx.append(xp);	
	# analisi della accelerazione
	Tna = r*pow(w, 2)*cmath.exp(complex(0,1)*valpha[i])+b*pow(betap,2)*cmath.exp(complex(0,1)*beta)-c*pow(gammap,2)*cmath.exp(complex(0,1)*gamma)
	xpp = np.dot(inv(np.array([[cbetap.real, cgammap.real], [cbetap.imag, cgammap.imag]])),np.array([Tna.real, Tna.imag]).transpose())
	matx1.append(-xpp);
	print(matx1)	
	betapp = -xpp[0]
	gammapp = -xpp[1] 
	# velocità e accelerazioni baricentri
	vG2 = complex(0,1)*w*r*cmath.exp(complex(0,1)*valpha[i])+complex(0,1)*betap*b/2*cmath.exp(complex(0,1)*beta);
	vG3 = complex(0,1)*gammap*c/2*cmath.exp(complex(0,1)*gamma)
	aG2 = -pow(w,2)*r*cmath.exp(complex(0,1)*valpha[i])+complex(0,1)*betapp*b/2*cmath.exp(complex(0,1)*beta)-pow(betap,2)*b/2*cmath.exp(complex(0,1)*beta)
	aG3 = complex(0,1)*gammapp*c/2*cmath.exp(complex(0,1)*gamma)-pow(gammap,2)*c/2*cmath.exp(complex(0,1)*gamma)
	g = 9.81
	P2 = np.array([0, -m2*g]) 
	P3 = np.array([0, -m3*g])
	vvG2 = np.array([vG2.real, vG2.imag]).transpose()
	vvG3 = np.array([vG3.real, vG3.imag]).transpose()
	vaG2 = np.array([aG2.real, aG2.imag])
	vaG3 = np.array([aG3.real, aG3.imag])	
	Mm.append((-np.dot(P2,vvG2)-np.dot(P3,vvG3)-(-m2*np.dot(vaG2,vvG2))-(-m3*np.dot(vaG3,vvG3))-(-J2*betapp*betap)-(-J3*gammapp*gammap))/w)	

X = [x[0] for x in matx]
X1 = [x[0] for x in matx1]
Y = [x[1] for x in matx]
Y1 = [x[1] for x in matx1]

plt.plot(valpha, Mm)
plt.show()
#mpl.style.use('seaborn')

#fig, (ax1, ax2) = plt.subplots(2);

#ax1.plot(valpha, X, 'C1', label='Forza motrice')
#ax2.plot(valpha, X1, 'C1', label='Forza motrice')
#plt.show()

#ax1.plot(valpha, Y, 'C', label='Forza motrice')
#ax2.plot(valpha, Y1, 'C', label='Forza motrice')
plt.show()
