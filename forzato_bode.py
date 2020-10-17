# Author: Mohsen Ghalbi 
# Date: 07/01/2020
# Description: Sistema vibrante Forzato

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab as pl
import cmath


a = np.arange(0,10,0.01);
xi = [0.001, 0.05, 0.1, 0.4, 1, 1.5, 2];


for i in np.arange(len(xi)):
	XXst = list(map(lambda x: 1/(1-pow(x,2)+2*xi[i]*complex(0,1)*x),a))
	XYO = list(map(lambda x: (1+2*xi[i]*complex(0,1)*x)/(1-pow(x,2)+2*xi[i]*complex(0,1)*x),a))
	plt.figure(1)
	ax1 = plt.subplot(211)
	ax1.set_xlim([0, 5])
	ax1.set_ylim(0,6)
	ax1.set_ylabel('$|kX/F_0|$')
	ax1.set_xlabel('$\omega/\omega_n$')	
	plt.legend(loc="upper right")
	ax1.grid()	
	plt.plot(a, np.abs(XXst), label=r'$\xi=$'+str(xi[i]))
	ax2 = plt.subplot(212)
	ax2.set_xlim([0, 5])
	ax2.set_ylabel('\u03B4[rad]')
	ax2.set_xlabel('$\omega/\omega_n$')
	ax2.grid()
	fig = plt.gcf()
	fig.canvas.set_window_title('Ampiezza e fase con forzante generica')
	plt.plot(a, -np.pi*(np.angle(XXst, deg = True)/360))	

	plt.figure(2)
	ax1 = plt.subplot(211)
	ax1.set_xlim([0, 5])
	ax1.set_ylim(0,6)
	ax1.set_ylabel('$|X/Y|=|F_TF_0|$')
	ax1.set_xlabel('$\omega/\omega_n$')	
	ax1.grid()
	plt.legend(loc="upper right")
	plt.plot(a, np.abs(XYO), label=r'$\xi=$'+str(xi[i]))
	ax2 = plt.subplot(212)
	ax2.set_xlim([0, 5])
	ax2.set_ylabel('\u03B4[rad]')
	ax2.set_xlabel('$\omega/\omega_n$')
	fig = plt.gcf()
	fig.canvas.set_window_title('Ampiezza e fase con spostamento di vincolo')
	ax2.grid()
	plt.plot(a, -np.pi*(np.angle(XYO, deg = True)/360))	

plt.show()
