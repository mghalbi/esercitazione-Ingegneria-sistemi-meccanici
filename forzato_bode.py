# Author: MohsenGhalbi 
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
	XXst = list(map(lambda x: 1/(1-pow(x,2)+2*xi[i]*complex(0,1)*x),a));
	XYO = list(map(lambda x: (1+2*xi[i]*complex(0,1)*x)/(1-pow(x,2)+2*xi[i]*complex(0,1)*x),a));
	plt.figure(1)
	plt.subplot(211)
	plt.plot(a, np.abs(XXst))
	plt.subplot(212)
	plt.plot(a, np.angle(XXst, deg = True))	

	plt.figure(2)
	plt.subplot(211)
	plt.plot(a, np.abs(XYO))
	plt.subplot(212)
	plt.plot(a, np.angle(XYO, deg = True))	

plt.show()
