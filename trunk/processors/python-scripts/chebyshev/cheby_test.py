#Quasi_tridiagonal system 

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from scipy.interpolate import splrep, splev


#coe = [3.0, 2.0, 2.0, -1.0]
#coe = [-0.005098767487, -0.005442552725, 0.000024452783, 0.001966230881, 0.001541011815, 0.001375185093, 0.003533302889, 0.002101136751]
#coe = [-0.005098767487, -0.005442552725, 0.000024452783, 0.001966230881, 0.0, 0.0, 0.0, 0.0]
coe_00 = [-0.008218693442, -0.006262365456, 0.003020900858, 0.002780311661, 0.001685770568, 0.001399485736, 0.003512022017, 0.002082568059]
real = []
for x in np.linspace(-1,1,100):
	real.append(np.polynomial.chebyshev.chebval(x,coe_00))

#print y[len(y)-1]
#print len(real)

x = np.linspace(-1,1,100)
tck = splrep(x, real)
Z = splev(x, tck, der = 1)

print "From derivative at 1 =  ",Z[-1]
print "From derivative at -1 =  ",Z[0]

result = 0.0
for i in range(0,8):
	#print i
	result += coe_00[i] * i**2
print "From sumation at 1= ", result

result2 = 0.0
for i in range(0,8):
	result2 += coe_00[i] * (i**2) * ((-1)**(i-1))

print "From sumation at -1= ", result2

'''
coe_01 = [-0.005098767487, -0.005442552725, 0.000024452783, 0.001966230881, 0.001541011815, 0.001375185093, 0.003533302889, 0.002101136751]
real = []
for x in np.linspace(-1,1,100):
	#print "x = ", x, "T(x) =", np.polynomial.chebyshev.chebval(x,coe_01)
	real.append(np.polynomial.chebyshev.chebval(x,coe_01))
#print real

#print len(real) 
x = np.linspace(-1,1,100)
tck = splrep(x, real)
Z = splev(x, tck, der = 1)

print "At -1.0 ", Z[0]
print "At 1.0 ", Z[99]
#plt.plot(x,x*x,'r')
#plt.plot(x,Z,'g')
#plt.show()
result = 0
for i in range(0,7):
	result += coe_01[i] * i**2
print result
'''



