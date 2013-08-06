#Quasi_tridiagonal system 

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


import scipy as sc
from scipy.interpolate import splrep, splev

Lambda = 17 + 1/0.01
x = np.linspace(-1,1,100)

coe =  [ 0.003332135514, -0.004931730643,  0.001376567737,  0.001412508310, -0.002472520720,  0.002265903610, -0.001598344871,  0.000955314732, -0.000508256581,  0.000245445379, -0.000109483763,  0.000045354576, -0.000017311924]
u = []
for x_i in x:
	#print "x_i = ", x_i, "T(x_i) =", np.polynomial.chebyshev.chebval(x_i,coe_01)
	u.append(np.polynomial.chebyshev.chebval(x_i,coe))
#print u

#Velocity
u = np.array(u)

#print len(u) 
tck = splrep(x, u)
Z = splev(x, tck, der = 2)

#Pressure
coe_p = [ 0.203557675451, -0.353599933525,  0.235684128474, -0.124832786853,  0.054120995507, -0.019767538741,  0.006202811572, -0.001705234158,  0.000415727245, -0.000091129103,  0.000018098872]
p = []
for x_i in x:
	#print "x = ", x, "T(x) =", np.polynomial.chebyshev.chebval(x,coe_01)
	p.append(np.polynomial.chebyshev.chebval(x_i,coe_p))

tck = splrep(x, p)
grad_p = splev(x, tck, der = 1)


# print len(Z)
# print len(u)
# print len(4*u)

g=-Z + Lambda*u + grad_p
print u[0], u[-1]
print g[0], g[-1]
plt.plot(x,g, 'b')
plt.plot(x,u, 'g')
plt.plot(x,p, 'r')

#print "At -1.0 ", Z[0]
#print "At 1.0 ", Z[99]
# plt.plot(x,x*x,'r')
##plt.plot(x,Z,'g')
plt.show()
#result = 0
#for i in range(0,2):
#	result += coe_01[i] * i**2
#print result




