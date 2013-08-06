#Chebyshev Transoformation


from __future__ import division
import numpy as np


x = np.linspace(0,np.pi,8)
y = np.cos(x)
print y
n=len(y)
f=np.zeros(n)

"""
#Forward Cosine
for k in range(1,n):
	f[k]=y[0]+(-1)**k * y[n-1]
	for j in range(1,n-1):
		f[k] = f[k] + 2.0 * y[j] * np.cos((np.pi/(n-1)) * k *j)
for j in range(0,n):
	f[0] = f[0] + y[j]

print f
"""

print "n = " , n

#Forward Chebyshev
"""
for k in range(1,n):
	f[k] = 2.0 * y[0] + 2.0 * (-1)**k * y[n-1]
	for j in range(1,n-1):
		f[k] = f[k] + 2.0 * y[j] * np.cos((np.pi/(n-1)) * k *j)
for j in range(0,n):
	f[0] = f[0] + y[j]
"""
"""
for k in range(0,n):
	for j in range(0,n):
		f[k] = f[k] + y[j] * np.cos((np.pi/(n-1)) * k *j)
"""
k=3
f=0
for j in range(0,n):
	f = f + np.cos(x[j]) * np.cos(k*x[j])

print f
