
#Quasi_tridiagonal system 

from __future__ import division
import numpy as np

A = np.matrix([ [1,1,1,1], [1,2,3,0], [0,1,2,1], [0,0,3,3] ])

a = [1, 1, 3]
b = [2, 2, 3]
c = [3, 1, 0]

d = np.matrix([4, 1, 2, 3]).T

n = len(d) - 1
X = np.zeros(n)
Y = np.zeros(n)
theta = np.zeros(n+1)
lamda = np.zeros(n+1)
THETA  = 0
LAMDA = 0
w = np.zeros(n+1)

#STEP 1
X[n-1] = - a[n-1]/b[n-1]
Y[n-1] =  d[n]/b[n-1]

for i in range(n-1,0,-1):
	X[i-1] = -a[i-1]/(b[i-1] + c[i-1] * X[i])
	Y[i-1] = ( d[i] - (c[i-1] * Y[i]) ) / (b[i-1] + c[i-1] * X[i])

#STEP 2
theta[0] = 1
lamda[0] = 0
n = len(A) - 1
for i in range(1,n+1):
	theta[i] = X[i-1] * theta[i-1]
	lamda[i] = ( X[i-1] * lamda[i-1] ) + Y[i-1]


#STEP 3
for i in range(0,n+1):
	THETA  = THETA + theta[i]
	LAMDA = LAMDA + lamda[i]
w[0] = (d[0] - LAMDA) / THETA


#STEP 4
for i in range (0,n):
	w[i+1] = X[i] * w[i] + Y[i]
w = np.matrix(w)
print "A x w = f\n\n", A, "\n\n x \n\n", w.T, "\n\n  = \n\n", A*w.T

	