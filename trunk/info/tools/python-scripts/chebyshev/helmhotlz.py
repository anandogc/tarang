#Quasi_tridiagonal system 

from __future__ import division
import numpy as np

A = np.matrix([ [1,1,1,1], [0.0,1.0,1.0,0], [0.0,65.542,-44.694,10.92], [0,0,5.461,-9.738] ])
f = np.matrix([ 0.5,0,0,0]).T

n = len(A) - 1
X = np.zeros(n)
Y = np.zeros(n)
theta = np.zeros(n+1)
lamda = np.zeros(n+1)
THETA  = 0
LAMDA = 0
w = np.zeros(n+1)

#STEP 1
X[n-1] = - A[n,n-1]/A[n,n]
Y[n-1] =  f[n]/A[n,n]

for i in range(n-1,0,-1):
	X[i-1] = -A[i,i-1]/(A[i,i] + A[i,i+1] * X[i])
	Y[i-1] = ( f[i] - (A[i,i+1] * Y[i]) ) / (A[i,i] + A[i,i+1] * X[i])

#STEP 2
theta[0] = 1
lamda[0] = 0
n = len(A) - 1
for i in range(1,n+1):
	theta[i] = X[i-1] * theta[i-1]
	lamda[i] = ( X[i-1] * lamda[i-1] ) + Y[i-1]


#STEP 3
for i in range(0,n+1):
	THETA  = THETA + A[0,:][0,i] * theta[i]
	LAMDA = LAMDA + A[0,:][0,i] * lamda[i]
w[0] = (f[0] - LAMDA) / THETA


#STEP 4
for i in range (0,n):
	w[i+1] = X[i] * w[i] + Y[i]
w = np.matrix(w)
print "A x w = f\n\n", A, "\n\n x \n\n", w.T, "\n\n  = \n\n", f

	
