#Chebyshev Transoformation


from __future__ import division
import numpy as np

k101 = np.array([1,0,1])
km101 = np.array([-1,0,1])
k002 = np.array([0,0,2])
k200 = np.array([2,0,0])


k011 = np.array([0,1,1])
k112 = np.array([1,1,2])
k213 = np.array([2,1,3])
k123 = np.array([1,2,3])
k1m10 = np.array([1,-1,0])
km110 = np.array([-1,1,0])

u101 = np.array([-1j,0,1j])
um101 = np.array([1j,0,1j])
u10m1 = np.conj(um101)

u011 = np.array([0,1,-1])
u112 = np.array([1,-1,0])
u213 = np.array([0,0,0])
u123 = np.array([0,0,0])

def N(k,up,uq):
    return 1j*(np.dot(k,uq)*up+np.dot(k,up)*uq)

def pressure(k,Nk):
    return 1j*np.dot(k,Nk)/np.dot(k,k)

N112 = N(k112, u101, u011)
N101 = N(k101, u112, u011)
N011 = N(k011, u112, u101)
N213 = N(k213, u101, u112)
N123 = N(k123, u011, u112)
N1m10 = N(k1m10, u101, u011)
Nm110 = N(km110, u011, u101)


print "u(101) = ", u101
print "u(-101) = ", um101
print "N(200) = ", N(k200, u101, um101)
print "N(002) = ", N(k002, u101, u10m1)

#
#print "N(112) = ", N(k112, u101, u011)
#
#print "N(101) = ", N(k101, u112, u011)
#
#print "N(011) = ", N(k011, u112, u101)
#
#print "N(213) = ", N(k213, u101, u112)
#
#print "N(123) = ", N(k123, u011, u112)
#
#print "N(1,-1,0) = ", N(k1m10, u101, u011)
#
#print "N(-1,1,0) = ", N(km110, u011, u101)
#
#u_dot_N = np.dot(N112,u112) + np.dot(N101,u101) + np.dot(N011,u011)
#
#print "pressure_101 = ", pressure(k101, N101)
#
#print "pressure_011 = ", pressure(k011, N011)
#
#print "pressure_112 = ",  pressure(k112, N112)
#
#print "pressure_213 = ", pressure(k213, N213)
#
#print "pressure_123 = ", pressure(k123, N123)
#
#print "pressure_1m10 = ", pressure(k1m10, N1m10)
#
#print "pressure_m110 = ", pressure(km110, Nm110)
#
#print "Nlin after additing pressure :"
#
#N101 = -N(k101, u112, u011) - 1j*k101*pressure(k101, N101)
#
#N011 = -N(k011, u112, u101) -1j*k011*pressure(k011, N011)
#
#N112 = -N(k112, u101, u011) -1j*k112*pressure(k112, N112)
#
#N213 = -N(k213, u101, u112) -1j*k213*pressure(k213, N213)
#
#N123 = -N(k123, u011, u112) -1j*k123*pressure(k123, N123)
#
#N1m10 = -N(k1m10, u101, u011) -1j*k1m10*pressure(k1m10, N1m10)
#
#Nm110 = -N(km110, u011, u101) -1j*km110*pressure(km110, Nm110)
#
#print "N(101) = ", -N(k101, u112, u011) - 1j*k101*pressure(k101, N101)
#
#print "N(011) = ", -N(k011, u112, u101) -1j*k011*pressure(k011, N011)
#
#print "N(112) = ", -N(k112, u101, u011) -1j*k112*pressure(k112, N112)
#
#print "N(213) = ", -N(k213, u101, u112) -1j*k213*pressure(k213, N213)
#
#print "N(123) = ", -N(k123, u011, u112) -1j*k123*pressure(k123, N123)
#
#print "N(1,-1,0) = ", -N(k1m10, u101, u011) -1j*k1m10*pressure(k1m10, N1m10)
#
#print "N(-1,1,0) = ", -N(km110, u011, u101) -1j*km110*pressure(km110, Nm110)
#
#
#print "sum(|N|^2) = ", np.dot(N101,np.conj(N101)) + np.dot(N011,np.conj(N011))+ np.dot(N112,np.conj(N112)) + np.dot(N213,np.conj(N213)) + np.dot(N123,np.conj(N123)) + np.dot(N1m10,np.conj(N1m10)) + np.dot(Nm110,np.conj(N1m10))
#








