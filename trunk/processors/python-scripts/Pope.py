"""To run the following script, do the following
1) python
2) %run corrsin.py
"""

from __future__ import division
import numpy as np
pi = np.pi
from scipy.integrate import quadrature
from scipy.integrate import quad
from scipy.interpolate import splrep, splev
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 3

#font = {'family' : 'serif', 'weight' : 'bold', 'size' : 24}
#rc('font', **font)


N = 256
kmax = N/2
nu = 1
epsilon = 1
eta = ((nu**3)/epsilon)**0.25

lower_limit = eta
upper_limit = 5
dns_lower_limit = eta
dns_upper_limit = eta*kmax
L = 2*np.pi

print "nu = ", nu
print "Grid size N = ", N
print "Epsilon (energy supply rate) = ",  epsilon
print "dns_upper_limit = ",  dns_upper_limit

print  "Based on input_energy supply rate: eta, kmaxeta = ", eta, kmax*eta
print
print

# Formula from Pope page 232
beta = 5.2
ceta = 0.4
cL = 6.78
p0 = 2
CKolm = 1.6


def fL(keta, cL, p0):
#    kL = keta*L/eta
    cLprime = cL*(eta/L)
    return (keta/np.sqrt(keta**2 + cLprime))**(5.0/3+p0)

def feta(keta, beta, ceta):
    return np.exp(-beta*(((keta)**4+ceta**4)**0.25-ceta))

def E(keta, cL, p0, beta, ceta):
	return keta**(-5.0/3)*fL(keta,cL, p0)*feta(keta, beta, ceta)

def kE(keta, cL, p0, beta, ceta):
	return keta*E(keta, cL, p0, beta, ceta)

def D(keta, cL, p0, beta, ceta):
	return keta**(1.0/3)*E(keta, cL, p0, beta, ceta)
	

tolerance = 1e-6
max_iter = 100

Etot = quad(E, lower_limit, upper_limit, args=(cL, p0, beta, ceta,))
Dtot = quad(D, lower_limit,  upper_limit, args=(cL, p0, beta, ceta,))

etastar = (2*CKolm*Dtot[0])**(3.0/4)
# eta = etastar(nu/epsilon**3)^(1/4)
print "eta_star = ", etastar

Etotal = CKolm*np.sqrt(epsilon*eta)*etastar**(2.0/3)*Etot[0]
# division by eta due integral wrt keta
print "Etotal = ",  Etotal
# print "Error_ Etot = ",  Etot[1]
print

Dtotal = 2*CKolm*epsilon*etastar**(-4.0/3)*Dtot[0]
print "Dtotal = ",  Dtotal
# print "Error_Dtot = ",  Dtot[1]/eta
print


dns_Etot = quad(E, dns_lower_limit, dns_upper_limit, args=(cL, p0, beta, ceta,))
dns_Dtot = quad(D, dns_lower_limit, dns_upper_limit, args=(cL, p0, beta, ceta,))
dns_kEtot = quad(kE, dns_lower_limit, dns_upper_limit, args=(cL, p0, beta, ceta,))

dns_Etotal = CKolm*np.sqrt(epsilon*eta)*etastar**(2.0/3)*dns_Etot[0]
# division by eta due integral wrt keta
print "dns_Etotal = ",  dns_Etotal
# print "Error_dns_Etot = ", dns_Etot[1]
print

dns_Dtotal = 2*CKolm*epsilon*etastar**(-4.0/3)*dns_Dtot[0]
print "dns_Dtotal = ",  dns_Dtotal
print

eta1 = ((nu**3)/dns_Dtotal)**0.25
print  "Based on DNS integral estimate: eta, kmaxeta = ", eta1, kmax*eta1


Urms = (2.0*Etotal/3.0)**0.5
lam = (15.0*nu*Urms**2/Etotal)**0.5
Re = Urms*L/nu
Rlambda = Urms*lam/nu
avg_k = dns_kEtot[0]/dns_Etot[0]

print "Urms = ", Urms
print "lambda = ", lam
print "Re = ", Re
print "Rlambda = ", Rlambda
print "avg_k = ", avg_k
dx = 2*np.pi/N
print "CFL_dt (Courant no=0.5) = ", 0.5*dx/Urms


keta =  np.arange(lower_limit, upper_limit, 0.1*eta)
k = keta/eta
ek = E(keta,cL, p0, beta, ceta)
fL_val = fL(keta, cL, p0)
feta_val = feta(keta, beta, ceta)
d_val = D(keta, cL, p0, beta, ceta)
#loglog(k, ek)
#show()
