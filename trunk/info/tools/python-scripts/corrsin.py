"""To run the following script, do the following
1) python
2) %run corrsin.py
"""

from scipy.integrate import quad
import numpy
import pylab

a = 1.0
q = 1.5
b = 1.0
nu = 0.1
energy_supply_rate = 1.0

lower_limit = 1
upper_limit = 1000
dns_upper_limit = 32

L = 2*numpy.pi

def E(k, a, q, b):
	return a*k**4*numpy.exp(-b*k**1.1)/(k**4 + q**4)**(1 + 2.8/12)

def kE(k, a, q, b):
	return k*E(k,a,q,b)

def D(k, a, q, b):
	return 2*nu*k**2*E(k,a,q,b)

Etot = quad(E, lower_limit, upper_limit, args=(a, q, b,))
Dtot = quad(D, lower_limit,  upper_limit, args=(a, q, b,))

print "Etot = ",  Etot[0]
# print "Error_ Etot = ",  Etot[1]
print

print "Dtot = ",  Dtot[0]
# print "Error_Dtot = ",  Dtot[1]
print


dns_Etot = quad(E, lower_limit, dns_upper_limit, args=(a, q, b,))
dns_Dtot = quad(D, lower_limit, dns_upper_limit, args=(a, q, b,))
dns_kEtot = quad(kE, lower_limit, dns_upper_limit, args=(a, q, b,))

print "dns_Etot = ", dns_Etot[0]
# print "Error_dns_Etot = ", dns_Etot[1]
print

print "dns_Dtot = ", dns_Dtot[0]
print

eta = ((nu**3)/Dtot[0])**0.25
kmaxeta = eta*dns_upper_limit
print  "eta, kmaxeta = ", eta, kmaxeta


eta = ((nu**3)/energy_supply_rate)**0.25
kmaxeta = eta*dns_upper_limit
print  "Based on input_energy supply rate: eta, kmaxeta = ", eta, kmaxeta
print
print

Urms = (2.0*Etot[0]/3.0)**0.5
lam = (15.0*nu*Urms**2/Etot[0])**0.5
Re = Urms*L/nu
Rlambda = Urms*lam/nu
avg_k = dns_kEtot[0]/dns_Etot[0]

print "Urms = ", Urms
print "lambda = ", lam
print "Re = ", Re
print "Rlambda = ", Rlambda
print "avg_k = ", avg_k


k =  numpy.arange(lower_limit, dns_upper_limit, 0.1)
pylab.loglog(k, E(k,a,q,b))
pylab.show()
