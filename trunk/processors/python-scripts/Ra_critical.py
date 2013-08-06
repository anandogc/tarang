from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math as m

c3=np.sqrt(3.0)

def q0(a,t):
	return a*np.sqrt(t-1)

def q1(a,t):
	return a*np.sqrt(0.5*np.sqrt(1+t+t**2) + 0.5*(1+t/2))

def q2(a,t):
	return a*np.sqrt(0.5*np.sqrt(1+t+t**2) - 0.5*(1+t/2))

a=3.117

t=np.linspace(3.5,4,100)
f = []
for i in t:
	f.append(q0(a,i)*m.tan(q0(a,i))+((q1(a,i)+q2(a,i)*c3)*m.sinh(2*q1(a,i))+(q1(a,i)*c3-q2(a,i))*m.sin(2*q2(a,i)))/(m.cosh(2*q1(a,i))+m.cos(2*q2(a,i))))
#	f.append(q0(a,i)*m.tan(q0(a,i)/2)+((q1(a,i)+q2(a,i)*c3)*m.sinh(q1(a,i))+(q1(a,i)*c3-q2(a,i))*m.sin(q2(a,i)))/(m.cosh(q1(a,i))+m.cos(q2(a,i))))

plt.plot(t,f)

# plt.ylim(-0.1,0.1)

tau=3.7

Ra=tau**3*a**4
print Ra
plt.show()
