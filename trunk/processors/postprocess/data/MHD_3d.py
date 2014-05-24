from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import tarang as t
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

pi = np.pi
t.init("/home/abhishek/Dropbox/my_tarang/sample/mhd/3d")
font = {'family' : 'serif', 'weight' : 'normal', 'size' : 15}
plt.rc('font', **font)

if t.para.program.basis_type == 'FFF':
	fff_slab = t.FFF_SLAB()
	max_radius_inside = fff_slab.Max_radius_inside(t.N)
elif t.para.program.basis_type == 'SFF':
	sff_slab = t.SFF_SLAB()
	max_radius_inside = sff_slab.Max_radius_inside(t.N)


########################################################################


################# glob #########################

"""
fig, axes = plt.subplots(figsize=(7,5))
axes.plot(t.time,t.U.cvf.total_energy,'r',lw = 3,label=r'$E_u$')
axes.plot(t.time,t.B.cvf.total_energy,'g',lw = 3,label=r'$E_{B}$')
axes.legend(loc=0,prop={'size':15})
axes.set_xlabel(r'$time$')
fig.tight_layout()
plt.show()
"""
#################### Spectrum #################################
"""
t.spectrum_showtime()
t.spectrum_read(3)


fig, axes = plt.subplots(figsize=(7,5))
axes.loglog(t.k, (t.U.cvf.Uek_shell_ek1[0]+t.U.cvf.Uek_shell_ek2[0]+t.U.cvf.Uek_shell_ek3[0]) ,'r',lw = 4,label=r'$E_u(k)$')
axes.loglog(t.k, (t.B.cvf.Uek_shell_ek1[0]+t.B.cvf.Uek_shell_ek2[0]+t.B.cvf.Uek_shell_ek3[0]) ,'b',lw = 4,label=r'$E_B(k)$')

axes.set_xlabel(r'$k$')
axes.legend(loc=0,prop={'size':15})
#axes.axis([2.5, 350, 10**-30, 10**1])
plt.show()
"""

#################### Flux ###############################
"""
t.flux_showtime()
t.flux_read(1)
kf = t.flux_radii(max_radius_inside)
print t.energytr.U2U[0][0:len(kf)]

#fig, axes = plt.subplots(figsize=(7,5))
#axes.loglog(kf, t.energytr.U2U[0][0:len(kf)] ,'r',lw = 3,label=r'$\Pi_u(k)$')
#axes.loglog(kf, t.energytr.T2T[0][0:len(kf)] ,'g',lw = 3,label=r'$\Pi_{\theta}(k)$')
#axes.set_xlabel(r'$k$')
#axes.legend(loc=0,prop={'size':15})
#plt.show()
"""
################## Fourier Probes ######################
"""
def Mag(real,imag):
	return ((real**2 + imag**2)**0.5)

t.field_k_read(1)
t.field_k_showmode()

mode1 = t.field_k_readmode(1,0,1)
time1 = np.array(mode1.Time)
ux_real1 = np.array(mode1.Ux_real)
ux_imag1 = np.array(mode1.Ux_imag)
uy_real1 = np.array(mode1.Uy_real)
uy_imag1 = np.array(mode1.Uy_imag)
uz_real1 = np.array(mode1.Uz_real)
uz_imag1 = np.array(mode1.Uz_imag)

Bx_real1 = np.array(mode1.Bx_real)
Bx_imag1 = np.array(mode1.Bx_imag)
By_real1 = np.array(mode1.By_real)
By_imag1 = np.array(mode1.By_imag)
Bz_real1 = np.array(mode1.Bz_real)
Bz_imag1 = np.array(mode1.Bz_imag)





mode2 = t.field_k_readmode(1,1,1)
time2 = np.array(mode2.Time)
ux_real2 = np.array(mode2.Ux_real)
ux_imag2 = np.array(mode2.Ux_imag)
uy_real2 = np.array(mode2.Uy_real)
uy_imag2 = np.array(mode2.Uy_imag)
uz_real2 = np.array(mode2.Uz_real)
uz_imag2 = np.array(mode2.Uz_imag)

Bx_real2 = np.array(mode2.Bx_real)
Bx_imag2 = np.array(mode2.Bx_imag)
By_real2 = np.array(mode2.By_real)
By_imag2 = np.array(mode2.By_imag)
Bz_real2 = np.array(mode2.Bz_real)
Bz_imag2 = np.array(mode2.Bz_imag)
"""



################### Real Probes ###########################
"""
t.field_r_read(1)
t.field_r_showmode()


mode1 = t.field_r_readmode(1,0,1)
time1 = np.array(mode1.Time)
ux1 = np.array(mode1.Ux)
uy1 = np.array(mode1.Uy)
uz1 = np.array(mode1.Uz)
Bx1 = np.array(mode1.Bx)
By1 = np.array(mode1.By)
Bz1 = np.array(mode1.Bz)
"""





###################### Shell-to-shell #########################
"""
t.shell_showtime()

u2u = t.shell_to_shell_read(1, "U2U")
u2u  = np.array(u2u)


(y,x) = u2u.shape


X = np.arange(0, x+1, 1)
Y = np.arange(0, y+1, 1)
n, m = np.meshgrid(X, Y)

font = {'family' : 'serif', 'weight' : 'normal', 'size' : 10}
plt.rc('font', **font)
fig, axes = plt.subplots(figsize=(5,5))
axes.axis([0, 19, 0, 19])
axes.xaxis.set_major_locator(MultipleLocator(5))
axes.yaxis.set_major_locator(MultipleLocator(5))
p = axes.pcolor(n,m,u2u)
axes.set_aspect(1)
axes.set_xlabel(r'$n$',fontsize = 28)
axes.set_ylabel(r'$m$',fontsize = 28)
fig.tight_layout()
fig.colorbar(p, ax = axes, shrink = 0.8)
plt.title('Shell to shell energy transfer for U2U')
plt.show()
"""

###################### Ring Spectrum ################################
"""
t.ring_spectrum_showtime()
data = t.ring_spectrum_read(6,'Uek')
data = np.array(data)
print data
"""
######################################################################
