from __future__ import division
import tarang as t
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

pi = np.pi
t.init("/Users/abhishek/Dropbox/my_tarang/sample/fluid/3d")
font = {'family' : 'serif', 'weight' : 'normal', 'size' : 15}
plt.rc('font', **font)


###################################################################
if t.para.program.basis_type == 'FFF':
	fff_slab = t.FFF_SLAB()
	max_radius_inside = fff_slab.Max_radius_inside(t.N)
elif t.para.program.basis_type == 'SFF':
	sff_slab = t.SFF_SLAB()
	max_radius_inside = sff_slab.Max_radius_inside(t.N)
###################################################################


################# glob #########################

"""
fig, axes = plt.subplots(figsize=(8,6))
axes.plot(t.time,t.U.cvf.total_energy,'r',lw=3,label=r'$E_u$')
axes.plot(t.time,t.U.cvf.total_dissipation,'g',lw=3,label=r'$\epsilon_u$')
#axes.plot(t.time,t.U.cvf.dissipation_coefficient,'b',label=r'$\nu$')
#axes.plot(t.time,t.Rlambda,'b',label=r'$Re_{\lambda}$')
#axes.plot(t.time,t.time_dt,'b',label=r'$dt$')
#axes.plot(t.time,t.U.cvf.total_E1,'g',label=r'$E_x$')
#axes.plot(t.time,t.U.cvf.total_E2,'b',label=r'$E_y$')
#axes.plot(t.time,t.U.cvf.total_E3,'c',label=r'$E_z$')
axes.set_xlabel(r'$time$',fontsize = 28)
axes.legend(loc=0,prop={'size':28})
fig.tight_layout()
plt.show()
"""
#################### Spectrum #################################

"""
t.spectrum_showtime()
t.spectrum_read(5,6)

fig, axes = plt.subplots(figsize=(8,6))
axes.loglog(t.k, (t.U.cvf.Uek_shell_ek1[0]+t.U.cvf.Uek_shell_ek2[0]+t.U.cvf.Uek_shell_ek3[0]), 'r',lw = 4)#,label=r'$E_u(k)$')
axes.loglog(t.k, (t.U.cvf.Uek_shell_ek1[1]+t.U.cvf.Uek_shell_ek2[1]+t.U.cvf.Uek_shell_ek3[1]), 'g',lw = 4)#,label=r'$E_u(k)$')

#axes.loglog(t.k, (t.U.cvf.UDk_shell_dissk1[0]+t.U.cvf.UDk_shell_dissk2[0]+t.U.cvf.UDk_shell_dissk3[0]),'b',lw=3, label=r'$D_u(k)$')
#axes.loglog(t.k, (t.U.cvf.Fv_v_shell_ek1[0]+t.U.cvf.Fv_v_shell_ek2[0]+t.U.cvf.Fv_v_shell_ek3[0]) ,'m',label=r'$F_u(k)$')
axes.set_xlabel(r'$k$',fontsize = 28)
axes.set_ylabel(r'$E_u(k)$',fontsize = 28)
axes.legend(loc=0,prop={'size':28})
fig.tight_layout()
plt.show()
"""

#################### Flux ###############################

"""
t.flux_showtime()
t.flux_read(5,6)
kf = t.flux_radii(max_radius_inside)
fig, axes = plt.subplots(figsize=(7,5))
axes.loglog(kf, t.energytr.U2U[0][0:len(kf)] ,'r',lw = 3)#, label=r'$\Pi_u(k)$')
axes.loglog(kf, t.energytr.U2U[1][0:len(kf)] ,'g',lw = 3)#, label=r'$\Pi_u(k)$')
axes.set_xlabel(r'$k$')
axes.set_ylabel(r'$\Pi_u(k)$')
#axes.legend(loc=0,prop={'size':35})
fig.tight_layout()
plt.show()
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

ux1 = Mag(ux_real1,ux_imag1)
uy1 = Mag(uy_real1,uy_imag1)
uz1 = Mag(uz_real1,uz_imag1)




mode2 = t.field_k_readmode(1,1,0)
time2 = np.array(mode2.Time)
ux_real2 = np.array(mode2.Ux_real)
ux_imag2 = np.array(mode2.Ux_imag)
uy_real2 = np.array(mode2.Uy_real)
uy_imag2 = np.array(mode2.Uy_imag)
uz_real2 = np.array(mode2.Uz_real)
uz_imag2 = np.array(mode2.Uz_imag)
theta_real2 = np.array(mode2.T_real)
theta_imag2 = np.array(mode2.T_imag)
ux2 = Mag(ux_real2,ux_imag2)
uy2 = Mag(uy_real2,uy_imag2)
uz2 = Mag(uz_real2,uz_imag2)



fig, axes = plt.subplots(figsize=(7,5))
axes.plot(time1, ux1 ,'r',label=r'$1 \, 0 \, 1$')
axes.plot(time2, ux2 ,'g',label=r'$1 \, 1 \, 0$')
axes.set_xlabel(r'$t$',fontsize = 20)
axes.set_ylabel(r'$|u_x(t)|$',fontsize = 20)
axes.legend(loc=0,prop={'size':15})
fig.tight_layout()
#plt.show()
"""
################### Real Probes ###########################
"""
t.field_r_read(1)
t.field_r_showmode()

mode1 = t.field_r_readmode(1,1,1)
time1 = np.array(mode1.Time)
ux1 = np.array(mode1.Ux)
uy1 = np.array(mode1.Uy)
uz1 = np.array(mode1.Uz)



mode2 = t.field_r_readmode(1,62,62)
time2 = np.array(mode2.Time)
ux2 = np.array(mode2.Ux)
uy2 = np.array(mode2.Uy)
uz2 = np.array(mode2.Uz)



fig, axes = plt.subplots(figsize=(7,5))
axes.plot(time1, ux1 ,'r',label=r'$1 \, 1 \, 1$')
axes.plot(time2, ux2 ,'g',label=r'$1 \, 62 \, 62$')

axes.set_xlabel(r'$t$',fontsize = 20)
axes.set_ylabel(r'$u_x(t)$',fontsize = 20)
axes.legend(loc=0,prop={'size':15})
fig.tight_layout()
plt.show()
"""
###################### Shell-to-shell #########################
"""
t.shell_showtime()
u2u = t.shell_to_shell_read(5, 'U2U')
u2u  = np.array(u2u)
(y,x) = u2u.shape


X = np.arange(0, x+1, 1)
Y = np.arange(0, y+1, 1)
n, m = np.meshgrid(X, Y)

font = {'family' : 'serif', 'weight' : 'normal', 'size' : 10}
plt.rc('font', **font)
fig, axes = plt.subplots(figsize=(8,8))
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
###############################################################

###################### Ring Spectrum ################################
"""
#t.ring_spectrum_showtime()
#data = t.ring_spectrum_read(5,'Uek')
#data = np.array(data)

data = np.loadtxt('/Users/abhishek/Dropbox/my_tarang/sample/fluid/3d/ring_specctrum_jd_80.d',comments='%%')
(m,n) = data.shape

r = np.matrix( (np.arange(1,m+1,1))/m )
theta = np.matrix( (0.5 * pi * np.arange(n-1,-1,-1))/(n-1))

X =  r.T * np.cos(theta)
Y =  r.T * np.sin(theta)

theta2 = (0.5 * pi *  np.arange(0,n+1,1) )/n

C = data
for i in np.arange(0,n,1):
	C[:,i] = C[:,i]/(np.abs( np.cos(theta2[i+1]) - np.cos(theta2[i]) ))



font = {'family' : 'serif', 'weight' : 'normal', 'size' : 10}
plt.rc('font', **font)

C = np.log10(np.abs(C))

fig, axes = plt.subplots(1,2,figsize=(7,5))
Z=np.log10(np.abs(C))
p1 = axes[0].pcolormesh(X,Y,Z)
fig.colorbar(p1, ax = axes[0], shrink = 0.5)
axes[0].set_aspect(1)
#p2 = axes[1].contour(X,Y,Z,15,lw=1.5)
axes[1].set_aspect(1)
#fig.colorbar(p2, ax = axes[1], shrink = 0.5)
#fig.tight_layout()
plt.show()
"""
#####################################################################



################## Profile ###########################################
"""
data = t.profile()

x = data[:,0]
v1_avg = data[:,1]
v2_avg = data[:,2]
v3_avg = data[:,3]

v1_rms = data[:,4]
v2_rms = data[:,5]
v3_rms = data[:,6]



fig, axes = plt.subplots(figsize=(7,7))
axes.plot(v2_avg,x,'g',lw = 3,label=r'$<U_y(z)>$')
axes.legend(loc=0,prop={'size':15})
axes.set_ylabel(r'$z$',fontsize = 28)
fig.tight_layout()
plt.show()
"""
####################################################################

################## Visualizations ##############################
################################################################
