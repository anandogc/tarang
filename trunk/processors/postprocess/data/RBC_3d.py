from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import tarang as t
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import mayavi.mlab as maya

pi = np.pi
t.init("/home/abhishek/Dropbox/my_tarang/sample/RBC/3d")
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
axes.plot(t.time,t.T.csf.total_energy,'g',lw = 3,label=r'$E_{\theta}$')
#axes.plot(t.time,t.U.cvf.total_dissipation,'r',label=r'$\epsilon_u$')
#axes.plot(t.time,t.T.csf.total_dissipation,'r',label=r'$\epsilon_{\theta}$')
#axes.plot(t.time,t.correlation.nusselt_no,'b',label=r'$Nu$')
#axes.plot(t.time,t.U.cvf.dissipation_coefficient,'r',label=r'$\nu$')
#axes.plot(t.time,t.T.csf.diffusion_coefficient,'b',label=r'$k$')
#axes.plot(t.time,t.Rlambda,'b',label=r'$Re_{\lambda}$')
#axes.plot(t.time,t.time_dt,'b',label=r'$dt$')
#axes.plot(t.time,t.U.cvf.total_E1,'g',label=r'$E_x$')
#axes.plot(t.time,t.U.cvf.total_E2,'b',label=r'$E_y$')
#axes.plot(t.time,t.U.cvf.total_E3,'c',label=r'$E_z$')
axes.legend(loc=0,prop={'size':15})
axes.set_xlabel(r'$time$')
fig.tight_layout()
plt.show()
"""
#################### Spectrum #################################
"""
t.spectrum_showtime()
t.spectrum_read(11)


fig, axes = plt.subplots(figsize=(7,5))
axes.loglog(t.k, (t.U.cvf.Uek_shell_ek1[0]+t.U.cvf.Uek_shell_ek2[0]+t.U.cvf.Uek_shell_ek3[0]) ,'r',lw = 4,label=r'$E_u(k)$')
axes.loglog(t.k, t.T.csf.Tek_shell_ek[0] ,'g',lw = 4,label=r'$E_{\theta}(k)$')
#axes.loglog(t.k, (t.U.cvf.UDk_shell_dissk1[0]+t.U.cvf.UDk_shell_dissk2[0]+t.U.cvf.UDk_shell_dissk3[0]) ,'b',label=r'$D_u(k)$')
#axes.loglog(t.k, (t.correlation.U_T_shell_ek1[0]+t.correlation.U_T_shell_ek2[0]+t.correlation.U_T_shell_ek3[0]) ,'c',label=r'$U(k)\theta(k)$')
#axes.loglog(t.k, (t.U.cvf.Fv_v_shell_ek1[0]+t.U.cvf.Fv_v_shell_ek2[0]+t.U.cvf.Fv_v_shell_ek3[0]) ,'m',label=r'$F_u(k)$')
#axes.loglog(t.k, t.T.csf.FT_T_shell_ek[0] ,'r--',label=r'$F_{\theta}(k)$')
axes.set_xlabel(r'$k$')
axes.legend(loc=0,prop={'size':15})
#axes.axis([2.5, 350, 10**-30, 10**1])
plt.show()
"""

#################### Flux ###############################


"""
t.flux_showtime()
t.flux_read(11)
kf = t.flux_radii(max_radius_inside)
fig, axes = plt.subplots(figsize=(7,5))
axes.loglog(kf, t.energytr.U2U[0][0:len(kf)] ,'r',lw = 3,label=r'$\Pi_u(k)$')
axes.loglog(kf, t.energytr.T2T[0][0:len(kf)] ,'g',lw = 3,label=r'$\Pi_{\theta}(k)$')
axes.set_xlabel(r'$k$')
axes.legend(loc=0,prop={'size':15})
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
theta_real1 = np.array(mode1.T_real)
theta_imag1 = np.array(mode1.T_imag)
ux1 = Mag(ux_real1,ux_imag1)
uy1 = Mag(uy_real1,uy_imag1)
uz1 = Mag(uz_real1,uz_imag1)
theta1 = Mag(theta_real1,theta_imag1)



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
theta2 = Mag(theta_real2,theta_imag2)



fig, axes = plt.subplots(figsize=(7,5))
axes.plot(time1, theta1 ,'r',label=r'$1 \, 0 \, 1$')
axes.plot(time2, ux2 ,'g',label=r'$1 \, 1 \, 0$')
axes.set_xlabel(r'$t$',fontsize = 20)
axes.set_ylabel(r'$|u_x(t)|$',fontsize = 20)
axes.legend(loc=0,prop={'size':15})
fig.tight_layout()
plt.show()
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
theta1 = np.array(mode1.T)


mode2 = t.field_r_readmode(1,1,32)
time2 = np.array(mode2.Time)
ux2 = np.array(mode2.Ux)
uy2 = np.array(mode2.Uy)
uz2 = np.array(mode2.Uz)
theta2 = np.array(mode2.T)




fig, axes = plt.subplots(figsize=(7,5))
axes.plot(time1, ux1 ,'r',label=r'$1 \, 1 \, 1$')
axes.plot(time2, ux2 ,'g',label=r'$1 \, 1 \, 32$')

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

T2T = t.shell_to_shell_read(5, 'T2T')
T2T  = np.array(T2T)
(y,x) = T2T.shape


X = np.arange(0, x+1, 1)
Y = np.arange(0, y+1, 1)
n, m = np.meshgrid(X, Y)

font = {'family' : 'serif', 'weight' : 'normal', 'size' : 10}
plt.rc('font', **font)
fig, axes = plt.subplots(figsize=(5,5))
axes.axis([0, 19, 0, 19])
axes.xaxis.set_major_locator(MultipleLocator(5))
axes.yaxis.set_major_locator(MultipleLocator(5))
p = axes.pcolor(n,m,T2T)
axes.set_aspect(1)
axes.set_xlabel(r'$n$',fontsize = 28)
axes.set_ylabel(r'$m$',fontsize = 28)
fig.tight_layout()
fig.colorbar(p, ax = axes, shrink = 0.8)
plt.title('Shell to shell energy transfer for T2T')
plt.show()
"""

################## Profile ###########################################
"""
data = t.profile()

x = data[:,0]
v1_avg = data[:,1]
v2_avg = data[:,2]
v3_avg = data[:,3]
theta_avg = data[:,4]
v1_rms = data[:,5]
v2_rms = data[:,6]
v3_rms = data[:,7]
theta_rms = data[:,8]
T = 1 - x + theta_avg

fig, axes = plt.subplots(figsize=(7,7))
axes.plot(theta_avg,x,'r',lw = 3,label=r'$<\theta(z)>$')
axes.plot(T,x,'g',lw = 3,label=r'$<T(z)>$')
axes.legend(loc=0,prop={'size':28})
axes.set_ylabel(r'$z$',fontsize = 28)
fig.tight_layout()
plt.show()
"""
####################################################################

################## Visualizations ##############################

X, Y, Z, V1, V2, V3, theta, T = t.visual()



#maya.figure(size=(800, 600)) # For Mac
maya.figure(bgcolor=(1.0, 1.0, 1.0),size=(800, 600))
maya.contour3d(T, colormap='jet', contours = 4, name = 'Temperature', vmax = 0.7, vmin = 0.3)
maya.colorbar(title='Temperature', orientation='vertical', nb_labels=3)



#src = maya.pipeline.scalar_field(T)
#maya.pipeline.surface(src, colormap='jet')
#maya.colorbar(title='Temperature', orientation='vertical', nb_labels=3)


#maya.quiver3d(X, Y, Z, V1, V2, V3, colormap = 'Oranges', line_width=0.2, scale_factor=0.5, mask_points=600)
#maya.colorbar(title='Velocity', nb_labels=3)

maya.outline()
maya.show()




################################################################

###################### Ring Spectrum ################################
"""
t.ring_spectrum_showtime()
data = t.ring_spectrum_read(6,'Uek')
data = np.array(data)
print data
"""
######################################################################
