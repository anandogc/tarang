from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import tarang as t

pi = np.pi
#t.init("/home/abhishek/Dropbox/shared/vijay_abhishek/sample/stratified")
t.init("/home/abhi/Dropbox/shared/vijay_abhishek/sample/RBC/2d")
font = {'family' : 'serif', 'weight' : 'normal', 'size' : 28}
plt.rc('font', **font)

#############  To be integrated in code ################################
N = t.para.field.N
#K = [pi,2*pi,2*pi]
K = [3.14159, 2.22144, 2.22144]

no_spheres = t.para.energy_transfer.flux.no_spheres

def Mag(real,imag):
	return ((real**2 + imag**2)**0.5)

def Energy_vector(A,B,C):
	return (A**2 + B**2 +C**2)/2.0
	
def Energy_scalar(A):
	return (A**2)/2.0
	
	
class FFF_SLAB:
	def Max_radius_inside(self):
		if N[1] > 1:
			Kmag = min( (N[0]/2) * K[0] , (N[1]/2) * K[1])
			Kmag = min(Kmag, (N[2]/2) * K[2])
		elif N[1] == 1:
			Kmag = min( (N[0]/2) * K[0] , (N[2]/2) * K[2])
		return int(Kmag)
		

class SFF_SLAB:
	def Max_radius_inside(self):
		if N[1] > 1:
			Kmag = min( (N[0]-1) * K[0] , (N[1]/2) * K[1])
			Kmag = min(Kmag, (N[2]/2) * K[2])
		elif N[1] == 1:
			Kmag = min( (N[0]-1) * K[0] , (N[2]/2) * K[2])
		return int(Kmag)


sff_slab = SFF_SLAB()
Max_radius_inside = sff_slab.Max_radius_inside()

def flux_radii():
	radii = np.zeros(no_spheres)
	radii[0] = 0.0
	radii[1] = 2.0
	radii[2] = 4.0
	radii[3] = 8.0
	radii[no_spheres - 2] = Max_radius_inside/2.0
	radii[no_spheres - 1] = Max_radius_inside
	
	if no_spheres > 6:
		s = np.log2(Max_radius_inside/16.0) / (no_spheres - 5)
		for i in np.arange(4,(no_spheres - 3) + 1, 1):
			radii[i] = 8 * 2**(s*(i-3))
	
	return radii
	
	
		


########################################################################


################# glob #########################
"""
#A = 60
#B = 100
fig, axes = plt.subplots(figsize=(15,10))
#axes.plot(t.time,t.U.cvf.total_energy,'r',label=r'$E_u$')
#axes.plot(t.time,t.T.csf.total_energy,'g',label=r'$E_{\theta}$')
#axes.plot(t.time,t.U.cvf.total_dissipation,'r',label=r'$\epsilon_u$')
#axes.plot(t.time,t.T.csf.total_dissipation,'r',label=r'$\epsilon_{\theta}$')
#axes.plot(t.time,t.correlation.nusselt_no,'b',label=r'$Nu$')
#axes.plot(t.time,t.U.cvf.dissipation_coefficient,'r',label=r'$\nu$')
#axes.plot(t.time,t.T.csf.diffusion_coefficient,'b',label=r'$k$')
#axes.plot(t.time,t.Rlambda,'b',label=r'$Re_{\lambda}$')
#axes.plot(t.time,t.time_dt,'b',label=r'$dt$')
axes.plot(t.time,t.U.cvf.total_E1,'r',label=r'$E_x$')
axes.plot(t.time,t.U.cvf.total_E2,'g',label=r'$E_y$')
axes.plot(t.time,t.U.cvf.total_E3,'b',label=r'$E_z$')
axes.legend(loc=0,prop={'size':35})
#print "kmax_eta_u = ",np.mean(t.U.cvf.kmax_eta[A:B])
#print "kmax_eta_theta = ",np.mean(t.T.csf.kmax_eta[A:B])
plt.show()
"""
#################### Spectrum #################################
"""
t.spectrum_showtime()
t.spectrum_read(11)

fig, axes = plt.subplots(figsize=(15,10))
axes.loglog(t.k, (t.U.cvf.Uek_shell_ek1[0]+t.U.cvf.Uek_shell_ek3[0]) ,'r',label=r'$E_u(k)$')
axes.loglog(t.k, t.T.csf.Tek_shell_ek[0] ,'g',label=r'$E_{\theta}(k)$')
#axes.loglog(t.k, (t.U.cvf.UDk_shell_dissk1[0]+t.U.cvf.UDk_shell_dissk2[0]+t.U.cvf.UDk_shell_dissk3[0]) ,'b',label=r'$D_u(k)$')
#axes.loglog(t.k, (t.correlation.U_T_shell_ek1[0]+t.correlation.U_T_shell_ek2[0]+t.correlation.U_T_shell_ek3[0]) ,'c',label=r'$U(k)\theta(k)$')
#axes.loglog(t.k, (t.U.cvf.Fv_v_shell_ek1[0]+t.U.cvf.Fv_v_shell_ek2[0]+t.U.cvf.Fv_v_shell_ek3[0]) ,'m',label=r'$F_u(k)$')
#axes.loglog(t.k, t.T.csf.FT_T_shell_ek[0] ,'r--',label=r'$F_{\theta}(k)$')
axes.set_xlabel(r'$k$')
axes.legend(loc=0,prop={'size':35})
plt.show()
"""

#################### Flux ###############################
"""
kf = flux_radii()
t.flux_showtime()
t.flux_read(11)

fig, axes = plt.subplots(figsize=(5,5))
axes.loglog(kf, t.energytr.U2U[0][0:len(kf)] ,'r',label=r'$\Pi_u(k)$')
axes.loglog(kf, t.energytr.T2T[0][0:len(kf)] ,'g',label=r'$\Pi_{\theta}(k)$')
axes.set_xlabel(r'$k$')
axes.legend(loc=0,prop={'size':35})
plt.show()
"""
################## Fourier Probes ######################
"""
def Mag(real,imag):
	return ((real**2 + imag**2)**0.5)

t.field_k_read()
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
E_u1 = Energy_vector(ux1,uy1,uz1)
E_theta1 = Energy_scalar(theta1)

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
E_u2 = Energy_vector(ux2,uy2,uz2)
E_theta2 = Energy_scalar(theta2)


fig, axes = plt.subplots(figsize=(15,10))
axes.plot(time1, E_u1 ,'r',label=r'$1 \, 0 \, 1$')
axes.plot(time2, E_u2 ,'g',label=r'$0 \, 1 \, 1$')
axes.set_xlabel(r'$t$')
axes.set_ylabel(r'$E_u(t)$')
axes.legend(loc=0,prop={'size':35})
plt.show()
"""
################### Real Probes ###########################
"""
t.field_r_read()
t.field_r_showmode()

mode1 = t.field_r_readmode(1,1,1)
time1 = np.array(mode1.Time)
ux1 = np.array(mode1.Ux)
uy1 = np.array(mode1.Uy)
uz1 = np.array(mode1.Uz)
theta1 = np.array(mode1.T)
E_u1 = Energy_vector(ux1,uy1,uz1)
E_theta1 = Energy_scalar(theta1)

mode2 = t.field_r_readmode(1,1,32)
time2 = np.array(mode2.Time)
ux2 = np.array(mode2.Ux)
uy2 = np.array(mode2.Uy)
uz2 = np.array(mode2.Uz)
theta2 = np.array(mode2.T)
E_u2 = Energy_vector(ux1,uy1,uz1)
E_theta2 = Energy_scalar(theta1)



fig, axes = plt.subplots(figsize=(15,10))
axes.plot(time1, E_u1 ,'r',label=r'$1 \, 1 \, 1$')
axes.plot(time2, E_u2 ,'g',label=r'$1 \, 1 \, 32$')
axes.set_xlabel(r'$t$')
axes.set_ylabel(r'$E_u(t)$')
axes.legend(loc=0,prop={'size':35})
plt.show()


"""
