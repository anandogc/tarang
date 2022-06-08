---
title: "Incompressible Fluid"
permalink: /fluid-incompress
---
{% include head.html %}

## Incompressible Fluid Solver

This solves the Navier-Stokes equation for hydrodynamic turbulence in cube geometry. 

\\[
\begin{aligned}
  \frac{\partial \bf u}{\partial t} + (\bf u \cdot \nabla)\bf u = &- \nabla p  + \bf F + \nu \nabla^2 \bf u
\end{aligned}
\\]

where,\
\\(\bf u\\) is the velocity field,\
\\(p\\) is the pressure,\
\\(\bf F\\) is the external Force,\
\\(\nu\\) is the kinamatic viscocity.\


program.kind: FLUID_INCOMPRESS

### Initial Conditions

#### 1. Complex field
Reads from a complex array of size \\(N_x \times N_y \times N_z/2+1\\) stored in an HDF5 format.\
\\(u_x\\) called as ```V1``` is stored in ```U.V1.hdf5/U.V1```\
\\(u_y\\) called as ```V2``` is stored in ```U.V1.hdf5/U.V2```\
\\(u_z\\) called as ```V3``` is stored in ```U.V1.hdf5/U.V3```

#### 2. Reduced complex field
Reads from a Complex array of lower grid resolution, file format is same as above.

#### 3. Real field
Reads from real arrays of size \\(N_x \times \y \times \N_z \\).\
```V1r``` is stored in ```U.V1r.h5/V1r```\
```V2r``` is stored in ```U.V2r.h5/V2r```\
```V3r``` is stored in ```U.V3r.h5/V3r```

#### 4. Modes
Activates the given modes only
<table>
	<thead>
		<col>
		<col pattern="\d+,\d+,\d+">
		<col pattern="(\d+,\d+)+">
	</thead>
	<tbody>
		<th>
			<td></td>
			<td>Mode</td>
			<td>Amplitude</td>
		</th>
		<tr>
			<td></td>
			<td></td>
			<td></td>
		</tr>
	</tbody>
</table>

#### 5. Energy helicity spectrum

#### 6. Taylor Green

#### 7. ABC
    
#### 8. Non helical to helical

#### 400. Vortex
    
#### 420. Channel flow
    
#### 501. User defined1
#### 502. User defined2

### Forcing
Following scenems are available:
1. decay
2. Carati_scheme
3. given_modes
4. Taylor Green
5. ABC
6. using_random_noise
7. Coriolis
8. Keplerian
9. Keplerian_SB
10. Liquid_metal
11. Kolmogorov_flow
12. Liquid_metal_const_energy_supply
13. Ekman_friction
14. Ekman_friction_const_energy_supply
15. pressure_grad
16. user_defined1
17. user_defined2

#### 1. Decay

No Forcing is applied

field_procedure: 0

#### 2. Carati Scheme

#### 3. Given modes
Only the given modes in the parameter file are exited.

field_procedure: 3

#### 4. Taylor Green

{% capture formula %}
{% raw %}
\begin{align*}
F_x &= amp*sin(k_0 \cdot x)cos(k_0 \cdot y) cos(k_0 \cdot z)\\
F_y &= -amp*cos(k_0 \cdot x)sin(k_0 \cdot y) cos(k_0 \cdot z)\\
F_z &= 0
\end{align*}
{% endraw %}
{% endcapture %}
<katex-block>{{ formula | escape }}</katex-block>

field_procedure: 4<br>
parameters: [\"\\(k_0\ amp\\)\"]


#### 5. ABC

{% capture formula %}
{% raw %}
\begin{align*}
F_x &= amp (B \cos(k_0 y) + C \sin(k_0 z))\\
F_y &= amp (A \sin(k_0 x) + C \cos(k_0 z))\\
F_z &= amp (A \cos(k_0 x) + C \cos(k_0 y))
\end{align*}
{% endraw %}
{% endcapture %}
<katex-block>{{ formula | escape }}</katex-block>

field_procedure: 5<br>
parameters: [\"\\(k_0\ A\ B\ C\ amp\\)\"]

#### 6. Random Noise

+ Initially, a 1D energy spectrum is created using Pre.py, typically using Pope/Pow model.
+ A 3D loop is made on \\(k_x\\), \\(k_y\\) and \\(k_z\\), \\(\|K\|\\) is found for each mode,
+ it's ceiling is taken to get the 1D energy spectrum index, such that \\(k_{i-1} < \|K\| \le k_i\\),
+ to find the modal energy, the shell energy id divided by number of moded in the shell, 
+ amplitude and phase are found in Craya-Herring basis then put to the Force array.

The random seed is updated every random_interval in eddy turn over time units.

{% capture formula %}
{% raw %}
\begin{align*}
Kmag &= \sqrt{kx^2 + ky^2 + kz^2}\\
index &= ceil(Kmag)\\
energy\_supply\_k &= U.energy\_supply\_spectrum(index)/ global.spectrum.shell.modes\_in\_shell(index);\\
amp\_plus &= sqrt((energy\_supply\_k+helicity\_supply\_k/Kmag)/random\_interval);\\
amp\_minus &= sqrt((energy\_supply\_k-helicity\_supply\_k/Kmag)/random\_interval);\\
phase\_plus &= 2*M\_PI * rand\_struct.random();\\
phase\_minus &= 2*M\_PI * rand\_struct.random();
\end{align*}
{% endraw %}
{% endcapture %}
<katex-block>{{ formula | escape }}</katex-block>

```c++
Put_or_add_force_vector(U, lx, ly, lz, amp_plus, amp_minus, phase_plus, phase_minus);
```

field_procedure: 6\
parameters: [inner_radius, outer_radius, random_interval]

#### 7. Coreolis Forcing

```c++
if (abs(two_omega1) > MYEPS){
	U.Force2 +=  two_omega1*(U.cvf.V3);
	U.Force3 += -two_omega1*(U.cvf.V2);
}

if(abs(two_omega2) > MYEPS){
	U.Force1 += -two_omega2*(U.cvf.V3);
	U.Force3 += two_omega2*(U.cvf.V1);

}
if(abs(two_omega3) > MYEPS){
	U.Force1 += two_omega3*(U.cvf.V2);
	U.Force2 += -two_omega3*(U.cvf.V1);
}
```

field_procedure: 21\
parameters: [\\(\omega_1\\), \\(\omega_2\\), \\(\omega_3\\)]\
Enable: <input type="checkbox" class="force coreolis">\
\\(\omega_1\\) = <input type="text" match="\d+\.\d+"/>\
\\(\omega_2\\) = <input type="text" match="\d+\.\d+"/>\
\\(\omega_3\\) = <input type="text" match="\d+\.\d+"/>