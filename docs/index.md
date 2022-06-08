---
title: "Tarang"
---
{% include head.html %}

## Welcome to Tarang

Tarang is a pseudo-spectra fluid solver in a regular cube geometry. Basically is solves this Navier-Stokes equation.

\\[
\begin{aligned}
  \frac{\partial \bf u}{\partial t} + (\bf u \cdot \nabla)\bf u = &- \nabla p  + \bf F + \nu \nabla^2 \bf u
\end{aligned}
\\]

where,\
\\(\bf u\\) is the velocity field,\
\\(p\\) is the pressure,\
\\(\bf F\\) is the external Force,\
\\(\nu\\) is the kinamatic viscocity.

For timestepping, *Euler*, *RK2*, and *RK4* integration schemes are available. Dealiasing is done using \\(2/3\\) dealiaing scheme. The viscous term is solved using exponential method.

The following solvers are available

+ [Incompressible Fluid](/fluid-incompress)
+ Incompressible Scalar
+ RBC


### Boundary Conditions

#### FFF
#### SFF
#### SSF
#### SSS

### Diagnostics

### Aliasing

### Integration Scheme

### Large Eddy Simulation

### FFT

### Field

### time
