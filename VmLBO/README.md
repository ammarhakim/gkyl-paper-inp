Input files used for the paper describing the implementation and benchmarking of the LBO collision operator in the Vlasov-Maxwell solver (VmLBO).

Relaxation tests (Problem 1)

- r1: Relaxation of step function to Maxwellian. p=1 case
- r2: Relaxation of step function to Maxwellian. p=2 case
- r3: Same as r1, except CFL factor of 2x smaller
- r4: Same as r1, except 2x more cells in V (entropy convergence comparison)
- r5: Relaxation of bi-Maxwellian to a Maxwellian

Sod-shock tests (Problem 2)

- s1: Sod-shock problem. MFP = 0.1*Lx
- s2: Sod-shock problem. MFP = 0.01*Lx
- s3: Sod-shock problem. MFP = 0.002*Lx
- s4: Exact solution to Euler Sod-shock

Sod-shock with sonic rarefaction (Problem 2)

- n1: Sod-shock problem. MFP = 0.01*Lx
- n2: Same as n2 except on a periodic domain to test momentum conservation
- n3: Same as n2, except p=1
- ne: Exact Euler solution to n1

Collisional Landau damping (Problem 3)

- d1: Collisional Landau damping, nu=0.0 (no collisions)
- d2: Collisional Landau damping, nu=0.05
- d3: Collisional Landau damping, nu=0.1

Magnetic pumping (Problem 4)

- mp9:  nu=0.0 (no collisions).
- mp10: nu=0.1*pi
- mp20: nu=1.5
