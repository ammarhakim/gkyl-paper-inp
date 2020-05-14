Input files used for PhD Thesis of James (Jimmy) Juno.
Python scripts also included to generate figures. Python scripts require at least Postgkyl 1.5.2.

Chapter 4 (Benchmarks)

- Vlasov--Fokker--Planck
  - Relaxation Tests
    - r1: Relaxation of step function to Maxwellian. p=1 case
    - r2: Relaxation of step function to Maxwellian. p=2 case
    - r3: Relaxation of two drifting Maxwellians to a single drifting Maxwellian
  - Kinetic Sod-Shock
    - s1: Sod-shock problem. mean-free-path = 0.1*Lx
    - s2: Sod-shock problem. mean-free-path = 0.01*Lx
    - s3: Sod-shock problem. mean-free-path = 0.002*Lx
    - s4: Exact solution to Euler Sod-Shock
    - n1: Sod-shock with sonic rarefaction, periodic domain, p=1, mean-free-path = 0.01*Lx
    - n2: Same as n1, p=2
- Vlasov--Maxwell
  - Conservation Tests (Asymmetric density gradient with background flow)
    - c1: p = 1, dx = 24 lambda_{D}, dv = 1 v_{th}, dt = 0.2 omega_{pe}^{-1}
    - c2: p = 1, dx = 24 lambda_{D}, dv = 1/2 v_{th}, dt = 0.2 omega_{pe}^{-1}
    - c3: p = 1, dx = 24 lambda_{D}, dv = 1/4 v_{th}, dt = 0.2 omega_{pe}^{-1}
    - c4: p = 2, dx = 24 lambda_{D}, dv = 1 v_{th}, dt = 0.4 omega_{pe}^{-1}
    - c5: p = 2, dx = 24 lambda_{D}, dv = 1 v_{th}, dt = 0.2 omega_{pe}^{-1}
    - c6: p = 2, dx = 24 lambda_{D}, dv = 1 v_{th}, dt = 0.1 omega_{pe}^{-1}
    - c7: p = 2, dx = 12 lambda_{D}, dv = 1 v_{th}, dt = 0.2 omega_{pe}^{-1}
    - c8: p = 2, dx = 6 lambda_{D}, dv = 1 v_{th}, dt = 0.1 omega_{pe}^{-1}
    - c9: p = 3, dx = 24 lambda_{D}, dv = 1 v_{th}, dt = 0.4 omega_{pe}^{-1}
    - c10: p = 3, dx = 24 lambda_{D}, dv = 1 v_{th}, dt = 0.2 omega_{pe}^{-1}
    - c11: p = 3, dx = 24 lambda_{D}, dv = 1 v_{th}, dt = 0.1 omega_{pe}^{-1}
    - c12: p = 3, dx = 12 lambda_{D}, dv = 1 v_{th}, dt = 0.2 omega_{pe}^{-1}
    - c13: p = 3, dx = 6 lambda_{D}, dv = 1 v_{th}, dt = 0.1 omega_{pe}^{-1}
    - c14: p = 2, dx = 12 lambda_{D}, dv = 1/2 v_{th}, dt = 0.2 omega_{pe}^{-1}
    - c15: p = 3, dx = 12 lambda_{D}, dv = 1/2 v_{th}, dt = 0.2 omega_{pe}^{-1}
    - c16: p = 2, dx = 6 lambda_{D}, dv = 1/4 v_{th}, dt = 0.1 omega_{pe}^{-1}
    - c17: p = 3, dx = 6 lambda_{D}, dv = 1/4 v_{th}, dt = 0.1 omega_{pe}^{-1}
    - c18: p = 2, dx = 3 lambda_{D}, dv = 1/4 v_{th}, dt = 0.05 omega_{pe}^{-1}
    - c19: p = 3, dx = 3 lambda_{D}, dv = 1/4 v_{th}, dt = 0.05 omega_{pe}^{-1}
    - c20: p = 2, dx = 3 lambda_{D}, dv = 1/8 v_{th}, dt = 0.05 omega_{pe}^{-1}
    - c21: p = 3, dx = 3 lambda_{D}, dv = 1/8 v_{th}, dt = 0.05 omega_{pe}^{-1}
    - c22: p = 2, dx = 1.5 lambda_{D}, dv = 1/4 v_{th}, dt = 0.025 omega_{pe}^{-1}
    - c23: p = 3, dx = 1.5 lambda_{D}, dv = 1/4 v_{th}, dt = 0.025 omega_{pe}^{-1}
    - c24: p = 2, dx = 1.5 lambda_{D}, dv = 1/8 v_{th}, dt = 0.025 omega_{pe}^{-1}
    - c25: p = 3, dx = 1.5 lambda_{D}, dv = 1/8 v_{th}, dt = 0.025 omega_{pe}^{-1}
    - c26: p = 2, dx = 0.75 lambda_{D}, dv = 1/4 v_{th}, dt = 0.0125 omega_{pe}^{-1}
    - c27: p = 3, dx = 0.75 lambda_{D}, dv = 1/4 v_{th}, dt = 0.0125 omega_{pe}^{-1}
    - c28: p = 2, dx = 0.75 lambda_{D}, dv = 1/8 v_{th}, dt = 0.0125 omega_{pe}^{-1}
    - c29: p = 3, dx = 0.75 lambda_{D}, dv = 1/8 v_{th}, dt = 0.0125 omega_{pe}^{-1}
    - c30: p = 2, dx = 0.375 lambda_{D}, dv = 1/4 v_{th}, dt = 0.00625 omega_{pe}^{-1}
    - c31: p = 3, dx = 0.375 lambda_{D}, dv = 1/4 v_{th}, dt = 0.00625 omega_{pe}^{-1}
    - c32: p = 2, dx = 0.375 lambda_{D}, dv = 1/8 v_{th}, dt = 0.00625 omega_{pe}^{-1}
    - c33: p = 3, dx = 0.375 lambda_{D}, dv = 1/8 v_{th}, dt = 0.00625 omega_{pe}^{-1}
  - Advection in Specified Electromagnetic Fields
    - a1: non-resonant solution, p = 2
    - a2: resonant solution, p = 2
    - a3: same as a1, but run for 1000 Omega_{c}^{-1}
    - a4: same as a3, but with p = 3
  - Landau Damping of Langmuir Waves
  - Three-Species Collisionless Shock
  - Lower Hybrid Drift Instability
  - Hybrid Two-Stream-Filamentation Mode
- Vlasov--Maxwell--Fokker--Planck
  - Collisional Landau Damping of Langmuir Waves
  - Magnetic Pumping

Chapter 5 (A Deep Dive into the Distribution Function)

- Perpendicular Collisionless Shock
- Phase Space Dynamics of Filamentation-Type Instabilities