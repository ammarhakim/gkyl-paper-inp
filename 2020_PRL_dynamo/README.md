This directory contains a collection of simulation inputs used in the paper [I. Pusztai et al. (2020) Phys. Rev. Lett., "Dynamo in weakly collisional non-magnetized plasmas impeded by Landau damping of magnetic fields", https://arxiv.org/abs/2001.11929]. References to figures below refer to this publication. 

These simulations are performed using the kinetic-Vlasov solver Gkeyll [version: cd65328c077f+ 2228+ default], for more information on the code visit https://gkyl.readthedocs.io/en/latest/index.html, or consult [J. Juno et al (2018) J. Comp. Phys 353, 110].

The input files are those with a `.lua` extension in each simulation directory below.

Contents:

* **Galloway-Proctor-flow_Fig1-kinetic-and-Fig2:** 
  Kinetic simulation of the Galloway-Proctor flow, corresponding to the solid lines in Fig. 1 and Fig. 2. 

* **Cnu-and-k-scan_Fig3-and-Fig4a:**
  This is a parameter scan in wavelength of the magnetic perturbations [ranging from L0 ("L0") to L0/8 ("L0per8"), with baseline domain size L0] and collision frequencies [ranging from 0.05 ("Cnu005") to 1 ("Cnu1") times the baseline values]. These results are presented in Fig. 3 and Fig. 4a.

* **Magnetization-scan_Fig4b:**
  Scan in magnetization shown in Fig. 4 b. The magnetic field varies between 1 and 100 T ["B1" and "B100", respectively].

* **Roberts-flow_Fig5:**
  Kinetic simulations of the Roberts flow, shown in Fig. 5. The collision frequency is scaled to 0.3 the physical value (dashed lines, "Roberts_Cnu03_Fig5"), and zero (solid lines, "Roberts_Cnu00_Fig5").   
