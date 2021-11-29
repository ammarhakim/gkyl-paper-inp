2021 cross-species LBO collisions paper

Input files used to generate the results presented in section 4.
The results for this publication were produced with revision [d884962](https://github.com/ammarhakim/gkyl/commit/d8849629fc5d3519ccda55910c47bf1fb7674779) of the gkyl code, although some were also produced with [ab9248845d02](https://github.com/ammarhakim/gkyl/commit/ab9248845d02a1005616938daf460f153b2c9e9d) since shared-memory worked then and this
allowed to produced some results more quickly.

1. vmConservation:
  - c0-vmConservation-1x1v-p1: LBO-EM in 1v.
  - c10-vmConservation-1x2v-p1: LBO-EM in 2v.
  - c20-vmConservation-1x3v-p1: LBO-EM in 3v.
  - c30-vmConservation-1x1v-p1: LBO-G in 1v.
  - c40-vmConservation-1x2v-p1: LBO-G in 2v.
  - c50-vmConservation-1x3v-p1: LBO-G in 3v.
  - c60-vmConservation-1x2v-p1: LBO-EM in 2v used for scanning mass ratio.
2. gkConservation:
  - g0-gkConservation-1x1v-p1: LBO-ET in 1v.
  - g10-gkConservation-1x2v-p1: LBO-ET in 2v.
  - g20-gkConservation-1x1v-p1: LBO-G in 1v.
  - g30-gkConservation-1x2v-p1: LBO-G in 2v.
3. langmuirLandau:
  - selfOnly: electron-electron collisions only
  - lboG: self-collisions and LBO-G cross-collisions.
  - lboEM: self-collisions and LBO-EM cross-collisions.
  - lboET: self-collisions and LBO-ET cross-collisions.
4. velocityTempRelax:
  - b17-lboCrossHager-1x2v-p1: LBO-G beta=0.
  - b18-lboCrossHager-1x2v-p1: LBO-EM.
  - b19-lboCrossHager-1x2v-p1: LBO-ET.
