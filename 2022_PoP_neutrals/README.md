2022 Kinetic modeling of neutral transport for a continuum GK code PoP paper

1. analytic-benchmarks:
   Input files used to generate the results presented in section III.A.
   The results for this section were produced with revision [d9eafc9f1234].
   (https://github.com/ammarhakim/gkyl/commit/d9eafc9f12343745e18514be006afac340e879bc) of the gkyl code.
   - 1x1v-cx-theory: test neutral CX interaction with static plasma and fixed recycling BCs.
   - 1x1v-iz-theory: test neutral ionization interaction with static plasma and fixed recycling BCs.

2. degas-benchmarks:
   Input files used to generate the results presented in section III.B.
   The results for this section were produced with revision [d2c5faba64ab].
   (https://github.com/ammarhakim/gkyl/commit/d2c5faba64abc2eccb6a381aad026c8eaa1dd0232) of the gkyl code.
   - 1x3v-lowDens-izCX: test neutral CX and IZ interactions with static background plasma and n0 = 1e18, Lz=40 m
   - 1x3v-lowDens-izOnly: test neutral IZ interaction with static background plasma and n0 = 1e18, Lz=40 m
   - 1x3v-lowDens-CXonly: test neutral CX interaction  with static background plasma and n0 = 1e18, Lz=40 m
   - 1x3v-hiDens-izCX: test neutral CX and IZ interactions with static background plasma and n0 = 1e19, Lz=10 m
   - 1x3v-hiDens-izOnly: test neutral IZ interactions with static background plasma and n0 = 1e19, Lz=10 m
   - 1x3v-hiDens-CXonly: test neutral CX interactions with static background plasma and n0 = 1e19, Lz=10 m

3. nstx-SOL-tests:
   Input files used to generate the results presented in section IV.
   The results for this section were produced with revision [416865089be6].
   (https://github.com/ammarhakim/gkyl/commit/223c7856f3fbf97787c768a62c115363362785f2) of the gkyl code
   - nstx-noNeut.lua: 3x2v sim with NSTX SOL parameters, no neutrals
   - nstx-neut.lua: 3x2v sim with NSTX SOL parameters and 3x3v neutrals, with CX, IZ, and recycling BCs
