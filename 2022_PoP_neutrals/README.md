2022 Kinetic modeling of neutral transport for a continuum GK code (PoP paper)

1. analytic-benchmarks:
   Input files used to generate the results presented in section III.A.
   The results for this section were produced with revision [c8b24bb98c25].
   (https://github.com/ammarhakim/gkyl/commit/c8b24bb98c2503bd83b9518f640a58be4e0f0a69) of the gkyl code.
   - 1x1v-cx-theory: test neutral CX interaction with static plasma and fixed recycling BCs.
   - 1x1v-iz-theory: test neutral ionization interaction with static plasma and fixed recycling BCs.

2. degas-benchmarks:
   Input files used to generate the results presented in section III.B.
   The results for this section were produced with revision [c8b24bb98c25].
   (https://github.com/ammarhakim/gkyl/commit/c8b24bb98c2503bd83b9518f640a58be4e0f0a69) of the gkyl code.
   - 1x3v-lowDens-izCX: test neutral CX and IZ interactions with static background plasma and n0 = 1e18, Lz=40 m
   - 1x3v-lowDens-izOnly: test neutral IZ interaction with static background plasma and n0 = 1e18, Lz=40 m
   - 1x3v-lowDens-CXonly: test neutral CX interaction  with static background plasma and n0 = 1e18, Lz=40 m
   - 1x3v-hiDens-izCX: test neutral CX and IZ interactions with static background plasma and n0 = 1e19, Lz=10 m
   - 1x3v-hiDens-izOnly: test neutral IZ interactions with static background plasma and n0 = 1e19, Lz=10 m
   - 1x3v-hiDens-CXonly: test neutral CX interactions with static background plasma and n0 = 1e19, Lz=10 m

3. nstx-SOL-tests:
   Input files used to generate the results presented in section IV.
   The results for this section were produced with revision [8a2b6a47522b].
   (https://github.com/ammarhakim/gkyl/commit/8a2b6a47522bdf4d215d174605d7dd4d976bbae7) of the gkyl code
   - nstx-noNeut.lua: 3x2v sim with NSTX SOL parameters, no neutrals
   - nstx-neut.lua: 3x2v sim with NSTX SOL parameters and 3x3v neutrals, with CX, IZ, and recycling BCs