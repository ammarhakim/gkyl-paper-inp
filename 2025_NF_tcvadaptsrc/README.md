2025 paper on adaptive source and NT/PT simulations of TCV discharges.
The results of the publications were produced with the version [902b3e9](https://github.com/ammarhakim/gkeyll/tree/902b3e9d104d43841465e8341978afd5004d107c)

Input files used to generate results presented in section 4:
- PT_coarse.c : coarse grid PT simulation (24x16x12x12x8), green in figure 3.
- PT_baseline.c : baseline grid PT simulation (48x32x16x12x8), orange in figure 3, red in figure 7,9,10 and figure 8 left.
- PT_hd.c : high resolution grid PT simulation (96x64x16x12x8), blue in figure 3, figure 5, and figure 6.
- NT_baseline.c : baseline grid NT simulation (48x32x16x12x8), blue in figure 7,9,10 and figure 8 right.

The provided Makefile will compile an input file named `gkyl.c` into an executable named `gkyl`. One must make sure to update line 26 i.e.
```bash
-include ${HOME}/gkyl_main/gkylsoft/gkylzero/share/config.mak
```
with the correct path to the `gkylsoft` installation on your machine.