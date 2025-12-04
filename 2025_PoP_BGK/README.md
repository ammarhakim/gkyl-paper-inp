2025 axisymmetric ASDEX SOL with implicit BGK paper

Input files used to generate results presented in section III.
In this section, the results for the benchmark tests 1&2 were produced with the version [31495ab](https://github.com/ammarhakim/gkylzero/commit/31495ab27330f4f988af5cce012bba0cf050a592), the results for the benchmark test 3 were produced with the version [3c8be52](https://github.com/ammarhakim/gkeyll/commit/3c8be5239cb82d520a8336f2779418c575a13903)
1. Conservation:
  - rt_gk_bgk_relax_1x2v_p1
2. Sod Shock:
  - nu1/rt_euler_periodic_sodshock:
  - nu1/rt_gk_bgk_im_periodic_sodshock_1x2v_p1
  - nu2/rt_gk_bgk_im_periodic_sodshock_1x2v_p1
  - nu2/rt_gk_lbo_periodic_sodshock_1x2v_p1
3. Cross Validation:
  - rt_gk_bgk_im_cross_relax_1x2v_p1
  - rt_gk_lbo_cross_relax_1x2v_p1

Input files used to generate results presented in section IV.
The results for the production runs in this section were produced with the version [9363a4b](https://github.com/ammarhakim/gkylzero/tree/9363a4b92191811ae54f4daa1e24d0d77f40ac34)
- rt_gk_bgk_im_asdex_IC2_2x2v_p1: low resolution BGK 
- rt_gk_bgk_im_asdex_IC2_Nz96_2x2v_p1: high resolution BGK
- rt_gk_lbo_asdex_IC2_Nz96_2x2v_p1: high resolution LBD
