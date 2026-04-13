mkdir -p python-plots

pgkyl "gk_lorentzian_mirror-field_[0-9]*.gkyl" interp sel --z0 0.0 col pl --title 'Midplane electrostatic potential' --xlabel 'Time (s)' --ylabel 'Electric Potential φ (V)' --saveas python-plots/field_potential_z0_vs_time.png --no-show &

pgkyl gk_lorentzian_mirror-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl pl --title 'Upper Boundary Flux' --xlabel 'Time (s)' --logy --scatter --saveas python-plots/ion_bflux_xupper_M0_vs_time.png --no-show --no-legend --subplot-ylabels 'M0,M1,M2Par,M2Perp' &

pgkyl gk_lorentzian_mirror-ion_integrated_moms.gkyl pl --title 'Integrated Ion Moments' --xlabel 'Time (s)' --scatter --saveas python-plots/ion_integrated_moments_vs_time.png --no-show --no-legend --subplot-ylabels 'M0,M1,M2Par,M2Perp' &

pgkyl "gk_lorentzian_mirror-ion_BiMaxwellianMoments_[0-9]*.gkyl" interp col ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' pl --xlabel 'Time (s)' --ylabel "Z, m" --saveas python-plots/ion_BiMaxwellian_moments_vs_time.png --no-show --no-legend --subplot-titles 'Density $m^3$, $U_\parallel m/s$, $T_\parallel$ eV, $T_\perp$ eV' &

frame=$(ls gk_lorentzian_mirror-ion_cflrate_*.gkyl | sed -E 's/.*_([0-9]+)\.gkyl/\1/' | sort -n | tail -1)
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp sel --z0 0.0 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=0.0 (Frame ${frame})" --saveas python-plots/ion_distf_z0.0_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp sel --z0 0.5 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=0.5 (Frame ${frame})" --saveas python-plots/ion_distf_z0.5_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp sel --z0 1.0 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=1.0 (Frame ${frame})" --saveas python-plots/ion_distf_z1.0_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp sel --z0 1.5 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=1.5 (Frame ${frame})" --saveas python-plots/ion_distf_z1.5_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp sel --z0 1.9 ev 'f abs' pl --xlabel 'Parallel Velocity $v_\parallel$ (m/s)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-25 --title "Ion Distribution Function at z=1.9 (Frame ${frame})" --saveas python-plots/ion_distf_z1.9_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp integ 2 ev 'f abs' pl --xlabel 'Axial Position z (m)' --ylabel 'Parallel Velocity $v_\parallel$ (m/s)' --logz --zmin 1e-35 --title "Ion Distribution Integrated over μ (Frame ${frame})" --saveas python-plots/ion_distf_integrated_mu_frame${frame}.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl interp integ 1 ev 'f abs' pl --xlabel 'Axial Position z (m)' --ylabel 'Magnetic Moment $\mu$ (J/T)' --logz --zmin 1e-10 --title "Ion Distribution Integrated over $v_\parallel$ (Frame ${frame})" --saveas python-plots/ion_distf_integrated_vpar_frame${frame}.png --no-show &

pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_source_0.gkyl gk_lorentzian_mirror-ion_0.gkyl interp sel --z2 0.5 pl -c --xlabel 'Axial Position z (m)' --ylabel 'Parallel Velocity $v_\parallel$ (m/s)' --title 'Ion Source vs Distribution at $z_2=0.5$ (Frame 0)' --logz --zmin 1e-25 --saveas python-plots/ion_source_vs_distf_z0.5_frame0.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_source_0.gkyl interp sel --z2 0.5 pl -c --xlabel 'Axial Position z (m)' --ylabel 'Parallel Velocity $v_\parallel$ (m/s)' --title 'Ion Source at $z_2=0.5$ (Frame 0)' --logz --zmin 1e-25 --saveas python-plots/ion_source_z0.5_frame0.png --no-show &
pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_0.gkyl interp sel --z2 0.5 pl -c --xlabel 'Axial Position z (m)' --ylabel 'Parallel Velocity $v_\parallel$ (m/s)' --title 'Ion Distribution at $z_2=0.5$ (Frame 0)' --logz --zmin 1e-25 --saveas python-plots/ion_distf_z0.5_frame0.png --no-show &

pgkyl gk_lorentzian_mirror-ion_source_HamiltonianMoments_0.gkyl interp pl --xlabel 'Axial Position z (m)' --title 'Source Hamiltonian Moments' --logz --zmin 1e-25 --saveas python-plots/ion_source_HamiltonianMoments_frame0.png --no-show &

# Plot the cfl rate at the last frame at z=0.0
pgkyl gk_lorentzian_mirror-ion-cflrate_${frame}.gkyl sel --z0 0.0 pl --title 'CFL Rate at z=0.0' --xlabel '$v_\parallel$ computational' --ylabel '$\mu$ computational' --logz --saveas python-plots/cfl_rate_z0.0_frame${frame}.png --no-show & 
# Plot cfl rate 5 frames ago
pgkyl gk_lorentzian_mirror-ion-cflrate_$(($frame - 5)).gkyl sel --z0 0.0 pl --title 'CFL Rate at z=0.0 (5 frames ago)' --xlabel '$v_\parallel$ computational' --ylabel '$\mu$ computational' --logz --saveas python-plots/cfl_rate_z0.0_frame$(($frame - 5)).png --no-show &

# Save data for the time trace of the field at z= 0.0, 0.5, 1.0, 1.5, and 1.9
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 0.0 col write -s -f field_time_trace_z0_eq_0 &
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 0.5 col write -s -f field_time_trace_z0_eq_0,5 &
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 1.0 col write -s -f field_time_trace_z0_eq_1 &
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 1.5 col write -s -f field_time_trace_z0_eq_1,5 &
pgkyl gk_lorentzian_mirror-field_[0-9]*.gkyl interp sel --z0 1.9 col write -s -f field_time_trace_z0_eq_1,9 &

# Plot bimaxwellian moments at the final frame
pgkyl gk_lorentzian_mirror-ion_BiMaxwellianMoments_${frame}.gkyl interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' pl --title 'Final Bi-Maxwellian Moments' --saveas python-plots/ion_BiMaxwellianMoments_frame${frame}.png --no-show --no-legend --subplot-ylabels 'Density $m^3$, $U_\parallel m/s$, $T_\parallel$ eV, $T_\perp$ eV' &

pgkyl gk_lorentzian_mirror-ion_570.gkyl gk_lorentzian_mirror-ion_575.gkyl interp ev 'f[1] f[0] - f[0] / abs' sel --z0 0.0 pl --title 'Fractional error between begining and end of 570 and 575' --saveas python-plots/ion_distf_diff_570_575.png --no-show &

# pgkyl gk_lorentzian_mirror-ion_BiMaxwellianMoments_12.gkyl gk_lorentzian_mirror-ion_BiMaxwellianMoments_0.gkyl ../initial-conditions/boltz-elc-288z-nu2000/gk_lorentzian_mirror-ion_BiMaxwellianMoments_1500.gkyl interp pl -f0

# pgkyl good-run-enhanced-nu-IC/gk_lorentzian_mirror-ion-cflrate_275.gkyl integ 1 pl --title 'cfl OAP integ 1 old' --logz &
# pgkyl good-run-enhanced-nu-IC/gk_lorentzian_mirror-ion-cflrate_275.gkyl integ 2 pl --title 'cfl OAP integ 2 old' --logz &
# pgkyl good-run-enhanced-nu-IC/gk_lorentzian_mirror-ion-cflrate_285.gkyl integ 1 pl --title 'cfl RDP integ 1' --logz &
# pgkyl good-run-enhanced-nu-IC/gk_lorentzian_mirror-ion-cflrate_285.gkyl integ 2 pl --title 'cfl RDP integ 2' --logz &


# pgkyl gk_lorentzian_mirror-ion-cflrate_7.gkyl integ 1 pl --title 'cfl OAP integ 1 relaxed' --logz &
# pgkyl gk_lorentzian_mirror-ion-cflrate_7.gkyl integ 2 pl --title 'cfl OAP integ 2 relaxed' --logz &

# pgkyl good-run-enhanced-nu-IC/gk_lorentzian_mirror-ion-cflrate_285.gkyl sel --z0 0.95 pl --title 'cfl OAP integ 1 old' --logz &
# pgkyl gk_lorentzian_mirror-ion-cflrate_1.gkyl sel --z0 0.95 pl --title 'cfl OAP integ 2 relaxed' --logz &

# pgkyl good-run-enhanced-nu-IC/gk_lorentzian_mirror-ion-cflrate_275.gkyl gk_lorentzian_mirror-ion-cflrate_1.gkyl sel --z0 0.0 ev 'f[0] f[1] -' pl --title 'difference' --logz &
# pgkyl good-run-enhanced-nu-IC/gk_lorentzian_mirror-ion-cflrate_275.gkyl gk_lorentzian_mirror-ion-cflrate_1.gkyl sel --z0 0.0 pl &

# pgkyl "gk_lorentzian_mirror-field_[0-9]*.gkyl" interp anim &
# pgkyl good-long-RDP-few-cycles/gk_lorentzian_mirror-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl sel -c0 pl --title 'bflux M0 moment xupper true loss cone' --xlabel 'time, s' --ylabel 'bflux M0 moment xupper' --logy --scatter &
# pgkyl "../stellar-wham1x-288z-restart-true-collisions/gk_lorentzian_mirror-ion_bflux_xupper_integrated_M0M1M2parM2perp.gkyl" sel -c0 pl --title 'bflux M0 moment xupper true loss cone' --xlabel 'time, s' --ylabel 'bflux M0 moment xupper' --logy --scatter --xmin '463e-6'&

# pgkyl gk_lorentzian_mirror-ion_BiMaxwellianMoments_0.gkyl gk_lorentzian_mirror-ion_BiMaxwellianMoments_790.gkyl "../stellar-wham1x-288z-restart-true-collisions/gk_lorentzian_mirror-ion_BiMaxwellianMoments_80.gkyl" interp pl --title 'BiMaxwellianMoments' -f0 --logy &

# pgkyl gk_lorentzian_mirror-ion_9.gkyl interp sel --z0 0.0 pl --title 'f at t=0' --logz --zmin 1e-25 &
# pgkyl gk_lorentzian_mirror-ion_19.gkyl interp sel --z0 0.0 pl --title 'f at t=19' --logz --zmin 1e-25 &
# pgkyl gk_lorentzian_mirror-ion_29.gkyl interp sel --z0 0.0 pl --title 'f at t=29' --logz --zmin 1e-25 &


# pgkyl gk_lorentzian_mirror-ion_11.gkyl gk_lorentzian_mirror-ion_19.gkyl interp sel --z2 0.5 ev 'f[0] f[1] -' pl --title 'df after OAP' &
# pgkyl gk_lorentzian_mirror-ion_20.gkyl gk_lorentzian_mirror-ion_21.gkyl interp sel --z2 0.2 ev 'f[0] f[1] -' pl --title 'df after RDP' &


# pgkyl gk_lorentzian_mirror-ion_20.gkyl gk_lorentzian_mirror-ion_21.gkyl interp integ 2 ev 'f[0] f[1] - abs' pl --title 'df after RDP' --logz --zmin 1e-25&

# pgkyl "gk_lorentzian_mirror-ion_BiMaxwellianMoments_[0-9]*.gkyl" interp sel -c0 anim --logy &