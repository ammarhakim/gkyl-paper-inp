#[ ........................................................... ]#
#[
#[ Check energy balance in a Gkeyll gyrokinetic simulation.
#[
#[ This script assumes the existence of the following files:
#[   - <sim_name>-<species_name>_fdot_integrated_moms.gkyl
#[   - <sim_name>-<species_name>_integrated_moms.gkyl
#[   - <sim_name>-<species_name>_source_integrated_moms.gkyl if using a source.
#[   - <sim_name>-<species_name>_bflux_<direction><side>_integrated_HamiltonianMoments.gkyl
#[       if non-periodic, non-zero-flux boundaries are used.
#[   - <sim_name>-field_energy.gkyl
#[   - <sim_name>-field_energy_dot.gkyl
#[   - <sim_name>-dt.gkyl
#[
#[ Manaure Francisquez.
#[
#[ ........................................................... ]#

import numpy as np
import postgkyl as pg
import matplotlib.pyplot as plt
import os

data_dir        = './gk_sheath_3x2v_p1/cfl0p5/' #[ Where Gkeyll data is located.
simulation_name = 'gk_sheath_3x2v_p1' #[ Name of Gkeyll simulation.
species_names   = ['elc','ion'] #[ Name of the particles species simulated.

plot_balance        = True #[ Balance of between various terms.
plot_relative_error = True #[ Conservation relative error.

save_fig_to_file = False #[ Output a figure file?.
fig_file_dir     = './' #[ Where to place figure written out.
fig_file_format  = '.png' #[ Can be .png, .pdf, .ps, .eps, .svg.


#[ ............... End of user inputs (MAYBE) ..................... ]#

#[ Some RGB colors. These are MATLAB-like.
defaultBlue    = [0, 0.4470, 0.7410]
defaultOrange  = [0.8500, 0.3250, 0.0980]
defaultGreen   = [0.4660, 0.6740, 0.1880]
defaultPurple  = [0.4940, 0.1840, 0.5560]
defaultRed     = [0.6350, 0.0780, 0.1840]
defaultSkyBlue = [0.3010, 0.7450, 0.9330]
grey           = [0.5, 0.5, 0.5]
#[ Colors in a single array.
defaultColors = [defaultBlue,defaultOrange,defaultGreen,defaultPurple,defaultRed,defaultSkyBlue,grey,'black']

#[ LineStyles in a single array.
lineStyles = ['-','--',':','-.','None','None','None','None']
markers    = ['None','None','None','None','o','d','s','+']

#[ Some fontsizes used in plots.
xyLabelFontSize       = 17
titleFontSize         = 17
colorBarLabelFontSize = 17
tickFontSize          = 14
legendFontSize        = 14

#.Set the font size of the ticks to a given size.
def setTickFontSize(axIn,fontSizeIn):
  axIn.tick_params(axis='both',labelsize=fontSizeIn)
  offset_txt = axIn.yaxis.get_offset_text() # Get the text object
  offset_txt.set_size(fontSizeIn) # # Set the size.
  offset_txt = axIn.xaxis.get_offset_text() # Get the text object
  offset_txt.set_size(fontSizeIn) # # Set the size.

#[ Check if a file exists............. ]#
def does_file_exist(fileIn):
  if os.path.exists(fileIn):
     return True
  else:
     return False
#[ .................................... ]#

#[ Read data and time stamps from a DynVector.............. ]#
def read_dyn_vector(dataFile):
  pgData = pg.GData(dataFile)  #[ Read data with pgkyl.
  time   = pgData.get_grid()  #[ Time stamps of the simulation.
  val    = pgData.get_values()  #[ Data values.
  return np.squeeze(time), np.squeeze(val)
#[ ......................................................... ]#

#[ Labels used to identify boundary flux files.
edges = ["lower","upper"]
dirs = ["x","y","z"]

#[ ............... End common utilities ..................... ]#

if plot_balance:
  #[ Read the Hamiltonian moment of df/dt, the source and the particle fluxes.
  
  data_path = data_dir + '/' + simulation_name + '-'
  
  #[ Plot each contribution.
  figProp1a = (7.5, 4.5)
  ax1aPos   = [0.09, 0.15, 0.87, 0.78]
  fig1a     = plt.figure(figsize=figProp1a)
  ax1a      = fig1a.add_axes(ax1aPos)
  
  hpl1a = list()
  hpl1a.append(ax1a.plot([-1.0,1.0], [0.0,0.0], color='grey', linestyle=':', linewidth=1))
  
  #[ Read particle energy data.
  for sI in range(len(species_names)):
    species = species_names[sI]
  
    time_fdot, fdot_s = read_dyn_vector(data_path + species + '_fdot_integrated_moms.gkyl')

    source_file = data_path + species + '_source_integrated_moms.gkyl'
    has_source = does_file_exist(source_file)
    if has_source:
      time_src, src_s = read_dyn_vector(source_file)

    nbflux = 0
    time_bflux, bflux_s = list(), list()
    has_bflux = False
    for d in dirs:
      for e in edges:
        bflux_file = data_path + species + '_bflux_' + d + e + '_integrated_HamiltonianMoments.gkyl'
        has_bflux_at_boundary = does_file_exist(bflux_file)
        if has_bflux_at_boundary:
          time_bflux_tmp, bflux_tmp = read_dyn_vector(bflux_file)
          time_bflux.append(time_bflux_tmp)
          bflux_s.append(bflux_tmp)
          has_bflux = has_bflux or has_bflux_at_boundary
          nbflux += 1

    #[ Select the Hamiltonian moment.
    fdot_s = fdot_s[:,2]
    if has_source:
      src_s = src_s[:,2]
      src_s[0] = 0.0 #[ Set source=0 at t=0 since we don't have fdot and bflux then.
    else:
      src_s = 0.0*fdot_s

    if has_bflux:
      for i in range(nbflux):
        bflux_s[i] = bflux_s[i][:,2]
    
      bflux_tot_s = bflux_s[0] #[ Total boundary flux loss.
      for i in range(1,nbflux):
        bflux_tot_s += bflux_s[i]
    else:
      bflux_tot_s = 0.0*fdot_s
  
    if sI == 0:
      fdot      = fdot_s     
      src       = src_s      
      bflux_tot = bflux_tot_s
    else:
      fdot      += fdot_s     
      src       += src_s      
      bflux_tot += bflux_tot_s
  
  #[ Read field energy data.
  field_file = data_path + 'field_energy_dot.gkyl'
  has_field = does_file_exist(field_file)
  if has_field:
    time_field_dot, field_dot = read_dyn_vector(field_file)
  else:
    field_dot = 0.0*fdot
  
  #[ Compute the error.
  mom_err = src - bflux_tot - (fdot - field_dot)
    
  #[ Plot each contribution.
  hpl1a.append(ax1a.plot(time_fdot, fdot, color=defaultColors[0], linestyle=lineStyles[0], linewidth=2, marker=markers[0]))
  legendStrings = [r'$\dot{f}$']

  if has_source:
    hpl1a.append(ax1a.plot(time_src, src, color=defaultColors[2], linestyle=lineStyles[2], linewidth=2, marker=markers[2]))
    legendStrings.append(r'$\mathcal{S}$')

  if has_bflux:
    hpl1a.append(ax1a.plot(time_bflux[0], -bflux_tot, color=defaultColors[1], linestyle=lineStyles[1], linewidth=2, marker=markers[1]))
    legendStrings.append(r'$-\int_{\partial \Omega}\mathrm{d}\mathbf{S}\cdot\mathbf{\dot{R}}f$')

  if has_field:
    hpl1a.append(ax1a.plot(time_field_dot, field_dot, color=defaultColors[4], linestyle=':', linewidth=2, marker='+',markevery=8))
    legendStrings.append(r'$\dot{\phi}$')

  hpl1a.append(ax1a.plot(time_fdot, mom_err, color=defaultColors[3], linestyle=lineStyles[3], linewidth=2, marker=markers[3]))
  legendStrings.append(r'$E_{\dot{\mathcal{E}}}=\mathcal{S}-\int_{\partial \Omega}\mathrm{d}\mathbf{S}\cdot\mathbf{\dot{R}}f-(\dot{f}-\dot{\phi})$')
  
  ax1a.set_xlabel(r'Time ($s$)',fontsize=xyLabelFontSize, labelpad=+4)
  ax1a.set_title(r'Energy balance',fontsize=titleFontSize)
  ax1a.set_xlim( time_fdot[0], time_fdot[-1] )
  ax1a.legend([hpl1a[i][0] for i in range(1,len(hpl1a))], legendStrings, fontsize=legendFontSize, frameon=False, loc='lower right')
  setTickFontSize(ax1a,tickFontSize)
  
  if save_fig_to_file:
    plt.savefig(fig_file_dir+simulation_name+'_energy_balance'+fig_file_format)
  else:
    plt.show()

#[ .......................................................... ]#

if plot_relative_error:
  #[ Plot the error normalized for different time steps.
  
  figProp2a = (7.5, 4.5)
  ax2aPos   = [0.12, 0.15, 0.87, 0.77]
  fig2a     = plt.figure(figsize=figProp2a)
  ax2a      = fig2a.add_axes(ax2aPos)
  
  hpl2a = list()
  hpl2a.append(ax2a.plot([-1.0,1.0], [0.0,0.0], color='grey', linestyle=':', linewidth=1))
  ylabelString = ""
  
  data_path = data_dir + simulation_name + '-'
  
  #[ Read the species energy data.
  for sI in range(len(species_names)):
    species = species_names[sI]

    time_fdot, fdot_s = read_dyn_vector(data_path + species + '_fdot_integrated_moms.gkyl')
    time_distf, distf_s = read_dyn_vector(data_path + species + '_integrated_moms.gkyl')

    source_file = data_path + species + '_source_integrated_moms.gkyl'
    has_source = does_file_exist(source_file)
    if has_source:
      time_src, src_s = read_dyn_vector(source_file)

    nbflux = 0
    time_bflux, bflux_s = list(), list()
    has_bflux = False
    for d in dirs:
      for e in edges:
        bflux_file = data_path + species + '_bflux_' + d + e + '_integrated_HamiltonianMoments.gkyl'
        has_bflux_at_boundary = does_file_exist(bflux_file)
        if has_bflux_at_boundary:
          time_bflux_tmp, bflux_tmp = read_dyn_vector(bflux_file)
          time_bflux.append(time_bflux_tmp)
          bflux_s.append(bflux_tmp)
          has_bflux = has_bflux or has_bflux_at_boundary
          nbflux += 1

    #[ Select the Hamiltonian moment and remove the t=0 data point.
    fdot_s = fdot_s[1:,2]
    distf_s = distf_s[1:,2]

    if has_source:
      src_s = src_s[1:,2]
    else:
      src_s = 0.0*fdot_s

    if has_bflux:
      for i in range(nbflux):
        bflux_s[i] = bflux_s[i][1:,2]
    
      bflux_tot_s = bflux_s[0] #[ Total boundary flux loss.
      for i in range(1,nbflux):
        bflux_tot_s += bflux_s[i]
    else:
      bflux_tot_s = 0.0*fdot_s

    if sI == 0:
      fdot      = fdot_s     
      distf     = distf_s      
      src       = src_s      
      bflux_tot = bflux_tot_s
    else:
      fdot      += fdot_s     
      distf     += distf_s      
      src       += src_s      
      bflux_tot += bflux_tot_s

  #[ Read time step.
  time_dt, dt = read_dyn_vector(data_path + 'dt.gkyl')
    
  #[ Read field energy data.
  field_file = data_path + 'field_energy_dot.gkyl'
  has_field = does_file_exist(field_file)
  if has_field:
    time_field_dot, field_dot = read_dyn_vector(field_file)
    time_field, field = read_dyn_vector(data_path + 'field_energy.gkyl')
    field_dot = field_dot[1:]
    field = field[1:]
  else:
    field_dot = 0.0*fdot
    field = 0.0*fdot

  #[ Compute the relative error.
  mom_err = src - bflux_tot - (fdot - field_dot)
  mom_err_norm = mom_err*dt/(distf-field)

  #[ Plot the relative error.
  hpl2a.append(ax2a.semilogy(time_dt, np.abs(mom_err_norm), color=defaultColors[0], linestyle=lineStyles[0], linewidth=2))
  
  ax2a.set_xlabel(r'Time ($s$)',fontsize=xyLabelFontSize, labelpad=+4)
  ax2a.set_ylabel(r'$|E_{\dot{\mathcal{E}}}~\Delta t/\mathcal{E}|$',fontsize=xyLabelFontSize, labelpad=0)
  ax2a.set_xlim( time_fdot[0], time_fdot[-1] )
  setTickFontSize(ax2a,tickFontSize)
  
  if save_fig_to_file:
    plt.savefig(fig_file_dir+simulation_name+'_energy_conservation_rel_error'+fig_file_format)
  else:
    plt.show()
  
