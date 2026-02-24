import numpy as np
import matplotlib.pyplot as plt
import tristanVis.isolde as isolde
from scipy.ndimage import gaussian_filter
import postgkyl as pg 
plt.style.use("postgkyl.mplstyle")

#import matplotlib as mpl
#mpl.rcParams["text.usetex"] = True

font_size = None
label_size = None
if font_size is None:
    font_size = plt.rcParams.get("font.size", 10.0)
if label_size is None:
    label_size = plt.rcParams.get("axes.labelsize", font_size)

filename='PPC10/usroutput'
f=open(filename,"r")
Lines = f.readlines()
f.close()

arr=Lines[-2].split()
e2_ppc10=[]
for i in range (len(arr)):
    e2_ppc10.append(float(arr[i][:-1]))

time10 = np.linspace(0,90,len(e2_ppc10))
    
filename='PPC100/usroutput'
f=open(filename,"r")
Lines = f.readlines()
f.close()

arr=Lines[-2].split()
e2_ppc100=[]
for i in range (len(arr)):
    e2_ppc100.append(float(arr[i][:-1]))

time100 = np.linspace(0,90,len(e2_ppc100))

filename='PPC10_cfl01/usroutput'
f=open(filename,"r")
Lines = f.readlines()
f.close()

arr=Lines[-2].split()
e2_ppc10_cfl01=[]
for i in range (len(arr)):
    e2_ppc10_cfl01.append(float(arr[i][:-1]))

time10cfl01 = np.linspace(0,90,len(e2_ppc10_cfl01))

step = 850
filename = 'PPC10/flds.tot.%05d'%step
fields = isolde.getFields(filename)

Ex = fields['ex'][0,0,:]
efft = np.fft.fft(Ex)
efft = np.fft.fftshift(efft)
e2fft = efft*np.conj(efft)
n = Ex.size
freq = np.fft.fftfreq(n, d=0.1)
freq = np.fft.fftshift(freq)

for step in range (851,901):
    filename = 'PPC10/flds.tot.%05d'%step
    fields = isolde.getFields(filename)

    Ex = fields['ex'][0,0,:]
    efft = np.fft.fft(Ex)
    efft = np.fft.fftshift(efft)
    e2fft = e2fft+efft*np.conj(efft)

e2fft_PPC10 = e2fft/50.0

step = 850
filename = 'PPC100/flds.tot.%05d'%step
fields = isolde.getFields(filename)

Ex = fields['ex'][0,0,:]
efft = np.fft.fft(Ex)
efft = np.fft.fftshift(efft)
e2fft = efft*np.conj(efft)
n = Ex.size
freq = np.fft.fftfreq(n, d=0.1)
freq = np.fft.fftshift(freq)

for step in range (851,901):
    filename = 'PPC100/flds.tot.%05d'%step
    fields = isolde.getFields(filename)

    Ex = fields['ex'][0,0,:]
    efft = np.fft.fft(Ex)
    efft = np.fft.fftshift(efft)
    e2fft = e2fft+efft*np.conj(efft)

e2fft_PPC100 = e2fft/50.0

step = 850
filename = 'PPC10_cfl01/flds.tot.%05d'%step
fields = isolde.getFields(filename)

Ex = fields['ex'][0,0,:]
efft = np.fft.fft(Ex)
efft = np.fft.fftshift(efft)
e2fft = efft*np.conj(efft)
n = Ex.size
freq = np.fft.fftfreq(n, d=0.1)
freq = np.fft.fftshift(freq)

print(freq)
for step in range (851,901):
    filename = 'PPC10_cfl01/flds.tot.%05d'%step
    fields = isolde.getFields(filename)

    Ex = fields['ex'][0,0,:]
    efft = np.fft.fft(Ex)
    efft = np.fft.fftshift(efft)
    e2fft = e2fft+efft*np.conj(efft)

e2fft_PPC10_cfl01 = e2fft/50.0

##########

directory='../Loading400pmin014'
baseName = "rt_vlasov_libby_pert_1x1v_p2"
e2 = pg.data.GData(directory+'/'+baseName+"-field-energy.gkyl")
interpGrid_ex2_400pmin014=e2.get_grid()
ex2_value_400pmin014=e2.get_values()
sigma400pmin014=5*len(ex2_value_400pmin014)/45000

step = 850
e = pg.data.GData(directory+'/'+baseName+'-field_'+str(step)+'.gkyl')
e_Interp = pg.data.GInterpModal(e, 2, 'ms')
interpGrid_ex_400pmin014, ex_value400pmin014 = e_Interp.interpolate(0)

Ex = ex_value400pmin014.reshape(len(ex_value400pmin014))

n = Ex.size
print("spacing ", interpGrid_ex_400pmin014[0][1]-interpGrid_ex_400pmin014[0][0])
freq_gkeyll = np.fft.fftfreq(n, d=interpGrid_ex_400pmin014[0][1]-interpGrid_ex_400pmin014[0][0])
freq_gkeyll = np.fft.fftshift(freq_gkeyll)

efft = np.fft.fft(Ex)
efft = np.fft.fftshift(efft)
e2fft = efft*np.conj(efft)

for step in range (851,901):
    e = pg.data.GData(directory+'/'+baseName+'-field_'+str(step)+'.gkyl')
    e_Interp = pg.data.GInterpModal(e, 2, 'ms')
    interpGrid_ex_400pmin014, ex_value400pmin014 = e_Interp.interpolate(0)
    Ex = ex_value400pmin014.reshape(len(ex_value400pmin014))
    
    efft = np.fft.fft(Ex)
    efft = np.fft.fftshift(efft)
    e2fft = e2fft+efft*np.conj(efft)
    
e2fft_400pmin014 = e2fft/50.0

directory='../Loading800pmin014'
baseName = "rt_vlasov_libby_pert_1x1v_p2"
e2 = pg.data.GData(directory+'/'+baseName+"-field-energy.gkyl")
interpGrid_ex2_800pmin014=e2.get_grid()
ex2_value_800pmin014=e2.get_values()
sigma800pmin014=5*len(ex2_value_800pmin014)/45000

step = 850
e = pg.data.GData(directory+'/'+baseName+'-field_'+str(step)+'.gkyl')
e_Interp = pg.data.GInterpModal(e, 2, 'ms')
interpGrid_ex_800pmin014, ex_value800pmin014 = e_Interp.interpolate(0)

Ex = ex_value800pmin014.reshape(len(ex_value800pmin014))

n = Ex.size

efft = np.fft.fft(Ex)
efft = np.fft.fftshift(efft)
e2fft = efft*np.conj(efft)

for step in range (851,901):
    e = pg.data.GData(directory+'/'+baseName+'-field_'+str(step)+'.gkyl')
    e_Interp = pg.data.GInterpModal(e, 2, 'ms')
    interpGrid_ex_800pmin014, ex_value800pmin014 = e_Interp.interpolate(0)
    Ex = ex_value800pmin014.reshape(len(ex_value800pmin014))
    
    efft = np.fft.fft(Ex)
    efft = np.fft.fftshift(efft)
    e2fft = e2fft+efft*np.conj(efft)
    
e2fft_800pmin014 = e2fft/50.0

directory='../Loading800pmin007'
baseName = "rt_vlasov_libby_pert_1x1v_p2"
e2 = pg.data.GData(directory+'/'+baseName+"-field-energy.gkyl")
interpGrid_ex2_800pmin007=e2.get_grid()
ex2_value_800pmin007=e2.get_values()
sigma800pmin007=5*len(ex2_value_800pmin007)/45000

step = 850
e = pg.data.GData(directory+'/'+baseName+'-field_'+str(step)+'.gkyl')
e_Interp = pg.data.GInterpModal(e, 2, 'ms')
interpGrid_ex_800pmin007, ex_value800pmin007 = e_Interp.interpolate(0)

Ex = ex_value800pmin007.reshape(len(ex_value800pmin007))

n = Ex.size

efft = np.fft.fft(Ex)
efft = np.fft.fftshift(efft)
e2fft = efft*np.conj(efft)

for step in range (851,901):
    e = pg.data.GData(directory+'/'+baseName+'-field_'+str(step)+'.gkyl')
    e_Interp = pg.data.GInterpModal(e, 2, 'ms')
    interpGrid_ex_800pmin007, ex_value800pmin007 = e_Interp.interpolate(0)
    Ex = ex_value800pmin007.reshape(len(ex_value800pmin007))
    
    efft = np.fft.fft(Ex)
    efft = np.fft.fftshift(efft)
    e2fft = e2fft+efft*np.conj(efft)
    
e2fft_800pmin007 = e2fft/50.0


elc = pg.data.GData(directory+'/'+baseName+'-elc_avg_'+str(step)+'.gkyl')
fp_elc = elc.get_values()
interpGrid_elc = elc.get_grid()

pos = pg.data.GData(directory+'/'+baseName+'-ion_avg_'+str(step)+'.gkyl')
fp_pos = pos.get_values()

data_vmap = pg.GData(directory+'/'+baseName+"-elc_vmap_avg.gkyl")
coords_vmap = data_vmap.get_values()

dist_no_jacob = fp_elc[:,:,0]+fp_pos[:,:,0]
dist_no_jacob_E = fp_elc[:,:,0]
dist_no_jacob_P = fp_pos[:,:,0]

print(freq_gkeyll[int(len(freq_gkeyll)/2)])
print(freq[int(len(freq)/2)])

fig1 = plt.figure(dpi=300, figsize=(6,4), facecolor='white')
ax1 = fig1.add_axes([0.12,0.12,0.83,0.82])
ax1.plot(freq[int(len(freq)/2)+1:],e2fft_PPC10[int(len(freq)/2)+1:]/np.sum(e2fft_PPC10[int(len(freq)/2)+1:]), linestyle='--', label=r'$N_{ppc}=10$')
ax1.plot(freq[int(len(freq)/2)+1:],e2fft_PPC100[int(len(freq)/2)+1:]/np.sum(e2fft_PPC100[int(len(freq)/2)+1:]), linestyle='--', label=r'$N_{ppc}=100$')
ax1.plot(freq[int(len(freq)/2)+1:],e2fft_PPC10_cfl01[int(len(freq)/2)+1:]/np.sum(e2fft_PPC10_cfl01[int(len(freq)/2)+1:]), linestyle='--', label=r'$N_{ppc}=10$, cfl=0.1')
ax1.plot(freq_gkeyll[int(len(freq_gkeyll)/2)+1:],e2fft_400pmin014[int(len(freq_gkeyll)/2)+1:]/np.sum(e2fft_400pmin014[int(len(freq_gkeyll)/2)+1:]),label=r'$N_u=400, u_{min}=0.14$')
ax1.plot(freq_gkeyll[int(len(freq_gkeyll)/2)+1:],e2fft_800pmin014[int(len(freq_gkeyll)/2)+1:]/np.sum(e2fft_800pmin014[int(len(freq_gkeyll)/2)+1:]),label=r'$N_u=800, u_{min}=0.14$')
ax1.plot(freq_gkeyll[int(len(freq_gkeyll)/2)+1:],e2fft_800pmin007[int(len(freq_gkeyll)/2)+1:]/np.sum(e2fft_800pmin014[int(len(freq_gkeyll)/2)+1:]),label=r'$N_u=800, u_{min}=0.07$')
ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_xlabel(r"$k d_e$", fontsize=label_size)
ax1.set_xlim((0.05, 3.0))
ax1.set_ylim((1.0e-7, 1.0))
ax1.set_yticks([1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0e0])
ax1.set_ylabel(r"$E_k^2$", fontsize=label_size)
ax1.tick_params(labelsize=font_size)
ax1.legend(framealpha=0, fontsize=13, loc='lower left')

fig1.savefig("FieldEnergy_kspace.png", dpi=300, bbox_inches='tight')
plt.close(fig1)

fig2 = plt.figure(dpi=300, figsize=(6,4), facecolor='white')
ax3 = fig2.add_axes([0.12,0.12,0.83,0.82])
ax3.plot(time10,gaussian_filter(e2_ppc10,sigma=5)/e2_ppc10[0], linestyle='--', label=r'$N_{ppc}=10$')
ax3.plot(time100,gaussian_filter(e2_ppc100,sigma=5)/e2_ppc100[0], linestyle='--', label=r'$N_{ppc}=100$')
ax3.plot(time10cfl01,gaussian_filter(e2_ppc10_cfl01,sigma=20)/e2_ppc10_cfl01[0], linestyle='--', label=r'$N_{ppc}=10$, cfl=0.1')
ax3.plot(interpGrid_ex2_400pmin014[0]/20,gaussian_filter(ex2_value_400pmin014[:,0],sigma=sigma400pmin014)/ex2_value_400pmin014[0,0],label=r'$N_u=400, u_{min}=0.14$')
ax3.plot(interpGrid_ex2_800pmin014[0]/20,gaussian_filter(ex2_value_800pmin014[:,0],sigma=sigma800pmin014)/ex2_value_800pmin014[0,0],label=r'$N_u=800, u_{min}=0.14$')
ax3.plot(interpGrid_ex2_800pmin007[0]/20,gaussian_filter(ex2_value_800pmin007[:,0],sigma=sigma800pmin007)/ex2_value_800pmin007[0,0],label=r'$N_u=400, u_{min}=0.07$')
#ax3.set_yscale("symlog",linthresh=1e-4)
ax3.set_ylim(1e-7,1e-3)
ax3.set_yticks([1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3])
ax3.set_xlim(0,90)
ax3.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
ax3.set_yscale("log")
ax3.set_xlabel(r"$t~[L_x/c]$", fontsize=label_size)
ax3.set_ylabel(r"$E_x^2(t)/E_x^2(0)$", fontsize=label_size)
ax3.tick_params(labelsize=font_size)
ax3.legend(
    loc='lower left',
    bbox_to_anchor=(0.0, -0.01),
    borderaxespad=0.05,
    framealpha=0,
    fontsize=13,
    ncol=2,
)

ax3_pos = ax3.get_position().bounds
ax3_pos = (ax3_pos[0] + 0.01, ax3_pos[1] + 0.015, ax3_pos[2], ax3_pos[3])
inset_width = 0.45 * ax3_pos[2]
inset_height = 0.45 * ax3_pos[3]
inset_pad_x = (0.02 / 0.6) * ax3_pos[2]
inset_pad_y = (0.02 / 0.4) * ax3_pos[3]
ax31 = fig2.add_axes([
    ax3_pos[0] + ax3_pos[2] - inset_pad_x - inset_width,
    ax3_pos[1] + ax3_pos[3] - inset_pad_y - inset_height,
    inset_width,
    inset_height,
])
ax31.plot(time10,gaussian_filter(e2_ppc10,sigma=5)/e2_ppc10[0], linestyle='--', label=r'$N_{ppc}=10$')
ax31.plot(time100,gaussian_filter(e2_ppc100,sigma=5)/e2_ppc100[0], linestyle='--', label=r'$N_{ppc}=100$')
ax31.plot(time10cfl01,gaussian_filter(e2_ppc10_cfl01,sigma=20)/e2_ppc10_cfl01[0], linestyle='--', label=r'$N_{ppc}=10$, cfl=0.1')
ax31.plot(interpGrid_ex2_400pmin014[0]/20,gaussian_filter(ex2_value_400pmin014[:,0],sigma=sigma400pmin014)/ex2_value_400pmin014[0,0],label=r'$N_u=400, u_{min}=0.14$')
ax31.plot(interpGrid_ex2_800pmin014[0]/20,gaussian_filter(ex2_value_800pmin014[:,0],sigma=sigma800pmin014)/ex2_value_800pmin014[0,0],label=r'$N_u=800, u_{min}=0.14$')
ax31.plot(interpGrid_ex2_800pmin007[0]/20,gaussian_filter(ex2_value_800pmin007[:,0],sigma=sigma800pmin007)/ex2_value_800pmin007[0,0],label=r'$N_u=800, u_{min}=0.07$')
#ax31.set_ylim(5e-7,1e-3)
ax31.set_xscale("log")
ax31.set_yscale("log")
ax31.set_xlim(0.09,90)
ax31.set_xlabel(r"$t~[L_x/c]$", fontsize=10)
ax31.set_ylabel(r"$E_x^2(t)/E_x^2(0)$", fontsize=10)
ax31.tick_params(labelsize=10)

fig2.savefig("FieldEnergy_time.png", dpi=300, bbox_inches='tight')
plt.close(fig2)

fig3 = plt.figure(dpi=300, figsize=(6,4), facecolor='white')
ax4 = fig3.add_axes([0.12,0.12,0.83,0.82])
ax4.plot(coords_vmap[:,0], np.mean(dist_no_jacob_E[:,:],axis=0),label='electrons')
ax4.plot(coords_vmap[:,0], np.mean(dist_no_jacob_P[:,:],axis=0),label='positrons')
ax4.set_yscale("log")
ax4.set_ylim(bottom=1e-3)
ax4.set_xscale("symlog", linthresh=1)
ax4.legend(framealpha=0, loc='upper left', fontsize=16)
ax4.set_xlabel(r"$u_x$", fontsize=label_size)
ax4.set_ylabel(r"$f(u_x)$", fontsize=label_size)
ax4.tick_params(labelsize=font_size)
ax4.set_ylim(top=1e3)
ax4.set_yticks([1.0e-3, 1.0e-2, 1.0e-1, 1.0e0, 1.0e1, 1.0e2, 1.0e3])

fig3.savefig("FieldEnergy_dist.png", dpi=300, bbox_inches='tight')
plt.close(fig3)
