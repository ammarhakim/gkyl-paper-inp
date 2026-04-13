import numpy as np
from scipy.integrate import quad
from scipy.special import erf
#[ Append postgkyl wrappers.
from scipy import special
from scipy import optimize
import postgkyl as pg


# Calculate the Pastukhov potential confinement time from an adiabatic electron simulation
# directory = '/scratch/gpfs/mr1884/scratch/gkylmax/stellar-lorentzian1x-exploration/stellar-lorentzian1x-orbit-average-nu-ie-maxwellian-odist-time-dilation-pos-extra-tdil20'
directory = '/home/mr1884/scratch/gkylmax/stellar-lorentzian1x-orbit-average'
sim_name = 'gk_lorentzian_mirror'
species = 'ion'
frame_num = 65
source_frame = 0
z_mirror_throat = 0.98
simulation_edge = 2.5


# Extract details about the magnetic field
bmag_data = pg.GData(	directory + f'/{sim_name}-bmag.gkyl',	mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
pg.GInterpModal(bmag_data).interpolate(0, overwrite=True)
_, bmag_sel = pg.data.select(bmag_data, z0='0.0')
Bmag0 = float(np.squeeze(bmag_sel))
_, bmag_sel = pg.data.select(bmag_data, z0=f'{z_mirror_throat}')
BmagThroat = float(np.squeeze(bmag_sel))
_, bmag_sel = pg.data.select(bmag_data, z0=f'{simulation_edge}')
BmagExpander = float(np.squeeze(bmag_sel))
R = BmagThroat / Bmag0  # Read from postgkyl

# Measure values of the bimaxwellian moments (midplane density, upar at the throat, upar at the wall)
bimax_data = pg.GData(	directory + f'/{sim_name}-{species}_BiMaxwellianMoments_{frame_num}.gkyl',	mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
pg.GInterpModal(bimax_data).interpolate(comp=slice(0, 3), overwrite=True)
_, upar_sel = pg.data.select(bimax_data, comp=1, z0=f'{z_mirror_throat}')
UparMirror = float(np.squeeze(upar_sel))
_, upar_sel = pg.data.select(bimax_data, comp=1, z0=f'{simulation_edge}')
UparExpander = float(np.squeeze(upar_sel))
_, n0_sel = pg.data.select(bimax_data, comp=0, z0='0.0')
n0 = float(np.squeeze(n0_sel))

# Find the integrated density inside the trap
m0_data = pg.GData(	directory + f'/{sim_name}-{species}_M0_{frame_num}.gkyl',	mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
pg.GInterpModal(m0_data).interpolate(comp=0, overwrite=True)
pg.data.select(m0_data, z0=f'{-z_mirror_throat}:{z_mirror_throat}', overwrite=True)
_, intM0dx = pg.tools.integrate(m0_data, 0, overwrite=True)
intM0dx = float(np.squeeze(intM0dx))

# Find the integrated source inside the trap
m0src_data = pg.GData(	directory + f'/{sim_name}-{species}_source_M0_{source_frame}.gkyl',	mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
pg.GInterpModal(m0src_data).interpolate(comp=0, overwrite=True)
pg.data.select(m0src_data, z0=f'{-z_mirror_throat}:{z_mirror_throat}', overwrite=True)
_, intM0Sdx = pg.tools.integrate(m0src_data, 0, overwrite=True)
intM0Sdx = float(np.squeeze(intM0Sdx))

eps0, mu0 = 8.8541878176204e-12, 1.2566370614359e-06
eV        = 1.602176487e-19
qe, qi    = -1.602176487e-19, 1.602176487e-19
me, mp    = 9.10938215e-31, 1.672621637e-27

mi        = 2.014*mp                         #[ Deuterium ion mass.
Te0       = 940*eV

Zpfl = 1.0 # Z in Najmabadi (1984). ee collisions would be 1.0, ee+ei collisions would be 1.5 (lambda_ei = lambda_ee, n_i = n_e)

# Electron-electron collision freq.
n0_cm3 = n0 * 1e-6
Te0_eV = Te0 / eV
logLambdaElc = 24 - np.log( n0_cm3**0.5 / Te0_eV ) # https://farside.ph.utexas.edu/teaching/plasma/Plasma/node39.html
nuElc = logLambdaElc * n0 * (4*np.pi / (2**1.5 * Te0**1.5 * me**0.5)) * (eV**2/(4*np.pi*eps0))**2 # Najmabadi (1984)

def tau_pe(x,R,frac):
  #[ Electron confinement time as a function of x=e*phi/Te and mirror ratio R.
  #[ This has a 1/4 (Cohen) instead of a 1/2 (Pastukhov). However, Cohen assumes
  #[ ee and ei collisions. If we want to consider only ee collisions use this
  #[ function with frac=0.5 (essentially turning the 1/4 into 1/2).
  def G(R):
    return ((2.*R+1)/(2.*R))*np.log(4.*R+2)
  def I_x(x):
    return 1.+0.5*np.sqrt(np.pi*x)*np.exp(1./x)*special.erfc(np.sqrt(1./x))
  return (np.sqrt(np.pi)/4.)*(1./(frac*nuElc))*G(R)*x*np.exp(x)/I_x(1./x)

# Compute analytical estimate from Najmabadi. It's more easily implemented from Post 1987
def Najmabadi_confinement_time(P, R, ZpFl=1, coeff = 0.84):
    w_term = np.sqrt(1 + 1/(R*(ZpFl - 1/(4*P)))) #Is this 1/4P or P/4
    u_eff  = P + np.log(w_term) #same as Najmabadi's a**2

    integrandNaj = lambda t: np.exp(-t) / t
    I_term, error = quad(integrandNaj, u_eff, np.inf)
    I_term = (ZpFl + 1/4)*u_eff*np.exp(u_eff)*I_term - 1/4 

    Loss_Najmabadi = 1/nuElc * \
        np.sqrt(np.pi)/4 * u_eff*np.exp(u_eff)/I_term \
        * (np.log((w_term+1)/(w_term-1)) - coeff)#0.84
    
    return Loss_Najmabadi


def Rosen_Dougherty_confinement_time(P, R, ZpFl, coeff = 0):
    w_term = np.sqrt(1 + 2*P/(R*ZpFl)) 
    a_term  = np.sqrt(P + np.log(w_term)) 

    Loss_Rosen = 1/nuElc * \
        1/(2*ZpFl / (np.log((w_term+1)/(w_term-1)) - coeff) * (1-erf(a_term)))
    
    return Loss_Rosen

def tau_pi(R):
  tau_pi = intM0dx / intM0Sdx
  return tau_pi

def rootEq_dough(x,R):
  return tau_pi(R)-Rosen_Dougherty_confinement_time(x,R,Zpfl,1.117)

def rootEq_past(x,R):
  return tau_pi(R)-tau_pe(x,R,1.0)

def rootEq_najmabadi(x,R):
  return tau_pi(R)-Najmabadi_confinement_time(x,R)

ephi_over_Te = optimize.ridder(rootEq_dough, Zpfl, 20, args=(R))

print(f"Dougherty e*phi/Te = {ephi_over_Te:.3f}")

ephi_over_Te = optimize.ridder(rootEq_past, Zpfl, 20, args=(R))
print(f"Pastukhov e*phi/Te = {ephi_over_Te:.3f}")

ephi_over_Te = optimize.ridder(rootEq_najmabadi, Zpfl, 20, args=(R))
print(f"Najmabadi e*phi/Te = {ephi_over_Te:.3f}")

expander_potential_drop = np.log(BmagThroat*UparExpander/BmagExpander/UparMirror)
print(f"Expander potential drop = {expander_potential_drop:.3f} Te")

print(f"tau_pi = {tau_pi(R):.3f} s")