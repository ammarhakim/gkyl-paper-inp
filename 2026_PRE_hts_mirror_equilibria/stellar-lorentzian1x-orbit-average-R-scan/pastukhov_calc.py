import numpy as np
from scipy.integrate import quad
from scipy.special import erf
#[ Append postgkyl wrappers.
from scipy import special
from scipy import optimize


# Calculate the Pastukhov potential confinement time from an adiabatic electron simulation
Bmag0 = 0.5273183
BmagThroat = 3.147247
BmagExpander = 0.1343069

UparMirror = 3.928870e+05
UparExpander = 1.695508e+06



R = BmagThroat / Bmag0
n0 = 9.640995e+18 # Read from the simulation at midplane

intM0Sdx = 1.711530e+20
intM0dx = 1.48e20

nu_ii = 2.972 # Read from the simulation at midplane

eps0, mu0 = 8.8541878176204e-12, 1.2566370614359e-06
eV        = 1.602176487e-19
qe, qi    = -1.602176487e-19, 1.602176487e-19
me, mp    = 9.10938215e-31, 1.672621637e-27

mi        = 2.014*mp                         #[ Deuterium ion mass.
Te0       = 940*eV

#[ Electron-electron collision freq.
logLambdaElc = 6.6 - 0.5*np.log(n0/1e20) + 1.5*np.log(Te0/eV)
nuElc        = logLambdaElc*(eV**4)*n0/(6*np.sqrt(2)*(np.pi**(3/2))*(eps0**2)*np.sqrt(me)*(Te0**(3/2)))

def Pastukhov_confinement_time(x,R,frac):
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
    print(Loss_Rosen)
    return Loss_Rosen

def tau_pi(R):
  #[ Ion confinement time as a function of mirror ratio R.
  # tau_i = 1/nu_ii
  # tau_pi = tau_i*np.log10(R) # Mana's paper
  # tau_pi = 2.6 * tau_i * np.log(R) #Novatron paper
  # tau_pi = 2.4 * tau_i * np.log10(R) # What I think Baldwin means
  # print(f"tau_pi = {tau_pi:.3f} s (from tau_i and log10(R))")
  # print(f"tau_pi = {tau_pi:.3f} s")
  # tau_pi = n0 / dnidt
  # print(f"tau_pi = {tau_pi:.3f} s (from dnidt)")
  # print(f"tau_pi = {tau_pi:.3f} s (from dnidt)")
  tau_pi = intM0dx / intM0Sdx
  return tau_pi

def rootEq_dough(x,R):
  return tau_pi(R)-Rosen_Dougherty_confinement_time(x,R,1,0.9)

def rootEq_past(x,R):
  return tau_pi(R)-Pastukhov_confinement_time(x,R,1.0)

def rootEq_najmabadi(x,R):
  return tau_pi(R)-Najmabadi_confinement_time(x,R)

print(rootEq_dough(1,R))
print(rootEq_dough(20,R))

ephi_over_Te = optimize.ridder(rootEq_dough, 1, 20, args=(R))

print(f"Dougherty e*phi/Te = {ephi_over_Te:.3f}")

ephi_over_Te = optimize.ridder(rootEq_past, 1, 20, args=(R))
print(f"Pastukhov e*phi/Te = {ephi_over_Te:.3f}")

ephi_over_Te = optimize.ridder(rootEq_najmabadi, 1, 20, args=(R))
print(f"Najmabadi e*phi/Te = {ephi_over_Te:.3f}")

expander_potential_drop = np.log(BmagThroat*UparExpander/BmagExpander/UparMirror)
print(f"Expander potential drop = {expander_potential_drop:.3f} Te")