import numpy as np

# Initial aerosol size distribution
# geometric standard deviation of the mode
sigma_g = 1.7
# geometric mean diameter of the mode [m]
mu = 150.0 * 1.0e-9
# total number of particles in the mode [m-3]
N_T = 1.0e8

# initial ambient conditions
# relative humidity [%]
RH0 = 95.0
# temperature [K]
T0 = 298.15
# pressure [Pa]
p0 = 101325.0
# updraft velocity for cloud parcel runs [m]
w = 0.
# integration time for the simulation
integration_time = 500
# Number of size bins
Nb = 50

# Available remapping schemes:
# 'hybrid-bin','quasi-stationary','moving-center','flat-top','full-moving'
remapping_scheme=['hybrid-bin','quasi-stationary','moving-center','flat-top','full-moving']
# flag for different model setups
# 1 = semivolatile organics condensation
# 2 = semivolatile organics condensation + nucleation
# 3 = water condensation
# 4 = water and semivolatile organics condensation
model_setup = 1

# ODE parameters
rtol=1.e-9           # Relative tolerance for VODE
atol=1e-30           # Absolute tolerance for VODE
dt=1.                # Time-step length (s)

# Physical constants
R_gas=8.31451        # Gas constant
g=9.8065             # Gravity constant
k_b=1.3806504e-23    # Boltzmann constant
c_p_water = 4185.8   # Heat capacity of water ( J/(kg K) )
c_p_air = 1004.67    # Heat capacity of air ( J/(kg K))
M_air = 0.02897      # Molar mass of dry air (kg / mol)
sigma = 0.07275      # surface tension of pure water[J m-2]@ 293 K
NA=6.0221367e+23   # Avogadro's constant [mol-1]  
N_low_limit=1.e-5    # Minimum number concentration in a bin (# / m3)


# Define species
# number of non-volatile species
num_nvbs_species=3
# number of VBS species
num_vbs_species=10
# List of all species
All_species=["H2O", "SO4"] + ['VBS{}'.format(i) for i in range(num_vbs_species)]

# List of Condensing species
Condensing_species=['H2O', 'SO4'] + ['VBS{}'.format(i) for i in range(num_vbs_species)]
#Condensing_species=['H2O', 'SO4']
# Compound indices
index ={
    "H2O"   : 1,
    "SO4"   : 2,
    }
# "VBS0","VBS1","VBS2","VBS3","VBS4","VBS5","VBS6","VBS7","VBS8","VBS9"
index.update({'VBS{}'.format(i) : i+3 for i in range(num_vbs_species)})

# molar weight [kg mol-1]
M={
    "H2O"   : 18.015e-3,
    "SO4"   : 98.08e-3,
    }
M.update({'VBS{}'.format(i) : 200e-3 for i in range(num_vbs_species)})

# density [kg m-3]
rho={
    "H2O"   : 995.,
    "SO4"   : 1830.,
    }
rho.update({'VBS{}'.format(i) : 2000. for i in range(num_vbs_species)})

# diffusivity of each species i
D_g={
    "H2O"   : 1.9e-4*(M['H2O']*1000.0)**(-2./3.),
    "SO4"   : 1.9e-4*(M['SO4']*1000.0)**(-2./3.),
    }
D_g.update({'VBS{}'.format(i) : 1.9e-4*(M['VBS{}'.format(i)]*1000.0)**(-2./3.) for i in range(num_vbs_species)})
    
# saturation vapor pressure of each species i in mol/m3
C_star={
    "H2O"   : lambda T : np.exp(20.386 - 5132 / T) / (R_gas * T) * 133.32,
    "SO4"   : lambda T : 0.0,
    }
C_star.update({
    'VBS{}'.format(i) : lambda T, i=i : 10**(i-6.0)/M['VBS{}'.format(i)] * 1.e-9 for i in range(num_vbs_species)})

# mass accommodation coefficient
alpha_d={
    "H2O"   : 1.0,
    "SO4"   : 1.0,
    }
alpha_d.update({'VBS{}'.format(i) : 1.0 for i in range(num_vbs_species)})

# hygroscopicity
kappa_i={
    "H2O"   : 0.0,
    "SO4"   : 0.57,
    }
kappa_i.update({'VBS{}'.format(i) : 0.1 for i in range(num_vbs_species)})

abundance={
    "SO4"  : 0.0001,
    "VBS0" : 0.1,
    "VBS1" : 0.1,
    "VBS2" : 0.15,
    "VBS3" : 0.22,
    "VBS4" : 0.36,
    "VBS5" : 0.47,
    "VBS6" : 0.58,
    "VBS7" : 0.69,
    "VBS8" : 0.84,
    "VBS9" : 1.0
}

# Latent heat of evaporation
L_e={
    "H2O"   : lambda T : 2.501e6-2370.0*(T-273.15)
    }

