# Import the relevant libraries, including an ODE solver from Scipy
# We only want to use the 'odeint' [Solve Initial Value Problems] from the
# scipy.integrate package
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt # Import Matplotlib so we can plot results
import pdb
import timeit # for timing
import numba as nb

# --- Physical constants --------------------------------------------
Lv_water_vapour=2.5e3 # Latent heat of vapourisation of water [J/g]
Rv=461.0 #Individual gas constant of water vapour [J/Kg.K]
Ra=287.0 #Gas constant for dry air [J/Kg.K]
R_gas=8.3144598 #Ideal gas constant [kg m2 s-2 K-1 mol-1]
R_gas_other=8.2057e-5 #Ideal gas constant [m3 atm K-1 mol-1]
GRAV=9.8; #Gravitational acceleration [(m/2)2]
cp=1005; #Specific heat capacity of air [J/Kg.K]
sigma=72.0e-3 # Assume surface tension of water (mN/m)
NA=6.0221409e+23 #Avogadros number
kb=1.380648E-23 #Boltzmanns constant
# ---------------------------------------------------------------------

Temp_K=298.15

# ------- Defining the gaseous and condensed phase components ---------
num_species = 10 # Number of condensing species from the gas phase
num_bins = 1 # Number of size bins

# The molecular weight of each condensing specie [g/mol]
mw_array=np.zeros((num_species), dtype=float)
mw_array[:]=200.0 # Assuming a constant value for all components

# The volatility of each specie, using the mass based C* convention
# Here we assuming a linear seperation in Log10 space
log_c_star = np.linspace(-6, 3, num_species)
Cstar = np.power(10.0,log_c_star)
# Convert C* to a pure component saturation vapour pressure [atm]
P_sat = (Cstar*R_gas_other*Temp_K)/(1.0e6*mw_array)

# Initialise an abundance of material in each volatility bin [micrograms/m3]
abundance = np.zeros((num_species), dtype = float)
abundance[0] = 0.1
abundance[1] = 0.1
abundance[2] = 0.15
abundance[3] = 0.22
abundance[4] = 0.36
abundance[5] = 0.47
abundance[6] = 0.58
abundance[7] = 0.69
abundance[8] = 0.84
abundance[9] = 1.0

# Unit conversion of gas abudance to molecules / cc
gas_abundance = ((abundance*1.0e-6)/(mw_array))*1.0e-6*NA

# Assumed accomodation coefficient
alpha_d_org=np.zeros((num_species), dtype=float)
alpha_d_org[:]=1.0
# Density of condensing species [kg/m3]
density_org=np.zeros((num_species), dtype=float)
density_org[:]=1400.0

# Molecular diffusion coefficient of each molecule in air. [cm2/s]
DStar_org = 1.9*np.power(mw_array,-2.0/3.0)
# Mean thermal velocity of each molecule [m/s] -
mean_them_vel=np.power((8.0*R_gas*Temp_K)/(np.pi*(mw_array*1.0e-3)),0.5)
# Mean free path for each molecule [m]
gamma_gas = ((3.0*DStar_org)/(mean_them_vel*1.0e2))*1.0e-2
# ---------------------------------------------------------------------

# ----------- Defining a monodisperse size distribution ---------------
# Create a monodisperse distribution with a specific number of particles
# Assume each particle starts with an involatile core of absorptive organic with
# a mass of 200 g/mol and density 1400 km/m3. Store this information in an array
# as 'core'. This will ensure, in this example, that we do not get 100%
# evaporative loss

# Define total number of particles [per cc]
N_total = 300.0

# We carry an inert and involatile in each size bin
core = np.zeros((num_bins), dtype=float)
core_abundance = np.zeros((num_bins), dtype=float)
density_core = np.zeros((num_bins), dtype=float)
core_mw = np.zeros((num_bins), dtype=float)
density_core[:] = 1400.0
core_mw[:] = 200.0

# The number of particles is only carried in one bin
N_per_bin = np.zeros((num_bins), dtype=float)
N_per_bin[0] = N_total

# Define the diameter of our particles [microns]
size_array = np.zeros((num_bins), dtype=float)
size_array[0] = 0.150 #microns

# Use the size to now calculate a concentration of a 'core' in molecules / cc
# This aligns with the units used for volatile species
core_abundance[0] = (N_per_bin[0])*((4.0/3.0)*np.pi*np.power(size_array[0]*0.5e-6,3.0)*density_core[0]*1.0e3)
core_abundance[0] = (core_abundance[0] / core_mw[0])*NA

# What is our existing dry mass in micrograms per cubic metre?
dry_mass = np.sum((core_abundance/NA)*core_mw)*(1.0e12)
print("Initial dry mass = ", dry_mass)

# ---------------------------------------------------------------------

# ------ Defining the gaseous and condensed concentration array -------
# Here we initialise the array used throughout the simulation and
# populate with initial concentrations in both the gas and condensed phase
array = np.zeros((num_species+num_species*num_bins), dtype=float)
# The first num_species cells hold the gas phase abundance
array[0:num_species] = gas_abundance
# The following cells hold the concentration of each specie in each size bin
array[num_species:num_species+num_species*num_bins] = 1.0e-20
# ---------------------------------------------------------------------

# --------------- Defining the Droplet Growth equation ----------------
# Here we use variable names as close to the equation syntax as possible
# Define the function name dy_dt and expected inputs

@nb.jit(nb.float64[:](nb.float64[:],nb.float64), nopython=True, cache=True)
def dy_dt(array,t):

    # Retrieve the gas phase concentrations
    Cg_i_m_t = array[0:num_species]
    # Retrieve the concentrations in our size bin as a slice of array
    temp_array=array[num_species:num_species+num_species*num_bins]
    # Sum the total molecules in the condensed phase, adding core material
    total_molecules=np.sum(temp_array)+core_abundance[0]
    # Calculate the mole fractions in the, assumed, liquid phase for use
    # in calculating the equilibrium pressure above the droplet
    mole_fractions=temp_array/total_molecules
    # Calculate the density of the assumed liquid phase
    density_array = np.zeros((num_species+num_bins), dtype=nb.float64)
    density_array[0:num_species]=density_org[0:num_species]
    density_array[num_species]=density_core[0]
    # Create an array that holds mass concentrations [g/cc], used for
    # calculation of solution density [kg/m3]
    mass_array = np.zeros((num_species+1), dtype=nb.float64)
    mass_array[0:num_species] = (temp_array/NA)*mw_array
    mass_array[num_species] = (core_abundance[0]/NA)*core_mw[0]
    total_mass=np.sum(mass_array)
    mass_fractions_array=mass_array/total_mass
    density=1.0/(np.sum(mass_fractions_array/density_array))

    # Now calculate the size [cm] of our particles according to the condensed
    # phase abundance of material
    # We need to remember that mass is in [g/cc] whilst density
    # is in [kg/m3]. Thus we convert mass and number concentrations to kg/m3 and /m3
    size=((3.0*((total_mass*1.0e3)/(N_per_bin[0]*1.0e6)))/(4.0*np.pi*density))**(1.0/3.0)

    # Calculate the Knudsen number for all condensing molecules based on this new size
    # This relies on mean free path for each species [cm] and particle radius [cm]
    Kn=gamma_gas/size
    # Non-continuum regime correction  \n')
    # Calculate a correction factor according to the continuum versus
    # non-continuum regimes
    Inverse_Kn=1.0/Kn
    Correction_part1=(1.33e0+0.71e0*Inverse_Kn)/(1.0e0+Inverse_Kn)
    Correction_part2=(4.0e0*(1.0e0-alpha_d_org))/(3.0e0*alpha_d_org)
    Correction_part3=1.0e0+(Correction_part1+Correction_part2)*Kn
    Correction=1.0/Correction_part3
    # Kelvin factor
    # Now calculate a kelvin factor for every semi-volatile compound in this
    # size bin
    kelvin_factor=np.exp((4.0E0*mw_array*1.0e-3*sigma)/(R_gas*Temp_K*size_array[0]*2.0e0*density))
    # Equilibrium pressure above droplets
    # Now calculate an equilibrium pressure [Pascals] for every compound in this
    # size bin
    #R_gas_other=8.2057e-5 #Ideal gas constant [m3 atm K-1 mol-1]
    Pressure_eq=kelvin_factor*mole_fractions*P_sat
    # Calculate the equilibrium concentration [molecules/cc] equivalent of this
    # pressure
    Cstar_i_m_t=Pressure_eq*(NA/(R_gas_other*1.0e6*Temp_K))

    # Implement the droplet growth equation
    # The following calculates the change in mass, per particle of this size.
    # The equation relies on the following parameters \n')
    # and units.
    # radius [m] \n')
    # DStar_org [m2/s] - pass in [cm2/s] so needs to be converted by *1.0E-4 \n')
    # Pressure_gas [Pascals] \n')
    # Pressure_eq [Pascals] \n')
    # R_gas [m3 Pascals /(K mol)] \n')
    # molw [g/mol] \n')
    # T [K] \n')
    # The units of the equation should therefore be g/s \n')
    k_i_m_t_part1 = DStar_org*Correction
    k_i_m_t=4.0e0*np.pi*size*1.0e2*N_per_bin[0]*k_i_m_t_part1
    dm_dt=k_i_m_t*(Cg_i_m_t-Cstar_i_m_t)

    # Here we store values of dm_dt in a matrix, ready for dealing with polydisperse
    # populations where a contribution/loss per size bin
    dy_dt_gas_matrix = np.zeros((num_species,num_bins), dtype=nb.float64)
    dy_dt_gas_matrix[0:num_species,0]=dm_dt

    # Now we initiase an array that passes back a set of derivatives for Each
    # compound in both the gaseous and condensed phases
    dy_dt_array = np.zeros((num_species+num_species*num_bins), dtype=nb.float64)
    dy_dt_array[num_species:num_species+num_species*num_bins]=dm_dt[0:num_species]
    dy_dt_array[0:num_species]=dy_dt_array[0:num_species]-np.sum(dy_dt_gas_matrix, axis=1)

    return dy_dt_array


# Define the time over which the simulation will take place
t = np.linspace(0, 10000, num=1000)

# Call the ODE solver with reference to our function, dy_dt, setting the
# absolute and relative tolerance.

# Do timing for numba: run once to compile and then time next 100 instances
def find_solution():
    solution = odeint(dy_dt, array, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)
find_solution
print("Time using numba: ", timeit.Timer(find_solution).timeit(number=100))
# ignore output since we continually integrating forward over t
# no plot required
