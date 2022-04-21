# Import the relevant libraries, including an ODE solver from Scipy
# We only want to use the 'odeint' [Solve Initial Value Problems] from the
# scipy.integrate package
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt # Import Matplotlib so we can plot results
import pdb
import numba as nb

# Define physical constants
const R_gas=8.31451 #Ideal gas constant [kg m2 s-2 K-1 mol-1]
const R_gas_other=8.20573e-5 #Ideal gas constant [m3 atm K-1 mol-1]
const sigma=72.0e-3 # Assume surface tension of water (mN/m)
const NA=6.0221409e+23 #Avogadros number

Temp_K=298.15

# Number of condensing species from the gas phase
num_species = 10

# The molecular weight of each condensing specie [g/mol]
# Assuming a constant value for all components
mw_array=np.zeros((num_species), dtype=float)
mw_array[:]=200.0 # Assuming a constant value for all components

# Define array of log10 C* values
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

# Set accomodation coefficient
alpha_d_org=np.zeros((num_species), dtype=float)
alpha_d_org[:]=1.0
# Set density of condensing species [kg/m3]
density_org=np.zeros((num_species), dtype=float)
density_org[:]=1400.0

# Molecular diffusion coefficient in air (Equation 2.22). [cm2/s]
DStar_org = 1.9*np.power(mw_array,-2.0/3.0)
# Mean thermal velocity of each molecule (Equation 2.21). [m/s]
mean_them_vel=np.power((8.0*R_gas*Temp_K)/(np.pi*(mw_array*1.0e-3)),0.5)
# Mean free path for each molecule (Equation 2.20). [m]
gamma_gas = ((3.0*DStar_org)/(mean_them_vel*1.0e2))*1.0e-2

# Additional array to include hypothetical activity coefficients
gamma_factor = np.zeros((num_species), dtype=float)
gamma_factor[0]= 3.0
gamma_factor[1]= 3.0
gamma_factor[2]= 3.0
gamma_factor[3]= 4.0
gamma_factor[4]= 5.5
gamma_factor[5]= 6.5
gamma_factor[6]= 7.0
gamma_factor[7]= 10.0
gamma_factor[8]= 15.0
gamma_factor[9]= 20.0

# Define a polydisperse size distribution
# Set the smallest size
d1 = 0.01
# Set the largest size
d_Nb = 1.0
# Number of size bins
num_bins = 8
# Volume ratio between bins (Equation 1.17)
V_rat = np.power((d_Nb/d1),3.0/(num_bins-1.0))
# Diameter ratio between bins  (Equation 1.18)
d_rat = V_rat**(1.0/3.0)

# Create an array of diameters
d_i=np.zeros((num_bins), dtype=float)
d_i[0]=d1
for step in range(num_bins):
    if step > 0:
       d_i[step]=d_i[step-1]*d_rat
log_di = np.log(d_i)

# Parameters of log-normal distribution
# Geometric standard deviation
sigmag1 = np.log(1.7)
# Mean particle diameter [150nm]
mean1 = np.log(0.15)
# Calculate the probability density distribution
distribution_1 = (np.exp(-(log_di - mean1)**2 / (2 * sigmag1**2)) / (sigmag1 * np.sqrt(2 * np.pi)))
# Seperate out the probability density function
d_width = d_i*np.power(2,1.0/3.0)*((np.power(V_rat,1.0/3.0)-1.0)/(np.power(1+V_rat,1.0/3.0)))
# Total number of particles [per cm-3]
N_total = 100.0
# Discrete number distribution
N_dist = N_total*(distribution_1*(d_width/d_i))

# Initialise a core abundance using the above size distribution
core = np.zeros((num_bins), dtype=float)
core_abundance = np.zeros((num_bins), dtype=float)
density_core = np.zeros((num_bins), dtype=float)
core_mw = np.zeros((num_bins), dtype=float)
density_core[:] = 1400.0
core_mw[:] = 200.0
N_per_bin = N_dist

# Define size array (radius )
size_array = d_i*0.5
# Use the size to now calculate a concentration of a 'core' in molecules / cc
core_abundance = (N_per_bin)*((4.0/3.0)*np.pi*np.power(size_array*1.0e-6,3.0)*density_core*1.0e3)
core_abundance = (core_abundance / core_mw)*NA
dry_mass = np.sum((core_abundance/NA)*core_mw)*(1.0e6)*(1.0e6)
print("Initial dry mass = ", dry_mass)

# New define an array that holds the molecular abundance of
# each gas and the concentration of each gas in each size bin
array = np.zeros((num_species+num_species*num_bins), dtype=float)
array[0:num_species] = gas_abundance
array[num_species:num_species+num_species*num_bins] = 1.0e-10 # assuming we start with nothing

# Define the RHS function (that includes droplet growth equation)
def dy_dt(array,t):

    # Retrieve the gas phase concentrations
    Cg_i_m_t = array[0:num_species]
    # We are working with 8 size bins, each of which has an involatile core
    size_array = np.zeros((num_bins), dtype=float)
    #Initialise empty dydt arrays
    dy_dt_array = np.zeros((num_species+num_species*num_bins), dtype=float)
    dy_dt_gas_matrix = np.zeros((num_species,num_bins), dtype=float)

    # Now cycle through each size bin
    for size_step in range(num_bins):

        # Select a slice of y that represents this size bin
        temp_array=array[num_species+size_step*num_species:num_species+num_species*(size_step+1)]
        # Sum the total molecules in the condensed phase, adding core material
        total_moles=np.sum(temp_array)+core_abundance[size_step]
        # Calculate the mole fractions in the, assumed, liquid phase for use
        # in calculating the equilibrium pressure above the droplet
        mole_fractions=temp_array/total_moles
        # Calculate the density of the assumed liquid phase
        density_array = np.zeros((num_species+1), dtype=float)
        density_array[0:num_species]=density_org[0:num_species]
        density_array[num_species]=density_core[size_step]
        # Create an array that holds mass concentrations [g/cc], used for
        # calculation of solution density [kg/m3]
        mass_array = np.zeros((num_species+1), dtype=float)
        mass_array[0:num_species] = (temp_array/NA)*mw_array
        mass_array[num_species] = (core_abundance[size_step]/NA)*core_mw[0]
        total_mass=np.sum(mass_array)
        mass_fractions_array=mass_array/total_mass
        density=1.0/(np.sum(mass_fractions_array/density_array))
        # Now calculate the size [cm] of our particles according
        # to the condensed phase abundance of material
        # We need to remember that mass is in [g/cc] whilst density
        # is in [kg/m3].
        # Thus we convert mass and number concentrations to kg/m3 and /m3
        size_array[size_step]=((3.0*((total_mass*1.0e3)/(N_per_bin[size_step]*1.0e6)))/(4.0*np.pi*density))**(1.0/3.0)
        # Calculate the Knudsen number for all condensing molecules based
        # on this new size
        # This relies on mean free path for each species [cm] and
        # particle radius [cm]
        Kn=gamma_gas/size_array[size_step]
        # Calculate Non-continuum regime correction (Equation 2.19)
        # Calculate a correction factor according to the continuum versus
        # non-continuum regimes
        Inverse_Kn=1.0/Kn
        Correction_part1=(1.33e0+0.71e0*Inverse_Kn)/(1.0e0+Inverse_Kn)
        Correction_part2=(4.0e0*(1.0e0-alpha_d_org))/(3.0e0*alpha_d_org)
        Correction_part3=1.0e0+(Correction_part1+Correction_part2)*Kn
        Correction=1.0/Correction_part3
        # Now calculate a kelvin factor for every semi-volatile compound in this
        # size bin (Equation 2.28)
        kelvin_factor=np.exp((4.0E0*mw_array*1.0e-3*sigma)/(R_gas*Temp_K*size_array[size_step]*2.0e0*density))
        # Calculate a hypothetical contribution from non-ideality
        activity_coeffs = 1.0+gamma_factor*mole_fractions
        # Calculate the equilibrium concentration [molecules/cc]
        Pressure_eq=kelvin_factor*mole_fractions*P_sat*activity_coeffs
        # Calculate the equilibrium concentration equivalent
        Cstar_i_m_t=Pressure_eq*(NA/(R_gas_other*1.0e6*Temp_K))
        # Implement the droplet growth equation (Equation 2.44)
        k_i_m_t_part1 = DStar_org*Correction
        k_i_m_t=4.0e0*np.pi*size_array[size_step]*1.0e2*N_per_bin[size_step]*k_i_m_t_part1
        dm_dt=k_i_m_t*(Cg_i_m_t-Cstar_i_m_t)
        # Now update the contribution to the ODEs being solved
        # Add contributory loss from the gas phase to particle phase
        dy_dt_gas_matrix[0:num_species,size_step]=dm_dt
        # Add a contributory gain to the particle phase from the gas phase
        dy_dt_array[num_species+size_step*num_species:num_species+num_species*(size_step+1)]=dm_dt[0:num_species]

    # Subtract the net condensational mass flux from the gas phase concentrations
    dy_dt_array[0:num_species]=dy_dt_array[0:num_species]-np.sum(dy_dt_gas_matrix, axis=1)

    return dy_dt_array

# Define the time over which the simulation will take place
t = np.linspace(0, 5000, num=1000)
# Call the ODE solver with reference to our function, dy_dt, setting the
# absolute and relative tolerance.
solution = odeint(dy_dt, array, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)

# Define an array that will store the size of our
# particles at each model output
size_array_matrix = np.zeros((1000,num_bins), dtype=float)

# Let's now try to plot the size variation as a function of time
for step in range(solution.shape[0]):

    # For each point in time, calculate the total SOA
    condensable_mass = 0

    for size_step in range(num_bins):

        temp_array=solution[step,num_species+size_step*num_species:num_species+num_species*(size_step+1)]
        total_moles=np.sum(temp_array)+core_abundance[size_step]
        mole_fractions=temp_array/total_moles
        density_array = np.zeros((num_species+1), dtype=float)
        density_array[0:num_species]=density_org[0:num_species]
        density_array[num_species]=density_core[size_step]
        mass_array = np.zeros((num_species+1), dtype=float)
        mass_array[0:num_species] = (temp_array/NA)*mw_array
        mass_array[num_species] = (core_abundance[size_step]/NA)*core_mw[size_step]
        total_mass=np.sum(mass_array)
        mass_fractions_array=mass_array/total_mass
        density=1.0/(np.sum(mass_fractions_array/density_array))
        size_array_matrix[step,size_step]=((3.0*((total_mass*1.0e3)/(N_per_bin[size_step]*1.0e6)))/(4.0*np.pi*density))**(1.0/3.0)

        condensable_mass+=np.sum(mass_array)*(1.0e6)*(1.0e6)

    condensable_mass=condensable_mass-dry_mass

print("Final condensed mass = ", condensable_mass)
