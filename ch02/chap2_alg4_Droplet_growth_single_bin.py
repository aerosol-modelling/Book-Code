# Import the relevant libraries, including an ODE solver from Scipy
# We only want to use the 'odeint' [Solve Initial Value Problems] from the
# scipy.integrate package
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt # Import Matplotlib so we can plot results

# Define physical constants
R_gas=8.31451 #Ideal gas constant [kg m2 s-2 K-1 mol-1]
R_gas_other=8.20573e-5 #Ideal gas constant [m3 atm K-1 mol-1]
sigma=72.0e-3 # Assume surface tension of water (mN/m)
NA=6.0221409e+23 #Avogadros number

# Specify the temperature
Temp_K=298.15

# Number of condensing species from the gas phase
num_species = 10
# Number of size bins
num_bins = 1

# The molecular weight of each condensing specie [g/mol]
# Assuming a constant value for all components
mw_array=np.zeros((num_species), dtype=float)
mw_array[:]=200.0

# Define array of log10 C* values
log_c_star = np.linspace(-6, 3, num_species)
Cstar = np.power(10.0,log_c_star)
# Convert C* to a pure component saturation vapour pressure [atm]
P_sat = (Cstar*R_gas_other*Temp_K)/(1.0e6*mw_array)

# Initialise abundance in each volatility bin [micrograms/m3]
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

# Define a monodisperse size distribution
# Assume each particle starts with an involatile core
# of absorptive organic with a mass of 200 g/mol and
# density 1400 km/m3. We store this information in an array
# as 'core'. This will ensure, in this example, that we do
# not get 100% evaporative loss

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

# Use the size to calculate concentration of 'core' [molecules / cc]
# This aligns with the units used for volatile component
core_abundance[0] = (N_per_bin[0])*((4.0/3.0)*np.pi*
    np.power(size_array[0]*0.5e-6,3.0)*density_core[0]*1.0e3)
core_abundance[0] = (core_abundance[0] / core_mw[0])*NA

# What is our existing dry mass in micrograms per cubic metre?
dry_mass = np.sum((core_abundance/NA)*core_mw)*(1.0e12)
print("Initial dry mass = ", dry_mass)


# Define the gas and condensed concentration arrays
# Here we initialise the array used throughout the simulation and
# populate with initial concentrations in gas and condensed phase
array = np.zeros((num_species+num_species*num_bins), dtype=float)
# The first num_species cells hold the gas phase abundance
array[0:num_species] = gas_abundance
# The following cells hold the concentration of each specie in each size bin
array[num_species:num_species+num_species*num_bins] = 1.0e-20

# Define the Droplet Growth equation
# Here we use variable names as close to the equation syntax as possible

def dy_dt(array,t):

    # Retrieve the gas phase concentrations
    Cg_i_m_t = array[0:num_species]
    # Retrieve concentrations in our size bin as a slice
    temp_array=array[num_species:num_species+num_species*num_bins]
    # Sum total molecules in the condensed phase, plus core
    total_molecules=np.sum(temp_array)+core_abundance[0]
    # Calculate mole fractions in the, assumed, liquid phase for
    # calculating the equilibrium pressure above the droplet
    mole_fractions=temp_array/total_molecules
    # Calculate the density of the assumed liquid phase
    density_array = np.zeros((num_species+num_bins), dtype=float)
    density_array[0:num_species]=density_org[0:num_species]
    density_array[num_species]=density_core[0]
    # Create array that holds mass concentrations [g/cc] for
    # calculation of solution density [kg/m3]
    mass_array = np.zeros((num_species+1), dtype=float)
    mass_array[0:num_species] = (temp_array/NA)*mw_array
    mass_array[num_species] = (core_abundance[0]/NA)*core_mw[0]
    total_mass=np.sum(mass_array)
    mass_fractions_array=mass_array/total_mass
    density=1.0/(np.sum(mass_fractions_array/density_array))

    # Now calculate the size [cm] of our particles according
    # to the condensed phase abundance of material
    # We need to remember that mass is in [g/cc] whilst density
    # is in [kg/m3].
    # We convert mass and number concentrations to kg/m3 and /m3
    size=((3.0*((total_mass*1.0e3)/(N_per_bin[0]*1.0e6)))/
        (4.0*np.pi*density))**(1.0/3.0)

    # Calculate the Knudsen number for all condensing
    # molecules based on this new size
    # This relies on mean free path for each species [cm]
    # and particle radius [cm]
    Kn=gamma_gas/size
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
    kelvin_factor=np.exp((4.0E0*mw_array*1.0e-3*sigma)/
        (R_gas*Temp_K*size_array[0]*2.0e0*density))
    # Calculate the equilibrium pressure at the surface
    Pressure_eq=kelvin_factor*mole_fractions*P_sat
    # Calculate the equilibrium concentration [molecules/cc]
    # equivalent of this pressure
    Cstar_i_m_t=Pressure_eq*(NA/(R_gas_other*1.0e6*Temp_K))

    # Implement the droplet growth equation (Equation 2.44)
    # The equation relies on the following parameters
    # radius [m]
    # DStar_org [m2/s] - pass in [cm2/s] so needs to be converted by *1.0E-4
    # Pressure_gas [Pascals]
    # Pressure_eq [Pascals]
    # R_gas [m3 Pascals /(K mol)]
    # molw [g/mol]
    # T [K]
    # The units of the equation should therefore be g/s
    k_i_m_t_part1 = DStar_org*Correction
    k_i_m_t=4.0e0*np.pi*size*1.0e2*N_per_bin[0]*k_i_m_t_part1
    dm_dt=k_i_m_t*(Cg_i_m_t-Cstar_i_m_t)

    # Store values of dm_dt in a matrix, ready for dealing with polydisperse
    # populations where a contribution/loss per size bin
    dy_dt_gas_matrix = np.zeros((num_species,num_bins), dtype=float)
    dy_dt_gas_matrix[0:num_species,0]=dm_dt

    #Initiase an array that passes back a set of derivatives for Each
    # compound in both the gaseous and condensed phases
    dy_dt_array = np.zeros((num_species+num_species*num_bins), dtype=float)
    dy_dt_array[num_species:num_species+num_species*num_bins]=dm_dt[0:num_species]
    dy_dt_array[0:num_species]=dy_dt_array[0:num_species]-np.sum(dy_dt_gas_matrix, axis=1)

    return dy_dt_array


# Define the time over which the simulation will take place
t = np.linspace(0, 10000, num=1000)

# Call ODE solver with reference to dy_dt, setting the
# absolute and relative tolerance.
solution = odeint(dy_dt, array, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)

# Define an array that will store the size of our
# particles at each model output
size_array = np.zeros((1000), dtype=float)

# Now plot the size variation as a function of time,
# converting the concentrations we get as output into particle size

for step in range(solution.shape[0]):

    # The ODE solver outputs a 2D array of values, where each row represents
    # the concentrations at each point in time
    temp_array=solution[step,num_species:num_species+num_species*num_bins]
    # The following code follows that within the dy_dt function
    total_moles=np.sum(temp_array)+core_abundance[0]
    mole_fractions=temp_array/total_moles
    density_array = np.zeros((num_species+num_bins), dtype=float)
    density_array[0:num_species]=density_org[0:num_species]
    density_array[num_species]=density_core[0]
    mass_array = np.zeros((num_species+1), dtype=float)
    mass_array[0:num_species] = (temp_array/NA)*mw_array
    mass_array[num_species] = (core_abundance[0]/NA)*core_mw[0]
    total_mass=np.sum(mass_array)
    mass_fractions_array=mass_array/total_mass
    density=1.0/(np.sum(mass_fractions_array/density_array))
    size_array[step]=((3.0*((total_mass*1.0e3)/(N_per_bin[0]*1.0e6)))
        /(4.0*np.pi*density))**(1.0/3.0)

print("Final mass = ", np.sum(mass_array)*(1.0e6)*(1.0e6) - dry_mass)

# Produce a figure of change in gas phase concentrations and particle size
fig = plt.figure(figsize=(10,10))
plt.subplot(1, 2, 1)
plt.plot(t, np.log10(solution[:,0:num_species]))
plt.grid()
plt.xlabel('time [seconds]')
plt.ylabel('Gas phase abundance [Log10/ molecules.cc]')
plt.title('Change in gas phase abundance')
plt.subplot(1, 2, 2)
plt.plot(t, size_array*2.0e6, 'r-', label='Mono-disperse population growth')
plt.grid()
plt.legend(loc='best')
plt.xlabel('time [seconds]')
plt.ylabel('Particle diameter [microns]')
plt.title('Particle growth under condensation')
plt.show()
