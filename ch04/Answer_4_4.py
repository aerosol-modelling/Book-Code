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

# Set the ambient temperature
Temp_K=298.15

# Number of condensing species from the gas phase
num_species = 10

# The molecular weight of each condensing specie [g/mol]
# Assuming a constant value for all components
mw_array=np.zeros((num_species+1), dtype=float)
mw_array[0:10]=200.0
mw_array[10]=47.998

# Define array of log10 C* values
log_c_star = np.linspace(-6, 3, num_species)
Cstar = np.power(10.0,log_c_star)
# Convert C* to a pure component saturation vapour pressure [atm]
P_sat = (Cstar*R_gas_other*Temp_K)/(1.0e6*mw_array[0:10])

# Initialise an abundance of material in each volatility bin [micrograms/m3]
concentration = np.zeros((num_species+1), dtype = float)
concentration[0] = 0.1
concentration[1] = 0.1
concentration[2] = 0.15
concentration[3] = 0.22
concentration[4] = 0.36
concentration[5] = 0.47
concentration[6] = 0.58
concentration[7] = 0.69
concentration[8] = 0.84
concentration[9] = 1.0

# Initialise abundance of a gas phase oxidant in ppb
ozone = 10.0
Cfactor= 2.46310e+10 #ppb-to-molecules/cc
# convert this to  molecules / cc
ozone_concentration = ozone*Cfactor
# Now convert to micrograms/m3
ozone_concentration = (ozone_concentration/(NA*1.0e-6))*(mw_array[10]*1.0e6)
concentration[10] = ozone_concentration # Ozone

# Define a rate coefficient between every volatility bin and ozone
rate_coefficient = 1.0e-15

# Unit conversion of gas abudance to molecules / cc
gas_concentration = ((concentration*1.0e-6)/(mw_array))*1.0e-6*NA

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

# Define a polydisperse size distribution
# Lowest size diameter
d1 = 0.01
# Diameter of largest bin [microns]
d_Nb = 1.0
# Number of size bins
num_bins = 8
# Volume ratio between bins (Equation 1.17)
V_rat = np.power((d_Nb/d1),3.0/(num_bins-1.0)) # Volume ratio between bins
# Diameter ratio between bins  (Equation 1.18)
d_rat = V_rat**(1.0/3.0) # Diameter ratio between bins

# Create an array of diameters
d_i=np.zeros((num_bins), dtype=float) # Diameter array
d_i[0]=d1
for step in range(num_bins):
    if step > 0:
       d_i[step]=d_i[step-1]*d_rat
log_di = np.log(d_i) # Log of Diameter array

# Define parameters of log-normal distribution
# Geometric standard deviation
sigmag1 = np.log(1.7)
 # Mean particle diameter [150nm]
mean1 = np.log(0.15)
# Calculate the probability density distribution
distribution_1 = (np.exp(-(log_di - mean1)**2 / (2 * sigmag1**2)) / (sigmag1 * np.sqrt(2 * np.pi))) # Probability density distribution
# Seperate out the probability density function
d_width = d_i*np.power(2,1.0/3.0)*((np.power(V_rat,1.0/3.0)-1.0)/(np.power(1+V_rat,1.0/3.0))) # Diameter width array of size bins
# Total number of particles [per cm-3]
N_total = 100.0
# Use pre-calculated distribution_1 to implement equation 1.16
N_dist = N_total*(distribution_1*(d_width/d_i))

# Initialise a core concentration using the above size distribution
core = np.zeros((num_bins), dtype=float)
core_concentration = np.zeros((num_bins), dtype=float)
density_core = np.zeros((num_bins), dtype=float)
core_mw = np.zeros((num_bins), dtype=float)
density_core[:] = 1400.0
core_mw[:] = 200.0
N_per_bin = N_dist

# Define size array (radius )
size_array = d_i*0.5
# Use the size to now calculate a concentration of a 'core' in molecules / cc
core_concentration = (N_per_bin)*((4.0/3.0)*np.pi*np.power(size_array*1.0e-6,3.0)*density_core*1.0e3)
core_concentration = (core_concentration / core_mw)*NA
dry_mass = np.sum((core_concentration/NA)*core_mw)*(1.0e6)*(1.0e6)
print("Initial dry mass = ", dry_mass)

# New define an array that holds the molecular abundance of
# each gas and the concentration of each gas in each size bin
array = np.zeros((num_species+num_species*num_bins+1), dtype=float)
array[0:num_species+1] = gas_concentration
array[num_species+1:num_species+num_species*num_bins+1] = 1.0e-10 # assuming we start with nothing

# Now let us convert the given emission rate of ozone into molecules / cc.second
emission_rate_ozone = ((10.0*1.0e-12)/(3600.0))*(NA/mw_array[10])

# We can also calculate the emission rate of the highest volatility bin
emission_rate_bin = ((1.5*1.0e-12)/(3600.0))*(NA/mw_array[9])

# RHS function
def dy_dt(array,t):

    # Retrieve the gas phase concentrations
    Cg_i_m_t = array[0:num_species]

    # We are working with 8 size bins, each of which has an involatile core
    size_array = np.zeros((num_bins), dtype=float)
    dy_dt_array = np.zeros((num_species+num_species*num_bins+1), dtype=float)
    dy_dt_gas_matrix = np.zeros((num_species,num_bins), dtype=float)

    # Calculate the rate of each reaction explicitly
    rate_1 = array[10]*array[1]*rate_coefficient
    rate_2 = array[10]*array[2]*rate_coefficient
    rate_3 = array[10]*array[3]*rate_coefficient
    rate_4 = array[10]*array[4]*rate_coefficient
    rate_5 = array[10]*array[5]*rate_coefficient
    rate_6 = array[10]*array[6]*rate_coefficient
    rate_7 = array[10]*array[7]*rate_coefficient
    rate_8 = array[10]*array[8]*rate_coefficient
    rate_9 = array[10]*array[9]*rate_coefficient

    # Now write out the RHS for each volatity bin according to our simple
    # mechanism
    dy_dt_array[0] = rate_1
    dy_dt_array[1] = rate_2 - rate_1
    dy_dt_array[2] = rate_3 - rate_2
    dy_dt_array[3] = rate_4 - rate_3
    dy_dt_array[4] = rate_5 - rate_4
    dy_dt_array[5] = rate_6 - rate_5
    dy_dt_array[6] = rate_7 - rate_6
    dy_dt_array[7] = rate_8 - rate_7
    dy_dt_array[8] = rate_9 - rate_8
    dy_dt_array[9] = -1.0*rate_9 + emission_rate_bin
    # Define the dydt for Ozone
    dy_dt_array[10] = emission_rate_ozone - (rate_1+rate_2+rate_3
        +rate_4+rate_5+rate_6+rate_7+rate_8+rate_9)

    # Now cycle through each size bin
    for size_step in range(num_bins):

        # Select a slice of y that represents this size bin
        temp_array=array[(num_species+size_step*num_species)+1:(num_species+num_species*(size_step+1))+1]
        # Sum the total molecules in the condensed phase, adding core material
        total_moles=np.sum(temp_array)+core_concentration[size_step]
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
        mass_array[0:num_species] = (temp_array/NA)*mw_array[0:10]
        mass_array[num_species] = (core_concentration[size_step]/NA)*core_mw[0]
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
        Kn=gamma_gas[0:num_species]/size_array[size_step]
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
        kelvin_factor=np.exp((4.0E0*mw_array[0:num_species]*1.0e-3*sigma)/(R_gas*Temp_K*size_array[size_step]*2.0e0*density))
        # Calculate the equilibrium concentration [molecules/cc]
        Pressure_eq=kelvin_factor*mole_fractions*P_sat
        # Calculate the equilibrium concentration equivalent
        Cstar_i_m_t=Pressure_eq*(NA/(R_gas_other*1.0e6*Temp_K))
        # Implement the droplet growth equation (Equation 2.44)
        k_i_m_t_part1 = DStar_org[0:num_species]*Correction
        k_i_m_t=4.0e0*np.pi*size_array[size_step]*1.0e2*N_per_bin[size_step]*k_i_m_t_part1
        dm_dt=k_i_m_t*(Cg_i_m_t-Cstar_i_m_t)
        # Now update the contribution to the ODEs being solved
        # Add contributory loss from the gas phase to particle phase
        dy_dt_gas_matrix[0:num_species,size_step]=dm_dt
        # Add a contributory gain to the particle phase from the gas phase
        dy_dt_array[(num_species+size_step*num_species)+1:(num_species+num_species*(size_step+1))+1]=dm_dt[0:num_species]

    dy_dt_array[0:num_species]=dy_dt_array[0:num_species]-np.sum(dy_dt_gas_matrix, axis=1)
    return dy_dt_array

# Define the time over which the simulation will take place
t = np.linspace(0, 7200, num=1000)
# Call the ODE solver with reference to our function, dy_dt, setting the
# absolute and relative tolerance.
solution_ozone = odeint(dy_dt, array, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)

# Now extract a size distribution for every point in the simulation
size_array_matrix_ozone = np.zeros((1000,num_bins), dtype=float)

for step in range(solution_ozone.shape[0]):
    # For each point in time, calculate the total SOA
    condensable_mass_ozone = 0
    for size_step in range(num_bins):
        temp_array=solution_ozone[step,(num_species+size_step*num_species)+1:(num_species+num_species*(size_step+1))+1]
        total_moles=np.sum(temp_array)+core_concentration[size_step]
        mole_fractions=temp_array/total_moles
        density_array = np.zeros((num_species+1), dtype=float)
        density_array[0:num_species]=density_org[0:num_species]
        density_array[num_species]=density_core[size_step]
        mass_array = np.zeros((num_species+1), dtype=float)
        mass_array[0:num_species] = (temp_array/NA)*mw_array[0:num_species]
        mass_array[num_species] = (core_concentration[size_step]/NA)*core_mw[size_step]
        total_mass=np.sum(mass_array)
        mass_fractions_array=mass_array/total_mass
        density=1.0/(np.sum(mass_fractions_array/density_array))
        size_array_matrix_ozone[step,size_step]=((3.0*((total_mass*1.0e3)/(N_per_bin[size_step]*1.0e6)))/(4.0*np.pi*density))**(1.0/3.0)
        condensable_mass_ozone+=np.sum(mass_array)*(1.0e6)*(1.0e6)
    condensable_mass_ozone=condensable_mass_ozone-dry_mass

print("Final condensed mass [with ozone reaction] = ", condensable_mass_ozone)

# Plot the decrease in gas phase abundance with the growth of the particle [both absolute and relative]

fig = plt.figure(figsize=(10,10))
plt.plot(t, size_array_matrix_ozone*2.0e6, 'r')
plt.grid()
plt.xlabel('time [seconds] ')
plt.ylabel('Particle diameter [microns]')
plt.title('Particle growth')
plt.show()
