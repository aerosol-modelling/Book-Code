# Import the relevant libraries, including an ODE solver from Scipy
# We only want to use the 'odeint' [Solve Initial Value Problems] from the
# scipy.integrate package
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt # Import Matplotlib so we can plot results

# --- Physical constants --------------------------------------------
Lv_water_vapour=2.5e3 # Latent heat of vapourisation of water [J/g]
Rv=461.0 #Individual gas constant of water vapour [J/Kg.K]
Ra=287.0 #Gas constant for dry air [J/Kg.K]
R_gas=8.3144598 #Ideal gas constant [kg m2 s-2 K-1 mol-1]
R_gas_other=8.2057e-5 #Ideal gas constant [m3 atm K-1 mol-1]
GRAV=9.8; #Gravitational acceleration [(m/2)2]
cp=1005; #Specific heat capacity of air [J/Kg.K]
sigma=72.0e-3 # Assume surface tension of water (mN/m)
Avo=6.0221409e+23 #Avogadros number
kb=1.380648E-23 #Boltzmanns constant
# ---------------------------------------------------------------------

Temp_K=298.15

# ------- Defining the gaseous and condensed phase components ---------
num_species = 10 # Number of condensing species from the gas phase

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
print( 'Initial vapor mass = ',np.sum(abundance))

# Initialise abundance of a gas phase oxidant in ppb
ozone = 10.0
Cfactor= 2.55e+10 #ppb-to-molecules/cc
# convert this to  molecules / cc
ozone_abundance = ozone*Cfactor
# Define a rate coefficient between every volatility bin and ozone
rate_coefficient = 1.0e-15

# Unit conversion of gas abudance to molecules / cc
gas_abundance = ((abundance*1.0e-6)/(mw_array))*1.0e-6*Avo

# Assumed accomodation coefficient
alpha_d_org=np.zeros((num_species), dtype=float)
alpha_d_org[:]=1.0
# Density of condensing species [kg/m3]
density_org=np.zeros((num_species), dtype=float)
density_org[:]=1400.0
rho_l = density_org * 1000.0 / mw_array * Avo  # [molec m-3]

# Molecular diffusion coefficient of each molecule in air. [cm2/s]
DStar_org = 1.9*np.power(mw_array,-2.0/3.0)
# Mean thermal velocity of each molecule [m/s] -
mean_them_vel = np.power((8.0*R_gas*Temp_K)/(np.pi*(mw_array*1.0e-3)),0.5)
# Mean free path for each molecule [m]
gamma_gas = ((3.0*DStar_org)/(mean_them_vel*1.0e2))*1.0e-2
# ---------------------------------------------------------------------

# -------- Create a modal distribution --------
# Parameters of log-normal distribution
sigmag1 = np.log(1.7) # Natural Log of Geometric standard deviation
Dp      = 150 * 1.0E-9 # Geometric Mean particle diameter [150nm]
N       = 100.0 # Total number of particles [per cm-3]

# Initialize all 3 moments of initial distribution
M0 = N * 1.0E6  # [particles per m-3]
M2 = M0 * Dp**2 * np.exp( 2 * sigmag1**2 )   # [m2 m-3]
M3 = M0 * Dp**3 * np.exp( 4.5 * sigmag1**2 ) # [m3 m-3]
# -------------------------------------------------------

# -------- Initialize Size Distribution of Particle Core abundance using the above size distribution -------
density_core = 1400.0  # kg m-3
core_mw      = 200.0   # g mol-1

dry_mass = np.pi/6.0 * density_core * M3 * 1.0E9 # ug m-3
print("Initial dry aerosol mass = ", dry_mass)
core_abundance = dry_mass * 1.0E-6 / core_mw * Avo * 1.0E-6
# [molec cm-3]
print("Initial core abundance = ", core_abundance )

# New define an array that holds the molecular abundance of each gas and its concentration
# in each aerosol mode [molecules / cc]
# Add 3 additional terms for each mode at the end representing M0, M2, and M3
array = np.zeros(( num_species + num_species + 3 ), dtype=float)
array[ 0:num_species ] = gas_abundance
array[ num_species:2*num_species ] = 1.0e-10 # assuming we start with no condensed 
                                             # material in the aerosol phase
nchem = num_species + num_species
array[nchem     :nchem+1 ] = M0
array[nchem + 1 :nchem+2 ] = M2
array[nchem + 2 :nchem+3 ] = M3

########################################################
####################################

# RHS function [with no ozone effects]
def dy_dt(array,t):

    # Load gas species concentrations in local array
    Cg_i_m_t = array[0:num_species]

    # We are working with 1 aerosol mode, which has an involatile core
    dy_dt_array = np.zeros(( nchem + 3 ), dtype='float64')
    dM0_dt = 0.
    dM2_dt = 0.
    dM3_dt = 0.

    # 1) Calculate total moles and mole fraction of each species in this mode,
    #    and assign modal parameters
    temp_array     = array[ num_species:nchem ]        # [molec/cc]
    total_moles    = np.sum(temp_array) + core_abundance # [molec/cc]
    mole_fractions = temp_array / total_moles            # [frac]

    M0 = array[ nchem   ]  # 0th moment = Number
    M2 = array[ nchem+1 ]  # 2nd moment prop. to surface area
    M3 = array[ nchem+2 ]  # 3rd moment prop. to volume

    # 2) Calculate density 
    density_array = np.zeros((num_species+1), dtype='float64')
    density_array[0:num_species] = density_org[0:num_species] # [kg m-3]
    density_array[num_species]   = density_core               # [kg m-3]

    mass_array = np.zeros((num_species+1), dtype='float64')
    mass_array[0:num_species] = (temp_array/Avo) * mw_array   # [g cm-3]
    mass_array[num_species] = (core_abundance/Avo)*core_mw    # [g cm-3]

    total_mass = np.sum(mass_array)                           # [g cm-3]
    mass_fractions_array = mass_array / total_mass            # [frac]
    density = 1.0 / (np.sum( mass_fractions_array / density_array ))  # [kg m-3]

    #3) Calculate Diameter and sigma from the three moments above. 
    
    sigmag1 = ( 1.0/3.0 * np.log( M0 ) - np.log( M2 ) + 2.0/3.0 * np.log( M3 ) ) **0.5
    Dp = ( M3 / ( M0 * np.exp( 4.5 * sigmag1**2 ) ) ) **(1.0/3.0)  # [m]
    M1 = M0 * Dp * np.exp( 0.5 * sigmag1**2 )   # [m m-3]

    #4) Free Molecular Regime Calculation for the 2nd and 3rd moments
    Ifm2 = np.pi * alpha_d_org * mean_them_vel / 4. * M1  # [ m2 m-3 s-1]
    Ifm3 = np.pi * alpha_d_org * mean_them_vel / 4. * M2  # [ m3 m-3 s-1]
    
    #5) Continuum Regime Calculation for the 2nd and 3rd Moments
    Ict2 = 2.0 * np.pi * DStar_org[:]*1.0e-4 * M0   # [m2 s-1 m-3]
    Ict3 = 2.0 * np.pi * DStar_org[:]*1.0e-4 * M1   # [m3 s-1 m-3]

    #6) Combine Ifm and Ict for moment- and species-dependent Beta correction terms
    Beta2 = ( Ifm2 * Ict2 ) / ( Ifm2 + Ict2 )  # [m2 m-3 s-1]
    Beta3 = ( Ifm3 * Ict3 ) / ( Ifm3 + Ict3 )  # [m3 m-3 s-1]

    #7) Modal growth equation
    Cstar_i_m_t = ((Cstar*1.0e-6)/(mw_array))*1.0e-6*Avo
    
    dm2_dt = 4.0 / np.pi * np.sum( (Cg_i_m_t - Cstar_i_m_t*mole_fractions) * 1.0e6 / rho_l * Beta2 )
    dm3_dt = 6.0 / np.pi * np.sum( (Cg_i_m_t - Cstar_i_m_t*mole_fractions) * 1.0e6 / rho_l * Beta3 )
    dy_dt_array[ num_species:nchem ]  = ( Cg_i_m_t[:] - Cstar_i_m_t*mole_fractions) * Beta3
    dy_dt_array[ 0:num_species ]  = -1.0 * dy_dt_array[num_species:2*num_species]
    dy_dt_array[ nchem     ] = 0.  # Particle number does not change. This is included in the ODE 
                                   # set for demonstration since it would change for other processes
                                   # like coagulation. It should be removed for operational use if 
                                   # only condensation is being considered.
    dy_dt_array[ nchem + 1 ] = dm2_dt
    dy_dt_array[ nchem + 2 ] = dm3_dt
   
       return dy_dt_array




def dy_dt_ozone(array,t):

    # Load gas species concentrations in local array
    Cg_i_m_t = array[0:num_species]

    # We are working with 1 aerosol mode, which has an involatile core
    dy_dt_array = np.zeros(( nchem + 3 ), dtype='float64')
    dM0_dt = 0.
    dM2_dt = 0.
    dM3_dt = 0.

    # Calculate the rate of each reaction explicitly
    # rate_1 = ozone_abundance*array[0]*rate_coefficient
    rate_1 = ozone_abundance*array[1]*rate_coefficient
    rate_2 = ozone_abundance*array[2]*rate_coefficient
    rate_3 = ozone_abundance*array[3]*rate_coefficient
    rate_4 = ozone_abundance*array[4]*rate_coefficient
    rate_5 = ozone_abundance*array[5]*rate_coefficient
    rate_6 = ozone_abundance*array[6]*rate_coefficient
    rate_7 = ozone_abundance*array[7]*rate_coefficient
    rate_8 = ozone_abundance*array[8]*rate_coefficient
    rate_9 = ozone_abundance*array[9]*rate_coefficient

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
    dy_dt_array[9] = -1.0*rate_9

    # 1) Calculate total moles and mole fraction of each species in this mode,
    #    and assign modal parameters
    temp_array     = array[ num_species:nchem ]          # [molec/cc]
    total_moles    = np.sum(temp_array) + core_abundance # [molec/cc]
    mole_fractions = temp_array / total_moles            # [frac]

    M0 = array[ nchem   ]  # 0th moment = Number
    M2 = array[ nchem+1 ]  # 2nd moment prop. to surface area
    M3 = array[ nchem+2 ]  # 3rd moment prop. to volume

    # 2) Calculate density 
    density_array = np.zeros((num_species+1), dtype='float64')
    density_array[0:num_species] = density_org[0:num_species] # [kg m-3]
    density_array[num_species]   = density_core               # [kg m-3]

    mass_array = np.zeros((num_species+1), dtype='float64')
    mass_array[0:num_species] = (temp_array/Avo) * mw_array   # [g cm-3]
    mass_array[num_species] = (core_abundance/Avo)*core_mw    # [g cm-3]

    total_mass = np.sum(mass_array)                           # [g cm-3]
    mass_fractions_array = mass_array / total_mass            # [frac]
    density = 1.0 / (np.sum( mass_fractions_array / density_array ))  # [kg m-3]

    #3) Calculate Diameter and sigma from the three moments above. 
    sigmag1 = ( 1.0/3.0 * np.log( M0 ) - np.log( M2 )
	  + 2.0/3.0 * np.log( M3 ) ) **0.5
    Dp = ( M3 / ( M0 * np.exp( 4.5 * sigmag1**2 ) ) )
	   **(1.0/3.0)  # [m]
    M1 = M0 * Dp * np.exp( 0.5 * sigmag1**2 )   # [m m-3]

    #4) Free Molecular Regime Calculation for the 2nd and 3rd moments
    Ifm2 = np.pi * alpha_d_org * mean_them_vel
	  / 4. * M1  # [ m2 m-3 s-1]
    Ifm3 = np.pi * alpha_d_org * mean_them_vel
	  / 4. * M2  # [ m3 m-3 s-1]
    
    #5) Continuum Regime Calculation for the 2nd and 3rd Moments
    Ict2 = 2.0 * np.pi * DStar_org[:]*1.0e-4 * M0   # [m2 s-1 m-3]
    Ict3 = 2.0 * np.pi * DStar_org[:]*1.0e-4 * M1   # [m3 s-1 m-3]

    #6) Combine Ifm and Ict for moment- and species-dependent Beta correction terms
    Beta2 = ( Ifm2 * Ict2 ) / ( Ifm2 + Ict2 )  # [m2 m-3 s-1]
    Beta3 = ( Ifm3 * Ict3 ) / ( Ifm3 + Ict3 )  # [m3 m-3 s-1]

    #7) Modal growth equation
    Cstar_i_m_t = ((Cstar*1.0e-6)/(mw_array))*1.0e-6*Avo
    
    dm2_dt = 4.0 / np.pi * np.sum( (Cg_i_m_t - Cstar_i_m_t*mole_fractions) * 1.0e6 / rho_l * Beta2 )
    dm3_dt = 6.0 / np.pi * np.sum( (Cg_i_m_t - Cstar_i_m_t*mole_fractions) * 1.0e6 / rho_l * Beta3 )
    dy_dt_array[ num_species:nchem ]  = ( Cg_i_m_t[:] - Cstar_i_m_t*mole_fractions) * Beta3
    dy_dt_array[ 0:num_species ]  = dy_dt_array[0:num_species] - dy_dt_array[num_species:2*num_species]
    dy_dt_array[ nchem     ] = 0.  # Particle number does not change. This is included in the ODE 
                                   # set for demonstration since it would change for other processes
                                   # like coagulation. It should be removed for operational use if 
                                   # only condensation is being considered.
    dy_dt_array[ nchem + 1 ] = dm2_dt
    dy_dt_array[ nchem + 2 ] = dm3_dt
   
    return dy_dt_array
########################################################

t = np.linspace(0, 10000, num=1000)

# Now test timing of the derivative function
dy_dt(array,0.0)

# Call the numba functions for first compilation
test=dy_dt(array,0.0)
solution_no_ozone = odeint(dy_dt, array, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)

# Now extract a size distribution for every point in the simulation
size_array_no_ozone = np.zeros((1000), dtype=float)
sigma_array_no_ozone = np.zeros((1000), dtype=float)

for step in range(solution_no_ozone.shape[0]):
    # For each point in time, calculate the total SOA
    condensable_mass = 0
    temp_array       = solution_no_ozone[step,num_species:nchem]
    total_moles      = np.sum(temp_array) + core_abundance
    mole_fractions   = temp_array / total_moles
    density_array    = np.zeros((num_species+1), dtype=float)
    density_array[0:num_species] = density_org[0:num_species]
    density_array[num_species]   = density_core
    mass_array                   = np.zeros((num_species+1), dtype=float)
    mass_array[0:num_species]    = (temp_array/Avo) * mw_array
    mass_array[num_species]      = (core_abundance/Avo) * core_mw
    total_mass                   = np.sum(mass_array)
    mass_fractions_array         = mass_array/total_mass

    # New Particle Mass [ug m-3]
    condensable_mass = np.sum(mass_array)*(1.0e6)*(1.0e6) - dry_mass
    
    # Size Distribution Properties
    density = 1.0 / (np.sum( mass_fractions_array / density_array ))
    M2 = solution_no_ozone[ step,nchem+1 ]
    M3 = solution_no_ozone[ step,nchem+2 ]
    sigma_array_no_ozone[step] = ( 1.0/3.0 * np.log( M0 ) - np.log( M2 ) + 2.0/3.0 * np.log( M3 ) ) **0.5
    size_array_no_ozone[step] = ( M3 / ( M0 * np.exp( 4.5 * sigma_array_no_ozone[step]**2 ) ) ) **(1.0/3.0)
    
    # Vapor Mass Remaining [ug m-3]
    vapor_mass = np.sum(solution_no_ozone[step,0:num_species] / Avo*mw_array*1.0e6*1.0e6)
 
print("Final condensed mass [no ozone] = ", condensable_mass)
 
# Plot the decrease in gas phase abundance with the growth of the particle [both absolute and relative]
fig = plt.figure(figsize=(10,10))

# Plot the Increase in Particle Size for both Cases
ax1 = plt.subplot(1,2,1)
ax1.plot(t, size_array_no_ozone*1.0e6, 'r-')
ax1.grid()
ax1.set_xlabel('time [seconds] ')
ax1.set_ylabel('Particle diameter [microns]')
ax1.set_title('Particle growth')

# Plot the Change in Mode Width for both cases
ax2 = plt.subplot(1,2,2)
ax2.plot(t, np.exp(sigma_array_no_ozone), 'r-',label='No Ozone')
ax2.grid()
ax2.set_xlabel('time [seconds] ')
ax2.set_ylabel(r'$\sigma_g$')
ax2.set_title('Mode Width')
ax2.legend()

plt.show(block=False) 






# Now test timing of the derivative function with ozone reaction
dy_dt_ozone(array,0.0)

# Call the numba functions for first compilation
test=dy_dt_ozone(array,0.0)
solution_ozone = odeint(dy_dt_ozone, array, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)

# Now extract a size distribution for every point in the simulation
size_array_ozone = np.zeros((1000), dtype=float)
sigma_array_ozone = np.zeros((1000), dtype=float)


for step in range(solution_ozone.shape[0]):
    # For each point in time, calculate the total SOA
    temp_array       = solution_ozone[step,num_species:nchem]
    total_moles      = np.sum(temp_array) + core_abundance
    mole_fractions   = temp_array / total_moles
    density_array    = np.zeros((num_species+1), dtype=float)
    density_array[0:num_species] = density_org[0:num_species]
    density_array[num_species]   = density_core
    mass_array                   = np.zeros((num_species+1), dtype=float)
    mass_array[0:num_species]    = (temp_array/Avo) * mw_array
    mass_array[num_species]      = (core_abundance/Avo) * core_mw
    total_mass                   = np.sum(mass_array)
    mass_fractions_array         = mass_array/total_mass

    # New Particle Mass [ug m-3]
    condensable_mass_ozone = np.sum(mass_array)*(1.0e6)*(1.0e6) - dry_mass
    
    # Size Distribution Properties
    density = 1.0 / (np.sum( mass_fractions_array / density_array ))
    M2 = solution_ozone[ step,nchem+1 ]
    M3 = solution_ozone[ step,nchem+2 ]
    sigma_array_ozone[step] = ( 1.0/3.0 * np.log( M0 ) - np.log( M2 ) + 2.0/3.0 * np.log( M3 ) ) **0.5
    size_array_ozone[step] = ( M3 / ( M0 * np.exp( 4.5 * sigma_array_ozone[step]**2 ) ) ) **(1.0/3.0)
    
    # Vapor Mass Remaining [ug m-3]
    vapor_mass_ozone = np.sum(solution_ozone[step,0:num_species] / Avo*mw_array*1.0e6*1.0e6)
 
print("Final condensed mass [ozone] = ", condensable_mass_ozone)


# Plot the decrease in gas phase abundance with the growth of the particle [both absolute and relative]
fig = plt.figure(figsize=(10,10))

# Plot the Increase in Particle Size for both Cases
ax1 = plt.subplot(1,2,1)
ax1.plot(t, size_array_ozone*1.0e6, 'r:')
ax1.grid()
ax1.set_xlabel('time [seconds] ')
ax1.set_ylabel('Particle diameter [microns]')
ax1.set_title('Particle growth')

# Plot the Change in Mode Width for both cases
ax2 = plt.subplot(1,2,2)
ax2.plot(t, np.exp(sigma_array_ozone), 'r:',label='With Ozone')
ax2.grid()
ax2.set_xlabel('time [seconds] ')
ax2.set_ylabel(r'$\sigma_g$')
ax2.set_title('Mode Width')
ax2.legend()

plt.show(block=False)


# Now Plot both together
# Plot the decrease in gas phase abundance with the growth of the particle [both absolute and relative]
fig = plt.figure(figsize=(10,10))

# Plot the Increase in Particle Size for both Cases
ax1 = plt.subplot(1,2,1)
ax1.plot(t, size_array_no_ozone*1.0e6, 'r-')
ax1.plot(t, size_array_ozone*1.0e6, 'r:')
ax1.grid()
ax1.set_xlabel('time [seconds] ')
ax1.set_ylabel('Particle diameter [microns]')
ax1.set_title('Particle growth')

# Plot the Change in Mode Width for both cases
ax2 = plt.subplot(1,2,2)
ax2.plot(t, np.exp(sigma_array_no_ozone), 'r-',label='No Ozone')
ax2.plot(t, np.exp(sigma_array_ozone), 'r:',label='With Ozone')
ax2.grid()
ax2.set_xlabel('time [seconds] ')
ax2.set_ylabel(r'$\sigma_g$')
ax2.set_title('Mode Width')
ax2.legend()

plt.show()

 
