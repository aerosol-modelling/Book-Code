# Import the relevant libraries, including an ODE solver from Scipy
# We only want to use the 'odeint' [Solve Initial Value Problems] from the
# scipy.integrate package
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt # Import Matplotlib so we can plot results

# Physical constants
R_gas=8.31451 #Ideal gas constant [kg m2 s-2 K-1 mol-1]
R_gas_other=8.20573e-5 #Ideal gas constant [m3 atm K-1 mol-1]# Avogadros number
NA=6.0221409e+23

# Set the ambient temperature
Temp_K=298.15

# Number of condensing species from the gas phase
num_species = 10

# The molecular weight of each condensing specie [g/mol]
# Assuming a constant value for all components
mw_array=np.zeros((num_species+1), dtype=float)
mw_array[0:10]=200.0 # Assuming a constant value for all components
mw_array[10]=47.998

# Define array of log10 C* values
log_c_star = np.linspace(-6, 3, num_species)
Cstar = np.power(10.0,log_c_star)
# Convert C* to a pure component saturation vapour pressure [atm]
P_sat = (Cstar*R_gas_other*Temp_K)/(1.0e6*mw_array[0:10])

# Initialise concentration of material in each volatility bin [micrograms/m3]
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

# Initialise concentration of a gas phase oxidant in ppb
ozone = 10.0
Cfactor= 2.46310e+10 #ppb-to-molecules/cc
ozone_concentration = ozone*Cfactor
# Now convert to micrograms/m3
ozone_concentration = (ozone_concentration/(NA*1.0e-6))*(mw_array[10]*1.0e6)
concentration[10] = ozone_concentration # Ozone

# Define a rate coefficient between every volatility bin and ozone
rate_coefficient = 1.0e-15

# Unit conversion of gas abudance to molecules / cc
gas_concentration = ((concentration*1.0e-6)/(mw_array))*1.0e-6*NA

# Now let us convert the given emission rate of ozone into molecules / cc.second
emission_rate_ozone = ((10.0*1.0e-12)/(3600.0))*(NA/mw_array[10])

# We can also calculate the emission rate of the highest volatility bin
emission_rate_bin = ((1.5*1.0e-12)/(3600.0))*(NA/mw_array[9])

# RHS function [with ozone as a global variable]
def dy_dt1(array,t):

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
    dy_dt_array = np.zeros((num_species+1), dtype=float)
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

    return dy_dt_array

# RHS function [with ozone as variable in our array]
def dy_dt2(array,t):

    # Calculate the rate of each reaction explicitly
    # We now track the concentration of ozone in our array
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
    dy_dt_array = np.zeros((num_species+1), dtype=float)
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

    return dy_dt_array


# Define the time over which the simulation will take place
t = np.linspace(0, 7200, num=1000)
# Call the ODE solver with reference to our function, dy_dt, setting the
# absolute and relative tolerance.
solution1 = odeint(dy_dt1, gas_concentration, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)
solution2 = odeint(dy_dt2, gas_concentration, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)

# Plot the decrease in gas phase abundance with the growth of the particle [both absolute and relative]
fig = plt.figure(figsize=(6,6))
plt.plot(t, solution1[:,9]/solution2[:,9])
plt.grid()
plt.xlabel('time [seconds] ')
plt.ylabel('Ratio')
#leg = plt.legend();
plt.legend(loc='lower right')
plt.show()
