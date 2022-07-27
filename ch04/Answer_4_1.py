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
mw_array=np.zeros((num_species), dtype=float)
mw_array[:]=200.0

# Define array of log10 C* values
log_c_star = np.linspace(-6, 3, num_species)
Cstar = np.power(10.0,log_c_star)
# Convert C* to a pure component saturation vapour pressure [atm]
P_sat = (Cstar*R_gas_other*Temp_K)/(1.0e6*mw_array)

# Initialise concentration of material in each volatility bin [micrograms/m3]
concentration = np.zeros((num_species), dtype = float)
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
# convert this to  molecules / cc
ozone_concentration = ozone*Cfactor
# Define a rate coefficient between every volatility bin and ozone
rate_coefficient = 1.0e-15

# Unit conversion of gas abudance to molecules / cc
gas_concentration = ((concentration*1.0e-6)/(mw_array))*1.0e-6*NA

# RHS function [with ozone as a global variable]
def dy_dt(array,t):

    # Calculate the rate of each reaction explicitly
    # rate_1 = ozone_abundance*array[0]*rate_coefficient
    rate_1 = ozone_concentration*array[1]*rate_coefficient
    rate_2 = ozone_concentration*array[2]*rate_coefficient
    rate_3 = ozone_concentration*array[3]*rate_coefficient
    rate_4 = ozone_concentration*array[4]*rate_coefficient
    rate_5 = ozone_concentration*array[5]*rate_coefficient
    rate_6 = ozone_concentration*array[6]*rate_coefficient
    rate_7 = ozone_concentration*array[7]*rate_coefficient
    rate_8 = ozone_concentration*array[8]*rate_coefficient
    rate_9 = ozone_concentration*array[9]*rate_coefficient

    # Now write out the RHS for each volatity bin according to our simple
    # mechanism
    dy_dt_array = np.zeros((num_species), dtype=float)
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


# Define the time over which the simulation will take place
t = np.linspace(0, 7200, num=1000)
# Call the ODE solver with reference to our function, dy_dt, setting the
# absolute and relative tolerance.
solution = odeint(dy_dt, gas_concentration, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)

# Plot the decrease in gas phase abundance with the growth of the particle [both absolute and relative]
fig = plt.figure(figsize=(6,6))
plt.plot(t, np.log10(solution[:,0]) , label='Bin 1')
plt.plot(t, np.log10(solution[:,1]) , label='Bin 2')
plt.plot(t, np.log10(solution[:,2]) , label='Bin 3')
plt.plot(t, np.log10(solution[:,3]) , label='Bin 4')
plt.plot(t, np.log10(solution[:,4]) , label='Bin 5')
plt.plot(t, np.log10(solution[:,5]) , label='Bin 6')
plt.plot(t, np.log10(solution[:,6]) , label='Bin 7')
plt.plot(t, np.log10(solution[:,7]) , label='Bin 8')
plt.plot(t, np.log10(solution[:,8]) , label='Bin 9')
plt.plot(t, np.log10(solution[:,9]) , label='Bin 10')
plt.grid()
plt.xlabel('time [seconds] ')
plt.ylabel('Gas phase concentration [Log10 (molecules.cc)]')
plt.legend(loc='lower right')
plt.show()
