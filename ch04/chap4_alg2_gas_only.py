# Import the relevant libraries, including an ODE solver from Scipy
# We only want to use the 'odeint' [Solve Initial Value Problems] from the
# scipy.integrate package
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt # Import Matplotlib so we can plot results

# Ideal gas constant [m3 atm K-1 mol-1]
R_gas_other=8.2057e-5
# Avogadro's number
NA=6.0221409e+23

# Set the ambient temperature
Temp_K=298.15

# Defining the gas phase components
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

# Initialise abundance of a gas phase oxidant in ppb
ozone = 10.0
Cfactor= 2.55e+10 #ppb-to-molecules/cc
# convert abundance of oxidant to molecules / cc
ozone_abundance = ozone*Cfactor
# Define a rate coefficient for second order reactions
rate_coefficient = 1.0e-15

# Unit conversion of gas abundance to molecules / cc
gas_abundance = ((abundance*1.0e-6)/(mw_array))*1.0e-6*NA

# New define an array that holds the molecular abundance
array = np.zeros((num_species), dtype=float)
array[0:num_species] = gas_abundance

###############################################################################
# Define RHS function
def dy_dt(array,t):

    dy_dt_array = np.zeros((num_species), dtype=float)

    # Calculate the rate of each reaction explicitly
    rate_1 = ozone_abundance*array[1]*rate_coefficient
    rate_2 = ozone_abundance*array[2]*rate_coefficient
    rate_3 = ozone_abundance*array[3]*rate_coefficient
    rate_4 = ozone_abundance*array[4]*rate_coefficient
    rate_5 = ozone_abundance*array[5]*rate_coefficient
    rate_6 = ozone_abundance*array[6]*rate_coefficient
    rate_7 = ozone_abundance*array[7]*rate_coefficient
    rate_8 = ozone_abundance*array[8]*rate_coefficient
    rate_9 = ozone_abundance*array[9]*rate_coefficient

    # Write RHS for each volatility bin according to our simple
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

    return dy_dt_array

t = np.linspace(0, 3000, num=1000)

solution_ozone = odeint(dy_dt, array, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)

# Plot the change in gas phase abundance
labels = ['-6', '-5', '-4', '-3', '-2', '-1', '0', '1', '2', '3']


fig = plt.figure(figsize=(10,10))

plt.subplot(1, 2, 1)
plt.plot(t, np.log10(solution_ozone[:,0]), label=r"$log10(C^{*})=-6$")
plt.plot(t, np.log10(solution_ozone[:,1]), label=r'$log10(C^{*})=-5$')
plt.plot(t, np.log10(solution_ozone[:,2]), label=r'$log10(C^{*})=-4$')
plt.plot(t, np.log10(solution_ozone[:,3]), label=r'$log10(C^{*})=-3$')
plt.plot(t, np.log10(solution_ozone[:,4]), label=r'$log10(C^{*})=-2$')
plt.plot(t, np.log10(solution_ozone[:,5]), label=r'$log10(C^{*})=-1$')
plt.plot(t, np.log10(solution_ozone[:,6]), label=r'$log10(C^{*})=0$')
plt.plot(t, np.log10(solution_ozone[:,7]), label=r'$log10(C^{*})=1$')
plt.plot(t, np.log10(solution_ozone[:,8]), label=r'$log10(C^{*})=2$')
plt.plot(t, np.log10(solution_ozone[:,9]), label=r'$log10(C^{*})=3$')
plt.grid()
plt.xlabel('time [seconds] ')
plt.ylabel('Concentration of each volatility bin [molecules/cc]')
plt.title('Evolving gas phase model')
plt.legend(loc=3)

plt.subplot(2, 2, 2)
plt.bar(labels, (solution_ozone[0,:]/(1.0e-12*NA))*(mw_array), color='green')
plt.xlabel(r'$log10(C^{*})$')
plt.ylabel(r'$\mu g.m^{-3}$')
plt.title("Initial conditions")
plt.xticks(labels, labels)

plt.subplot(2, 2, 4)
plt.bar(labels, (solution_ozone[-1,:]/(1.0e-12*NA))*(mw_array), color='red')
plt.xlabel(r'$log10(C^{*})$')
plt.ylabel(r'$\mu g.m^{-3}$')
plt.title("Final conditions")
plt.xticks(labels, labels)
plt.tight_layout()
plt.show()
