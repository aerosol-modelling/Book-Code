import numpy as np
import matplotlib.pyplot as plt

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

core=20.0
COA_first_guess = 1.0 + core

def partitioning_coefficient(Cstar,COA_total):
    # In the first instance we need to calculate mass that condenses from the core alone
    # Partitioning coefficient
    epsilon = np.power(1.0+(Cstar/(COA_total)),-1.0)
    return epsilon


###############################################################################
# Define a linearly seperated volatility array
log_c_star = np.linspace(-6, 3, 10) # Array of log10 C* values
Cstar = np.power(10.0,log_c_star) # Array of C* values

# loop through values of total condensed mass, COA_total, and record the
# values of the partitioning coefficient
log_COA = np.linspace(-8,5,100)
# Define a 2D Numpy array to store values of the partitioning coefficient
partitioning_coefficients = np.zeros((100,10),dtype=float)

step=0
for entry in log_COA:
    epsilon = partitioning_coefficient(Cstar,np.power(10.0,entry))
    partitioning_coefficients[step,:]=epsilon
    step+=1

for i in range(10):
    plt.plot(log_COA, partitioning_coefficients[:,i].T,'k')
plt.vlines(log_c_star, 0,1,color='k', linestyle='dashed')
plt.hlines(0.5, log_COA[0], log_COA[99], colors='k', linestyles='dashed')
plt.xlabel(r"Log10($C_{OA}$)", size=14)
plt.ylabel(r'Partitioning coefficient $\varepsilon _{i}$', size=14)
plt.show()
