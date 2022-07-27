# Calculate the equilibrium mass based on initialising
# a monodisperse population
import numpy as np

#Avogadros number
NA=6.0221409e+23

# Implement equation 2.8
def partitioning(Cstar,abundance,COA_c,core):
    # Partitioning coefficient
    epsilon = np.power(1.0+(Cstar/(COA_c+core)),-1.0)
    # Partitionined mass
    COA_c = np.sum(epsilon*abundance)
    return COA_c

# Implement equation 2.13
def partitioning_dash(Cstar,abundance,COA_c,core):
    epsilon = np.power(1.0+(Cstar/(COA_c+core)),-2.0)*(Cstar/((COA_c+core)**2.0))
    COA_dash = np.sum(epsilon*abundance)
    return COA_dash

# Implement Newtons method
def Newtons_method(Cstar,abundance,COA_init,core):

    COA_c = partitioning(Cstar,abundance,COA_init,core)
    f = COA_c - (partitioning(Cstar,abundance,COA_c,core))
    f_dash = 1.0 - partitioning_dash(Cstar,abundance,COA_c,core)
    COA_new = COA_c - f/f_dash

    while (abs((COA_new-COA_c)/COA_new) > 1.0e-3):

        COA_c = COA_new
        f = COA_c - (partitioning(Cstar,abundance,COA_c,core))
        f_dash = 1.0 - partitioning_dash(Cstar,abundance,COA_c,core)
        COA_new = COA_c - f/f_dash

    return COA_c

# Define the gaseous and condensed phase components
# Number of condensing species from the gas phase
num_species = 10
# Number of size bins
num_bins = 1
# The volatility of each specie, using the mass based C* convention
# Here we assuming a linear seperation in Log10 space
log_c_star = np.linspace(-6, 3, num_species)
Cstar = np.power(10.0,log_c_star)

# Populate volatility basis set with gas phase abundance
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
density_core[:] = 1400.0 #kg.m-3

# The number of particles is only carried in one bin
N_per_bin = np.zeros((num_bins), dtype=float)
N_per_bin[0] = N_total

# Define the diameter of our particles [microns]
size_array = np.zeros((num_bins), dtype=float)
size_array[0] = 0.150 #microns

# Use the size to now calculate a concentration of a 'core' in molecules / cc
# This aligns with the units used for volatile species
core_abundance[0] = (N_per_bin[0])*((4.0/3.0)*
                     np.pi*np.power(size_array[0]*0.5e-6,3.0)*
                     density_core[0]*1.0e3)

# What is our existing dry mass in micrograms per cubic metre?
dry_mass = np.sum((core_abundance)*(1.0e12))

print("Initial mass = ", dry_mass)

# Initialise a value for COA
COA_first_guess = 1.0
# Now use Newtons method for iterating to a final mass
COA_final = Newtons_method(Cstar,abundance,COA_first_guess,dry_mass)

print("Secondary mass = ", COA_final)

# We can also calculate the partitioning coefficients for all volatility bins
epsilon = np.power(1.0+(Cstar/(COA_final+dry_mass)),-1.0)
print("Partitioning coefficients = ", epsilon)
