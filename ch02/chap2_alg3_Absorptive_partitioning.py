import numpy as np

#Avogadros number
NA=6.02214076e+23

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

    # Iterate estimates of COA until tolerance met
    while (abs((COA_new-COA_c)/COA_new) > 1.0e-3):

        COA_c = COA_new
        f = COA_c - (partitioning(Cstar,abundance,COA_c,core))
        f_dash = 1.0 - partitioning_dash(Cstar,abundance,COA_c,core)
        COA_new = COA_c - f/f_dash

    return COA_c

# Set an initial bulk core mass and condensed secondary mass
core=5.0
COA_first_guess = 1.0
print("Core absorptive mass = ", core)

# Populate volatility basis set with gas phase abundance
abundance = np.zeros((10), dtype = float)
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
print("Available secondary mass = ", np.sum(abundance))

# Define array of log10 C* values
log_c_star = np.linspace(-6, 3, 10)
Cstar = np.power(10.0,log_c_star)

# Call our function 'Newtons_method'
COA_final = Newtons_method(Cstar,abundance,COA_first_guess,core)
print("Secondary mass = ", COA_final)
