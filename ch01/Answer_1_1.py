import numpy as np
import matplotlib.pyplot as plt

# Volume ratio discrete distribution
d1 = 0.01 # Lowest size diameter
d_Nb = 1.0 # Diameter of largest bin [microns]
Nb = 30 # Number of size bins

di = np.exp(np.linspace(np.log(d1), np.log(d_Nb), num=30)) # values for x-axis

V_rat = np.power((d_Nb/d1),3.0/(Nb-1.0)) # Volume ratio between bins
d_rat = V_rat**(1.0/3.0) # Diameter ratio between bins

# Use the volume ratio to create an array of diameters as follows
# Create an empty diameter array
d_i=np.zeros((Nb), dtype=float)
d_i[0]=d1
for step in range(Nb):
    if step > 0:
       d_i[step]=d_i[step-1]*d_rat

# Log of Diameter array
log_di = np.log(d_i)

# Create an array of lower bin boundaries
vi = (4.0/3.0)*np.pi*np.power(d_i/2.0,3.0)
vi_low = 2.0*vi/(1.0+V_rat)
di_low = 2.0*(np.power(vi_low/((4.0/3.0)*np.pi),(1.0/3.0)))

d_width = d_i*np.power(2,1.0/3.0)*((np.power(V_rat,1.0/3.0)-1.0)/
          (np.power(1+V_rat,1.0/3.0))) # Diameter width array of size bins

# Parameters of log-normal distribution
sigmag1 = np.log(1.7) # Geometric standard deviation
mean1 = np.log(0.15) # Mean particle size [150nm]
# Seperate out the probability density function
distribution_1 = (np.exp(-(log_di - mean1)**2 / (2 * sigmag1**2)) /
     (sigmag1 * np.sqrt(2 * np.pi)))

Ntot = 600.0 # Total number of particles [per cm-3]
# Use pre-calculated distribution_1 to implement equation 1.16
N_dist = Ntot*(distribution_1*(d_width/d_i)) # Discrete number distribution

plt.bar(log_di, N_dist)
plt.ylabel(r'dNdLogDp')
plt.xlabel(r'LogDp')
plt.title(r'Discrete size distribution')
plt.show()
plt.show()
