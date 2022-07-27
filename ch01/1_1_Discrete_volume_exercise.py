import numpy as np
import matplotlib.pyplot as plt

# Volume ratio discrete distribution
d1 = 0.01 # Lowest size diameter
d_Nb = 1.0 # Diameter of largest bin [microns]
Nb = 30 # Number of size bins

V_rat =  # Volume ratio between bins
d_rat =  # Diameter ratio between bins

# Use the volume ratio to create an array of diameters as follows
# Create an empty diameter array
d_i=np.zeros((Nb), dtype=float) # Diameter array
d_i[0]=d1
# Populate values in array d_i
for step in range(Nb):
    if step > 0:
       d_i[step]=

# Log of Diameter array
log_di = np.log(d_i) # Log of Diameter array

# Create an array of lower bin boundaries
vi = (4.0/3.0)*np.pi*np.power(d_i/2.0,3.0)
vi_low = 2.0*vi/(1.0+V_rat)
di_low = 2.0*(np.power(vi_low/((4.0/3.0)*np.pi),(1.0/3.0)))

# Create a diameter width array
d_width =

# Set parameters of log-normal distribution
sigmag1 =  # Geometric standard deviation
mean1 =  # Mean particle size [150nm]

Ntot = 600.0 # Total number of particles [per cm-3]
# Implement equation 1.16
N_dist =
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

plt.bar(log_di, N_dist)
plt.ylabel(r'dNdLogDp')
plt.xlabel(r'LogDp')
plt.title(r'Discrete size distribution')
plt.show()
plt.show()
