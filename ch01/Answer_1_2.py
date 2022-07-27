# Import the relevant libraries, including an ODE solver from Scipy
# We only want to use the 'odeint' [Solve Initial Value Problems] from the
# scipy.integrate package
import numpy as np
import matplotlib.pyplot as plt # Import Matplotlib so we can plot results
import math
# -------------------------------------------------------

# ----- Exercise 2 -----
# Initialize Number, Surface Area, and Volume Concentrations
N = 2300. * 1.0E6  #particles m-3
S = 2.144e-4       #m2 m-3
M = 27.2          #ug m-3

# Calculate Moments of Modal Distribution
V  = M*1.0e-9 / 1400 # m3 m-3
M3 = V / np.pi * 6 # m3 m-3
M2 = S / np.pi     # m2 m-3

# Calculate and print properties of log-normal mode
sigma = np.exp( np.sqrt( 1./3. *np.log(N) \
               -np.log(M2) + 2./3. * np.log(M3) ) )
mu = ( M3 / N / \
       np.exp( 9./2. * (np.log(sigma))**2 ) ) **(1./3.)
print( 'mu=',mu,' sigma=',sigma )

# ----- Exercise 3 -----
# Define  lower  and  upper  limit of the  size  distribution [microns]
Dlo = 50e-9
Dhi = 400e-9

# Calculate Number Concentration between Dlo and Dhi
mug = mu
Nslice = N * 0.5 * \
        ( math.erf((np.log(Dhi) - np.log(mug))/  \
                   (np.sqrt(2)*np.log(sigma)))   \
         -math.erf((np.log(Dlo) - np.log(mug))/  \
                   (np.sqrt(2)*np.log(sigma))) )

# Calculate Surface Area Concentration between Dlo and Dhi
mug = np.exp( np.log(mu) + 2 * np.log(sigma)**2 )
Sslice = S * 0.5 * \
        ( math.erf((np.log(Dhi) - np.log(mug))/  \
                    (np.sqrt(2)*np.log(sigma)))  \
         -math.erf((np.log(Dlo) - np.log(mug))/  \
                    (np.sqrt(2)*np.log(sigma))) )

# Calculate Mass Concentration between Dlo and Dhi
mug = np.exp( np.log(mu) + 3 * np.log(sigma)**2 )
Vslice = V * 0.5 * \
        ( math.erf((np.log(Dhi) - np.log(mug))/  \
                   (np.sqrt(2)*np.log(sigma)))   \
         -math.erf((np.log(Dlo) - np.log(mug))/  \
                   (np.sqrt(2)*np.log(sigma))) )
Mslice = 1400 * Vslice * 1.0e9 # ug m-3

print( 'N(50-400nm)=',Nslice*1.0e-6, ' S(50-400nm)=',Sslice, ' M(50-400nm)=',Mslice)

