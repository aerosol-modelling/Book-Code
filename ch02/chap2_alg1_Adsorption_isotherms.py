# Script to apply 3 different adsorption isotherms
# Import the libraries we will be using
import numpy as np
import matplotlib.pyplot as plt
import time

# Isotherm parameters for water on NaCl
A = 0.91
B = 0.67
c = 1.68

# Define functions for both adsorption theories
def BET(c,S):
    theta = c*S/((1.0-S)*(1.0-S+c*S))
    return theta

def FHH(A,B,S):
    theta = np.power((A/(-1.0*np.log(S))),1.0/B)
    return theta

# Initialise a list of S values
# Arrays for 1st option
S_array1 = []
BET_theta_array1 = []
FHH_theta_array1 = []

# Arrays for 2nd option
S_array2 = []
BET_theta_array2 = []
FHH_theta_array2 = []

# Arrays for 3rd option
S_array3 = np.linspace(0.01, 0.9, 90) # Create a Numpy array of S
BET_theta_array3 = []
FHH_theta_array3 = []

###############################################################################
# 1) Loop through values of S, the saturation ratio
for step in range(1,90):
    S = step / 100.0
    S_array1.append(S)
    BET_theta =  c*S/((1.0-S)*(1.0-S+c*S))
    FHH_theta = np.power((A/(-1.0*np.log(S))),1.0/B)
    BET_theta_array1.append(BET_theta)
    FHH_theta_array1.append(FHH_theta)

# Plot results from the first option
fig, ax = plt.subplots()
ax.plot(S_array1, BET_theta_array1, "-b", label="BET")
ax.plot(S_array1, FHH_theta_array1, "-r", label="FHH")
plt.legend(loc="upper left")
#plt.ylim(-1.5, 2.0)
ax.set(xlabel='Saturation ratio, S', ylabel=r'$\theta$',
       title='BET and FHH adsorption isotherms')
ax.grid()
plt.show()

# 2) Loop through values of S by calling the dedicated functions
for step in range(1,90):
    S = step / 100.0
    S_array2.append(S)
    BET_theta = BET(c,S)
    FHH_theta = FHH(A,B,S)
    BET_theta_array2.append(BET_theta)
    FHH_theta_array2.append(FHH_theta)

# Plot results from the first option
fig, ax = plt.subplots()
ax.plot(S_array2, BET_theta_array2, "-b", label="BET")
ax.plot(S_array2, FHH_theta_array2, "-r", label="FHH")
plt.legend(loc="upper left")
#plt.ylim(-1.5, 2.0)
ax.set(xlabel='Saturation ratio, S', ylabel=r'$\theta$',
       title='BET and FHH adsorption isotherms')
ax.grid()
plt.show()

# 3) Pass S as a Numpy array
BET_theta_array3 = BET(c,S_array3)
FHH_theta_array3 = FHH(A,B,S_array3)

# Plot results from the first option
fig, ax = plt.subplots()
ax.plot(S_array3, BET_theta_array3, "-b", label="BET")
ax.plot(S_array3, FHH_theta_array3, "-r", label="FHH")
plt.legend(loc="upper left")
#plt.ylim(-1.5, 2.0)
ax.set(xlabel='Saturation ratio, S', ylabel=r'$\theta$',
       title='BET and FHH adsorption isotherms')
ax.grid()
plt.show()
###############################################################################

###############################################################################
# Profile the time spent on each option through some basic time differences
tic = time.perf_counter()
for step1 in range(1,100000):
    for step in range(1,90):
        S = step / 100.0
        BET_theta =  c*S/((1.0-S)*(1.0-S+c*S))
        FHH_theta = np.power((A/(-1.0*np.log(S))),1.0/B)
toc = time.perf_counter()
print(f"Ran through the calculations for option 1 in {toc - tic:0.4f} seconds")

tic = time.perf_counter()
for step1 in range(1,100000):
    for step in range(1,90):
        S = step / 100.0
        BET_theta = BET(c,S)
        FHH_theta = FHH(A,B,S)
toc = time.perf_counter()
print(f"Ran through the calculations for option 2 in {toc - tic:0.4f} seconds")

tic = time.perf_counter()
for step1 in range(1,100000):
    BET_theta_array3 = BET(c,S_array3)
    FHH_theta_array3 = FHH(A,B,S_array3)
toc = time.perf_counter()
print(f"Ran through the calculations for option 3 in {toc - tic:0.4f} seconds")
###############################################################################
