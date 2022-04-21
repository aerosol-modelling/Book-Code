# Import the relevant libraries
import numpy as np
import matplotlib.pyplot as plt

# Define lower and upper limits
lower_size = -10.0
upper_size = 10.0

# Create 400 linearly seperated values for our 'x' variable
x = np.linspace(lower_size, upper_size, num=400) # values for x-axis

# Define the mean and standard deviation for each normal distribution
sigma_1 = 0.3
mean_1 = 0.0

sigma_2 = 1.0
mean_2 = 0.0

sigma_3 = 1.3
mean_3 = 0.0

# Implement equation 1.5
distribution_1 = (np.exp(-(x - mean_1)**2 / (2 * sigma_1**2)) / 
    (sigma_1 * np.sqrt(2 * np.pi)))
distribution_2 = (np.exp(-(x - mean_2)**2 / (2 * sigma_2**2)) / 
    (sigma_2 * np.sqrt(2 * np.pi)))
distribution_3 = (np.exp(-(x - mean_3)**2 / (2 * sigma_3**2)) / 
    (sigma_3 * np.sqrt(2 * np.pi)))

# Plot the results
line1 = plt.plot(x, distribution_1 , linewidth=2, color='r',
    label='$\sigma$ = 0.3')
line2 = plt.plot(x, distribution_2 , linewidth=2, color='b',
    label='$\sigma$ = 1.0')
line3 = plt.plot(x, distribution_3 , linewidth=2, color='k',
    label='$\sigma$ = 1.3')
ax = plt.gca()
ax.legend()
ax.set_xlabel('x')
ax.set_ylabel('p(x)')
plt.show()
