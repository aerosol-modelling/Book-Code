# Here we plot 3 different log-normal distributions with a mean
# diameter of 800nm and geometric standard deviation of 1.1, 1.3 and 2.0
# We assume our smallest size is 10nm and largest is 20 microns

# Import the relevant libraries
import numpy as np
import matplotlib.pyplot as plt

# Define lower and upper limit of the size distribution [microns]
lower_size = 0.01
upper_size = 20

# Create an array of values in log space
x_log = np.linspace(np.log(lower_size), np.log(upper_size), num=400)

# Define geomatric standard deviations and mean values
sigmag_1 = 0.3
mean_1 = np.log(0.8)
sigmag_2 = 0.6
mean_2 = np.log(0.8)
sigmag_3 = 1.1
mean_3 = np.log(0.8)

# Implement equation 1.7
distribution_1 = (np.exp(-(x_log - mean_1)**2 / (2 * sigmag_1**2)) / (sigmag_1 * np.sqrt(2 * np.pi)))
distribution_2 = (np.exp(-(x_log - mean_2)**2 / (2 * sigmag_2**2)) / (sigmag_2 * np.sqrt(2 * np.pi)))
distribution_3 = (np.exp(-(x_log - mean_3)**2 / (2 * sigmag_3**2)) / (sigmag_3 * np.sqrt(2 * np.pi)))

# Plot the results
plt.plot(x_log, distribution_1 , linewidth=2, color='r')
plt.plot(x_log, distribution_2 , linewidth=2, color='b')
plt.plot(x_log, distribution_3 , linewidth=2, color='k')
ax = plt.gca()
ax.legend()
ax.set_xlabel('x')
ax.set_ylabel('p(x)')
plt.show()
