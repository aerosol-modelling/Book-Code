import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.stats
import kernel_helper

mono = False
t_final = 12*3600

def get_lognormal(x_axis, d_mean, sigma, n_total):
    """
    Inputs
    x_axis: array of diameters (m)
    d_mean: geometric mean diameter (m)
    sigma: geometric standard deviation (1)
    n_total: total number concentration
    -------
    Output: total number concentration, units of n_total
    """
    return (n_total / (np.sqrt(2*np.pi)*np.log10(sigma))) * \
      np.exp(-(np.log10(x_axis) - np.log10(d_mean))**2 / (2*np.log10(sigma)**2))

def integrate_lognormal(d1, d2, d_mean, sigma):
    """
    Inputs
    d1: lower diameter range (m)
    d2: upper diameter (m)
    d_mean: geometric mean diameter (m)
    sigma: geometric standard deviation (1)
    ------
    Output: integrate number concentration between d1 and d2.
    """
    total =  .5 * (math.erf((np.log(d2) - np.log(d_mean))/(np.sqrt(2)*np.log(sigma))) - \
        math.erf((np.log(d1) - np.log(d_mean))/(np.sqrt(2)*np.log(sigma))))
    return(total) 

### Setup the sectional grid
dmin = 1e-9
dmax = 1e-6
V_ratio = 1.3
nbins = 1 + int(np.log((dmax/dmin)**3)/np.log(V_ratio))
print(nbins)
diams_sectional = np.logspace(np.log10(dmin),np.log10(dmax),nbins)

vol = np.zeros(nbins)
for i in range(nbins):
    vol[i] = (1.0/6.0)*np.pi*diams_sectional[i]**3

f = np.zeros((nbins,nbins,nbins))
for i in range(nbins):
    for j in range(nbins):
        total_volume = vol[i] + vol[j]
        for k in range(nbins):
            if (k < nbins-1):
                if (total_volume >= vol[k] and total_volume < vol[k+1]):
                    f[i,j,k] = ((vol[k+1] - total_volume)/(vol[k+1]-vol[k])) * (vol[k] / total_volume)
                if (k > 0):
                    if (total_volume > vol[k-1] and total_volume < vol[k]):
                        f[i,j,k] = 1.0 - f[i,j,k-1]
            elif (k == nbins-1):
                if (total_volume >= vol[k]):
                    f[i,j,k] = 1.0

# Kernel is constant in time since the environment is constant
B = np.zeros((nbins,nbins))
for i in range(nbins):
    for j in range(nbins):
        B[i,j] = kernel_helper.GetKernel(vol[i],vol[j])

C = np.zeros(nbins)

if (mono):
    C[0] = 1e12
else:
    num_conc = [3.2e9,2.9e9]
    log10_std_dev_radius = [.161, .217]
    char_radius = [2e-8/2,1.16e-7/2]
    for i_mode in range(2):
        for i_bin in range(nbins-1):
            C[i_bin] += num_conc[i_mode] * integrate_lognormal(diams_sectional[i_bin],diams_sectional[i_bin+1],
                                    2*char_radius[i_mode],10**log10_std_dev_radius[i_mode])
C_init = C.copy()
delta_t = 60.0
ntimes = int((t_final)/delta_t)
C_new = np.zeros(nbins)
print('total number:', np.sum(C))
for t in range(ntimes):
    for k in range(nbins):
        sum_gain = 0.0
        sum_loss = 0.0
        for j in range(k+1):
            for i in range(k):
                sum_gain += f[i,j,k] * B[i,j] * vol[i] * C_new[i] * C[j]
        for j in range(nbins):
            sum_loss += (1.0 - f[k,j,k]) * B[k,j] * C[j]
        C_new[k] = (vol[k]*C[k] + ((delta_t * sum_gain))) / (1.0 + delta_t * sum_loss )/ vol[k]
    C = C_new

diam_centers = np.zeros(nbins)
for i_bin in range(nbins):
    if (i_bin == nbins-1):
        diam_centers[i_bin] = diams_sectional[i_bin]
    else:
        diam_centers[i_bin] = ((diams_sectional[i_bin]**-1 + diams_sectional[i_bin+1]**-1)/2)**-1


fig = plt.figure()
axes = fig.add_subplot(1,1,1)
    
axes.plot(diam_centers,C/np.log10(diams_sectional[1]/diams_sectional[0]),'*--')
axes.plot(diam_centers,C_init/np.log10(diams_sectional[1]/diams_sectional[0]),'*--')
axes.set_xscale('log')
#plt.ylim([1e0,1e8])
axes.set_xlabel('Particle diameter (m)')
axes.grid(True)
axes.set_ylabel('dN/dlogD cm$^{-3}$')
axes.set_xlim([1e-9,1e-6])
fig.savefig('sectional_method.pdf')
