import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib

def get_lognormal(x_axis, d_mean, sigma, n_total):
    """
    Inputs
    x_axis: array of diameters (m)
    d_mean: geometric mean diameter (m)
    sigma: geometric standard deviation (1)
    n_total: total number concentration
    Output: total number concentration, units of n_total
    """

    return (n_total / (np.sqrt(2*np.pi)*np.log10(sigma))) * \
      np.exp(-(np.log10(x_axis) - np.log10(d_mean))**2 / (2*np.log10(sigma)**2))


def sample_lognormal(num_conc, log10_std, char_radius, n_part):
    radius = np.zeros(n_part)
    counter = 0
    V = n_part / np.sum(num_conc)
    n_actual = 0
    for i_mode in range(len(char_radius)):
        n_part_mode = int(n_part * num_conc[i_mode] / np.sum(num_conc))
        x_mean_prime = np.log10(char_radius[i_mode])
        for i_part in range(n_part_mode):
            radius[counter] = 10**np.random.normal(x_mean_prime, log10_std_dev_radius[i_mode])
            counter +=1
        n_actual += n_part_mode
    return(radius[0:n_actual]*2, V)

produce_data = True

if (produce_data):
    num_conc = [3e9,2e9]
    log10_std_dev_radius = [np.log10(1.3),np.log10(1.5)]
    char_radius = [2e-7/2,1e-6]
    diams,V = sample_lognormal(num_conc, log10_std_dev_radius,char_radius,100000)
    np.savetxt('particle_samples.txt', diams)

diams = np.loadtxt('particle_samples.txt')
diams *= 1e6 # Convert to micrometers

n_bins = 200
d_min = 1e-9 * 1e6
d_max = 1e-4 * 1e6
diam_linear = np.linspace(d_min,d_max,n_bins+1)
diam_linear = np.logspace(np.log10(d_min),np.log10(d_max),n_bins+1)
delta_linear = diam_linear[1:]-diam_linear[0:-1]

diam_log = np.logspace(np.log10(d_min),np.log10(d_max),n_bins+1)
delta_log = np.log(diam_log[1]/diam_log[0])

centers_linear = (diam_linear[0:-1] + diam_linear[1:])/2
n_n_linear, bin_edges, binnumber = scipy.stats.binned_statistic(
    diams, diams, statistic='count', bins=diam_linear)

centers_log = ((diam_log[0:-1]**-1 + diam_log[1:]**-1)/2)**-1
n_n_log, bin_edges, binnumber = scipy.stats.binned_statistic(
    diams, diams,statistic='count', bins=diam_log)

surface_area = np.pi*(diams)**2 # Units: um^2
n_s_linear, bin_edges, binnumber = scipy.stats.binned_statistic(
    diams, surface_area, statistic='sum', bins=diam_linear)
n_s_log, bin_edges, binnumber = scipy.stats.binned_statistic(
    diams, surface_area, statistic='sum', bins=diam_log)

volumes = (np.pi/6)*(diams)**3 # Units: um^3
n_v_linear, bin_edges, binnumber = scipy.stats.binned_statistic(
    diams, volumes, statistic='sum', bins=diam_linear)
n_v_log, bin_edges, binnumber = scipy.stats.binned_statistic(
    diams, volumes, statistic='sum', bins=diam_log)

fig = plt.figure(figsize=(12.5,10))
plt.subplots_adjust(wspace=.35,hspace=.3)

plt.subplot(321)
plt.plot(centers_linear,n_n_linear/delta_linear)
plt.ylabel(r'$n(d_{\rm p})$ ($\mu$m$^{-1}$ m$^{-3}$)')
plt.subplot(322)
plt.plot(centers_log,n_n_log/delta_log)
plt.ylabel(r'$\tilde{n}(\log d_{\rm p})$ (m$^{-3}$)')
plt.subplot(323)
plt.plot(centers_linear,n_s_linear/delta_linear)
plt.plot(centers_linear,np.pi*centers_linear**2*n_n_linear/delta_linear,'k')
plt.ylabel(r'$S(d_{\rm p})$ ($\mu$m$^{1}$ m$^{-3}$)')
plt.subplot(324)
plt.plot(centers_log,n_s_log/delta_log)
plt.plot(centers_log,np.pi*centers_log**2*n_n_log/delta_log,'k')
plt.ylabel(r'$\tilde{S}(\log d_{\rm p})$ ($\mu$m$^{2}$ m$^{-3}$)')
plt.subplot(325)
plt.plot(centers_linear,n_v_linear/delta_linear)
plt.plot(centers_linear,np.pi/6*centers_linear**3*n_n_linear/delta_linear,'k')
plt.ylabel(r'$V(d_{\rm p})$ ($\mu$m$^{2}$ m$^{-3}$)')
plt.subplot(326)
plt.plot(centers_log,n_v_log/delta_log)
plt.plot(centers_log,np.pi/6*centers_log**3*n_n_log/delta_log,'k')
plt.ylabel(r'$\tilde{V}(\log d_{\rm p})$ ($\mu$m$^{3}$ m$^{-3}$)')

allaxes = fig.get_axes()
for i_axes in range(len(allaxes)):
    allaxes[i_axes].set_xlabel(r"Particle diameter $d_{\rm p}$($\mu$m)")
    allaxes[i_axes].set_xscale('log')
    allaxes[i_axes].set_ylim(bottom=0)
    allaxes[i_axes].set_xlim([1e-3,1e2])
    
plt.savefig('ch05_q01a.pdf')


# ## 1(c) solution

# Generate data
fig = plt.figure(figsize=(5,5))
num_conc = [3e9]
sigma = [1.3]
char_radius = [2e-7/2]
d,V = sample_lognormal(num_conc, np.log10(sigma), char_radius, 100000)
d *= 1e6
# Plot data
diam_axis = np.logspace(-3,1,1000)
bin_values, bin_edges, binnumber = scipy.stats.binned_statistic(
    d, d, statistic='count', bins=diam_axis)
centers = ((diam_axis[0:-1]**-1 + diam_axis[1:]**-1)/2)**-1
plt.plot(centers,bin_values/np.log10(diam_axis[1]/diam_axis[0])/V, label='sampled data')
plt.xscale('log')

vals = get_lognormal(diam_axis*1e-6, char_radius[0]*2, sigma[0], num_conc[0])
plt.plot(diam_axis,vals,'k',label='log normal (exact)')
plt.legend()
plt.ylabel(r'$\tilde{n}(\log d_{\rm p})$ (m$^{-3}$)')
plt.xlabel(r'Particle diameter $d_{\rm p}$ ($\mu$m)');
plt.savefig('ch05_q01c.pdf')
