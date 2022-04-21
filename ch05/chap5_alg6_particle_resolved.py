import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.stats
import kernel_helper

mono = False
nbins = 50
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

def sample_lognormal(n_part):
    radius = np.zeros(n_part)
    num_conc = [3.2e9,2.9e9]
    log10_std_dev_radius = [.161, .217]
    char_radius = [2e-8/2,1.16e-7/2]
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

t = 0
delta_t = 60.0
diam_edges = np.logspace(-10,-5,nbins+1)
vol_edges = np.zeros(nbins+1)
for i in range(nbins+1):
    vol_edges[i] = (1.0/6.0)*np.pi*diam_edges[i]**3

## Compute the max kernel values per bin pair
MaxKernel = np.zeros((nbins,nbins))
nsample = 3
for i in range(nbins):
    for j in range(nbins):
        k_max = 0.0
        i_sub = np.linspace(vol_edges[i],vol_edges[i+1],nsample)
        j_sub = np.linspace(vol_edges[j],vol_edges[j+1],nsample)
        for ii in range(nsample):
            for jj in range(nsample):
                B_sub = kernel_helper.GetKernel(i_sub[ii],j_sub[jj])
                k_max = max(k_max, B_sub)
        MaxKernel[i,j] = k_max

# Make some initial particles
N_part = int(1e4)
N = np.zeros(nbins,dtype=int)

# Do nested list
M = [[]]
for i in range(nbins):
    M.append([])

#M = np.zeros((nbins,N_part+1))
if (mono):
    for i_part in range(N_part):
        vol = vol_edges[0] #(vol_edges[0] + vol_edges[1])*.5
        i_bin = np.where(vol >= vol_edges)[0][-1]
        M[i_bin].append(vol)
        N[i_bin] += 1
    V = N_part / 1e12
else:
    diams, V = sample_lognormal(N_part)
    N_part = len(diams)
    for i_part in range(N_part):
        vol = diams[i_part]**3 * (np.pi/6) #(vol_edges[0] + vol_edges[1])*.5
        i_bin = np.where(vol >= vol_edges)[0][-1]
        M[i_bin].append(vol)
        N[i_bin] += 1
print(N_part / V)
N_init = N.copy()
rng = np.random.default_rng(12345)
while t < t_final:
    for k in range(nbins):
        for l in range(k+1):
            Kmax = MaxKernel[k, l]  
            if (k != l):
                N_pairs = N[k]*N[l]
            else:
                N_pairs = .5*N[k]*(N[l]-1)
            N_events_exact = N_pairs * Kmax / V
            N_events = np.random.poisson(N_events_exact*delta_t)
            for ll in range(N_events):
                if (N[k] > 0 and N[l] > 0):
                    r_i = rng.integers(low=0, high=N[k]-1, size=1,endpoint=True)[0]
                    r_j = rng.integers(low=0, high=N[l]-1, size=1,endpoint=True)[0]
                    #if (M[k,r_i] == 0 or M[l,r_j] == 0):
                    #    print('problem', N_events,M[k,r_i], M[l,r_j], k ,l, N[k], N[l])
                    #B = kernel(M[k,r_i],M[l,r_j])
                    B = kernel_helper.GetKernel(M[k][r_i],M[l][r_j])
                    r = rng.random()
                    if (r < B/Kmax):
                        new_vol = M[k][r_i] + M[l][r_j]
                        dest_bin = np.where(vol_edges[k:] <= new_vol)[0][-1] + k - 1
                        M[dest_bin].append(new_vol)
                        N[dest_bin] += 1
                        if (r_i > r_j):
                            M[k].pop(r_i)
                            N[k] -= 1
                            M[l].pop(r_j)
                            N[l] -= 1     
                        else:
                            M[l].pop(r_j)
                            N[l] -= 1 
                            M[k].pop(r_i)
                            N[k] -= 1
    t += delta_t
    # When the number of particles decreases too much, we will duplicate and increase volume
    if (np.sum(N) < N_part /2):
        for i_bin in range(nbins):
            M[i_bin] = M[i_bin] + M[i_bin]
            N[i_bin] = 2* N[i_bin]
        V *= 2


fig = plt.figure()
axes = fig.add_subplot(1,1,1)

diams_part = .5*(diam_edges[0:-1] + diam_edges[1:])
#plt.plot(diams_part,(N_init)/V/np.log10(diam_edges[1]/diam_edges[0]),'o--',label='particle')
axes.plot(diams_part,(N)/V/np.log10(diam_edges[1]/diam_edges[0]),'o--',label='particle (final)')
axes.set_xlabel('Particle diameter (m)')
axes.set_ylabel('dN/dlogD')
#plt.yscale('log')
axes.set_xscale('log')
axes.set_ylim([0,1e10])
axes.set_xlim([1e-9,1e-6])

A = np.loadtxt('file.txt')
M0i = A[0]
M3i = A[1]
M0j = A[2]
M3j = A[3]

log_sigma = [.161, .217]
dgi = (M3i/(M0i*np.exp((9/2)*np.log(10**log_sigma[0])**2)))**(1.0/3)
dgj = (M3j/(M0j*np.exp((9/2)*np.log(10**log_sigma[1])**2)))**(1.0/3)

# Lets rebin the particles somehow as the actual grid is coarse.
particle_diams = np.zeros(np.sum(N))
counter = 0
for i_bin in range(nbins):
    for i_part in range(len(M[i_bin])):
        particle_diams[counter] = ((6/np.pi)*M[i_bin][i_part])**(1/3)
        counter += 1

## Maybe we want to rebin the particles onto the sectional grid?
edges = np.logspace(-9,-6,100)
bin_values, bin_edges, binnumber = scipy.stats.binned_statistic(
    particle_diams, particle_diams,statistic='count',bins=edges)
centers = ((edges[0:-1]**-1 + edges[1:]**-1)/2)**-1

i_time = 0
x_axis = np.logspace(-9,-6,250)
vali = get_lognormal(x_axis, dgi[i_time], 10**log_sigma[0], M0i[i_time])
valj = get_lognormal(x_axis, dgj[i_time], 10**log_sigma[1], M0j[i_time])
plt.plot(x_axis,valj+vali,'r',label='modal (init)')
i_time = len(M0i)  -2
vali = get_lognormal(x_axis, dgi[i_time], 10**log_sigma[0], M0i[i_time])
valj = get_lognormal(x_axis, dgj[i_time], 10**log_sigma[1], M0j[i_time])
print(M0i[0]+M0j[0])
print(M0i[i_time] + M0j[i_time])
plt.plot(x_axis,valj+vali,'k',label='modal (final)')
plt.ylabel('$dN/d \log D_p$')
plt.legend();
plt.savefig('size_distribution.pdf')


fig.savefig('particle_resolved.pdf')
