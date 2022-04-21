import numpy as np
import kernel_helper
import matplotlib.pyplot as plt
import time


(mb, v_factor, t_final) = kernel_helper.GetParameters()

def main():

    t = 0.0
    N_part = 100
    N = np.zeros(N_part)
    N_bin = len(N)
    V = (N_part / 100) * v_factor
    print('Initial number concentration %e' %(N_part/V))
    max_size = N_part 
    hist = np.zeros((N_part,N_bin),dtype=int)
    N[0] = 1.0*N_part
    bin_size = np.arange(1,N_bin+1)
    times = [0.0]
    hist[0,:] = N 
    index = 1 
    t1 = time.time()
    while t < t_final and np.sum(N) > 1:
        tau = np.full((N_bin,N_bin),np.Inf, dtype=float)
        for ind_k, k in enumerate(bin_size):
            for ind_l, l in enumerate(bin_size[0:ind_k+1]):
                K12 = kernel_helper.GetKernel(k*mb,l*mb)
                if (k == l):
                    N_pairs = (N[ind_k]*(N[ind_l]-1)) / 2
                else:
                    N_pairs = N[ind_k]*N[ind_l]
                lamb = N_pairs*K12/V
                u = np.random.rand()
                if (lamb > 0):
                    tau[ind_k,ind_l] = -np.log(u)/lamb
        # Find the index of the min tau
        index_array = np.argmin(tau) # argmin returns index into flattened array
        (p1,p2) = np.unravel_index(index_array,shape=tau.shape) # convert flat index to index tuple
        dt = tau[p1,p2]
        N[p1] -= 1
        N[p2] -= 1
        N[p1+p2+1] += 1
        t += dt
        times.append(t)
        hist[index,:] = N
        index += 1

    t2 = time.time()
    print(t2-t1)

    return(times, hist[0:len(times),:])

fig = plt.figure(figsize=(4.5,3.5))
axes = fig.add_subplot(1,1,1)
max_size_plot = 8 
times, histogram = main()
for i in range(max_size_plot):
    axes.plot(times, histogram[:,i], label='$N_{%i}$' %(i+1))
axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0,t_final,num=7))
axes.set_ylim([0,100])
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')
fig.savefig('fig_method2.pdf',bbox_inches='tight')
### Do an ensemble

from scipy import interpolate
import scipy.stats

n_interp = 3600
x = np.linspace(0,t_final,n_interp)
fig = plt.figure(figsize=(4.5,3.5))
axes = fig.add_subplot(1,1,1)
n_runs = 20
max_size = 8 
values = np.zeros((n_interp,n_runs,max_size))
for i_run in range(n_runs):
    times, hist = main()
    for i_size in range(max_size):
        f = interpolate.interp1d(times, hist[:,i_size],fill_value="extrapolate")
        values[:,i_run,i_size] = f(x)

for i_size in range(max_size):
    mean_vals = np.mean(values[:,:,i_size],axis=1)
    std_val = np.std(values[:,:,i_size],axis=1)
    error = scipy.stats.t.ppf(.975, n_runs)*std_val #/ np.sqrt(n_runs)
    axes.plot(x,mean_vals, label='$N_{%i}$' %(i_size+1))
    axes.fill_between(x, mean_vals-error, mean_vals+error,color=axes.lines[-1].get_color(), alpha=.3,lw = 0)
#    axes.errorbar(x, mean_vals, yerr=error)
axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0,t_final,num=7))
axes.set_ylim([0,100])
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')
fig.savefig('fig_method2_ensemble.pdf',bbox_inches='tight')
