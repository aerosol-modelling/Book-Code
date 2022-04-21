import numpy as np
import kernel_helper
import matplotlib.pyplot as plt
import time
import scipy.stats


(mb, v_factor, t_final) = kernel_helper.GetParameters()
nt = 60

def main():

    t = 0.0
    N_part = 100
    N = np.zeros(N_part)
    N_bin = len(N)
    V = (N_part / 100) * v_factor 
    print('Initial number concentration %e' %(N_part/V))
    max_size = N_part 
    index = 0
    N[0] = 1.0*N_part
    bin_size = np.arange(1,101)
    times = np.linspace(0,t_final,nt+1)
    dt = t_final / nt
    hist = np.zeros((len(times),N_part))
    hist[0,:] = N 
    index = 1
    t1 = time.time()
    while t < t_final and np.sum(N) > 1:
        for ind_k, k in enumerate(bin_size):
            for ind_l, l in enumerate(bin_size[0:ind_k+1]):
                K12 = kernel_helper.GetKernel((k)*mb,(l)*mb)
                if (k == l):
                    N_pairs = (N[ind_k]*(N[ind_l]-1)) / 2
                else:
                    N_pairs = N[ind_k]*N[ind_l]
                if (N_pairs < 0):
                    print(N_pairs,k,l, N[ind_k],N[ind_l])
                lamb = N_pairs*K12/V
                N_event = np.random.poisson(lamb*dt)
                N[ind_k] -= N_event
                N[ind_l] -= N_event
                if (N[ind_k] < 0):
                    N[ind_k] = 0
                if (N[ind_l] < 0):
                    N[ind_l] = 0
                dest = min(k + l, N_bin)
                N[dest-1] += N_event
        t += dt
        hist[index,:] = N 
        index += 1

    t2 = time.time()
    print(t2-t1)

    return(times, hist)

times, hist = main()
fig = plt.figure(figsize=(4.5,3.5))
axes = fig.add_subplot(1,1,1)
max_size_plot = 8 
for i in range(max_size_plot):
    axes.plot(times, hist[:,i], label='$N_{%i}$' %(i+1))
axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0,t_final,num=7))
axes.set_ylim([0,100])
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')

print(axes.get_position())
fig.savefig('fig_method3.pdf',bbox_inches='tight')

### Do an ensemble

from scipy import interpolate
fig = plt.figure(figsize=(4.5,3.5))
axes = fig.add_subplot(1,1,1)
n_runs = 20
max_size = 100
values = np.zeros((nt+1,n_runs,max_size))
for i_run in range(n_runs):
    times, hist = main()
    values[:,i_run,:] = hist 

for i_size in range(max_size_plot):
    mean_vals = np.mean(values[:,:,i_size],axis=1)
    std_val = np.std(values[:,:,i_size],axis=1)
    error = scipy.stats.t.ppf(.975, n_runs)*std_val #/ np.sqrt(n_runs)
    axes.plot(times,mean_vals,label='$N_{%i}$' %(i_size+1))
    axes.fill_between(times, mean_vals-error, mean_vals+error,color=axes.lines[-1].get_color(), alpha=.33, lw =0)

axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0,t_final,num=7))
axes.set_ylim([0,100])
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')
fig.savefig('fig_method3_ensemble.pdf',bbox_inches='tight')
