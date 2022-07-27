import numpy as np
import kernel_helper
import matplotlib.pyplot as plt
import time


(mb, v_factor, t_final) = kernel_helper.GetParameters()
memo = {}

def kernel_memo(p1, p2):
  global memo
  if (p1, p2) in memo:
      return memo[(p1, p2)]
  else:
      memo[(p1,p2)] = kernel_helper.GetKernel(p1*mb, p2*mb)
  return memo[(p1, p2)]

def main(save_tau=True):

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
    ###### BEGIN MODIFICATION ######
    first_time = True
    while t < t_final and np.sum(N) > 1:
        if(first_time):
            tau = np.full((N_bin,N_bin),np.Inf, dtype=float)
            for ind_k, k in enumerate(bin_size):
                for ind_l, l in enumerate(bin_size[0:ind_k+1]):
                    K12 = kernel_memo(k,l) #Q[i], Q[j])
                    #K12 = kernel_helper.GetKernel(k*mb,l*mb)
                    if (k == l):
                        N_pairs = (N[ind_k]*(N[ind_l]-1)) / 2
                    else:
                        N_pairs = N[ind_k]*N[ind_l]
                    lamb = N_pairs*K12/V
                    u = np.random.rand()
                    if (lamb > 0):
                        tau[ind_k,ind_l] = -np.log(u)/lamb
            if (save_tau):
                first_time = False 
        else:
            # Loop over all entries in a few rows
            for ind_k in ([p1,p2,p1+p2+1]):
                k = bin_size[ind_k]
                for ind_l, l in enumerate(bin_size[0:ind_k+1]):
                   K12 = kernel_memo(k,l)
#                   K12 = kernel_helper.GetKernel(k*mb,l*mb)
                   if (k == l):
                       N_pairs = (N[ind_k]*(N[ind_l]-1)) / 2
                   else:
                       N_pairs = N[ind_k]*N[ind_l]
                   lamb = N_pairs*K12/V
                   u = np.random.rand()
                   if (lamb > 0):
                       tau[ind_k,ind_l] = -np.log(u)/lamb
            # Loop over all rows to change a few entries
            for ind_k, k in enumerate(bin_size):
               for ind_l in ([p1,p2,p1+p2+1]):
                 l = bin_size[ind_l]
                 if (l <= k): # Need to be careful on where to sample
                   K12 = kernel_memo(k,l)
#                   K12 = kernel_helper.GetKernel(k*mb,l*mb)
                   if (k == l):
                       N_pairs = (N[ind_k]*(N[ind_l]-1)) / 2
                   else:
                       N_pairs = N[ind_k]*N[ind_l]
                   lamb = N_pairs*K12/V
                   u = np.random.rand()
                   if (lamb > 0):
                       tau[ind_k,ind_l] = -np.log(u)/lamb

    ###### END MODIFICATION #######
        # Find the index of the min tau
        index_array = np.argmin(tau) # argmin returns index into flattened array
        (p1,p2) = np.unravel_index(index_array,shape=tau.shape) # convert flat index to index tuple
        dt = tau[p1,p2]
        N[p1] -= 1
        N[p2] -= 1
        N[p1+p2+1] += 1
        # These are the tau to reset?
        tau[p1,:] = np.inf
        tau[p2,:] = np.inf
        tau[p1+p2+1,:] = np.inf
        tau[:,p1] = np.inf
        tau[:,p2] = np.inf
        tau[:,p1+p2+1] = np.inf
        t += dt
        times.append(t)
        hist[index,:] = N
        index += 1
    t2 = time.time()
    print(t2-t1)

    return(times, hist[0:len(times),:])

### Do an ensemble

from scipy import interpolate
import scipy.stats

n_interp = 360
x = np.linspace(0,t_final,n_interp)
fig = plt.figure(figsize=(4.5,3.5))
axes = fig.add_subplot(1,1,1)
n_runs = 100
max_size = 3 
values = np.zeros((n_interp,n_runs,max_size))
for i_run in range(n_runs):
    times, hist = main()
    for i_size in range(max_size):
        f = interpolate.interp1d(times, hist[:,i_size],fill_value="extrapolate")
        values[:,i_run,i_size] = f(x)

for i_size in range(max_size):
    mean_vals = np.mean(values[:,:,i_size],axis=1)
    std_val = np.std(values[:,:,i_size],axis=1)
    error = scipy.stats.t.ppf(.975, n_runs)*std_val
    axes.plot(x,mean_vals, label='$N_{%i}$ optimized' %(i_size+1))
    axes.fill_between(x, mean_vals-error, mean_vals+error,color=axes.lines[-1].get_color(),
         alpha=.3,lw = 0)
values = np.zeros((n_interp,n_runs,max_size))
for i_run in range(n_runs):
    times, hist = main(save_tau=False)
    for i_size in range(max_size):
        f = interpolate.interp1d(times, hist[:,i_size],fill_value="extrapolate")
        values[:,i_run,i_size] = f(x)

for i_size in range(max_size):
    mean_vals = np.mean(values[:,:,i_size],axis=1)
    std_val = np.std(values[:,:,i_size],axis=1)
    error = scipy.stats.t.ppf(.975, n_runs)*std_val
    axes.plot(x,mean_vals, label='$N_{%i}$ (original)' %(i_size+1))
    axes.fill_between(x, mean_vals-error, mean_vals+error,color=axes.lines[-1].get_color(),
         alpha=.3,lw = 0)
axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0,t_final,num=7))
axes.set_ylim([0,100])
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')
fig.savefig('ch05_q09.pdf',bbox_inches='tight')
