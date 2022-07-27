import numpy as np
import kernel_helper
import matplotlib.pyplot as plt
import time
import scipy.stats


(mb, v_factor, t_final) = kernel_helper.GetParameters()
v_factor = 1.0
t_final = 3600.0

def main(dt):
    t = 0.0
    N_part = int(1e12)
    N_bin = 2
    N = np.zeros(N_bin)
    V = 1.0
    print('Initial number concentration %e' %(N_part/V))
    index = 0
    N[0] = 1.0*N_part
    bin_size = np.arange(1,101)
    nt = int(t_final / dt)
    times = np.linspace(0,t_final,nt+1)
    dt = t_final / nt
    hist = np.zeros((len(times),1))
    hist[0,:] = N[0]
    index = 1
    t1 = time.time()
    K12 = kernel_helper.GetKernel(mb,mb)
    ind_k = 0
    k = 1
    ind_l = 0
    l = 1
    dest = min(k + l, N_bin)
    while t < t_final and np.sum(N) > 1:
        N_pairs = (N[ind_k]*(N[ind_l]-1)) / 2
        lamb = N_pairs*K12/V
        N_event = np.random.poisson(lamb*dt)
        N[ind_k] -= N_event
        N[ind_l] -= N_event
        N[dest-1] += N_event
        t += dt
        hist[index,:] = N[0] 
        index += 1

    t2 = time.time()
    print(t2-t1)
    return(times, hist)

fig = plt.figure(figsize=(4.5,3.5))
axes = fig.add_subplot(1,1,1)
max_size_plot = 1 
dts = [1,10,100]
n_runs = 1
for i_dt in (dts):
    n_times = int(t_final / i_dt) + 1
    number_conc = np.zeros((n_runs,n_times))
    for i_run in range(n_runs):
       times, hist = main(i_dt)
       number_conc[i_run,:] = hist[:,0]
#    for i_run in range(n_runs):
#       axes.plot(times, number_conc[i_run,:])
    axes.plot(times, np.mean(number_conc,axis=0), label='$\Delta t =$ %i' %(i_dt))
axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0,t_final,num=7))
axes.set_yscale('log')
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')

print(axes.get_position())
fig.savefig('ch05_q10.pdf',bbox_inches='tight')
