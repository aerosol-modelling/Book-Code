import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

N = 100000 # Number of particles
tau = np.zeros(N, dtype=float)
r = 6.0
vals = np.zeros(100)
for i_repeat in range(100):
    u = np.random.rand(N)
    tau = -np.log(u)/r
    vals[i_repeat] = np.mean(tau)

print(np.mean(vals),np.std(vals),np.mean(vals)-np.std(vals),np.mean(vals)+np.std(vals))
print(vals)

bin_count, bin_edges, binnumber = scipy.stats.binned_statistic(tau,tau,statistic='count', bins=100)
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2
fig = plt.figure()
axes = fig.add_subplot(1,1,1)
axes.plot(bin_centers, bin_count,label='samples')
print(np.mean(tau))
axes.plot([np.mean(tau),np.mean(tau)],[0,np.max(bin_count)],'k--',label='mean')
axes.plot(bin_centers,bin_count[0]*np.exp(-6*bin_centers),'k',label='analytical')
axes.set_xlabel(r'Time $\tau$ (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.set_yscale('log')
axes.set_xlim([0,3]) #bin_edges[-1]])
axes.legend()
fig.savefig('ch05_q05.pdf', bbox_inches='tight')
