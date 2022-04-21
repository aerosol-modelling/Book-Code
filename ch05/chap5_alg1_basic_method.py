import numpy as np
import kernel_helper
import matplotlib.pyplot as plt
import time

(mb, v_factor, t_final) = kernel_helper.GetParameters()

t = 0.0
N_part = 100 
Q = np.ones(N_part, dtype='int')
V = (N_part / 100) * v_factor
print('Initial number concentration %e' %(N_part/V))

max_size_track = N_part # Max size that we are tracking in histogram
hist = np.zeros((N_part, max_size_track), dtype=int)
hist[0,0] = N_part
i_time = 1 
times = [0.0]
t1 = time.time()
Q_hist = []
Q_hist.append(Q.copy())
pairs = np.zeros((N_part,2), dtype=int)
while t < t_final and len(Q) > 1:
    N_particles = len(Q) # Number of remaining particles
#    tau = np.full((N_particles,N_particles),np.Inf, dtype=float)
    tau = np.zeros((N_particles,N_particles))
    for i in range(N_particles):
        for j in range(i):
            K12 = kernel_helper.GetKernel(Q[i]*mb, Q[j]*mb)
            lamb = K12/V
            u = np.random.rand()
            tau[i,j] = -np.log(u)/lamb
    # Find the index of the min tau
    index_array = np.argmin(tau)
    (p1,p2) = np.unravel_index(index_array, shape=tau.shape) 
    pairs[i_time,0] = Q[p1]
    pairs[i_time,1] = Q[p2] 
    Q = np.append(Q,Q[p1]+Q[p2]) # Append new particle
    Q = np.delete(Q,p1) # Remove particle 1
    Q = np.delete(Q,p2) # Remove particle 2
    dt = tau[p1,p2]
    print(dt)
    t += dt
    times.append(t)
    # Save result by histogramming the list of particles
    hist[i_time,:], edges = np.histogram(Q, bins=max_size_track,
                                        range=(.5,max_size_track+.5))
    i_time += 1
    Q_hist.append(Q.copy())
t2 = time.time()
print(N_part,t2-t1)
# Plot results for different monomer sizes
fig = plt.figure()
axes = fig.add_subplot(1,1,1)
max_size_plot = 10
for i in range(max_size_plot):
    axes.plot(times, hist[0:len(times),i], label='$m_{%i}$' %(i+1))
axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0, t_final, num=7))
axes.set_ylim([0,100])
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')
fig.savefig('fig_method1.pdf', bbox_inches='tight')
# Make the log version
fig = plt.figure()
axes = fig.add_subplot(1,1,1)
#hist = np.ma.masked_less(hist,1)
for i in range(max_size_plot):
    axes.step(times, hist[0:len(times),i], where='post', label='$m_{%i}$' %(i+1))

values = plt.rcParams['axes.prop_cycle'].by_key()['color']
ms = 8 
lw = 1 
for i in range(1,len(times)):
    p1 = pairs[i,0] - 1
    p2 = pairs[i,1] - 1
    if (p1 < max_size_plot):
        axes.plot(times[i],hist[i,p1], 'X',
                  color=values[np.mod(p1, len(values))],ms=ms,lw=lw)
    if (p2 < max_size_plot):
        axes.plot(times[i],hist[i,p2], 'X',
                  color=values[np.mod(p2,len(values))],ms=ms,lw=lw)
    if (p1 + p2 < max_size_plot):
        axes.plot(times[i],hist[i,p1+p2+1], 'P',
                  color=values[np.mod(p1+p2+1,len(values))], ms=ms,lw=lw)
axes.set_ylim([.5,100])
axes.set_yscale('log')
axes.set_xlim([0,t_final])
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.set_xticks(np.linspace(0, t_final, num=7))
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')
fig.savefig('fig_method1_log.pdf', bbox_inches='tight')

np.savetxt('method1_results.txt', hist)
np.savetxt('method1_times.txt', times)

fig = plt.figure(figsize=(18/2.54,4.5))
axes = fig.add_axes([.05,.125,.9,.85])
for i_time in range(len(times)):
    parts = len(Q_hist[i_time])
    sizes = np.array(Q_hist[i_time])**2
    axes.scatter(np.ones(parts)*times[i_time], np.linspace(1,parts,parts),sizes,
                 clip_on=False)
axes.set_xlabel('Time (s)')
axes.set_ylabel("Particles")
axes.set_xlim([0,3600])
axes.set_ylim([.5,N_part+.5])
axes.set_xticks([0,600,1200,1800,2400,3000,3600])
axes.spines['right'].set_visible(False)
axes.spines['left'].set_visible(False)
axes.spines['top'].set_visible(False)
axes.get_yaxis().set_ticks([])
fig.savefig('fig_method1_particles.pdf')
