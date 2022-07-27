import numpy as np
import kernel_helper
import matplotlib.pyplot as plt
import time

(mb, v_factor, t_final) = kernel_helper.GetParameters()

t = 0.0
N_part = 100 
Q = [1] * N_part
N = len(Q)
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

K12 = np.zeros((max_size_track,max_size_track))
for i in range(max_size_track):
    for j in range(i+1):
       K12[i,j] = kernel_helper.GetKernel((i+1)*mb, (j+1)*mb)

while t < t_final and len(Q) > 1:
    N_particles = len(Q) # Number of remaining particles
    tau = np.full((N_particles,N_particles),np.Inf, dtype=float)
    for i in range(N_particles):
        for j in range(i):
            r = K12[i,j]/V
            u = np.random.rand()
            tau[i,j] = -np.log(u)/r
    # Find the index of the min tau
    index_array = np.argmin(tau)
    (p1,p2) = np.unravel_index(index_array, shape=tau.shape) 
    pairs[i_time,0] = Q[p1]
    pairs[i_time,1] = Q[p2] 
    Q.append(Q[p1]+Q[p2]) # Append new particle
    Q.pop(p1) # Remove particle 1
    Q.pop(p2) # Remove particle 2
    dt = tau[p1,p2]
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
