import numpy as np
import matplotlib.pyplot as plt
import time

#np.random.seed(seed=10)

decay = 1e-13

t = 0.0
t_final = 3600.0 
N_part =  100
M = [1] * N_part
N = len(M)
V = (N_part / 100) * 2e-10 
hist = np.zeros(N_part, dtype=int)
hist[0] = N_part
i_time = 1 
times = [0.0]
t1 = time.time()
while t < t_final and len(M) > 1:
    N = len(M) # Number of remaining particles
    tau = np.full((N),np.Inf, dtype=float)
    for i in range(N):
        lamb = decay/V
        u = np.random.rand()
        tau[i] = -np.log(u)/lamb
    p1 = np.argmin(tau)
    M.pop(p1) # Remove particle 1
    dt = tau[p1]
    t += dt
    times.append(t)
    hist[i_time] = len(M) 
    i_time += 1
t2 = time.time()
print(N_part,t2-t1)

fig = plt.figure()
axes = fig.add_subplot(1,1,1)
axes.plot(times, hist[0:len(times)])
axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0, t_final, num=7))
axes.set_ylim([0,100])
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
#axes.legend(bbox_to_anchor=(1, 1), loc='upper right')
fig.savefig('ch05_q11_alg1.pdf', bbox_inches='tight')

np.savetxt('alg1.txt', [times,hist[0:len(times)]])
