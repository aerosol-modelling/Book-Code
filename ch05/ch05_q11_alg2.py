import numpy as np
import matplotlib.pyplot as plt
import time

decay = 1e-13

d0 = 1e-8 # Monomer diameter
M0 = (1/6)*np.pi*d0**3 # Monomer volume

t = 0.0
t_final = 3600.0
N_part = 100
M = np.zeros(N_part)
N = len(M)
V = (N_part / 100) * 2e-10
print('Initial number concentration %e' %(N_part/V))
max_size = N_part 
hist = np.zeros(N_part,dtype=int)
M = 1.0*N_part
bin_size = np.arange(1,100)
x = np.arange(1,N+1)
times = [0.0]
hist[0] = M
index = 1 
t1 = time.time()
while t < t_final and np.sum(M) > 1:
    tau = np.full((N),np.Inf, dtype=float)
    K12 = decay
    N_pairs = M
    lamb = N_pairs*K12/V
    u = np.random.rand()
    if (lamb > 0):
        tau = -np.log(u)/lamb
    # Find the index of the min tau
    index_array = np.argmin(tau)
    p1 = index_array #np.unravel_index(index_array,shape=tau.shape) 
    dt = tau
    M -= 1
    t += dt
    times.append(t)
    hist[index] = M
    index += 1

t2 = time.time()
print(t2-t1)

fig = plt.figure()
axes = fig.add_subplot(1,1,1)
axes.plot(times, hist[0:len(times)]) #, label='N$_{%i}$' %(i+1))
axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0,t_final,num=7))
axes.set_ylim([0,100])
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')
fig.savefig('ch05_q11_alg2.pdf',bbox_inches='tight')

np.savetxt('alg2.txt', [times,hist[0:len(times)]])
