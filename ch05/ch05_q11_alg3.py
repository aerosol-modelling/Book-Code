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
N_bin = len(M)
V = (N_part / 100) * 2e-10
print('Initial number concentration %e' %(N_part/V))
max_size = N_part 
index = 0
M = 1.0*N_part
bin_size = np.arange(1,101)
nt = 60
times = np.linspace(0,t_final,nt+1)
dt = t_final / nt
hist = np.zeros(len(times))
hist[0] = M
print(times)
index = 1
t1 = time.time()
while t < t_final and np.sum(M) > 1:
    K12 = decay
    lamb = M*K12/V
    N_event = np.random.poisson(lamb*dt)
    M -= N_event
    t += dt
    hist[index] = M
    index += 1

t2 = time.time()
print(t2-t1)
fig = plt.figure()
axes = fig.add_subplot(1,1,1)
max_size_plot = 10
axes.plot(times, hist) #, label='N$_{%i}$' %(i+1))
axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0,t_final,num=7))
axes.set_ylim([0,100])
axes.set_xlabel('Time (s)')
axes.set_ylabel('Number of particles $N_k$')
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')
fig.savefig('ch05_q11_alg3.pdf',bbox_inches='tight')

np.savetxt('alg3.txt', [times,hist[0:len(times)]])
