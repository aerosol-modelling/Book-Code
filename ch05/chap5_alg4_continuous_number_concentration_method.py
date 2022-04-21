import numpy as np
import kernel_helper
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import time

(mb, v_factor, t_final) = kernel_helper.GetParameters()

n_init = 100.0 / v_factor 

def CoagulationRates(C):
    N_bin = len(C)
    r = np.zeros(N_bin)
    bin_masses = np.arange(1,N_bin+1)
    for ind_k,k in enumerate(bin_sizes):
        # Destruction
        for ind_l,l in enumerate(bin_masses):
            K = kernel_helper.GetKernel((k)*mb,(l)*mb)
            r[ind_k] -= C[ind_k]*C[ind_l]*K
        # Production
        for ind_l,l in enumerate(bin_masses[0:(ind_k)]):
            K = kernel_helper.GetKernel(l*mb,(k-l)*mb)
            r[ind_k] += .5*K*C[ind_l]*C[ind_k-ind_l-1]

    return(r)

max_size = 8 
fig = plt.figure(figsize=(4.5,3.5))
axes = fig.add_subplot(1,1,1)
pos = axes.get_position()
print(pos)
pos.y0 = .15
pos.y1 = pos.y1 + .04
axes.set_position(pos)
print(pos)
C = np.zeros(max_size)
C[0] = n_init 
dt =  10
nt = int(3600/ dt)
values = np.zeros((nt+1,max_size))
values[0,:] = C 
t1 = time.time()
for i in range(nt):
    rates = CoagulationRates(C)
    C += rates*dt
    values[i+1,:] = C 
t2 = time.time()
print(t2-t1)
for ii in range(8):
   axes.plot(np.linspace(0,3600,nt+1), values[:,ii],label='$C_{%i}$' %(ii+1))
axes.set_xlim([0,t_final])
axes.ticklabel_format(useOffset=False)
axes.set_xticks(np.linspace(0,t_final,num=7))
axes.set_xlabel('Time (s)')
axes.set_ylim([0,n_init])
axes.set_ylabel('Number concentration $C_k$ (# m$^{-3}$)')
axes.legend(bbox_to_anchor=(1, 1), loc='upper right')
#axes.yaxis.get_offset_text().set_size(0)
fig.savefig('fig_method4.pdf') #,bbox=[0,-.1,1,1])
