import numpy as np
import matplotlib.pyplot as plt
import time

decay = 1e-13 * 5e9
t_final = 3600.0
ts = np.linspace(0, t_final, 3600)
P0 =  100 
Nt = P0 * np.exp(-decay*ts)
fig = plt.figure()
axes = fig.add_subplot(1,1,1)
axes.plot(ts,Nt) 
axes.set_xlim([0,t_final])
axes.set_xticks(np.linspace(0,t_final,num=7))
axes.set_xlabel('Time (s)')
axes.set_ylim([0,P0])
axes.set_ylabel('Number concentration # m$^{-3}$')
axes.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
fig.savefig('ch05_q11_alg4.pdf',bbox_inches='tight')

np.savetxt('alg4.txt', [ts,Nt])
