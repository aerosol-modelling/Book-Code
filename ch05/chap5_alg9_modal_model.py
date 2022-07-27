import numpy as np
from scipy.integrate import odeint
from scipy import constants
import matplotlib.pyplot as plt

# Define parameters
lambd_air = 0.0686*1.0e-6 # mean free path (seinfeld & pandas, page 602, Figure 13.5) 
b10 = 0.8
b13 = 0.9
kb = constants.k  # Boltzmann constant = 1.38064852 × 1e-23 m2 kg s-2 K-1
mu_SI = 1.83*1e-5 # viscosity
T = 25+273.15
rho_p = 1.0e3

Knc = (2*kb*T)/(3*mu_SI)

Kfm = np.sqrt((3*kb*T)/rho_p)

# define the initial parameters 
sigma_g1 = 10**0.161
sigma_g2 = 10**0.217

def get_R_0_nc(dg1,dg2,sigma_g1,sigma_g2):
    Kng1 = 2.0*lambd_air/dg1
    Kng2 = 2.0*lambd_air/dg2
    A1 = 1.392 * (Kng1)**0.0783
    A2 = 1.392 * (Kng2)**0.0783
    return(Knc *(2+A1*Kng1*(np.exp((4/8)*np.log(sigma_g1)**2)
           + (dg2/dg1)*np.exp((16/8)*np.log(sigma_g1)**2)
           * np.exp((4/8)*np.log(sigma_g2)**2))
           + A2*Kng2*(np.exp((4/8)*np.log(sigma_g2)**2)
           + (dg1/dg2)*np.exp((16/8)*np.log(sigma_g2)**2)*np.exp((4/8)*np.log(sigma_g1)**2))
           + (dg1/dg2+dg2/dg1)*(np.exp((4/8)*np.log(sigma_g2)**2))
           * (np.exp((4/8)*np.log(sigma_g1)**2))))

def get_R_3_nc(dg1,dg2,sigma_g1,sigma_g2):
    Kng1 = 2.0*lambd_air/dg1
    Kng2 = 2.0*lambd_air/dg2
    A1 = 1.392 * (Kng1)**0.0783
    A2 = 1.392 * (Kng2)**0.0783
    return(Knc*(dg1**3)*(2*np.exp((36/8)*np.log(sigma_g1)**2)*A1*Kng1
           * (np.exp((16/8)*np.log(sigma_g1)**2)
           + (dg2/dg1)*np.exp((4/8)*np.log(sigma_g1)**2)*np.exp((4/8)*np.log(sigma_g2)**2))
           + A2*Kng2*(np.exp((36/8)*np.log(sigma_g1)**2)*np.exp((4/8)*np.log(sigma_g2)**2)\
           + (dg1/dg2)*np.exp((64/8)*np.log(sigma_g1)**2)*np.exp((16/8)*np.log(sigma_g2)**2))\
           + (dg2/dg1)*np.exp((16/8)*np.log(sigma_g1)**2)*np.exp((4/8)*np.log(sigma_g2)**2)\
           + (dg1/dg2)*np.exp((64/8)*np.log(sigma_g1)**2)*np.exp((4/8)*np.log(sigma_g2)**2)))

def get_R_0_fm(dg1,dg2,sigma_g1,sigma_g2):
    Kng1 = 2.0*lambd_air/dg1
    Kng2 = 2.0*lambd_air/dg2
    A1 = 1.392 * (Kng1)**0.0783
    A2 = 1.392 * (Kng2)**0.0783

    return(Kfm*b10*np.sqrt(dg1)*(np.exp((1/8)*np.log(sigma_g1)**2)
           + np.sqrt(dg2/dg1)*np.exp((1/8)*np.log(sigma_g2)**2)
           + 2*(dg2/dg1)*np.exp((1/8)*np.log(sigma_g1)**2)
           * np.exp((4/8)*np.log(sigma_g2)**2)
           + (dg2**2/dg1**2)*np.exp((9/8)*np.log(sigma_g1)**2)
           * np.exp((16/8)*np.log(sigma_g2)**2)
           + (np.sqrt(dg1/dg2)**3)*np.exp((16/8)*np.log(sigma_g1)**2)
           * np.exp((9/8)*np.log(sigma_g2)**2)
           + 2*np.sqrt(dg1/dg2)*np.exp((4/8)*np.log(sigma_g1)**2)
           * np.exp((1/8)*np.log(sigma_g2)**2)))

def get_R_3_fm(dg1,dg2,sigma_g1,sigma_g2):
    Kng1 = 2.0*lambd_air/dg1
    Kng2 = 2.0*lambd_air/dg2
    A1 = 1.392 * (Kng1)**0.0783
    A2 = 1.392 * (Kng2)**0.0783
    return (Kfm*b13*dg1**(7/2)*(np.exp((49/8)*np.log(sigma_g1)**2)
            + np.sqrt(dg2/dg1)*np.exp((36/8)*np.log(sigma_g1)**2)
            * np.exp((1/8)*np.log(sigma_g2)**2)
            + 2*(dg2/dg1)*np.exp((25/8)*np.log(sigma_g1)**2)
            * np.exp((4/8)*np.log(sigma_g2)**2)
            + (dg2**2/dg1**2)*np.exp((9/8)*np.log(sigma_g1)**2)
            * np.exp((16/8)*np.log(sigma_g2)**2)
            + (np.sqrt(dg1/dg2)**3)*np.exp((100/8)*np.log(sigma_g1)**2)
            * np.exp((9/8)*np.log(sigma_g2)**2)
            + 2*np.sqrt(dg1/dg2)*np.exp((64/8)*np.log(sigma_g1)**2)
            * np.exp((1/8)*np.log(sigma_g2)**2)))

def get_R(dg1,dg2,sigma_g1,sigma_g2):
    Kng1 = 2.0*lambd_air/dg1
    Kng2 = 2.0*lambd_air/dg2
    A1 = 1.392 * (Kng1)**0.0783
    A2 = 1.392 * (Kng2)**0.0783
    
    R_0_ij_nc = get_R_0_nc(dg1,dg2,sigma_g1,sigma_g2)
    R_3_ij_nc = get_R_3_nc(dg1,dg2,sigma_g1,sigma_g2)
    R_0_ij_fm = get_R_0_fm(dg1,dg2,sigma_g1,sigma_g2)
    R_3_ij_fm = get_R_3_fm(dg1,dg2,sigma_g1,sigma_g2)
    
    R_0_ii_nc = get_R_0_nc(dg1,dg1,sigma_g1,sigma_g1)
    R_0_jj_nc = get_R_0_nc(dg2,dg2,sigma_g2,sigma_g2)
    R_0_ii_fm = get_R_0_fm(dg1,dg1,sigma_g1,sigma_g1)
    R_0_jj_fm = get_R_0_fm(dg2,dg2,sigma_g2,sigma_g2)
    
    return(R_0_ij_nc, R_3_ij_nc, R_0_ij_fm, R_3_ij_fm, R_0_ii_nc, R_0_jj_nc,
           R_0_ii_fm, R_0_jj_fm)

def model(z,t):
    # z = [M_0i, M_0j, M_3i, M_3j] 
    dg1 = (z[2]/(z[0]*np.exp((9/2)*np.log(sigma_g1)**2)))**(1.0/3)
    dg2 = (z[3]/(z[1]*np.exp((9/2)*np.log(sigma_g2)**2)))**(1.0/3)
    
    R_0_ij_nc, R_3_ij_nc, R_0_ij_fm, R_3_ij_fm, R_0_ii_nc, R_0_jj_nc, R_0_ii_fm, R_0_jj_fm = \
         get_R(dg1,dg2,sigma_g1,sigma_g2)
    
    Ca_0_ij_nc = z[0]*z[1]*R_0_ij_nc
    Ca_3_ij_nc = z[0]*z[1]*R_3_ij_nc
    Ca_0_ij_fm = z[0]*z[1]*R_0_ij_fm
    Ca_3_ij_fm = z[0]*z[1]*R_3_ij_fm
    
    Ca_0_ii_nc = z[0]*z[0]*R_0_ii_nc
    Ca_0_jj_nc = z[1]*z[1]*R_0_jj_nc
    
    Ca_0_ii_fm = z[0]*z[0]*R_0_ii_fm
    Ca_0_jj_fm = z[1]*z[1]*R_0_jj_fm
 
    Ca_0_ii = Ca_0_ii_nc * Ca_0_ii_fm / (Ca_0_ii_nc + Ca_0_ii_fm)
    Ca_0_jj = Ca_0_jj_nc * Ca_0_jj_fm / (Ca_0_jj_nc + Ca_0_jj_fm)
    Ca_0_ij = Ca_0_ij_nc * Ca_0_ij_fm / (Ca_0_ij_nc + Ca_0_ij_fm)
    Ca_3_ij = Ca_3_ij_nc * Ca_3_ij_fm / (Ca_3_ij_nc + Ca_3_ij_fm)   

    dM0idt = -Ca_0_ii - Ca_0_ij
    dM0jdt = -Ca_0_jj
    dM3idt = -Ca_3_ij
    dM3jdt = Ca_3_ij
    
    return ([dM0idt, dM0jdt, dM3idt, dM3jdt])


# ## solve the model, and plot the results
Ni = 3.2e9
Nj = 2.9e9
dg1 = 2.0e-8
dg2 = 1.16e-7

M0i_init = Ni
M0j_init = Nj
M3i_init = (dg1**3)*(Ni*np.exp((9/2)*np.log(sigma_g1)**2))
M3j_init = (dg2**3)*(Nj*np.exp((9/2)*np.log(sigma_g2)**2))
z0 = [M0i_init, M0j_init, M3i_init, M3j_init]

# number of time points
n = 60*12
t_final = 12*3600
t = np.linspace(0, t_final+1, n)

# store solution
M0i = np.empty_like(t)
M0j = np.empty_like(t)
M3i = np.empty_like(t)
M3j = np.empty_like(t)

# record initial conditions
M0i[0] = z0[0]
M0j[0] = z0[1]
M3i[0] = z0[2]
M3j[0] = z0[3]

for i in (range(1,n)):
    tspan = [t[i-1],t[i]]
    z = odeint(model, z0, tspan)
    M0i[i] = z[1][0]
    M0j[i] = z[1][1]
    M3i[i] = z[1][2]
    M3j[i] = z[1][3]
    z0 = z[1]
    
fig = plt.figure(figsize=(12,3))
# M_0,i and M_0,j
ax = fig.add_subplot(1,2,1)
ax.plot(t,M0i,'b-',label=r'$M_\mathrm{0,i}$')
ax.plot(t,M0j,'r--',label=r'$M_\mathrm{0,j}$')
ax.set_title("number conc.")
ax.legend(loc='best')
# M_3,i and M_3,j
ax = fig.add_subplot(1,2,2)
ax.plot(t,M3i,'b-',label=r'$M_\mathrm{3,i}$')
ax.plot(t,M3j,'r--',label=r'$M_\mathrm{3,j}$')
ax.set_yscale('log')
ax.set_title("mass conc.")
ax.legend(loc='best')
plt.savefig('modal_result.pdf')

np.savetxt('file.txt',[M0i,M3i,M0j,M3j])

