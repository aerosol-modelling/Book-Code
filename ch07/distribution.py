import numpy as np

from parameters import *

# distribution.py calculates a normal and
# log-normal # distribution for:
# mu - geometric mean diameter of the mode 
# sigma - deviation
# N_l, total number concentration of a mode 

def normal_distribution(mu, sigma, N_T):
    # Define lower and upper limit of the size distribution [m]
    d_1  = 0.
    d_Nb = mu*2. 
    
    # diameters of each size bin [m]
    d = np.linspace(d_1, d_Nb, num=Nb)
    
    # width of the size bin [m}
    delta_d=d[1]-d[0]
    
    # number concentration in each size class
    n = N_T*delta_d/(np.sqrt(2.0*np.pi)*sigma)*np.exp(-(d-mu)**2/(2.0*sigma)**2)
    
    return n, d

def log_normal_distribution(mu, sigma_g, N_T):

    # Define lower and upper limit of the size distribution [m]
    d_1 = 10.0e-9
    d_Nb = 1.0e-6
    
    # Volume ratio between bins
    v_rat = (d_Nb/d_1)**(3.0/(Nb-1.0))
    # Diameter ratio between bins
    d_rat = v_rat**(1.0/3.0)

    # diameters of each size bin [m]
    d=np.zeros((Nb), dtype=float)
    d[0] = d_1

    for step in range(Nb):
        if step > 0:
            d[step]=d[step-1]*d_rat
            #d = np.exp(np.linspace(np.log(d_1), np.log(d_Nb),
            #                       num=Nb))
            
    # volumes of each size bin [m3]
    v = 4.0 / 3.0 * np.pi *(d / 2.0)**3

    # lower limit of each size bin [m]
    v_lo = 2 * v / (1.0 + v_rat)
    
    # upper limit of each size bin [m]
    v_hi = v_rat * v_lo
    
    # width of the size bin [m}
    delta_d=d*2**(1.0/3.0)*(v_rat**(1.0/3.0) - 1)/\
        (1 + v_rat)**(1.0/3.0)
    
    # log normal distribution p(x)
    n = N_T*delta_d/(d*np.sqrt(2.0*np.pi)*np.log(sigma_g))*\
        np.exp(-(np.log(d)-np.log(mu))**2/(2.0*np.log(sigma_g)**2))

    return n, d, v, v_lo, v_hi
