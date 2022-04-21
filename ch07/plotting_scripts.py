import matplotlib.pyplot as plt
import numpy as np
from auxiliary_functions import initialize_variables, Kelvin_effect
from parameters import *
from sectional_methods import quasi_stationary

def plot_dN_dlogDp(y, v, v_lo, v_hi, sectional_method, linestyle):

    _, R_wet, _, N_m, _, _, _, _, _, _, _, _ = initialize_variables(y)

    # initialize dlogDp
    dlogDp=np.zeros(Nb)
    
    if sectional_method == 'full-moving':

        # for full-moving system, we assume that dlogDp is
        # equal to the difference of logarithms of the
        # mean diameter of adjacent bins
        dlogDp[1:]=np.log(2.0*R_wet[1:])-np.log(2.0*R_wet[0:-1])
        dlogDp[0]=dlogDp[1]
        
    else:

        # find empty bins 
        empty_bins = np.where(N_m < N_low_limit)
        # set mean size of empty bins the the center of the bin
        R_wet[empty_bins] = (v[empty_bins]*3.0/(4.0*np.pi))**(1./3.)
        # calculate dlogDp
        dlogDp=np.log(v_hi**(1./3.)/v_lo**(1./3.))

    # plot the size distribution
    plt.plot(R_wet*2.0e6, N_m/dlogDp,linestyle)

    # set x-axis to log scale
    plt.xscale('log')

    # set the limits for the x-axis
    plt.xlim(5.0e-2, 1.0)

    # add labels for x and y
    plt.ylabel('$dN/d\mathrm{log}D_\mathrm{p}(\mathrm{m}^{-3})$')
    plt.xlabel('$D_\mathrm{p}(\mu \mathrm{m})$')

    # add the legend
    plt.legend(remapping_scheme)

def plot_N_R_t(y, d, n, v, sectional_method, linestyle):

    R_dry_t = d / 2.
    R_wet_t = d / 2.
    N_m_t = n
    
    for i in range(integration_time):

        if sectional_method == 'full-moving':
            # remap particles to a fixed grid
            # using the quasi-stationary method
            # (this is only for plotting purposes)
            array = quasi_stationary(v, y[i,:])

        else:

            array = y[i,:]

        # retrieve the wet radii and
        # number concentrations from the ODE array
        _, R_wet, _, N_m, _, _, _, _, _, _, _, _ = initialize_variables(array)

        # combine the time steps to one array
        R_wet_t = np.vstack([R_wet_t, R_wet])
        N_m_t = np.vstack([N_m_t, N_m])

    # plot the size distribution as a function of time
    x=range(integration_time)
    y=np.log10(d)
    z=np.maximum(N_m_t[1:,:], N_low_limit)
    yy, xx = np.meshgrid(y,x)
    plt.pcolormesh(xx,yy, np.log(z), shading='auto')

    # add labels to x and y
    plt.xlabel('$t(\mathrm{s})$')
    plt.ylabel('$\mathrm{log_{10}}(R_\mathrm{wet})(\mathrm{m})$')

    # set limits to y-axis
    plt.ylim(-8,-6)

    plt.show()
    plt.close()
    
def plot_R_t(y, d, n, sectional_method, linestyle):

    R_dry_t = d / 2.
    R_wet_t = d / 2.
    N_m_t = n
    
    for i in range(integration_time):

        array = y[i,:]

        # retrieve the wet radii and
        # number concentrations from the ODE array
        _, R_wet, _, N_m, _, _, _, _, _, _, _, _ = initialize_variables(array)

        # combine the time steps to one array
        R_wet_t = np.vstack([R_wet_t, R_wet])
        N_m_t = np.vstack([N_m_t, N_m])

    # plot the size distribution as a function of time
    for i in range(Nb):

        plt.plot(range(integration_time), R_wet_t[1:,i],'-k')    

    # set y-axis to log scale
    plt.yscale('log')

    # set limits to y-axis
    #plt.ylim(10.0e-9,1.0e-6)

    # add labels to x and y
    plt.xlabel('$t(\mathrm{s})$')
    plt.ylabel('$R(\mathrm{m})$')
    
    plt.show()
    plt.close()

def plot_variables(y, smallest_activated):

    S_w_t = np.zeros(Nb)
    S_t = np.zeros(integration_time)
    T_t = np.zeros(integration_time)
    
    for i in range(integration_time):

        R_dry, R_wet, T, N_m, mol_tot, k_i_m, C_g, X, dydt, Vs, Vw, kappa\
            = initialize_variables(y[i][:])

        # calculate the saturation ratio S_w of water
        # over a liquid droplet:
        # Kelvin effect of water on the droplet surface
        K = Kelvin_effect(R_wet, T, 'H2O')
        # saturation ratio of water on the droplet surface
        S_w = (R_wet**3 - R_dry**3)/\
            np.maximum((R_wet**3-R_dry**3*(1.0-kappa)), 1.e-30) * K

        # calculate the saturation ratio S of water
        # in gas phase:
        # find the index of the concentration of the
        # gas phase water in the array y
        gas = Nb+Nb*len(index)+index['H2O']-1
        # gas phase concentration of water
        C_h2o = y[i][gas]
        # saturation ratio of water in the gas phase
        S_t[i] = C_h2o/C_star['H2O'](T)
        # Air temperature
        T_t[i] = T
        
        # combine the time steps for S_w to one array
        S_w_t = np.vstack([S_w_t, S_w])

    # plot the supersaturation ratio of water in air
    plt.plot(range(integration_time), S_t,'-k')
    # plot the smallest activated bin
    plt.plot(range(integration_time), S_w_t[1:,smallest_activated],'--r')
    # plot the largest non-activated bin
    plt.plot(range(integration_time), S_w_t[1:,smallest_activated-1],'--b')
    # plot the smallest bin
    plt.plot(range(integration_time), S_w_t[1:,0],'--m')
    # plot the largest bin
    plt.plot(range(integration_time), S_w_t[1:,Nb-1],'--g')
    plt.xlim(1000,2500)
    plt.ylim(0.985,1.01)
    plt.legend(['$S$',
                '$S_\mathrm{w, smallest\,activated}$',
                '$S_\mathrm{w, largest\,non-activated}$',
                '$S_\mathrm{w, smallest\,bin}$',
                '$S_\mathrm{w, largest\,bin}$'])                
    # add labels to x and y
    plt.ylabel('$S$')
    plt.xlabel('$t$(s)')
 
    plt.show()
    plt.close()

    # plot the air temperature
    plt.plot(range(integration_time), T_t,'-k')
    
    # add labels to x and y
    plt.ylabel('$T$(K)')
    plt.xlabel('$t$(s)')

    plt.show()
    plt.close()

