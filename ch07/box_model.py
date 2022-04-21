#!/usr/bin/python3
# import physical constants and the properties
# of chemical species
from parameters import *

# import remapping schemes for the sectional scheme
from sectional_methods import *

# import the module for creating a log-normal
# distribution
from distribution import log_normal_distribution

# import module which defines ODEs
from differential_equations import dy_dt

# import matplotlib for plotting
import matplotlib.pyplot as plt

# import the ODE solver module from SciPy
from scipy.integrate import odeint

# import griddata module
from scipy.interpolate import griddata

# import plotting scripts
from plotting_scripts import plot_dN_dlogDp, plot_N_R_t,\
    plot_R_t, plot_variables

# import extra functions
from auxiliary_functions import find_activated

# initialize number size distribution
# using the module log_normal_distribution
n, d, v, v_lo, v_hi = log_normal_distribution(mu, sigma_g, N_T)

# initialize the vector aerosol composition size distribution
c = np.zeros((2+num_vbs_species,n.shape[0]))

# concentration of sulfate [mol/m3]
c[1,:] = n * np.pi/6.0*d**3 * rho['SO4'] / M['SO4']

# set the gas phase concentrations for dry aerosol cases
if model_setup <= 2:
    
    C_gases=[
        0.0,                                 # water
        0.0,                                 # sulfuric acid
        abundance['VBS0']/M['VBS0']*1.e-09,  # VBS species
        abundance['VBS1']/M['VBS1']*1.e-09,
        abundance['VBS2']/M['VBS2']*1.e-09,
        abundance['VBS3']/M['VBS3']*1.e-09,
        abundance['VBS4']/M['VBS4']*1.e-09,
        abundance['VBS5']/M['VBS5']*1.e-09,
        abundance['VBS6']/M['VBS6']*1.e-09,
        abundance['VBS7']/M['VBS7']*1.e-09,
        abundance['VBS8']/M['VBS8']*1.e-09,
        abundance['VBS9']/M['VBS9']*1.e-09,
    ]

# setup for cloud parcel cases
if model_setup >= 3:
    
    # calculate the volume of water [m3]
    # based on Kappa Köhler
    a_w = RH0 / 100.0
    v_w = a_w / (1 - a_w) * kappa_i['SO4'] * np.pi/6.0*d**3

    # convert to molar concentration [mol/m3]
    c[0,:] = n * v_w * rho['H2O'] / M['H2O']
    
    C_gases=[
        a_w * C_star['H2O'](T0),               # water
        0.0,                                   # sulfuric acid
        0., 0., 0., 0., 0., 0., 0., 0., 0.,0., # VBS species
    ] 

if model_setup == 4:
    C_gases[2:-1]=[
        abundance['VBS0']/M['VBS0']*0.1e-09,  # VBS species
        abundance['VBS1']/M['VBS1']*0.1e-09,
        abundance['VBS2']/M['VBS2']*0.1e-09,
        abundance['VBS3']/M['VBS3']*0.1e-09,
        abundance['VBS4']/M['VBS4']*0.1e-09,
        abundance['VBS5']/M['VBS5']*0.1e-09,
        abundance['VBS6']/M['VBS6']*0.1e-09,
        abundance['VBS7']/M['VBS7']*0.1e-09,
        abundance['VBS8']/M['VBS8']*0.1e-09,
        abundance['VBS9']/M['VBS9']*0.1e-09,
    ]

# initialize aerosol phase concentration of condensing organic species as zero
c[2:2+num_vbs_species,:] = 0.0

for j in remapping_scheme:

    print('simulating using the '+j+' method') 
    # Determine the initial values for ODEs
    array=n.tolist()                          # number concentration
    array.extend(c[0,:].tolist())             # water
    array.extend(c[1,:].tolist())             # sulfate
    for i in range(2,2+num_vbs_species):      # condensing organics
        array.extend(c[i,:].tolist())
    array.extend(C_gases)                     # gas phase species
    array.extend([T0, p0])                    # temperature and pressure

    # plot the initial size distribution
    #plot_dN_dlogDp(np.array(array),v, v_lo, v_hi, j, '-r')
    
    if j == 'full-moving':

        # set the boundary conditions for the time step
        t = np.linspace(1, integration_time, num=integration_time)
        
        solution = odeint(dy_dt, array, t,
                          rtol=rtol, atol=atol,
                          tcrit=None, mxstep=5000)

        # full-moving needs no remapping no remapping
        # the solution is the last line of array
        array = solution[-1,:]
        linestyle = 'k--'
            
    elif j != 'full-moving' and model_setup <= 2: 

        # initialize the solution array
        solution = np.zeros((integration_time,len(array)))
        # split the time interval for remapping schemes
        for k in range(integration_time):

            # set the boundary conditions for a 1 sec time step 
            t = [float(k), float(k+1)]

            # solve ODE's over one time step
            solution[range(k-1,k+1),:] = odeint(dy_dt, array, t,
                              rtol=rtol, atol=atol,
                              tcrit=None, mxstep=5000)

            # remap the size distribution
            if j == 'quasi-stationary':
                
                # call the quasi-stationary scheme
                array = quasi_stationary(v, solution[k,:])
                linestyle = 'g-'
                
            if j == 'moving-center':
                # call the moving center scheme
                array = moving_center(v_hi, solution[k,:])
                linestyle = 'b-'

            if j == 'hybrid-bin':
                # call the hybrid bin scheme
                array = hybrid_bin(v, v_lo, v_hi, solution[k,:])
                linestyle = 'r-'
                
            if j == 'flat-top':
                # call the hybrid bin flat top
                array = flat_top(v, v_lo, v_hi, solution[k,:])
                linestyle = '-m'

    else:

        print('Only full-moving scheme can used with the parcel setup')

    # print the elapsed time
    print('Elapsed time =',integration_time,
          '/',integration_time,'s',end='\r')

    # plot the final size distribution at the end of the simulation
    plot_dN_dlogDp(solution[-1][:],v, v_lo, v_hi, j, linestyle)
    # plot the size distribution as a function of time
    #plot_N_R_t(solution, d, n, v, j, linestyle)
    # plot radii as a function of time
    #plot_R_t(solution, d, n, j, linestyle)

    # find the index of the smallest activated bin
    # at the end of the simulation
    #smallest_activated=find_activated(solution[-1][:])
    # plot different parameters
    #plot_variables(solution, smallest_activated)
    
plt.show()
