import numpy as np
from parameters import *

def diffusion_coefficient(N_m, R_wet, T, S_w, species):

    # Mean thermal velocity of each molecule [m/s]
    mean_them_vel=np.sqrt((8.0*R_gas*T)/(np.pi*M[species]))

    # Mean free path for each molecule [m]
    gamma_gas = 3.0*D_g[species]/mean_them_vel

    # Knudsen number
    # Calculate the Knudsen number for all condensing molecules based on this new size \n')
    Kn = gamma_gas / np.maximum(2.0 * R_wet, 1.e-10)

    # Non-continuum regime correction
    # Calculate a correction factor according to the continuum versus non-continuum regimes
    Inverse_Kn=1.0/Kn
    Correction_part1=(1.33e0+0.71e0*Inverse_Kn)/(1.0e0+Inverse_Kn)
    Correction_part2=(4.0e0*(1.0e0-alpha_d[species]))/(3.0e0*alpha_d[species])
    Correction_part3=1.0e0+(Correction_part1+Correction_part2)*Kn
    Correction=1.0/Correction_part3
    
    D_eff = D_g[species] * Correction
    
    # Corrected thermal conductivity of moist air
    kappa_air = Correction * 0.023807 + 7.1128e-5 * (T - 273.15)
    
    if species == 'H2O':

        # Diffusion coefficient
        D_eff = D_eff \
            / (1.0 + S_w * C_star[species](T) *(D_eff * M[species]**2 * L_e[species](T)**2) / (kappa_air * R_gas * T))

    return D_eff

def Kelvin_effect(R_wet, T, species):

    Kelvin = np.exp((4.0*(M[species])*sigma)/(R_gas*T*rho[species]*np.maximum(R_wet,1.e-10)*2))

    return Kelvin

def calculate_kappa(v_si, v_s):

    kappa=np.zeros(Nb)
    
    for species in All_species:

         if species != 'H2O':

             # calculate the kappa in each size bin
             kappa=kappa+v_si[species]/np.maximum(v_s, 1.e-30)*kappa_i[species]

    return kappa

def initialize_variables(y):   

    # initialize values
    v_T=np.zeros(Nb)
    mol_tot=np.zeros(Nb)
    k_i_m=dict()
    C_g=dict()
    X=dict()
    dydt=np.zeros(len(y))
    v_si=dict()
    v_s=np.zeros(Nb)
    v_w=np.zeros(Nb)
    
    # retrieve values for temperature in y-vector
    # y[-1-1] is the second to last value in y
    T = y[-1-1]
    
    # retrieve values for number concentrations
    N_m = y[range(Nb)]

    # Sum up the volumes of all species to obtain
    # the total volume of the droplet
    for species in All_species:

        # find the indices of each species
        bins = range(Nb*index[species],Nb*index[species]+Nb)

        # volume of one soluble species in each bin
        v_si[species]=y[bins]*M[species]/rho[species]/\
            np.maximum(N_m, N_low_limit)

        # total volume of particles in each bin
        v_T+=v_si[species]

        if species == 'H2O':
            # volume of water in droplets
            v_w+=v_si[species]

        else:
            # volume os the solid part of droplets
            v_s+=v_si[species]
        
        # total number of moles in individual bins
        mol_tot+=y[bins]

    # Calculate the droplet radius
    R_wet = (3.0/(4.0*np.pi)*v_T)**(1./3.)

    # Calculate the dry radius
    R_dry = (3.0/(4.0*np.pi)*v_s)**(1./3.)

    # Calculate kappa
    kappa = calculate_kappa(v_si, v_s)

    return R_dry, R_wet, T, N_m, mol_tot, k_i_m, C_g, X, dydt, v_w, v_s, kappa

def find_activated(array):

    _, R_wet, _, _, _, _, _, _, _, _, _, _ = initialize_variables(array)

    # find the index for the smallest activated bin
    smallest_activated = np.argmax(R_wet[1:]/R_wet[:-1]) + 1

    return smallest_activated
