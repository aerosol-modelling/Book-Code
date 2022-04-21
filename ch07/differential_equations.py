from parameters import *
import numpy as np
from auxiliary_functions import diffusion_coefficient, Kelvin_effect, initialize_variables

def dy_dt(y, t):

    y = np.maximum(y, 0.0)
    
    R_dry, R_wet, T, N_m, mol_tot, k_i_m, C_g, X, dydt, Vs, Vw, kappa\
        = initialize_variables(y)

    for species in Condensing_species:

        # Calculate the Kelvin effect
        K = Kelvin_effect(R_wet, T, species)

        # find the indices of each species in Array y
        # particle phase
        bins = range(Nb*index[species],Nb*index[species]+Nb)
        # gas phase
        gas = Nb+len(bins)*len(index)+index[species]-1

        # mole fraction of an individual species in a bin
        X[species] = y[bins]/np.maximum(mol_tot, 1.e-30)

        # gas phase concentration of a species
        C_g[species]=y[gas]

        # gas phase concentration above the droplet surface
        if species == 'H2O':
            # for water:
        
            # Saturation ratio of water
            S_w = (R_wet**3 - R_dry**3)/np.maximum((R_wet**3-R_dry**3*(1.0-kappa)), 1.e-30) * K
            # Equilibrium concentration on the surface
            C_surf = S_w * C_star[species](T)

        else:
            # for other species 
            C_surf = K * X[species]*C_star[species](T)
            
        # Calculate the diffusion coefficient for each individual species
        D_eff = diffusion_coefficient(N_m, R_wet, T, S_w, species)

        # condensation equation for individual particle species
        dydt[bins] = N_m * 4.0 * np.pi * R_wet * D_eff * (C_g[species]-C_surf)

        # condensation equation for gases
        dydt[gas] = - sum(dydt[bins])

        # temperature change due to
        # condensation /evaporation of water
        if species == 'H2O':

            # pressure
            p = y[-1]
            # moles of water in all size bins per m3 air
            C_water = sum(y[bins])
            # density of the air parcel
            rho_air = p * M_air / (R_gas * T)
            # mass density of the air parcel
            m_tot = rho_air + (C_water + C_g[species]) * \
                M[species]
            # heat capacity of the air parcel
            c_tot = (rho_air * c_p_air +                 \
                     (C_water + C_g[species]) *          \
                     M[species] * c_p_water)/m_tot
            # temperature change
            dydt[-2] = 1.0 / (rho_air * c_tot) *         \
                L_e[species](T) * (-dydt[gas])

        # Change in the number concentration and volume concentration
        # of the smallest size bin due to nucleation
        if model_setup == 2 and species == 'SO4':

            if t < 50.:
                # inject sulfuric acid in the gas phase for 50 seconds
                dydt[gas] = dydt[gas] + abundance['SO4']/M['SO4']*1.e-9
                
            # formation rate of new particles
            dydt[0] = dydt[0] + 1.0e20*C_g[species]

            # volume of the new particles
            v_NPF = np.pi / 6.0 * 10.0e-9**3

            # moles of sulfate in the new particles
            C_NPF = v_NPF * rho[species] / M[species]

            # change in molar concentration
            dydt[bins[0]] = dydt[bins[0]] + dydt[0] * C_NPF

            # change in gas phase concentration
            dydt[gas] = dydt[gas] - dydt[0] * C_NPF

    # Ordinary differential equations for an ascending air parcel
    if model_setup >= 3:
        # temperature change due to adiabatic
        # expansion of the air parcel
        dydt[-2] = dydt[-2] - g * w / c_tot

        # pressure change due to adiabatic
        # expansion of the air parcel
        dydt[-1] = -g * p * M_air * w / (R_gas * T)
        # term gamma accounting for adiabatic expansion
        gamma = dydt[-1] /p - dydt[-2] / T

        for species in Condensing_species:
            
            # gas phase
            gas = Nb+len(bins)*len(index)+index[species]-1
            # add the term for air parcel volume change
            # for gas phase concentrations
            dydt[gas] = dydt[gas] + y[gas] * gamma
            
        for species in All_species:
            
            # find the indices of each species in Array y
            # particle phase
            bins = range(Nb*index[species],Nb*index[species]+Nb)
            # add the term for air parcel volume change
            # for liquid phase concentrations
            dydt[bins] = dydt[bins] + y[bins] * gamma
            
        # for number concentrations,
        dydt[0:Nb] = dydt[0:Nb] + y[0:Nb] * gamma

    # print the elapsed time
    print('Elapsed time =',np.int(np.modf(t)[1]),
          '/',integration_time,'s',end='\r')

    return dydt
