import numpy as np

def GetParameters():
    d0 = 1e-8 # Monomer diameter
    rho = 1000.0
    mb = rho*(1/6)*np.pi*d0**3 # Monomer volume
    v_factor = 1e-10
    t_final = 3600.0
    return (mb,v_factor,t_final)

def GetKernel(m1,m2):
    # Assume that m1, m2 are actually volumes for now
    den_i = 1000.0 #000.0 # kg m^-3
    den_j = 1000.0 #000.0 # kg m^-3

    boltz = 1.3806505e-23
    pressure = 1e5
    univ_gas_const =  8.314472e0
    temp = 273.15 + 25.0 # K 
    MW_air =  2.89644e-2  # kg mole^-1
    rhoair = MW_air * pressure / (univ_gas_const * temp) # kg m^-3 
    viscosd = 1.8325e-05 * (416.16 / (temp + 120)) * \
         (temp / 296.16)**1.5
    avagadro = 6.02214179e23

    viscosk = viscosd / rhoair
    air_molec_weight = 2.89e-2
    gasspeed = np.sqrt((8.0 * boltz * temp * avagadro) / \
         (np.pi * air_molec_weight))

    gasfreepath = 2.0 * viscosk / gasspeed
 
    vol_i = m1 / den_i
    rad_i = ((3.0/(4.0 * np.pi)) * vol_i)**(1.0/3.0)
    Rme_i = rad_i
    knud      = gasfreepath/Rme_i
    cunning   = 1.0 + knud*(1.249 + 0.42*np.exp(-0.87/knud))
    diffus_i  = (boltz * temp * cunning) / \
         (6.0 * np.pi * Rme_i * viscosd)
    speedsq_i = 8.0 * boltz * temp / (np.pi * den_i * vol_i)
    freepath  = 8.0*diffus_i/(np.pi*np.sqrt(speedsq_i))
    tmp1      = (2.0*Rme_i + freepath)**3
    tmp2      = (4.0*Rme_i*Rme_i + freepath*freepath)**1.5
    deltasq_i = ( (tmp1-tmp2)/(6.0*Rme_i*freepath) - 2.0*Rme_i )**2

    vol_j = m2 / den_j 
    rad_j = ((3.0/(4.0 * np.pi)) * vol_j)**(1.0/3.0)
    Rme_j = rad_j
    knud      = gasfreepath/Rme_j
    cunning   = 1.0 + knud*(1.249 + 0.42*np.exp(-0.87/knud))
    diffus_j  = (boltz * temp * cunning) / \
        (6.0 * np.pi * Rme_j * viscosd)
    speedsq_j = 8.0 * boltz * temp / (np.pi * den_j * vol_j)
    freepath  = 8.0*diffus_j/(np.pi*np.sqrt(speedsq_j))
    tmp1      = (2.0*Rme_j + freepath)**3
    tmp2      = (4.0*Rme_j*Rme_j + freepath*freepath)**1.5
    deltasq_j = ( (tmp1-tmp2)/(6.0*Rme_j*freepath) - 2.0*Rme_j )**2

    rad_sum    = rad_i + rad_j
    diffus_sum = diffus_i + diffus_j
    tmp1       = rad_sum/(rad_sum + np.sqrt(deltasq_i + deltasq_j))
    tmp2       = 4.0*diffus_sum/(rad_sum*np.sqrt(speedsq_i + speedsq_j))
    bckernel  = 4.0*np.pi*rad_sum*diffus_sum/(tmp1 + tmp2)

    return(bckernel)
