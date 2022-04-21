import numpy as np

def Kohler_curve(Dp,Temp_K,moles_solute,diss_num,molar_vol_solute,
    molar_vol_water,density_solute,density_water,
    molar_weight_water,surf_tens,R_gas):

    # Convert Droplet size to water content and mole fraction
    # First calculate the number of moles of water
    volume_total = (4.0/3.0)*np.pi*np.power((Dp*0.5),3.0)
    volume_water = volume_total - (moles_solute*molar_vol_solute)
    moles_water = volume_water/molar_vol_water
    # Then calculate a mole fraction, and density, assuming additive mixing rule
    X_w = moles_water / (moles_water+moles_solute*diss_num)
    # Mass fractions. Create a numpy array of mass loadings
    total_mass = moles_solute*molar_vol_solute+moles_water*molar_vol_water
    mass_fractions = np.zeros((2),dtype=float)
    mass_fractions[0] = (moles_solute*molar_vol_solute)/total_mass
    mass_fractions[1] = (moles_water*molar_vol_water)/total_mass
    # Now create a numpy array of densities for the solute and water
    density_pure = np.zeros((2),dtype=float)
    density_pure[0] = density_solute
    density_pure[1] = density_water
    # Now calculate the density of the solution
    density_solution = 1.0/(np.sum(mass_fractions/density_pure))
    # Calculate the Kelvin effect
    Kelvin = np.exp((4.0*(molar_weight_water*1.0e-3)*surf_tens)/
             (R_gas*Temp_K*density_water*Dp))
    # Calculate the equilibrium saturation ratio
    Sw = Kelvin*X_w

    return Sw

def Kohler_curve_min(Dp,Temp_K,moles_solute,diss_num,molar_vol_solute,
    molar_vol_water,density_solute,density_water,
    molar_weight_water,surf_tens,R_gas):

    Dp=Dp*1.0e-9
    # Convert Droplet size to water content and mole fraction
    # First calculate the number of moles of water
    volume_total = (4.0/3.0)*np.pi*np.power((Dp*0.5),3.0)
    volume_water = volume_total - (moles_solute*molar_vol_solute)
    moles_water = volume_water/molar_vol_water
    # Then calculate a mole fraction, and density, assuming additive mixing rule
    X_w = moles_water / (moles_water+moles_solute*diss_num)
    # Mass fractions. Create a numpy array of mass loadings
    total_mass = moles_solute*molar_vol_solute+moles_water*molar_vol_water
    mass_fractions = np.zeros((2),dtype=float)
    mass_fractions[0] = (moles_solute*molar_vol_solute)/total_mass
    mass_fractions[1] = (moles_water*molar_vol_water)/total_mass
    # Now create a numpy array of densities for the solute and water
    density_pure = np.zeros((2),dtype=float)
    density_pure[0] = density_solute
    density_pure[1] = density_water
    # Now calculate the density of the solution
    density_solution = 1.0/(np.sum(mass_fractions/density_pure))
    # Calculate the Kelvin effect
    Kelvin = np.exp((4.0*(molar_weight_water*1.0e-3)*surf_tens)/
             (R_gas*Temp_K*density_water*Dp))

    # Calculate the equilibrium saturation ratio
    Sw = Kelvin*X_w
    return -1.0*Sw


def Crit_2param(Dp,Temp_K,moles_solute,surf_tens,diss_num,
    molar_vol_water,density_water,molar_weight_water,R_gas):

    # Calculate the parameters A and B
    A = (4.0*(molar_weight_water*1.0e-3)*surf_tens)/\
        (R_gas*Temp_K*density_water)
    B = (6.0*diss_num*moles_solute*(molar_weight_water*1.0e-3))/\
        (np.pi*density_water)
    ln_Sc = np.power((4.0*np.power(A,3.0))/(27.0*B),0.5)
    Sc = np.exp(ln_Sc)
    Dpc = np.power(((3.0*B)/A),0.5)

    return Sc, Dpc

def Kappa_kohler(Dp,D_dry,Temp_K,solubility_C,Kappa_array,
    surf_tens,vol_frac,density_water,molar_weight_water,R_gas):

    # Calculate the growth factor
    g=Dp/D_dry
    # Calculate the dissolved volume fraction
    xi=((g**3.0)-1.0)*(solubility_C/vol_frac)
    H_xi = np.zeros(len(vol_frac),dtype=float)
    step=0
    for entry in xi:
        if entry < 1.0:
            H_xi[step]=xi[step]
        else:
            H_xi[step]=1.0
        step+=1

    # Calculate the Kappa parameter
    Kappa = np.sum(vol_frac*Kappa_array*H_xi)

    # Calculate the Kelvin effect
    Kelvin = np.exp((4.0*(molar_weight_water*1.0e-3)*surf_tens)/
             (R_gas*Temp_K*density_water*Dp))

    # Calculate the equilibrium saturation ratio
    Sw = Kelvin*((Dp**3.0-D_dry**3.0)/(Dp**3.0-(D_dry**3.0)*(1.0-Kappa)))
    return Sw
