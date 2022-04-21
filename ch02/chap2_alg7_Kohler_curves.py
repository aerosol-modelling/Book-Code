# Import the libraries we will be using
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import Kohler_theory_modules
# --- Physical constants --------------------------------------------
Lv_water_vapour=2.5e3 # Latent heat of vapourisation of water [J/g]
Rv=461.0 #Individual gas constant of water vapour [J/Kg.K]
Ra=287.0 #Gas constant for dry air [J/Kg.K]
R_gas=8.3144598 #Ideal gas constant [kg m2 s-2 K-1 mol-1]
R_gas_other=8.2057e-5 #Ideal gas constant [m3 atm K-1 mol-1]
GRAV=9.8; #Gravitational acceleration [(m/2)2]
cp=1005; #Specific heat capacity of air [J/Kg.K]
sigma=72.0e-3 # Assume surface tension of water (mN/m)
NA=6.0221409e+23 #Avogadros number
kb=1.380648E-23 #Boltzmanns constant
# ---------------------------------------------------------------------

Temp_K=298.15

surf_tens = 72.0e-3 #N.m
density_water = 1000.0 #kg/m3
density_solute = 1770.0 #kg/m3
molar_weight_water = 18.01528 #g/mol
molar_weight_solute = 132.14 #g/mol
molar_vol_solute = (molar_weight_solute*1.0e-3) / density_solute
molar_vol_water = molar_weight_water*1.0e-3 / density_water
diss_num = 3.0

option = 3

if option == 1:
    #1) For 3 sizes, calculate and plot the equilibrium curve to visualise the critical
    # droplet diameter of a single component aerosol particle
    dry_size_list = [50.0,100.0,250.0] #nm
    # Initialise a figure canvas to plot on
    plt.figure()
    for size in dry_size_list:

        mass_particle = (4.0/3.0)*np.pi*\
        np.power((0.5*size*1.0e-9),3.0)*density_solute #kg
        moles_solute = (mass_particle*1.0e3)/molar_weight_solute

        #pdb.set_trace()

        # Initialise a list of saturation ratios, Sw
        Sw_list = []
        droplet_size_list =[]
        # Now iterative through droplet sizes. To do this we add 1nm to the
        # initial dry particle size until we reach 10'000nm.
        droplet_size = size+1.0

        while droplet_size <= 10000.0:

            Dp = droplet_size*1.0e-9
            droplet_size_list.append(Dp)

            # Call the Kohler_curve function
            Sw=Kohler_theory_modules.Kohler_curve(\
                Dp,Temp_K,moles_solute,\
                diss_num,molar_vol_solute,\
                molar_vol_water,density_solute,\
                density_water,molar_weight_water,\
                surf_tens,R_gas)
            Sw_list.append(Sw)

            droplet_size=droplet_size+1.0

        # Plot the results on the figure
        # Note we are plotting supersaturation levels
        plt.plot(np.array(droplet_size_list),
        (np.array(Sw_list)-1.0)*100.0,label='%s nm' % size)
    # Constrain the y axis to show supersaturation
    plt.xscale('log')
    plt.ylim(-0.3, 0.5)
    ax = plt.gca()
    ax.axhline(0, color='black',lw=0.5)
    ax.set_ylabel('Supersaturation (%)')
    ax.set_xlabel('Droplet diamater (m)')
    plt.legend()
    plt.show()

elif option == 2:

    # For 10 sizes calculate the critical saturation ratio
    dry_size_list=[10,50,100,200,500]
    iterative_Sc=[]
    analytical_Sc=[]
    iterative_Dc=[]
    analytical_Dc=[]
    upper_limit= 50000 #nm

    for size in dry_size_list:

        Dp = size*1.0e-9 # Convert from nm to m
        mass_particle = (4.0/3.0)*np.pi*\
           np.power((0.5*Dp),3.0)*density_solute #kg
        moles_solute = (mass_particle*1.0e3)/\
           molar_weight_solute #moles

        # Calculate the critical saturation ratio by Brents method in minimize_scalar
        result = minimize_scalar(\
           Kohler_theory_modules.Kohler_curve_min,
           bounds=(size, upper_limit), \
           args=(Temp_K,moles_solute,stoich_coeff,\
           molar_vol_solute,molar_vol_water,\
           density_solute,density_water,\
           molar_weight_water,surf_tens,R_gas),\
           method='bounded')
        iterative_Sc.append(((result.fun*-1.0)-1.0)*100.0)
        iterative_Dc.append(result.x)

        # Calculate the critical saturation ratio using the analytical 2 parameter solution
        Sc_it, Dpc_it = Kohler_theory_modules.Crit_2param(Dp,\
           Temp_K,moles_solute,surf_tens,diss_num,\
           molar_vol_water,density_water,\
           molar_weight_water,R_gas)
        analytical_Sc.append((Sc_it-1.0)*100.0)
        analytical_Dc.append(Dpc_it*1.0e9)

    #Plot 2 subplots, one with supersaturation and the other critical diameter
    # Creates two subplots and unpacks the output array immediately
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(dry_size_list, iterative_Sc,label='Minimise Kohler')
    ax1.scatter(dry_size_list, analytical_Sc,label='Analytical solution')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylabel('Supersaturation %')
    ax1.set_xlabel('Dry diamater (nm)')
    ax2.plot(dry_size_list, iterative_Dc,label='Minimise Kohler')
    ax2.scatter(dry_size_list, analytical_Dc,label='Analytical solution')
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_ylabel('Critical Diameter (nm)')
    ax2.set_xlabel('Dry diamater (nm)')
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.legend()
    plt.show()

elif option == 3:

    dry_diameter = 80.0
    volume_frac_list = [0.002,0.004,0.02] #nm
    solubility_C=[1.6e-1,5.6e-2]
    Kappa_array=[1.28,0.23]
    # Initialise a figure canvas to plot on
    plt.figure()

    for vol_frac in volume_frac_list:

        # Initialise a list of saturation ratios, Sw
        Sw_list = []
        droplet_size_list =[]
        # Now iterative through droplet sizes. To do this we add 1nm to the
        # initial dry particle size until we reach 10'000nm.
        droplet_size = dry_diameter+1.0

        vol_frac_array=np.zeros((2),dtype=float)
        vol_frac_array[0]=vol_frac
        vol_frac_array[1]=1.0-vol_frac

        while droplet_size <= 10000.0:

            Dp = droplet_size*1.0e-9
            droplet_size_list.append(Dp)

            # Call the Kohler_curve function
            Sw=Kohler_theory_modules.Kappa_kohler(\
               Dp,dry_diameter*1.0e-9,Temp_K,\
               solubility_C,Kappa_array,surf_tens,\
               vol_frac_array,
            density_water,molar_weight_water,R_gas)
            Sw_list.append(Sw)

            droplet_size=droplet_size+1.0

        # Plot the results on the figure
        # Note we are plotting supersaturation levels
        plt.plot(np.array(droplet_size_list),
        (np.array(Sw_list)-1.0)*100.0,label=r'%s' % vol_frac)
    # Constrain the y axis to show supersaturation
    plt.xscale('log')
    plt.ylim(-1.0, 0.7)
    ax = plt.gca()
    ax.axhline(0, color='black',lw=0.5)
    ax.set_ylabel('Supersaturation (%)')
    ax.set_xlabel('Droplet diamater (m)')
    plt.legend()
    plt.show()
