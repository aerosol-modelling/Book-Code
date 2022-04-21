# Import the relevant libraries
using DifferentialEquations
using DiffEqBase
using DiffEqCallbacks
using DiffEqOperators
using Sundials
using Plots

# Define physical constants
const R_gas=8.31451 #Ideal gas constant [kg m2 s-2 K-1 mol-1]
const R_gas_other=8.20573e-5 #Ideal gas constant [m3 atm K-1 mol-1]
const sigma=72.0e-3 # Assume surface tension of water (mN/m)
const NA=6.0221409e+23 #Avogadros number

Temp_K=298.15

# Number of condensing species from the gas phase
const num_species = 10

# The molecular weight of each condensing specie [g/mol]
# Assuming a constant value for all components
mw_array=zeros(num_species)
mw_array.=200.0

# Define array of log10 C* values
log_c_star = collect(range(-6, 3, length = num_species))
Cstar = broadcast(^, 10.0, log_c_star)
# Convert C* to a pure component saturation vapour pressure [atm]
P_sat = (Cstar*R_gas_other*Temp_K)./(1.0e6*mw_array)

# Initialise abundance in each volatility bin [micrograms/m3]
concentration = zeros(num_species)
concentration[1] = 0.1
concentration[2] = 0.1
concentration[3] = 0.15
concentration[4] = 0.22
concentration[5] = 0.36
concentration[6] = 0.47
concentration[7] = 0.58
concentration[8] = 0.69
concentration[9] = 0.84
concentration[10] = 1.0

# Unit conversion of gas abudance to molecules / cc
gas_concentration = ((concentration*1.0e-6)./(mw_array))*1.0e-6*NA

# Set accomodation coefficient
alpha_d_org=zeros(num_species)
alpha_d_org.=1.0
# Set density of condensing species [kg/m3]
density_org=zeros(num_species)
density_org.=1400.0

# Molecular diffusion coefficient in air (Equation 2.22). [cm2/s]
DStar_org = 1.9*(mw_array.^(-2.0/3.0))
# Mean thermal velocity of each molecule (Equation 2.21). [m/s]
mean_them_vel=((8.0*R_gas*Temp_K)./(pi*(mw_array*1.0e-3))).^0.5
# Mean free path for each molecule (Equation 2.20). [m]
gamma_gas = ((3.0*DStar_org)./(mean_them_vel*1.0e2))*1.0e-2

# Define a polydisperse size distribution
# Set the smallest size
d1 = 0.01
# Set the largest size
d_Nb = 1.0
# Number of size bins
num_bins = 8

# Volume ratio between bins (Equation 1.17)
V_rat =(d_Nb/d1).^(3.0/(num_bins-1.0))
# Diameter ratio between bins  (Equation 1.18)
d_rat = V_rat.^(1.0/3.0)

# Create an array of diameters
d_i=zeros(num_bins)
d_i[1]=d1
for step in range(1,num_bins,step=1)
    if step > 1
       d_i[step]=d_i[step-1]*d_rat
   end
end
# Log of Diameter array
log_di = log.(d_i)

# Define parameters of log-normal distribution
# Geometric standard deviation
sigmag1 = log.(1.7)
# Mean particle diameter [150nm]
mean1 = log.(0.15)
# Calculate the probability density distribution
distribution_1 = (exp.((-1.0.*(log_di .- mean1).^2) ./ (2 * sigmag1^2)) ./ (sigmag1 * sqrt.(2 * pi)))
# Seperate out the probability density function
d_width = d_i.*(2^(1.0/3.0))*((V_rat.^(1.0/3.0)-1.0)/((1+V_rat).^(1.0/3.0)))
# Total number of particles [per cm-3]
N_total = 100.0
# Discrete number distribution
N_dist = N_total.*(distribution_1.*(d_width./d_i))

# Initialise a core abundance using the above size distribution
core = zeros(num_bins)
core_concentration =zeros(num_bins)
density_core = zeros(num_bins)
core_mw = zeros(num_bins)
density_core.= 1400.0
core_mw.= 200.0
N_per_bin = N_dist

# Define size array (radius )
size_array = d_i.*0.5
# Use the size to now calculate a concentration of a 'core' in molecules / cc
core_concentration = (N_per_bin).*((4.0/3.0)*pi*((size_array.*1.0e-6).^3.0).*density_core.*1.0e3)
core_concentration = (core_concentration ./ core_mw).*NA
dry_mass = sum((core_concentration./NA).*core_mw).*(1.0e6)*(1.0e6)
print("Initial dry mass = ", dry_mass)

# New define an array that holds the molecular abundance of
# each gas and the concentration of each gas in each size bin
array = zeros(num_species+num_species*num_bins)
array[1:num_species].= gas_concentration
array[num_species+1:num_species+num_species*num_bins].= 1.0e-10 # assuming we start with nothing

# Define the RHS function (that includes droplet growth equation)
function dy_dt!(u,p,t)

    # Retrieve abundance
    array=u
    # Access parameters through parameter dictionary
    num_species,N_per_bin,num_bins=[p[i] for i in ["num_species","N_per_bin","num_bins"]]
    core_concentration,density_org,density_core,mw_array,core_mw=[p[i] for i in ["core_concentration","density_org","density_core","mw_array","core_mw"]]
    alpha_d_org,gamma_gas,DStar_org=[p[i] for i in ["alpha_d_org","gamma_gas","DStar_org"]]
    sigma,R_gas,Temp_K,P_sat,NA,R_gas_other=[p[i] for i in ["sigma","R_gas","Temp_K","P_sat","NA","R_gas_other"]]
    # Retrieve gas phase abundance
    Cg_i_m_t = array[1:num_species]

    # We are working with 8 size bins, each of which has an involatile core
    size_array = zeros(num_bins)
    #Initialise empty dydt arrays
    dy_dt_array = zeros(num_species+num_species*num_bins)
    dy_dt_gas_matrix = zeros(num_species,num_bins)

    # Now cycle through each size bin
    for size_step = 1:num_bins

        # Select a slice of y that represents this size bin
        temp_array=array[(num_species+(size_step-1)*num_species)+1:num_species+num_species*(size_step)]
        # Sum the total molecules in the condensed phase, adding core material
        total_moles=sum(temp_array)+core_concentration[size_step]
        # Calculate the mole fractions in the, assumed, liquid phase for use
        # in calculating the equilibrium pressure above the droplet
        mole_fractions=temp_array./total_moles
        # Calculate the density of the assumed liquid phase
        density_array = zeros(num_species+1)
        density_array[1:num_species].=density_org[1:num_species]
        density_array[num_species+1]=density_core[size_step]
        # Create an array that holds mass concentrations [g/cc], used for
        # calculation of solution density [kg/m3]
        mass_array = zeros(num_species+1)
        mass_array[1:num_species].= (temp_array./NA).*mw_array
        mass_array[num_species+1]= (core_concentration[size_step]/NA)*core_mw[size_step]
        total_mass=sum(mass_array)
        mass_fractions_array=mass_array./total_mass
        density=1.0./(sum(mass_fractions_array./density_array))
        # Now calculate the size [cm] of our particles according
        # to the condensed phase abundance of material
        # We need to remember that mass is in [g/cc] whilst density
        # is in [kg/m3].
        # Thus we convert mass and number concentrations to kg/m3 and /m3
        size_array[size_step]=((3.0*((total_mass*1.0e3)/(N_per_bin[size_step]*1.0e6)))/(4.0*pi*density))^(1.0/3.0)
        # Calculate the Knudsen number for all condensing molecules based
        # on this new size
        # This relies on mean free path for each species [cm] and
        # particle radius [cm]
        Kn=gamma_gas./size_array[size_step]
        # Calculate Non-continuum regime correction (Equation 2.19)
        # Calculate a correction factor according to the continuum versus
        # non-continuum regimes
        Inverse_Kn=1.0./Kn
        Correction_part1=(1.33e0.+0.71e0.*Inverse_Kn)./(1.0e0.+Inverse_Kn)
        Correction_part2=(4.0e0.*(1.0e0.-alpha_d_org))./(3.0e0.*alpha_d_org)
        Correction_part3=1.0e0.+(Correction_part1.+Correction_part2).*Kn
        Correction=1.0./Correction_part3
        # Now calculate a kelvin factor for every semi-volatile compound in this
        # size bin (Equation 2.28)
        kelvin_factor=exp.((4.0E0.*mw_array.*1.0e-3.*sigma)/(R_gas.*Temp_K.*size_array[size_step].*2.0e0.*density))
        # Calculate the equilibrium concentration [molecules/cc]
        Pressure_eq=kelvin_factor.*mole_fractions.*P_sat
        # Calculate the equilibrium concentration equivalent
        Cstar_i_m_t=Pressure_eq.*(NA./(R_gas_other.*1.0e6.*Temp_K))
        # Implement the droplet growth equation (Equation 2.44)
        k_i_m_t_part1 = DStar_org.*Correction
        k_i_m_t=4.0e0.*pi.*size_array[size_step].*1.0e2.*N_per_bin[size_step].*k_i_m_t_part1
        dm_dt=k_i_m_t.*(Cg_i_m_t.-Cstar_i_m_t)
        # Now update the contribution to the ODEs being solved
        # Add contributory loss from the gas phase to particle phase
        dy_dt_gas_matrix[1:num_species,size_step].=dm_dt[1:num_species]
        # Add a contributory gain to the particle phase from the gas phase
        dy_dt_array[(num_species+(size_step-1)*num_species)+1:num_species+num_species*(size_step)].=dm_dt[1:num_species]

    end

    # Subtract the net condensational mass flux from the gas phase concentrations
    dy_dt_array[1:num_species]=dy_dt_array[1:num_species].-sum(dy_dt_gas_matrix, dims = 2)[:,1]

    du=dy_dt_array

    return du

end
# Create a dictionary that contains properties of condensing gases and
# the liquid mixture
param_dict=Dict("num_species"=>num_species,"N_per_bin"=>N_per_bin,"num_bins"=>num_bins,#"dydt"=>dydt,
                    "core_concentration"=>core_concentration,"density_org"=>density_org,"density_core"=>density_core,
                    "mw_array"=>mw_array,"core_mw"=>core_mw,"alpha_d_org"=>alpha_d_org,
                    "gamma_gas"=>gamma_gas,"DStar_org"=>DStar_org,"sigma"=>sigma,
                    "R_gas"=>R_gas,"Temp_K"=>Temp_K,"P_sat"=>P_sat,
                    "NA"=>NA,"R_gas_other"=>R_gas_other)
# Define the time over which the simulation will take place
tspan = (0.0, 10000.0)
# Define the problem to be solved [using DifferentialEquations]
prob = ODEProblem(dy_dt!,array,tspan,param_dict)
saved_values = SavedValues(Float64, Vector{Float64})
cb = SavingCallback((u,t,integrator)->integrator(t,Val{1}), saved_values)
# Call the ODE solver with reference to our function, dy_dt, setting the
# absolute and relative tolerance.
sol = solve(prob, alg=CVODE_BDF(), abstol = 1e-4, reltol = 1e-6, dt=1e-6, dtmax= 10.0, callback=cb)

# Now prouce a plot of the size bin growth over time
# Define an array that will store the size of our
# particles at each model output
size_array_matrix = zeros((size(sol.t)[1],num_bins))

# Let's now try to plot the size variation as a function of time
for step_var = 1:size(sol.t)[1]
   # For each point in time, calculate the total SOA
   for size_step = 1:num_bins#
       temp_array=sol[(num_species+(size_step-1)*num_species)+1:num_species+num_species*(size_step),step_var]
       total_moles=sum(temp_array)+core_concentration[size_step]
       mole_fractions=temp_array./total_moles
       density_array = zeros(num_species+1)
       density_array[1:num_species].=density_org[1:num_species]
       density_array[num_species+1]=density_core[size_step]
       mass_array = zeros(num_species+1)
       mass_array[1:num_species].= (temp_array./NA).*mw_array
       mass_array[num_species+1]= (core_concentration[size_step]/NA)*core_mw[size_step]
       total_mass=sum(mass_array)
       mass_fractions_array=mass_array./total_mass
       density=1.0./(sum(mass_fractions_array./density_array))
       size_array_matrix[step_var,size_step]=((3.0*((total_mass*1.0e3)/(N_per_bin[size_step]*1.0e6)))/(4.0*pi*density))^(1.0/3.0)

   end
end

# Create a matrix that holds the relative change in size
size_change_matrix = zeros((size(sol.t)[1],num_bins))
for step_var = 1:size(sol.t)[1]
    size_change_matrix[step_var,:]=size_array_matrix[step_var,:]./size_array_matrix[1,:]
end

# Plot the decrease in gas phase abundance with the growth of the particle [both absolute and relative]

p1 = plot(sol.t, log10.(sol[1:num_species,:]'), ylabel = "Gas phase abundance [Log10/ molecules.cc]", label="",lw = 3, size = (600, 700),linecolor = :red) # Make a line plot
xlabel!("time [seconds]")
p2 = plot(sol.t, size_array_matrix.*2.0e6, ylabel = "Particle diameter [microns]", label="",lw = 3, size = (600, 700),linecolor = :red)
xlabel!("time [seconds]")
p3 = plot(sol.t, log10.(size_change_matrix), ylabel = "Relative change in particle size [log10]", label="",lw = 3, size = (600, 700),linecolor = :red)
xlabel!("time [seconds]")
# Four histograms each with 10 points? Why not!
plot(p1, p2, p3, layout = (1, 3), legend = false)
savefig("Julia_multibin.png")
