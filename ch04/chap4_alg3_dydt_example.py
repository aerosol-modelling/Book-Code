# A demonstration of two different approaches to solving
# the ODES for a gas phase simulation.
import numpy as np
# The function that is called by the ODE solver and accepts
# time and concentration arrays as inputs as follows:
# input:
# t - time variable [internal to solver]
# y - concentrations of all compounds in both phases [molecules/cc]
# output:
# dydt - the dy_dt of each compound in each phase [molecules/cc.sec]

def dydt_func(t,y):

    # make sure the y array is not a list.
    y_asnumpy=np.array(y)

    # Initialise the output array
    dy_dt=np.zeros((num_species,1),)

    # Calculate time of day
    time_of_day_seconds=start_time+t

    # Calculate the concentration of RO2 species, using an index file created
    # during parsing
    RO2=np.sum(y_asnumpy[RO2_indices])

    # Call a separate function to evaluate reaction coefficients as a function
    # of temperature and other variables.
    rates=evaluate_rates(time_of_day_seconds,RO2,H2O,temp)

    # Call a separate function to calculate product of all reactants and
    # stoichiometry for each reaction [A^a*B^b etc]
    reactants=reactant_product(y_asnumpy)

    # Multiply product of reactants with rate coefficient to get reaction rate
    reactions = np.multiply(reactants,rates)

    # Option (1):
    # Use an explicit set of procedures
    dydt_gas=dydt_eval(np.zeros((num_species)),reactions)
    dydt_gas[0]=-1.0*reactions[0]-reactions[1]\
       -reactions[2]-reactions[3]
    dydt_gas[1]=-1.0*reactions[7]+reactions[0]\
       +reactions[41]
    dydt_gas[2]=-1.0*reactions[12]-reactions[13]
    dydt_gas[4]=-1.0*reactions[2]-reactions[3]\
       +reactions[178]+reactions[593]
    dydt_gas[5]=-1.0*reactions[18]-reactions[19]\
       +reactions[2]
    # .......

    # Option (2):
    # Use generic function that multiplies stoichiometric matrix by the
    # variables in reactions array
    dydt_gas=np.dot(stoich,reactions)
    # .......

    dy_dt[0:num_species,0]=dydt_gas

    return dy_dt
