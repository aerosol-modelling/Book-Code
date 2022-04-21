import numpy as np
from parameters import *

def quasi_stationary(v, y):

    # initialize variables
    delta_n, delta_V, n, V, v_j, bins = initialize_remapping(y)

    for j in range(Nb):

        # find the index of bins between which the size of each bin
        # has evolved: k an l in Equation (7.1)
        l = np.argmax(v_j[j] < v)
        k = l - 1

        # if the index was found and the bin is not empty
        if l > 0 and n[j] > N_low_limit:

            # calculate the fraction of n_j remapped to n_l (Equation 7.3)
            x_l = (v_j[j] - v[k])/(v[l] - v[k])

            x_k = 1.0 - x_l

            # fraction of n_j to be remapped to n_l (Eq (7.3))
            delta_n[l] += x_l*n[j]
            # fraction of n_j to be remapped to n_k (Eq (7.4))
            delta_n[k] += x_k*n[j] 
            
            # loop over individual species
            for i in All_species:
                # fraction of V_i,j to be remapped to V_i,k (Eq (7.5))
                delta_V[i][l] += x_l*V[i][j]*v[l]/v_j[j] 
                # fraction of V_i,j to be remapped to V_i,l (Eq (7.5))
                delta_V[i][k] += x_k*V[i][j]*v[k]/v_j[j]

        # if the size bin does not fall between any bins
        else:
            
            delta_n[j] += n[j]

            # loop over individual species
            for i in All_species:
                delta_V[i][j] += V[i][j]

    # substitute the remapped values to number concentrations
    n = delta_n

    # substitute these values to array y
    y[range(Nb)] = n

    # volume concentrations of individual species
    for i in All_species:
    
        V[i] = delta_V[i]

        # convert volume concentration back to molar concentration
        # and substitute the values in array y
        y[bins[i]] = V[i]*rho[i]/M[i]

    return y

def moving_center(v_hi, y):

    delta_n, delta_V, n, V, v_j, bins = initialize_remapping(y)

    for j in range(Nb):

        # find the index of the size bin
        # where v_lo < v_j < v_hi
        l = np.argmax(v_j[j] < v_hi)

        # check if the index was found, the size class has grown/evaporated
        # out of its initial size bin, and the bin is not empty
        if l != 0 and j != l and n[j] > N_low_limit:

            # particle numbers are transferred from bin j to bin l
            delta_n[l] += n[j]
            # the same amount is subtracted from bin j
            delta_n[j] -= n[j] 
            
            # loop over individual species
            for i in All_species:
                # all volume of compound i in bin j is transferred to size bin l
                delta_V[i][l] += V[i][j] 
                # this amount is subtracted from bin j
                delta_V[i][j] -= V[i][j]

    # transfer the remapped values to number concentrations
    n = n + delta_n

    # substitute these values to array y
    y[range(Nb)] = n

    # volume concentrations of individual species
    for i in All_species:
    
        V[i] = V[i] + delta_V[i]

        # convert volume concentration back to molar concentration
        # and substitute the values in array y
        y[bins[i]] = V[i]*rho[i]/M[i]

    return y

def hybrid_bin(v, v_lo, v_hi, y):

    # initialize variables for remapping
    delta_n, delta_V, n, V, v_t_dt, bins = initialize_remapping(y)

    # volume concentration of all species
    # in the grown/evaporated size bins
    V_t_dt = v_t_dt * np.maximum(n, N_low_limit)
    
    for j in range(Nb-1):

        # If bin is not empty 
        if n[j] > N_low_limit:

            # Calculate the change in volume
            # during a time step (Eq (7.13))
            delta_v = v_t_dt[j] - v[j]

            # Volumes of the bounds
            # of size bins (Eq (7.12))
            # assuming pre-growth linear method
            v_1 = v_lo[j] + delta_v
            v_2 = v_hi[j] + delta_v

            # calculate n_0, k, (Eq (7.10))
            # and x_0 for integrating Eqs (7.16) and (7.17)
            n_0 = n[j] / (v_2-v_1)
            v_0 = (v_2+v_1) / 2.0
            k = 12.0 * (V_t_dt[j]-v_0*n[j])/(v_2-v_1)**3
            x_0 = v_0
            
            # if the particles have grown and
            # we are not remapping the largest bin
            if delta_v > 0 and j != Nb:

                # set the limits a and b
                # for integrating Eqs (7.16) and (7.17)
                a = v_hi[j]
                b = v_2
                # the index of the bin to which
                # particles are transferred to
                m = j+1
                
                # positivity check
                # y(v_1) in Eq (7.21)            
                if k * (v_1 - v_0) + n_0 < 0:

                    # Equation (7.21)
                    # Notice that the volume concentrations
                    # of individual species are summed
                    v_ast = 3.0*V_t_dt[j]/n[j] - 2.0*v_2
                    n_0 = 0.0
                    x_0 = v_ast
                    k = k_ast = 2.0*n[j]/(v_2-v_0)**2
                    a = v_hi
                    b = v_2

                
                # v(v_2) in Eq (7.22)
                if k * (v_2 - v_0) + n_0 < 0:
                    # Equation (7.22)
                    v_ast = 3.0*V_t_dt[j]/n[j] - 2.0*v_1
                    n_0 = 0.0
                    x_0 = v_ast
                    k = k_ast = -2.0*n[j]/(v_1-v_0)**2
                    a = v_2
                    b = v_ast
                    
            # if the particles have shrunk in size
            elif delta_v < 0 and j != 0:
                
                # Equation (7.19)
                a = v_1
                b = v_lo[j]
                # the index of the bin to which
                # particles are tranferred to
                m = j-1

                # positivity check
                # y(v_1) in Eq (7.21)            
                if k * (v_1 - v_0) + n_0 < 0:

                    # Equation (7.21)
                    # Notice that the volume concentrations
                    # of individual species are summed
                    v_ast = 3.0*V_t_dt[j]/n[j] - 2.0*v_2
                    n_0 = 0.0
                    x_0 = v_ast
                    k = k_ast = 2.0*n[j]/(v_2-v_0)**2
                    a = v_ast
                    b = v_lo[j]

                
                # v(v_2) in Eq (7.22)
                if k * (v_2 - v_0) + n_0 < 0:

                    # Equation (7.22)
                    v_ast = 3.0*V_t_dt[j]/n[j] - 2.0*v_1
                    n_0 = 0.0
                    x_0 = v_ast
                    k = k_ast = -2.0*n[j]/(v_1-v_0)**2
                    a = v_1
                    b = v_lo[j]
                
                
            else:
                # nothing is remapped
                k = 0.0
                a = v[j]
                b = v[j]
                m = j
                
            # if a or b are outside the fixed
            # bin boundaries, the bin is remapped

            if b > v_hi[j] or a < v_lo[j]:
                # Equation (7.16)
                delta_n[m] += (b-a) * (n_0-k*(x_0-(a+b)/2.0))
                delta_n[j] -= (b-a) * (n_0-k*(x_0-(a+b)/2.0))

                # Equation (7.17)
                delta_V_all = n_0 * (b**2-a**2)/2.0 \
                    + k*((b**3-a**3)/3.0 - x_0*((b**2-a**2)/2.0))

                # loop over individual species
                for i in All_species:
                    # V[i][j]/np.sum(V[:][j]) is the volume fraction
                    # of species i in size bin j
                    delta_V[i][m] += delta_V_all * V[i][j]/(v_t_dt[j]*n[j])
                    delta_V[i][j] -= delta_V_all * V[i][j]/(v_t_dt[j]*n[j])

    # substitute the remapped values to number concentrations
    n = np.maximum(n + delta_n, 0.0)

    # substitute these values to array y
    y[range(Nb)] = n

    # volume concentrations of individual species
    for i in All_species:
    
        V[i] = np.maximum(V[i] + delta_V[i], 0.0)

        # convert volume concentration back to molar concentration
        # and substitute the values in array y
        y[bins[i]] = V[i]*rho[i]/M[i]

    return y

def flat_top(v, v_lo, v_hi, y):

    # initialize variables for remapping
    delta_n, delta_V, n, V, v_t_dt, bins = initialize_remapping(y)

    # volume concentration of all species
    # in the grown/evaporated size bins
    V_t_dt = v_t_dt * np.maximum(n, N_low_limit)
    
    for j in range(Nb-1):

        # If bin is not empty 
        if n[j] > N_low_limit:

            # Calculate the change in volume
            # during a time step (Eq (7.13))
            delta_v = v_t_dt[j] - v[j]

            # Volumes of the bounds
            # of size bins
            v_1 = v_lo[j] + delta_v
            v_2 = v_hi[j] + delta_v

            # Equation (7.10)
            n_0 = n[j] / (v_2 - v_1)

            # if the particles have grown and
            # we are not remapping the largest bin
            if delta_v > 0 and j != Nb:

                b = v_2
                a = v_hi[j]
                m = j + 1

            # if the particles have shrunk in size               
            elif delta_v <= 0 and j != 0:

                b = v_lo[j]
                a = v_1
                m = j - 1

            else:
                # nothing is remapped
                k = 0.0
                a = v[j]
                b = v[j]
                m = j

            delta_n[j] -= np.minimum(n[j], (b-a)*n_0)
            delta_n[m] += np.minimum(n[j], (b-a)*n_0)

            # loop over individual species
            for i in All_species:
                
                # V[i][j]/np.sum(V[:][j]) is the volume fraction
                # of species i in size bin j
                delta_V[i][j] -= np.minimum(V[i][j],n_0 * (b**2-a**2)/2.0 * V[i][j]/V_t_dt[j])
                delta_V[i][m] += np.minimum(V[i][j],n_0 * (b**2-a**2)/2.0 * V[i][j]/V_t_dt[j])

    # substitute the remapped values to number concentrations
    n = n + delta_n

    # substitute these values to array y
    y[range(Nb)] = n

    # volume concentrations of individual species
    for i in All_species:
    
        V[i] = np.maximum(V[i] + delta_V[i], 0.0)

        # convert volume concentration back to molar concentration
        # and substitute the values in array y
        y[bins[i]] = V[i]*rho[i]/M[i]

    return y

def initialize_remapping(y):
    
    # initialize variables
    v_j = np.zeros(Nb)     # Volumes of grown/evaporated bins
    delta_n = np.zeros(Nb) # fraction of remapped number concentration in size bins
    delta_V = dict()       # fraction of remapped volume concentrations in size bins
    for i in All_species:
        delta_V[i] = np.zeros(Nb)
    bins = dict()

    V = dict()             # volume concentrations 
    
    # number concentration in each bin
    n = y[range(Nb)]
    
    # Calculate the total volume of each bin
    for i in All_species:

        # find the indices of each species i
        bins[i] = range(Nb*index[i],Nb*index[i]+Nb)

        # convert molar concentration of individual i
        # to volume concentrations [m3/m3]
        V[i] = y[bins[i]]*M[i]/rho[i]
        
        # volume of individual particles in grown/evaporated bins
        v_j+=V[i]/np.maximum(n, N_low_limit)

    return delta_n, delta_V, n, V, v_j, bins

