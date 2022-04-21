# Import modules
import numpy as np
from scipy.integrate import odeint


# --  define some initial conditions --
# .....
# -------------------------------------

def dydt(input):

    # -- perform some calculations --
    # .....
    # -------------------------------

    return dydt_array

# Define the time over which the simulation will take place
t = np.linspace(0, 10000, num=1000)

# Call the ODE solver with reference to our function, dydt, setting the
# absolute and relative tolerance, atol and rtol respectively.
solution = odeint(dy_dt, array, t, rtol=1.0e-6, atol=1.0e-4, tcrit=None)

 # Do something with the output contained within variable 'solution'..
