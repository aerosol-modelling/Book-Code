import numpy as np
import pandas as pd
from scipy import constants
import matplotlib.pyplot as plt
from q03_helper import *

# generate data
Dp1 = 1
# create a range of Dp2
Dp_ls = np.logspace(-3,2,(2-(-3))*200+1)
        
df_dict = {"D_p2":[], 
           "Brownian":[], 
           "Sedimentation":[]}

for Dp2 in Dp_ls:
    df_dict["D_p2"].append(Dp2)
    df_dict["Brownian"].append(get_K12(Dp1,Dp2))
    df_dict["Sedimentation"].append(get_K12_sed(Dp1,Dp2))        
df = pd.DataFrame(df_dict)

############ plot figure ############
plt.plot(df["D_p2"],df["Brownian"],label="Brownian")
plt.plot(df["D_p2"],df["Sedimentation"],label="Sedimentation")
plt.xscale("log")
plt.xlabel(r"$d_\mathrm{p2}$, $\mathrm{\mu m}$")
plt.yscale("log")
plt.ylabel(r"$K_\mathrm{12}$, $\mathrm{cm^3\ s^{-1}}$")
plt.legend()
plt.tight_layout()
plt.savefig("./ch05_q03_b.pdf") 

# find the relative importance
df["diff"] = df["Sedimentation"]-df["Brownian"] 

print(r"when Sedimentation is more important, the diameter should be greater than")
print(df.loc[[df.loc[df['diff'] >= 0, "diff"].idxmin()]])

print("\n")
print("when Brownian is more important, the diameter should be smaller than")
print(df.loc[[df.loc[df['diff'] <= 0, "diff"].idxmax()]])