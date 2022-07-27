import numpy as np
import pandas as pd
from scipy import constants
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.image as mpimg

rc={'axes.labelsize':15.0, 
    'legend.fontsize':15.0, 
    'axes.titlesize':15.0,
    'xtick.labelsize': 12.0,
    'ytick.labelsize': 12.0}
plt.rcParams.update(**rc)      
    
def get_vt(Dpi):
    """
    parameters
    Dpi: diameter, unit: um
    -------
    return:
    vt: terminal velocity, unit: m s-1
    """  
    # convert Dpi from "um" to "m"
    Dpi_SI = Dpi*1.0e-6 # m
    g = constants.g # m s-2
    rho = 1.0*1e3 # kg m-3
    mu = 1.83*1e-5 # viscosity, unit: kg m-1 s-1
    vt = Dpi_SI**2 * g *rho / (18.0 * mu)
    return vt

def get_K12_sed(Dp1,Dp2):        
    """
    parameters
    Dp1: diameter, unit: um
    Dp1: diameter, unit: um
    -------
    return:
    K12: coagulation kernel, unit: m3 s-1
    """    
    E = 1
    Dp1_SI = Dp1*1.0e-6 
    Dp2_SI = Dp2*1.0e-6
    vt1 = get_vt(Dp1) # m s-1
    vt2 = get_vt(Dp2) # m s-1
    K12_SI = E*np.pi/4.0*((Dp1_SI+Dp2_SI)**2)*np.abs(vt1-vt2) 
    # convert K12_SI from "m3 s-1" to "cm3 s-1"
    K12 = K12_SI * 1e6
    return K12

# generate data
Dp_ls = []
order_ls = range(-3,2)
for order in order_ls:
    for i in range(1,10):
        Dp_ls.append(i*10**order)

df_dict = {"D_p1":[], 
           "D_p2":[], 
           "K_12":[]}

for Dp1 in Dp_ls:
    for Dp2 in Dp_ls:
        df_dict["D_p1"].append(Dp1)
        df_dict["D_p2"].append(Dp2)
        df_dict["K_12"].append(get_K12_sed(Dp1,Dp2))        
df = pd.DataFrame(df_dict)

############ left figure ############
# plot the left figure 
fig = plt.figure(figsize=(5.2,5))
ax = fig.add_subplot(111)
D_p2_ls = [0.01,0.1,1.0,10]
color_ls = ["#0c2c84","#225ea8","#1d91c0","#41b6c4"]
for color, D_p2 in zip(color_ls,D_p2_ls):
    df_temp = df[df["D_p2"]==D_p2]
    df_temp = df_temp[df_temp["D_p2"]>=df_temp["D_p1"]]
    ax.plot(df_temp["D_p1"],
            df_temp["K_12"],
            label = str(D_p2)+ r" $\mu$m",
            color = color,
            lw = 2)  

# set up the scale and label
ax.set_yscale('log')
ax.set_ylim(10**(-22),10**(-1))
ax.set_ylabel(r"$K_\mathrm{12}$, $\mathrm{cm^3\ s^{-1}}$")
ax.set_xscale('log')
ax.set_xlim(0.001,20)
ax.set_xlabel(r"$d_\mathrm{p1}$, $\mathrm{\mu m}$")
plt.grid()
plt.legend(ncol=2)
plt.tight_layout()
plt.savefig("./ch05_q03_a1.pdf") 

############ right figure ############
# plot the right figure 
x=10**np.linspace(-3,2,200)
y=10**np.linspace(-3,2,200)
xx, yy = np.meshgrid(x, y)
z = get_K12_sed(xx,yy)
fig, ax = plt.subplots(figsize=(6.6,5))
lev_exp = np.arange(-23,
                    np.ceil(np.log10(z.max())+1))
levs = np.power(10, lev_exp)
im = ax.contourf(xx,yy, z, 
                 levs, norm=LogNorm(),
                 cmap="RdYlBu_r")
# set up the scale and label
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r"$d_\mathrm{p1}$, $\mathrm {\mu m}$")
ax.set_ylabel(r"$d_\mathrm{p2}$, $\mathrm {\mu m}$")
ax.set_xlim(1e-3,20)
ax.set_ylim(1e-3,20)
cbar = plt.colorbar(im)
cbar.set_label(r"$K_{12}$, $\mathrm{cm^3\ s^{-1}}$")
plt.tight_layout()
plt.savefig("./ch05_q03_a2.pdf") 

