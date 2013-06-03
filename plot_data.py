import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import asciitable

PROTOSTELLARDISK = True
SUPERSONICTURBULENCE = False
SNSHOCKWAVE = False

def plot_data_sife():
    # load data
    data = asciitable.read("DATA/sih_data.csv",delimiter=',')
    x = data['feh']
    y = data['sih']
    ulindex = np.where(data['lsi'] == '<')[0]
    
    # create figure
    fig = plt.figure(1, figsize=(7.5,7.5))
    
    ax = plt.subplot(111)
    p = ax.scatter(x, y, c='k', s=15.0, facecolors='none')
    
    # plot arrows
    arrowdx = 0.0
    arrowdy = -0.1
    arrowhw = .04
    arrowhl = .08
    for ix in ulindex:
        ax.arrow(x[ix],y[ix],arrowdx,arrowdy,fc='r',ec='r',head_width=arrowhw,head_length=arrowhl)
    
    # plot upper limit bars for HE0557 and HE1327
    barwidth = .02
    ax.add_patch(plt.Rectangle((-4.81-barwidth/2.,-5.4),barwidth,-4.84+5.4,color='orange'))
    ax.arrow(-4.81,-5.4,arrowdx,arrowdy,fc='r',ec='r',head_width=arrowhw,head_length=arrowhl)
    ax.add_patch(plt.Rectangle((-5.76-barwidth/2.,-5.4),barwidth,-5.2+5.4,color='orange'))
    ax.arrow(-5.76,-5.4,arrowdx,arrowdy,fc='r',ec='r',head_width=arrowhw,head_length=arrowhl)

    return fig,ax

#def plot_data_dtranssi():
    

if PROTOSTELLARDISK:    
    fig,ax = plot_data_sife()
    # plot [Si/H]crit lines
    sihcrit = -5.33
    sihcriterr = .2
    
    ax.plot([-6.0,-1.0],[sihcrit+sihcriterr,sihcrit+sihcriterr],'b:')
    ax.plot([-6.0,-1.0],[sihcrit,sihcrit],'b--')
    ax.plot([-6.0,-1.0],[sihcrit-sihcriterr,sihcrit-sihcriterr],'b:')
    # std size distribution
    ax.plot([-6.0,-1.0],[-4.53,-4.53],'g--')
    ax.plot([-6.0,-1.0],[-4.58,-4.58],'g--')
    ax.plot([-6.0,-1.0],[-4.63,-4.63],'g--')
    ax.plot([-6.0,-1.0],[-4.45,-4.45],'g--')
    # x5.5 size distribution
    ax.plot([-6.0,-1.0],[-5.21,-5.21],'r--')
    ax.plot([-6.0,-1.0],[-5.34,-5.34],'r--')
    ax.plot([-6.0,-1.0],[-5.38,-5.38],'r--')
    ax.plot([-6.0,-1.0],[-5.25,-5.25],'r--')

    ax.set_xlabel(r'[Fe/H]')
    ax.set_ylabel(r'[Si/H]')
    ax.set_xlim(-6.0,-1.0)
    ax.set_ylim(-6.0,-1.0)
    ax.set_title(r"Critical [Si/H] Protostellar Disk  $n=10^{12}$ $T=1000$ (20 $M_\odot$)")
    plt.savefig("PLOTS/yongplus_abundances_psdisk_m20.pdf")

## # std size distribution
## ax.plot([-6.0,-1.0],[-4.53,-4.53],'g--')
## ax.plot([-6.0,-1.0],[-4.58,-4.58],'g--')
## ax.plot([-6.0,-1.0],[-4.63,-4.63],'g--')
## ax.plot([-6.0,-1.0],[-4.45,-4.45],'g--')
## # x5.5 size distribution
## ax.plot([-6.0,-1.0],[-5.21,-5.21],'r--')
## ax.plot([-6.0,-1.0],[-5.34,-5.34],'r--')
## ax.plot([-6.0,-1.0],[-5.38,-5.38],'r--')
## ax.plot([-6.0,-1.0],[-5.25,-5.25],'r--')
## # std size distribution
## ax.plot([-6.0,-1.0],[-4.53,-4.53],'g--')
## ax.plot([-6.0,-1.0],[-4.58,-4.58],'g--')
## ax.plot([-6.0,-1.0],[-4.63,-4.63],'g--')
## ax.plot([-6.0,-1.0],[-4.45,-4.45],'g--')
## # x5.5 size distribution
## ax.plot([-6.0,-1.0],[-5.21,-5.21],'r--')
## ax.plot([-6.0,-1.0],[-5.34,-5.34],'r--')
## ax.plot([-6.0,-1.0],[-5.38,-5.38],'r--')
## ax.plot([-6.0,-1.0],[-5.25,-5.25],'r--')

#ax.set_title(r"Critical [Si/H] $n=10^{12}$ $T=1000$ (20 $M_\odot$)")
#plt.savefig("PLOTS/yongplus_abundances_n12_T1000.pdf")
#ax.set_title(r"Critical [Si/H] $n=10^{12}$ $T=1000$ (20 $M_\odot$)")
#plt.savefig("PLOTS/yongplus_abundances_n12_T1000.pdf")
