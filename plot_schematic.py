import numpy as np
import pylab as plt

##These are logT as a function of logn
def dustcooling(logn, logS, logD):
    return 9.24 - logn - 2*logS - 2*logD
def jeans(logn):
    return logn/3.0 - 1.57
def opacity(logn):
    return -3.0*logn + 51.1 + 3

if __name__ == "__main__":
    Tsubl = 1500
    Tcmb  = 50
    
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    
    lognarr = np.linspace(0,20)
    
    ax.plot(lognarr,dustcooling(lognarr,5.5,-8),'b')
    ax.plot(lognarr,dustcooling(lognarr,5.5,-7),'b--')
    ax.plot(lognarr,jeans(lognarr),'r')
    ax.plot(lognarr,2/3.+jeans(lognarr),'r--')
    ax.plot(lognarr,opacity(lognarr),'g')
    ax.plot(lognarr,opacity(lognarr)+7,'g--')
    ax.plot(lognarr,np.zeros(len(lognarr))+np.log10(Tcmb),'m')
    ax.plot(lognarr,np.zeros(len(lognarr))+np.log10(Tsubl),'m')
    
    ax.scatter([8],[3],marker='*', s=36,c='blue', edgecolor='blue')
    ax.scatter([10],[3],marker='*', s=36,c='green', edgecolor='green')
    ax.scatter([12],[3],marker='*', s=36,c='yellow', edgecolor='orange')
    
    plt.xlim((7,20))
    plt.ylim((1.5,3.5))
    plt.xlabel(r'log n [cm$^{-3}$]')
    plt.ylabel(r'log T [K]')
    
    plt.savefig("PLOTS/schematic.pdf")
