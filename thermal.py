import numpy as np
import pylab as plt
from scipy.integrate import odeint

from dust_calculations import get_wavelen, dopckekappa, powerkappa
from cooling import get_kappafn
from cooling import find_Td, find_Td_D

if __name__ == "__main__":
    n0 = 10**10 #cm^(-3)
    T0 = 500 #K
    klambda,Sgeom = np.load("DATA/UM-ND-20_klambda.npy")
    klambda2,Sgeom2 = np.load("DATA/SOIF06-PISN_klambda.npy")
    print "Heating: adiabatic contraction heating"
    print "Cooling: dust with different models (hard coded and changing)"
    
    wavelen, logwavelen = get_wavelen()
    narr = 10**np.linspace(np.log10(n0),16)
    kappafn = get_kappafn(wavelen,klambda,Tsub=1500)
    #kappafn2 = get_kappafn(wavelen,klambda2,Tsub=2000)
    #kappafn2 = dopckekappa
    #kappafn2 = lambda Td: dopckekappa(Td,k0=400)
    kappafn2 = lambda Td: powerkappa(Td,k0=400)

    ## Physical Constants
    kB = 1.38 * 10**(-16)
    mu = 1. + 4.0*0.083 #1 + 4yHe where yHe = nHe/nH
    mP = 1.67e-24 #g
    G  = 6.67e-8  #cm^3/(g s^2)
    
    ## analytic solution to dn/dt = n/t_ff
    #c0 = 2.0/np.sqrt(n0) #integration constant
    A = 6.156e-16 # cm^(3/2)/sec
    #def n_ff(t):
    #    return np.log10(4)-2*np.log10(np.abs(A*t - c0))
    
    plotstyle = ['b','g','r','c','m']
    def Gamma(T,n): #only adiabatic heating
        return 1.5*n*kB*T/(np.sqrt(3*np.pi/(8*G*mu*mP*n)))
    for i,D in enumerate([1e-5,1e-6,1e-7,1e-8,1e-9]):
        def Lambda(T,n):
            #Td, Hd, l = find_Td(n,T,Sgeom,kappafn)
            Td, Hd, l = find_Td_D(n,T,Sgeom,kappafn,D)
            l = D*l
            return l
        def dTdn(T,n):
            return (Gamma(T,n)-Lambda(T,n))/(1.5*A*kB*n**2.5)

        Tarr = odeint(dTdn, T0, narr)
        Larr = np.array([Lambda(Tarr[j],narr[j]) for j in xrange(len(narr))])
        plt.subplot(211)
        plt.plot(narr,Tarr,plotstyle[i],lw=.5)
        plt.subplot(212)
        plt.plot(narr,Larr,plotstyle[i],lw=.5)
    
    plotstyle = ['b:','g:','r:','c:','m:']
    for i,D in enumerate([1e-5,1e-6,1e-7,1e-8,1e-9]):
        def Lambda(T,n):
            #Td, Hd, l = find_Td(n,T,Sgeom2,kappafn2)
            Td, Hd, l = find_Td_D(n,T,Sgeom2,kappafn2,D)
            l = D*l
            return l
        def dTdn(T,n):
            return (Gamma(T,n)-Lambda(T,n))/(1.5*A*kB*n**2.5)

        Tarr = odeint(dTdn, T0, narr)
        Larr = np.array([Lambda(Tarr[j],narr[j]) for j in xrange(len(narr))])
        plt.subplot(211)
        plt.plot(narr,Tarr,plotstyle[i],lw=2.0)
        plt.subplot(212)
        plt.plot(narr,Larr,plotstyle[i],lw=2.0)

    plt.subplot(211)
    plt.plot(narr,1500+np.zeros(len(narr)),"k:")
    #plt.plot(narr,2000+np.zeros(len(narr)),"k:")
    plt.gca().loglog()
    plt.xlabel(r'n [cm$^{-3}$]')
    plt.ylabel(r'T [K]')
    plt.legend(('-5','-6','-7','-8','-9'),loc="upper left")

    plt.subplot(212)
    plt.gca().loglog()
    plt.xlabel(r'n [cm$^{-3}$]')
    plt.ylabel(r'$\Lambda$ [erg s$^{-1}$ cm$^{-3}$]')
    
    plt.gcf().set_size_inches(8,11)
    plt.savefig("PLOTS/thermal.pdf")
