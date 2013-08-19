from calcq import calcq
import numpy as np
import pylab as plt
from scipy.integrate import simps
from scipy.integrate import trapz
from scipy.optimize import brentq
from scipy.interpolate import UnivariateSpline
from math import copysign
from dust_calculations import Blambda, find_normalization, std_size_distr
from dust_calculations import Ikappa, kappaPlanck, dopckekappa
from dust_calculations import get_wavelen, Sgeometric
from kappa_gas import kappa_gas
from cooling import find_Td, find_Td_CMB_D, get_kappafn
from paper_comparedustmodels import Gamma_ad

if __name__ == "__main__":
    wavelen, logwavelen = get_wavelen()
    klambda,S = np.load("DATA/x5.5_UM-ND-20_klambda.npy")
    kappafn = get_kappafn(wavelen,klambda,Tsub=1500)

    D = 10.**(-7)

    #Tarr = np.array([10, 50, 100, 200, 500, 1000, 1500, 2000])
    Tarr = np.array([100, 200, 500, 1000, 1500, 2000])
    narr = 10**np.linspace(4,20,160)
    
    Tdarr = np.zeros((len(narr),len(Tarr)))
    Hdarr = np.zeros((len(narr),len(Tarr)))
    Ldarr = np.zeros((len(narr),len(Tarr)))

    Tdarr2 = np.zeros((len(narr),len(Tarr)))
    Hdarr2 = np.zeros((len(narr),len(Tarr)))
    Ldarr2 = np.zeros((len(narr),len(Tarr)))

    fig = plt.figure(1,figsize=(7.5,7.5))
    plt.plot(narr,Gamma_ad(narr,1000),'k',lw=1.5)
    plotstyle=['b','g','r','c','m','y','k','b']
    plotstyle1=[x+'' for x in plotstyle]
    plotstyle2=[x+':' for x in plotstyle]
    linelist = ['' for x in xrange(len(Tarr))]
    for i,T in enumerate(Tarr):
        for j,n in enumerate(narr):
            Tdarr[j,i],Hdarr[j,i],Ldarr[j,i] = find_Td_CMB_D(n,T,S,kappafn,50,1500,D)
            Tdarr2[j,i],Hdarr2[j,i],Ldarr2[j,i] = find_Td(n,T,S,kappafn)

    for i,T in enumerate(Tarr):
        linelist[i], = plt.plot(narr, D*Ldarr[:,i],plotstyle1[i],lw=2)
        plt.plot(narr, D*Ldarr2[:,i],plotstyle2[i],lw=2)

    plt.xlabel(r'n [cm $^{-3}$]',fontsize=16)
    plt.ylabel(r'$\Lambda_d$ [ergs s$^{-1}$ cm$^{-3}$]',fontsize=16)
    plt.gca().loglog()
    plt.gca().tick_params(axis='both',which='major',labelsize=16)
    plt.legend(linelist, ['T='+str(x)+'K' for x in Tarr],loc="upper left",prop={'size':18},frameon=False)
    #plt.gcf().set_size_inches(6,12)
    plt.savefig("PLOTS/paper_optthickthin.eps",bbox_inches="tight")
