from calcq import calcq
import numpy as np
import pylab as plt
from scipy.integrate import simps
from scipy.integrate import trapz
from scipy.optimize import brentq
from scipy.interpolate import UnivariateSpline
from dust_calculations import Blambda, find_normalization, std_size_distr
from dust_calculations import Ikappa, kappaPlanck, dopckekappa
from dust_calculations import get_wavelen, Sgeometric
from kappa_gas import kappa_gas

from math import copysign

import os.path
import sys

def get_kappafn(wavelen,klambda,Tsub=1500,kfloor=1e-5):
    wavelen_cm = wavelen/10**4
    Tarr = 10**np.linspace(0,4,200)
    kParr = np.zeros(len(Tarr))
    for i,T in enumerate(Tarr):
        kParr[i] = kappaPlanck(T, wavelen_cm, klambda)
    kappaspl = UnivariateSpline(Tarr,kParr,s=0)
    def kappafn(x):
        if x > Tsub:
            return 0#kfloor
        else:
            return kappaspl(x)
    return kappafn

def plot_kappafn(wavelen, klambda):
    wavelen_cm = wavelen/10**4
    kappafn = get_kappafn(wavelen,klambda)
    Tarr = 10**np.linspace(0,4,200)
    kParr = np.zeros(len(Tarr))
    for i,T in enumerate(Tarr):
        kParr[i] = kappaPlanck(T, wavelen_cm, klambda)
    plt.plot(Tarr,kParr,'o',markersize=2)
    plt.plot(Tarr,[kappafn(x) for x in Tarr],'r-')
    plt.xlabel('T')
    plt.ylabel(r'\kappa_P')
    plt.gca().loglog()
    plt.savefig("PLOTS/check_kappafn.pdf")

def dust_cooling(D,n,T,kappafn,S):
    Td, Hd_over_D, Ld_over_D = find_Td(n, T, S, kappafn)
    return D*Ld_over_D #_dust_cooling(n,T,kappafn,S)

def _Hd(Td, n, T, S):
    CH = np.sqrt(8*1.38e-16*1.67e-24/np.pi) * (1+4*.083) * (.3536 + .5*.083) * 2 * 1.38e-16
    return CH * S * n**2 * np.sqrt(T) * (T-Td)
def _beta(kappa, n, T):
    return 1

def _Lambda(Td, n, T, kappafn):
    CL = 4*5.67e-5*(1+4*.083)*1.67e-24
    kappa = kappafn(Td)
    return CL * n * Td**4 * kappa * _beta(kappa, n, T)
def _Lambda_CMB(Td, n, T, kappafn, Tcmb):
    CL = 4*5.67e-5*(1+4*.083)*1.67e-24
    kappa = kappafn(Td)
    return CL * n * (Td**4-Tcmb**4) * kappa * _beta(kappa, n, T)# * copysign(1,Td-Tcmb)

def find_Td(n, T, S, kappafn):
    myHd = lambda Td: _Hd(Td, n, T, S)
    myLambda = lambda Td: _Lambda(Td, n, T, kappafn)
    Td = brentq(lambda x: myHd(x)-myLambda(x), 0, T+100)
    return Td, myHd(Td), myLambda(Td)

def _beta_all(kdust, n, T, D):
    kgas = kappa_gas(n,T)
    kappa = (D*kdust + kgas)
    if kappa <= 0:
        return 1
    return min(1,(6.67e-8)/(kappa**2 * np.pi * 1.38e-16 * n * T))
##     kappa = kdust
##     kappacrit = np.sqrt(6.67e-8/(np.pi*1.38e-16*n*T))
##     if kappa < kappacrit:
##         return 1
##     else:
##         return (kappacrit/kappa)**2
def _Lambda_D(Td,n,T,kappafn,D):
    CL = 4*5.67e-5*(1+4*.083)*1.67e-24
    kappa = kappafn(Td)
    return CL * n * Td**4 * kappa * _beta_all(kappa, n, T, D)
def find_Td_D(n,T,S,kappafn,D):
    myHd = lambda Td: _Hd(Td, n, T, S)
    myLambda = lambda Td: _Lambda_D(Td, n, T, kappafn,D)
    Td = brentq(lambda x: myHd(x)-myLambda(x), 0, T+100)
    return Td, myHd(Td), myLambda(Td)
def _Lambda_CMB_D(Td, n, T, kappafn, Tcmb, D):
    CL = 4*5.67e-5*(1+4*.083)*1.67e-24
    kappa = kappafn(Td)
    return CL * n * (Td**4-Tcmb**4) * kappa * _beta_all(kappa, n, T, D)
    #return CL * n * (Td-Tcmb)**4 * kappa * _beta_all(kappa, n, T, D) * copysign(1,Td-Tcmb)
def find_Td_CMB_D(n,T,S,kappafn,Tcmb,Tsubl,D):
    myHd = lambda Td: _Hd(Td, n, T, S)
    myLambda = lambda Td: _Lambda_CMB_D(Td, n, T, kappafn,Tcmb,D)
    if myHd(Tsubl) > myLambda(Tsubl):
        return np.nan, np.nan, np.nan
    if T < Tcmb:
        return Tcmb, myHd(Tcmb), myLambda(Tcmb)
    try:
        Td = brentq(lambda x: myHd(x)-myLambda(x), Tcmb, Tsubl)
        return Td, myHd(Td), myLambda(Td)
    except ValueError:
        return np.nan, np.nan, np.nan
##     try:
##         #Td = brentq(lambda x: myHd(x)-myLambda(x), 0, T+100)
##         #Td = brentq(lambda x: myHd(x)-myLambda(x), Tcmb, T-1)
##         Td = brentq(lambda x: myHd(x)-myLambda(x), Tcmb, Tsubl)
##     except ValueError:
##         if T < Tcmb:
##             return Tcmb,myHd(Tcmb),myLambda(Tcmb)
##         else:
##             return np.nan,np.nan,np.nan
##     return Td, myHd(Td), myLambda(Td)

if __name__ == "__main__":
    wavelen, logwavelen = get_wavelen()
    klambda,S = np.load("DATA/SOIF06-PISN_klambda.npy")
    kappafn = get_kappafn(wavelen,klambda,Tsub=2000)

    Dcrit = 4.4*10**-9.
    Tdarr = 10**np.linspace(0,4,400)
    Ldarr100 = [Dcrit*_Lambda_D(Td,10**4,100,kappafn,Dcrit) for Td in Tdarr]
    Ldarr200 = [Dcrit*_Lambda_D(Td,10**4,200,kappafn,Dcrit) for Td in Tdarr]
    Ldarr500 = [Dcrit*_Lambda_D(Td,10**4,500,kappafn,Dcrit) for Td in Tdarr]
    Ldarr1000= [Dcrit*_Lambda_D(Td,10**4,1000,kappafn,Dcrit) for Td in Tdarr]
    Hdarr100 = [Dcrit*_Hd(Td,10**4,100,S) for Td in Tdarr]
    Hdarr200 = [Dcrit*_Hd(Td,10**4,200,S) for Td in Tdarr]
    Hdarr500 = [Dcrit*_Hd(Td,10**4,500,S) for Td in Tdarr]
    Hdarr1000= [Dcrit*_Hd(Td,10**4,1000,S) for Td in Tdarr]

    linearr = [0,0,0,0]
    plt.subplot(211)
    linearr[0], = plt.plot(Tdarr, Ldarr100,'b')
    plt.plot(Tdarr, Hdarr100,'b:')
    linearr[1], = plt.plot(Tdarr, Ldarr200,'g')
    plt.plot(Tdarr, Hdarr200,'g:')
    linearr[2], = plt.plot(Tdarr, Ldarr500,'r')
    plt.plot(Tdarr, Hdarr500,'r:')
    linearr[3], = plt.plot(Tdarr, Ldarr1000,'c')
    plt.plot(Tdarr, Hdarr1000,'c:')
    plt.legend(linearr,['100','200','500','1000'],loc="lower right")
    #plt.xlabel(r'T$_d$ [K]')
    plt.ylabel(r'$\Lambda_d$, $H_d$')
    plt.gca().loglog()

    Ldarr4 = [Dcrit*_Lambda_D(Td,10**4,500,kappafn,Dcrit) for Td in Tdarr]
    Ldarr6 = [Dcrit*_Lambda_D(Td,10**6,500,kappafn,Dcrit) for Td in Tdarr]
    Ldarr8 = [Dcrit*_Lambda_D(Td,10**8,500,kappafn,Dcrit) for Td in Tdarr]
    Ldarr10= [Dcrit*_Lambda_D(Td,10**10,500,kappafn,Dcrit) for Td in Tdarr]
    Hdarr4 = [Dcrit*_Hd(Td,10**4,500,S) for Td in Tdarr]
    Hdarr6 = [Dcrit*_Hd(Td,10**6,500,S) for Td in Tdarr]
    Hdarr8 = [Dcrit*_Hd(Td,10**8,500,S) for Td in Tdarr]
    Hdarr10= [Dcrit*_Hd(Td,10**10,500,S) for Td in Tdarr]
    plt.subplot(212)
    linearr[0], = plt.plot(Tdarr, Ldarr4,'b')
    plt.plot(Tdarr, Hdarr4,'b:')
    linearr[1], = plt.plot(Tdarr, Ldarr6,'g')
    plt.plot(Tdarr, Hdarr6,'g:')
    linearr[2], = plt.plot(Tdarr, Ldarr8,'r')
    plt.plot(Tdarr, Hdarr8,'r:')
    linearr[3], = plt.plot(Tdarr, Ldarr10,'c')
    plt.plot(Tdarr, Hdarr10,'c:')
    plt.legend(linearr,['4','6','8','10'],loc="lower right")
    plt.xlabel(r'T$_d$ [K]')
    plt.ylabel(r'$\Lambda_d$, $H_d$')
    plt.gca().loglog()
    plt.savefig("PLOTS/check_intersection.pdf")
    
    Tarr = np.array([10, 50, 100, 200, 500, 1000, 1500, 2000])
    #Tarr = np.array([100, 200, 500, 1000, 1500, 2000])
    narr = 10**np.linspace(0,20,200)
    
    Tdarr = np.zeros((len(narr),len(Tarr)))
    Hdarr = np.zeros((len(narr),len(Tarr)))
    Ldarr = np.zeros((len(narr),len(Tarr)))

    Tdarr2 = np.zeros((len(narr),len(Tarr)))
    Hdarr2 = np.zeros((len(narr),len(Tarr)))
    Ldarr2 = np.zeros((len(narr),len(Tarr)))

    plotstyle=['b','g','r','c','m','y','k','b']
    plotstyle2=[x+':' for x in plotstyle]
    linelist = ['' for x in xrange(len(Tarr))]
    plt.clf()
    plt.subplot(211)
    for i,T in enumerate(Tarr):
        for j,n in enumerate(narr):
            Tdarr[j,i],Hdarr[j,i],Ldarr[j,i] = find_Td_D(n,T,S,kappafn,Dcrit)
            #Tdarr[j,i],Hdarr[j,i],Ldarr[j,i] = find_Td_D(n,T,S,dopckekappa,Dcrit)
            Tdarr2[j,i],Hdarr2[j,i],Ldarr2[j,i] = find_Td(n,T,S,kappafn)
            #Tdarr2[j,i],Hdarr2[j,i],Ldarr2[j,i] = find_Td_D(n,T,S,kappafn,Dcrit*100.)
        linelist[i], = plt.plot(narr, Tdarr[:,i],plotstyle[i],lw=1)
        plt.plot(narr,Tdarr2[:,i],plotstyle2[i],lw=2)
        #plt.plot(narr, T+np.zeros(len(narr)), plotstyle2[i])

    plt.xlabel(r'n [cm $^{-3}$]')
    plt.ylabel(r'T$_d$ [K]')
    plt.title(r'S+06 PISN, Numerical $T_d$')
    plt.gca().loglog()

    plt.subplot(212)
    for i,T in enumerate(Tarr):
        plt.plot(narr, Dcrit*Ldarr[:,i],plotstyle[i],lw=1)
        plt.plot(narr, Dcrit*Ldarr2[:,i],plotstyle2[i],lw=2)
    plt.xlabel(r'n [cm $^{-3}$]')
    plt.ylabel(r'H$_d$,L$_d$ [ergs s$^{-1}$ cm$^{-3}$]')
    plt.title(r'S+06 PISN, Dust Cooling (Schneider+12 Dcrit, thin)')
    plt.gca().loglog()
    plt.legend(linelist, ['T='+str(x) for x in Tarr],loc="upper left",ncol=2,prop={'size':8},frameon=False)
    plt.gcf().set_size_inches(6,12)
    plt.savefig("PLOTS/check_Td_LdHd.pdf")
