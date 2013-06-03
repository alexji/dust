import numpy as np
import pylab as plt

from dust_calculations import get_wavelen, dopckekappa
from cooling import get_kappafn, find_Td, find_Td_D

from scipy.optimize import curve_fit

def _powerkappafit(Td,Tflat,Tsubl,k0):
    if Td < Tflat:
        return k0*(Td/Tflat)**2
    elif Td < Tsubl:
        return k0
    else:
        return k0*(Td/Tsubl)**(-12) #0
def powerkappafit(Td,Tflat,Tsubl,k0):
    try:
        return [_powerkappafit(T,Tflat,Tsubl,k0) for T in Td]
    except TypeError:
        return _powerkappafit(Td,Tflat,Tsubl,k0)
def logpowerkappafit(Td,Tflat,Tsubl,k0):
    try:
        return [np.log10(_powerkappafit(T,Tflat,Tsubl,k0)) for T in Td]
    except TypeError:
        return np.log10(_powerkappafit(Td,Tflat,Tsubl,k0))

if __name__ == "__main__":
    wavelen, logwavelen = get_wavelen()
    dustmodels = ['SOIF06-PISN','UM-D-20','UM-ND-20','UM-D-170','UM-ND-170']
    Tsubarr = [2000,1500,1500,1500,1500]
    fig1 = plt.figure()
    fig1.subplots_adjust(wspace=0.3, hspace=0.3)
    fig2 = plt.figure()
    fig2.subplots_adjust(wspace=0.3, hspace=0.3)
    for idust,dustmodel in enumerate(dustmodels):
        print dustmodels[idust]
        ##Load 'true' kappa
        klambda,S = np.load("DATA/"+dustmodel+"_klambda.npy")
        kappatrue = get_kappafn(wavelen,klambda,Tsub=Tsubarr[idust])
        ##Fit powerkappa
        fitTarr = 10**np.linspace(0,4,200)
        fitkParr = np.array([kappatrue(Td) for Td in fitTarr])
        fitTarrcut  = fitTarr[fitkParr > 0]
        fitkParrcut = fitkParr[fitkParr > 0]
        ## Fit on log kappa, but not log T (because that shouldn't matter too much)
        popt, pcov = curve_fit(logpowerkappafit, fitTarrcut, np.log10(fitkParrcut),p0=[200,1500,500])
        kappapow = lambda Td: powerkappafit(Td, popt[0], popt[1], popt[2])
        paramstr=r'$T_{flat}='+str(int(round(popt[0])))+r'$, $T_{subl}='+str(int(round(popt[1])))+r'$, $\kappa_0='+str(int(round(popt[2])))+'$'
        print popt

        ##Plot Kappa
        ax = fig1.add_subplot(3,2,idust)
        ax.plot(fitTarr,fitkParr,'b')
        ax.plot(fitTarr,kappapow(fitTarr),'b:')
        ax.set_xlabel(r'$T_d$ [K]')
        ax.set_ylabel(r'$\kappa_P$ [cm$^2$ g$^{-1}$]')
        ax.set_ylim((10**-2,10**3))
        ax.set_title(dustmodel)
        ax.loglog()
        ax.text(10,.03,paramstr,fontsize=8)

        ##Plot cooling
        ax = fig2.add_subplot(3,2,idust)
        Tarr = np.array([10, 50, 100, 200, 500, 1000, 1500, 2000])
        narr = 10**np.linspace(0,20,200)
        trueTdarr = np.zeros(len(narr))
        trueHdarr = np.zeros(len(narr))
        trueLdarr = np.zeros(len(narr))
        powTdarr  = np.zeros(len(narr))
        powHdarr  = np.zeros(len(narr))
        powLdarr  = np.zeros(len(narr))
        plotstyle = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
        plotstyle2 =[x+'--' for x in plotstyle]
        linelist = ['' for x in xrange(len(Tarr))]
        for i,T in enumerate(Tarr):
            for j,n in enumerate(narr):
                trueTdarr[j],trueHdarr[j],trueLdarr[j] = find_Td(n,T,S,kappatrue)
                powTdarr[j], powHdarr[j], powLdarr[j]  = find_Td(n,T,S,kappapow)
            linelist[i], = ax.plot(narr, trueLdarr,plotstyle[i],lw=.6)
            ax.plot(narr, powLdarr,plotstyle2[i],lw=1.2)
        ax.set_xlabel(r'n [cm $^{-3}$]')
        ax.set_ylabel(r'$L_d/D$ [ergs s$^{-1}$ cm$^{-3}$]')
        ax.set_title(dustmodel)
        ax.loglog()
    fig2.legend(linelist, ['T='+str(x) for x in Tarr],loc="lower left",ncol=2,prop={'size':8},frameon=False)
    fig1.set_size_inches(10,12)
    fig1.savefig("PLOTS/powerlaw-kappa.pdf")
    fig2.set_size_inches(10,12)
    fig2.savefig("PLOTS/powerlaw-lambda.pdf")
##     ## Load 'true' kappafn
##     wavelen,logwavelen = get_wavelen()
##     klambda,S = np.load("DATA/SOIF06-PISN_klambda.npy")
##     kappatrue = get_kappafn(wavelen,klambda,Tsub=2000)

##     ## Fit powerkappa
##     fitTarr = 10**np.linspace(0,4,200)
##     fitkParr = np.array([kappatrue(Td) for Td in fitTarr])
##     fitTarr  = fitTarr[fitkParr > 0]
##     fitkParr = fitkParr[fitkParr > 0]

##     popt, pcov = curve_fit(logpowerkappafit, fitTarr, np.log10(fitkParr),p0=[200,1500,500])
##     kappapow = lambda Td: powerkappafit(Td, popt[0], popt[1], popt[2])
##     #print popt
##     #print pcov

##     plt.subplot(211)
##     plt.plot(fitTarr,fitkParr,'b')
##     plt.plot(fitTarr,kappapow(fitTarr),'b:')
##     plt.xlabel(r'$T_d$ [K]')
##     plt.ylabel(r'$\kappa_P$ [cm$^2$ g$^{-1}$]')
##     plt.gca().loglog()

##     ## Calculate Lambda(n,T)
##     plt.subplot(212)
##     Tarr = np.array([10, 50, 100, 200, 500, 1000, 1500, 2000])
##     narr = 10**np.linspace(0,20,200)

##     trueTdarr = np.zeros((len(narr),len(Tarr)))
##     trueHdarr = np.zeros((len(narr),len(Tarr)))
##     trueLdarr = np.zeros((len(narr),len(Tarr)))

##     powTdarr = np.zeros((len(narr),len(Tarr)))
##     powHdarr = np.zeros((len(narr),len(Tarr)))
##     powLdarr = np.zeros((len(narr),len(Tarr)))

##     plotstyle = ['b','g','r','c','m','y','k','b']
##     plotstyle2 =[x+'--' for x in plotstyle]
##     linelist = ['' for x in xrange(len(Tarr))]
##     for i,T in enumerate(Tarr):
##         for j,n in enumerate(narr):
##             trueTdarr[j,i],trueHdarr[j,i],trueLdarr[j,i] = find_Td(n,T,S,kappatrue)
##             powTdarr[j,i], powHdarr[j,i], powLdarr[j,i]  = find_Td(n,T,S,kappapow)
##         linelist[i], = plt.plot(narr, trueLdarr[:,i],plotstyle[i],lw=.6)
##         plt.plot(narr, powLdarr[:,i],plotstyle2[i],lw=1.2)

##     plt.xlabel(r'n [cm $^{-3}$]')
##     plt.ylabel(r'$L_d/D$ [ergs s$^{-1}$ cm$^{-3}$]')
##     plt.gca().loglog()
##     plt.ylim((10**-18,10**6))
##     plt.legend(linelist, ['T='+str(x) for x in Tarr],loc="upper left",ncol=2,prop={'size':8},frameon=False)
##     plt.gcf().set_size_inches(6,12)
##     plt.savefig("PLOTS/powerlaw.pdf")
