import numpy as np
import pylab as plt

#from cooling import find_Td_D

from dust_calculations import get_wavelen
from cooling import get_kappafn
from cooling import find_Td_CMB_D
from environmentcooling import Gamma_ad
from matplotlib import cm

from plot_schematic import dustcooling,jeans,opacity

## Global vars: the density and temperature range
_narr = 10.**np.linspace(0,20,200)
_Tarr = 10.**np.linspace(1.5,3.5,200)

## def yes_dust_opacity(narr,Tarr,dustmodel,kappafn,S,Tcmb,Tsubl):
##     logDarr = np.linspace(-10,-3,71)

##     logDcritarr = np.zeros((len(narr),len(Tarr)))*np.nan
##     logTdarr = np.zeros((len(narr),len(Tarr)))*np.nan
##     for i,n in enumerate(narr):
##         for j,T in enumerate(Tarr):
##             heating = Gamma_ad(n,T)
##             #Go from high Dcrit to low Dcrit
##             for k,logD in enumerate(logDarr):
##                 D = 10.**logD
##                 Td,Hd,Ld = find_Td_CMB_D(n,T,S,kappafn,Tcmb,Tsubl,D)
##                 if Ld*D >= heating:
##                     logDcritarr[i,j] = logD
##                     logTdarr[i,j] = np.log10(Td)
##                     break
##     np.save("DATA/"+dustmodel+"_ntplane_dustopacity.npy",(logDcritarr,logTdarr))

## def no_dust_opacity(narr,Tarr,dustmodel,kappafn,S,Tcmb,Tsubl):
##     logDcritarr = np.zeros((len(narr),len(Tarr)))
##     logTdarr = np.zeros((len(narr),len(Tarr)))
##     for i,n in enumerate(narr):
##         for j,T in enumerate(Tarr):
##             heating = Gamma_ad(n,T)
##             ## ignore dust self-opacity in finding Td
##             Td,Hd,Ld = find_Td_CMB_D(n,T,S,kappafn,Tcmb,Tsubl,0.0)
##             if Ld <= 0: ## e.g. if dust is sublimated
##                 logDcritarr[i,j] = 15 #very large
##                 logTdarr[i,j] = 0
##             else:
##                 logDcritarr[i,j] = np.log10(heating/Ld)
##                 logTdarr[i,j] = np.log10(Td)
##     np.save("DATA/"+dustmodel+"_ntplane.npy",(logDcritarr,logTdarr))

def plot_ntplane(dustmodel,dust_opacity=False):
    if dust_opacity:
        Dcrit,Td = np.load("DATA/"+dustmodel+"_ntplane_dustopacity.npy")
        #mask = Dcrit > -5
        #Dcrit[mask] = np.nan
        #Td[mask] = np.nan
    else:
        Dcrit,Td = np.load("DATA/"+dustmodel+"_ntplane.npy")
        mask = Dcrit > -3
        Dcrit[mask] = np.nan
        Td[mask] = np.nan

    fig = plt.figure(figsize=(8.5,7.5))
    img = plt.imshow(-1*Dcrit.transpose(),extent=(0,20,1.5,3.5),origin='lower',aspect=8.0)
    img.set_clim(3,9)
    cbar = plt.colorbar(shrink=0.75,ticks=np.arange(3,10))#,ticks=[-15,-12,-9,-6,-3,0,3,6,8])
    cbar.ax.set_yticklabels([r'$10^{-'+str(x)+'}$' for x in np.arange(3,10)])
    cbar.ax.tick_params(axis='y',which='major',labelsize=16)
    plt.plot([0,20],[np.log10(50),np.log10(50)],'k--')
    plt.plot([0,20],[np.log10(1500),np.log10(1500)],'k--')
    plt.plot([0,20],2./3+jeans(np.array([0.,20.])),'k--')
    plt.plot([0,20],jeans(np.array([0.,20.])),'k--')
    plt.plot([0,20],opacity(np.array([0.,20.])),'k--')
    plt.plot([0,20],dustcooling(np.array([0.,20.]),np.log10(S),-6),'k--')
    plt.xlabel(r'log n [cm$^{-3}$]',fontsize=16)
    plt.ylabel(r'log T [K]',fontsize=16)
    plt.gca().tick_params(axis='both',which='major',labelsize=16)
    #plt.title(r'-log $\mathcal{D}_{crit}$ for '+dustmodel)
    plt.xlim((4,18))
    plt.ylim((1.5,3.5))
    if dust_opacity:
        plt.savefig("PLOTS/paper_"+dustmodel+"_Dcritplane_dustopacity.pdf",bbox_inches='tight')
    else:
        plt.savefig("PLOTS/paper_"+dustmodel+"_Dcritplane.pdf",bbox_inches='tight')

    fig = plt.figure(figsize=(8.5,7.5))
    plt.imshow(Td.transpose(),extent=(0,20,1.5,3.5),origin='lower',aspect=8.0)
    plt.colorbar(shrink=0.75)
    plt.plot([0,20],[np.log10(50),np.log10(50)],'k--')
    plt.plot([0,20],[np.log10(1500),np.log10(1500)],'k--')
    plt.plot([0,20],2./3+jeans(np.array([0.,20.])),'k--')
    plt.plot([0,20],jeans(np.array([0.,20.])),'k--')
    plt.plot([0,20],opacity(np.array([0.,20.])),'k--')
    plt.plot([0,20],dustcooling(np.array([0.,20.]),np.log10(S),-6),'k--')
    plt.xlabel(r'log $n$ [cm$^{-3}$]')
    plt.ylabel(r'log $T$ [K]')
    plt.title(r'log $T_d$ for '+dustmodel)
    plt.xlim((0,20))
    plt.ylim((1.5,3.5))
    if dust_opacity:
        plt.savefig("PLOTS/paper_"+dustmodel+"_Tdplane_dustopacity.pdf",bbox_inches='tight')
    else:
        plt.savefig("PLOTS/paper_"+dustmodel+"_Tdplane.pdf",bbox_inches='tight')

if __name__ == "__main__":
    """
    Two ways of proceeding:
    a) ignore dust opacity when calculating optically thick gas
    b) pick multiple D, find when you break adiabatic cooling
    """
    dustmodels = ['x5.5_UM-ND-20', 'std_UM-ND-20',
                  'x5.5_UM-D-20',  'std_UM-D-20',
                  'x5.5_UM-ND-170','std_UM-ND-170',
                  'x5.5_UM-D-170', 'std_UM-D-170',
                  'x5.5_M-D-20',   'std_M-D-20',
                  'x5.5_M-ND-20',  'std_M-ND-20',
                  'x5.5_M-D-170',  'std_M-D-170',
                  'x5.5_M-ND-170', 'std_M-ND-170',
                  'x5.5_SOIF06-CCSN','std_SOIF06-CCSN',
                  'x5.5_SOIF06-PISN','std_SOIF06-PISN',
                  'x5.5_Caff20','std_Caff20',
                  'x5.5_Caff35','std_Caff35']

    wavelen,lw = get_wavelen()
    Tcmb = 50
    Tsubl = 1500
    for dustmodel in dustmodels:                 
        print dustmodel
        klambda,S = np.load("DATA/"+dustmodel+"_klambda.npy")
        kappafn = get_kappafn(wavelen,klambda,Tsub=Tsubl)
        
        #no_dust_opacity(_narr,_Tarr,dustmodel,kappafn,S,Tcmb,Tsubl)
        #plot_ntplane(dustmodel,dust_opacity=False)
        
        #yes_dust_opacity(_narr,_Tarr,dustmodel,kappafn,S,Tcmb,Tsubl)
        plot_ntplane(dustmodel,dust_opacity=True)
