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
#_narr = 10.**np.linspace(0,20,200)
#_Tarr = 10.**np.linspace(1.5,3.5,200)
_narrzoom = 10.**np.linspace(10,14,100)
_Tarrzoom = 10.**np.linspace(1.5,3.5,100)

def yes_dust_opacity(narr,Tarr,dustmodel,kappafn,S,Tcmb,Tsubl):
    logDarr = np.linspace(-10,-4,61)

    logDcritarr = np.zeros((len(narr),len(Tarr)))*np.nan
    logTdarr = np.zeros((len(narr),len(Tarr)))*np.nan
    for i,n in enumerate(narr):
        for j,T in enumerate(Tarr):
            heating = Gamma_ad(n,T)
            #Go from high Dcrit to low Dcrit
            for k,logD in enumerate(logDarr):
                D = 10.**logD
                Td,Hd,Ld = find_Td_CMB_D(n,T,S,kappafn,Tcmb,Tsubl,D)
                if Ld*D >= heating:
                    logDcritarr[i,j] = logD
                    logTdarr[i,j] = np.log10(Td)
                    break
    np.save("DATA/"+dustmodel+"_ntplane_dustopacityzoom.npy",(logDcritarr,logTdarr))

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

## def plot_ntplane(dustmodel,dust_opacity=False,make_colorbar=True,MsiMdust=None):
##     if dust_opacity:
##         Dcrit,Td = np.load("DATA/"+dustmodel+"_ntplane_dustopacityzoom.npy")
##         #mask = Dcrit > -5
##         #Dcrit[mask] = np.nan
##         #Td[mask] = np.nan
##     else:
##         Dcrit,Td = np.load("DATA/"+dustmodel+"_ntplane.npy")
##         mask = Dcrit > -3
##         Dcrit[mask] = np.nan
##         Td[mask] = np.nan

##     if make_colorbar:
##         fig = plt.figure(figsize=(8.5,7.5))
##     else:
##         fig = plt.figure(figsize=(7.5,7.5))
    
##     if MsiMdust == None:
##         img = plt.imshow(-1*Dcrit.transpose(),extent=(10,14,1.5,3.5),origin='lower',aspect=8.0)
##         img.set_clim(3,9)
##         if make_colorbar:
##             cbar = plt.colorbar(shrink=0.75,ticks=np.arange(3,10))
##             cbar.ax.set_yticklabels([r'$10^{-'+str(x)+'}$' for x in np.arange(3,10)])
##             cbar.ax.tick_params(axis='y',which='major',labelsize=16)
##     else:
##         Sicrit = Dcrit + np.log10(4./(3*28.1)) + np.log10(MsiMdust) +12-7.51
##         img = plt.imshow(Sicrit.transpose(),extent=(10,14,1.5,3.5),origin='lower',aspect=8.0)
##         #img.set_clim(3,9)
##         if make_colorbar:
##             cbar = plt.colorbar(shrink=0.75)#,ticks=np.arange(3,10))
##             #cbar.ax.set_yticklabels([r'$10^{-'+str(x)+'}$' for x in np.arange(3,10)])
##             cbar.ax.tick_params(axis='y',which='major',labelsize=16)
##     plt.plot([10,14],[np.log10(50),np.log10(50)],'k--')
##     plt.plot([10,14],[np.log10(1500),np.log10(1500)],'k--')
##     plt.plot([10,14],2./3+jeans(np.array([10.,14.])),'k--',lw=2)
##     plt.plot([10,14],jeans(np.array([10.,14.])),'k--',lw=2)
##     plt.plot([10,14],opacity(np.array([10.,14.])),'k-.')
##     #plt.plot([0,20],dustcooling(np.array([0.,20.]),np.log10(S),-6),'k:')
##     plt.plot([10,14],dustcooling(np.array([10.,14.]),np.log10(S),-7),'k:')
##     #plt.plot([0,20],[-1.7,4.0],'k--')
##     plt.xlabel(r'log n [cm$^{-3}$]',fontsize=16)
##     plt.ylabel(r'log T [K]',fontsize=16)
##     plt.gca().tick_params(axis='both',which='major',labelsize=16)
##     #plt.title(r'-log $\mathcal{D}_{crit}$ for '+dustmodel)

##     if MsiMdust == None:
##         plotname = "Dcritplane"
##         plt.xlim((4,18))
##         plt.ylim((1.5,3.5))
##         if dustmodel[0:4]=="std_":
##             plt.gca().text(14,1.6,"standard",fontsize=16)
##         elif dustmodel[0:5]=="x5.5_":
##             plt.gca().text(14,1.6,"shock",fontsize=16)
##     else:
##         plotname = "Sicritplane"
##         #plt.xlim((10,14))
##         plt.ylim((1.5,3.5))
##         if dustmodel[0:4]=="std_":
##             plt.gca().text(13,1.7,"standard",fontsize=16)
##         elif dustmodel[0:5]=="x5.5_":
##             plt.gca().text(13,1.7,"shock",fontsize=16)
##     if dust_opacity:
##         plt.savefig("PLOTS/paper_"+dustmodel+"_"+plotname+"_dustopacityzoom.pdf",bbox_inches='tight')
##     else:
##         plt.savefig("PLOTS/paper_"+dustmodel+"_"+plotname+"zoom.pdf",bbox_inches='tight')

##     fig = plt.figure(figsize=(8.5,7.5))
##     plt.imshow(Td.transpose(),extent=(0,20,1.5,3.5),origin='lower',aspect=8.0)
##     if make_colorbar:
##         plt.colorbar(shrink=0.75)
##     plt.plot([0,20],[np.log10(50),np.log10(50)],'k--')
##     plt.plot([0,20],[np.log10(1500),np.log10(1500)],'k--')
##     plt.plot([0,20],2./3+jeans(np.array([0.,20.])),'k--',lw=2)
##     plt.plot([0,20],jeans(np.array([0.,20.])),'k--',lw=2)
##     plt.plot([0,20],opacity(np.array([0.,20.])),'k-.')
##     #plt.plot([0,20],dustcooling(np.array([0.,20.]),np.log10(S),-6),'k:')
##     plt.plot([0,20],dustcooling(np.array([0.,20.]),np.log10(S),-7),'k:')
##     plt.xlabel(r'log $n$ [cm$^{-3}$]')
##     plt.ylabel(r'log $T$ [K]')
##     plt.title(r'log $T_d$ for '+dustmodel)
##     plt.xlim((12,20))
##     plt.ylim((1.7,3.5))

##     if dust_opacity:
##         plt.savefig("PLOTS/paper_"+dustmodel+"_Tdplane_dustopacity.eps",bbox_inches='tight')
##     else:
##         plt.savefig("PLOTS/paper_"+dustmodel+"_Tdplane.eps",bbox_inches='tight')

def plot_paperfig():
    dustmodels = ['x5.5_UM-ND-20', 'std_UM-ND-20']
    MsiMdust = .469
    fig = plt.figure(figsize=(7.5,7.5))

    Dcrit,Td = np.load("DATA/"+dustmodels[1]+"_ntplane_dustopacityzoom.npy")
    Sicrit = Dcrit + np.log10(4./(3*28.1)) + np.log10(MsiMdust) +12-7.51
    plt.subplot(121)
    img = plt.imshow(Sicrit.transpose(),extent=(10,14,1.5,3.5),origin='lower',aspect=4.0)
    img.set_clim(-6,-1)
    plt.plot([10,14],2./3+jeans(np.array([10.,14.])),'k--',lw=2)
    plt.plot([10,14],jeans(np.array([10.,14.])),'k--',lw=2)
    plt.plot([10,14],[np.log10(50),np.log10(50)],'k--')
    plt.plot([10,14],[np.log10(1500),np.log10(1500)],'k--')
    plt.gca().text(12,1.6,"standard",fontsize=16,ha='center')
    plt.xlabel(r'log n [cm$^{-3}$]',fontsize=16)
    plt.ylabel(r'log T [K]',fontsize=16)
    plt.xticks([10,12])
    plt.ylim([1.5,3.5])
    plt.gca().tick_params(axis='both',which='major',labelsize=16)

    Dcrit,Td = np.load("DATA/"+dustmodels[0]+"_ntplane_dustopacityzoom.npy")
    Sicrit = Dcrit + np.log10(4./(3*28.1)) + np.log10(MsiMdust) +12-7.51
    plt.subplot(122)
    img = plt.imshow(Sicrit.transpose(),extent=(10,14,1.5,3.5),origin='lower',aspect=4.0)
    img.set_clim(-6,-1)
    plt.plot([10,14],2./3+jeans(np.array([10.,14.])),'k--',lw=2)
    plt.plot([10,14],jeans(np.array([10.,14.])),'k--',lw=2)
    plt.plot([10,14],[np.log10(50),np.log10(50)],'k--')
    plt.plot([10,14],[np.log10(1500),np.log10(1500)],'k--')
    plt.gca().text(12,1.6,"shock",fontsize=16,ha='center')
    plt.xlabel(r'log n [cm$^{-3}$]',fontsize=16)
    plt.xticks([10,12,14])
    plt.ylim([1.5,3.5])
    plt.gca().tick_params(axis='x',which='major',labelsize=16)
    plt.gca().set_yticklabels(tuple(['','','','','']))

    fig.subplots_adjust(wspace=0.0,right=0.8)
    cbar_ax = fig.add_axes([.85,.15,.05,.7])
    fig.colorbar(img,cax=cbar_ax)
    plt.savefig("PLOTS/paper_ntplane_sih.pdf",bbox_inches='tight')

    #################################################

    Dcrit,Td = np.load("DATA/"+dustmodels[1]+"_ntplane_dustopacityzoom.npy")
    plt.subplot(121)
    img = plt.imshow(Dcrit.transpose(),extent=(10,14,1.5,3.5),origin='lower',aspect=4.0)
    img.set_clim(-9,-4)
    plt.plot([10,14],2./3+jeans(np.array([10.,14.])),'k--',lw=2)
    plt.plot([10,14],jeans(np.array([10.,14.])),'k--',lw=2)
    plt.plot([10,14],[np.log10(50),np.log10(50)],'k--')
    plt.plot([10,14],[np.log10(1500),np.log10(1500)],'k--')
    #plt.plot([10,14],opacity(np.array([10.,14.])),'k-.')
    #plt.plot([10,14],dustcooling(np.array([10.,14.]),np.log10(S),-7),'k:')
    plt.gca().text(12,1.6,"standard",fontsize=16,ha='center')
    plt.xlabel(r'log n [cm$^{-3}$]',fontsize=16)
    plt.ylabel(r'log T [K]',fontsize=16)
    plt.xticks([10,12])
    plt.ylim([1.5,3.5])
    plt.gca().tick_params(axis='both',which='major',labelsize=16)

    Dcrit,Td = np.load("DATA/"+dustmodels[0]+"_ntplane_dustopacityzoom.npy")
    plt.subplot(122)
    img = plt.imshow(Dcrit.transpose(),extent=(10,14,1.5,3.5),origin='lower',aspect=4.0)
    img.set_clim(-9,-4)
    plt.plot([10,14],2./3+jeans(np.array([10.,14.])),'k--',lw=2)
    plt.plot([10,14],jeans(np.array([10.,14.])),'k--',lw=2)
    plt.plot([10,14],[np.log10(50),np.log10(50)],'k--')
    plt.plot([10,14],[np.log10(1500),np.log10(1500)],'k--')
    #plt.plot([10,14],opacity(np.array([10.,14.])),'k-.')
    #plt.plot([10,14],dustcooling(np.array([10.,14.]),np.log10(S),-7),'k:')
    plt.gca().text(12,1.6,"shock",fontsize=16,ha='center')
    plt.xlabel(r'log n [cm$^{-3}$]',fontsize=16)
    plt.xticks([10,12,14])
    plt.ylim([1.5,3.5])
    plt.gca().tick_params(axis='x',which='major',labelsize=16)
    plt.gca().set_yticklabels(tuple(['','','','','']))

    fig.subplots_adjust(wspace=0.0,right=0.8)
    cbar_ax = fig.add_axes([.85,.15,.05,.7])
    fig.colorbar(img,cax=cbar_ax)
    plt.savefig("PLOTS/paper_ntplane_Dcrit.pdf",bbox_inches='tight')

    #cbar = plt.colorbar(shrink=0.75)#,ticks=np.arange(3,10))
    #cbar.ax.set_yticklabels([r'$10^{-'+str(x)+'}$' for x in np.arange(3,10)])
    #cbar.ax.tick_params(axis='y',which='major',labelsize=16)
    #plt.plot([10,14],[np.log10(50),np.log10(50)],'k--')
    #plt.plot([10,14],[np.log10(1500),np.log10(1500)],'k--')
    #plt.plot([10,14],2./3+jeans(np.array([10.,14.])),'k--',lw=2)
    #plt.plot([10,14],jeans(np.array([10.,14.])),'k--',lw=2)
    #plt.plot([10,14],opacity(np.array([10.,14.])),'k-.')
    #plt.plot([0,20],dustcooling(np.array([0.,20.]),np.log10(S),-6),'k:')
    #plt.plot([10,14],dustcooling(np.array([10.,14.]),np.log10(S),-7),'k:')
    #plt.plot([0,20],[-1.7,4.0],'k--')
    #plt.xlabel(r'log n [cm$^{-3}$]',fontsize=16)
    #plt.ylabel(r'log T [K]',fontsize=16)
    #plt.gca().tick_params(axis='both',which='major',labelsize=16)
    #plt.title(r'-log $\mathcal{D}_{crit}$ for '+dustmodel)

    #plotname = "Sicritplane"
    #plt.xlim((10,14))
    #plt.ylim((1.5,3.5))
    #if dustmodel[0:4]=="std_":
    #    plt.gca().text(13,1.7,"standard",fontsize=16)
    #elif dustmodel[0:5]=="x5.5_":
    #    plt.gca().text(13,1.7,"shock",fontsize=16)
    #if dust_opacity:
    #    plt.savefig("PLOTS/paper_"+dustmodel+"_"+plotname+"_dustopacityzoom.pdf",bbox_inches='tight')
    #else:
    #    plt.savefig("PLOTS/paper_"+dustmodel+"_"+plotname+"zoom.pdf",bbox_inches='tight')

if __name__ == "__main__":
    """
    Two ways of proceeding:
    a) ignore dust opacity when calculating optically thick gas
    b) pick multiple D, find when you break adiabatic cooling
    """
##     dustmodels = ['x5.5_UM-ND-20', 'std_UM-ND-20',
##                   'x5.5_UM-D-20',  'std_UM-D-20',
##                   'x5.5_UM-ND-170','std_UM-ND-170',
##                   'x5.5_UM-D-170', 'std_UM-D-170',
##                   'x5.5_M-D-20',   'std_M-D-20',
##                   'x5.5_M-ND-20',  'std_M-ND-20',
##                   'x5.5_M-D-170',  'std_M-D-170',
##                   'x5.5_M-ND-170', 'std_M-ND-170']
##                   'x5.5_SOIF06-CCSN','std_SOIF06-CCSN',
##                   'x5.5_SOIF06-PISN','std_SOIF06-PISN',
##                   'x5.5_Caff20','std_Caff20',
##                   'x5.5_Caff35','std_Caff35']
    dustmodels = ['x5.5_UM-ND-20', 'std_UM-ND-20']
    MsiMdust = .469

    wavelen,lw = get_wavelen()
    Tcmb = 50
    Tsubl = 1500
    for dustmodel in dustmodels:
        print dustmodel
        #klambda,S = np.load("DATA/"+dustmodel+"_klambda.npy")
        #kappafn = get_kappafn(wavelen,klambda,Tsub=Tsubl)
        
        #no_dust_opacity(_narr,_Tarr,dustmodel,kappafn,S,Tcmb,Tsubl)
        #plot_ntplane(dustmodel,dust_opacity=False)
        
        #yes_dust_opacity(_narrzoom,_Tarrzoom,dustmodel,kappafn,S,Tcmb,Tsubl)
        #plot_ntplane(dustmodel,dust_opacity=True,MsiMdust=MsiMdust)
    plot_paperfig()
