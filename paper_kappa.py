import numpy as np
import pylab as plt
from dust_calculations import calc_klambda, Blambda, std_size_distr, shock_size_distr
from dust_calculations import get_wavelen, get_densities, get_mass_fractions
from dust_calculations import kappaPlanck, Sgeometric
import csv

if __name__ == "__main__":
    wavelen,lw = get_wavelen()
    wavelen_cm = wavelen/10**4
    def makeplot(prefix):
        plt.clf()
        r = csv.reader(open("DATA/schneider_ccsn_kappa.txt",'r'))
        header = r.next()
        schneiderccsn = np.array([(float(x),float(y)) for x,y in r])
        r = csv.reader(open("DATA/schneider2012_norev_kappa.txt",'r'))
        header = r.next()
        schneidernorev = np.array([(float(x),float(y)) for x,y in r])

        klambdaUMD20,SgeomUMD20 = np.load("DATA/"+prefix+"UM-D-20_klambda.npy")
        klambdaUMND20,SgeomUMND20= np.load("DATA/"+prefix+"UM-ND-20_klambda.npy")
        klambdaUMD170,SgeomUMD170 = np.load("DATA/"+prefix+"UM-D-170_klambda.npy")
        klambdaUMND170,SgeomUMND170= np.load("DATA/"+prefix+"UM-ND-170_klambda.npy")
        klambdaMD20,SgeomMD20 = np.load("DATA/"+prefix+"M-D-20_klambda.npy")
        klambdaMND20,SgeomMND20= np.load("DATA/"+prefix+"M-ND-20_klambda.npy")
        klambdaMD170,SgeomMD170 = np.load("DATA/"+prefix+"M-D-170_klambda.npy")
        klambdaMND170,SgeomMND170= np.load("DATA/"+prefix+"M-ND-170_klambda.npy")
##     klambdaSOIFPISN,SgeomSOIFPISN= np.load("DATA/"+prefix+"SOIF06-PISN_klambda.npy")
##         klambdaSOIFCCSN,SgeomSOIFCCSN= np.load("DATA/"+prefix+"SOIF06-CCSN_klambda.npy")
##     klambdaCaff20,SgeomCaff20 = np.load("DATA/"+prefix+"Caff20_klambda.npy")
##     klambdaCaff35,SgeomCaff35 = np.load("DATA/"+prefix+"Caff35_klambda.npy")

        Tarr = 10**np.linspace(1,4,150)
        kPlanckUMD20_arr = np.zeros(len(Tarr))
        kPlanckUMND20_arr= np.zeros(len(Tarr))
        kPlanckUMD170_arr = np.zeros(len(Tarr))
        kPlanckUMND170_arr= np.zeros(len(Tarr))
        kPlanckMD20_arr = np.zeros(len(Tarr))
        kPlanckMND20_arr= np.zeros(len(Tarr))
        kPlanckMD170_arr = np.zeros(len(Tarr))
        kPlanckMND170_arr= np.zeros(len(Tarr))
##     kPlanckSOIFPISN_arr= np.zeros(len(Tarr))
        kPlanckSOIFCCSN_arr= np.zeros(len(Tarr))
##     kPlanckCaff20_arr= np.zeros(len(Tarr))
##     kPlanckCaff35_arr= np.zeros(len(Tarr))
        for i,T in enumerate(Tarr):
            kPlanckUMD20_arr[i] = kappaPlanck(T, wavelen_cm, klambdaUMD20)
            kPlanckUMND20_arr[i]= kappaPlanck(T, wavelen_cm, klambdaUMND20)
            kPlanckUMD170_arr[i] = kappaPlanck(T, wavelen_cm, klambdaUMD170)
            kPlanckUMND170_arr[i]= kappaPlanck(T, wavelen_cm, klambdaUMND170)
            kPlanckMD20_arr[i] = kappaPlanck(T, wavelen_cm, klambdaMD20)
            kPlanckMND20_arr[i]= kappaPlanck(T, wavelen_cm, klambdaMND20)
            kPlanckMD170_arr[i] = kappaPlanck(T, wavelen_cm, klambdaMD170)
            kPlanckMND170_arr[i]= kappaPlanck(T, wavelen_cm, klambdaMND170)
##         kPlanckSOIFPISN_arr[i] = kappaPlanck(T, wavelen_cm, klambdaSOIFPISN)
##         kPlanckSOIFCCSN_arr[i] = kappaPlanck(T, wavelen_cm, klambdaSOIFCCSN)
##         kPlanckCaff20_arr[i] = kappaPlanck(T, wavelen_cm, klambdaCaff20)
##         kPlanckCaff35_arr[i] = kappaPlanck(T, wavelen_cm, klambdaCaff35)
        
        plt.plot(Tarr,kPlanckUMD20_arr)
        plt.plot(Tarr,kPlanckUMND20_arr)
        plt.plot(Tarr,kPlanckMD20_arr)
        plt.plot(Tarr,kPlanckMND20_arr)
        plt.plot(Tarr,kPlanckUMD170_arr)
        plt.plot(Tarr,kPlanckUMND170_arr)
        plt.plot(Tarr,kPlanckMD170_arr)
        plt.plot(Tarr,kPlanckMND170_arr)
##         plt.plot(Tarr,kPlanckSOIFPISN_arr)
##         plt.plot(Tarr,kPlanckSOIFCCSN_arr)
        plt.plot(schneiderccsn[:,0], schneiderccsn[:,1],'k:')
        plt.plot(schneidernorev[:,0],schneidernorev[:,1],'k-.')
        plt.plot([1500,1500],[.1,10**5],'k--')
        plt.xlim((10,10**4))
        plt.ylim((.1,10**4))
        plt.xlabel(r'$T_{dust}$ [K]',fontsize=16)
        plt.ylabel(r'$\kappa_P$ [cm$^2$/g]',fontsize=16)
        #plt.legend([r'UM D 20M$_\odot$',r'UM ND 20M$_\odot$',r'M D 20M$_\odot$',r'M ND 20M$_\odot$',r'S+06 CCSN',r'Caff20',r'Caff35'],loc="lower right")
        plt.gca().loglog()
        plt.gca().tick_params(axis='both',which='major',labelsize=16)
        plt.gcf().set_size_inches(7.5,7.5)
        
        if prefix=="std_":
            plt.text(13,4000,"standard",fontsize=20)
        if prefix=="x5.5_":
            plt.text(13,4000,"shock",fontsize=20)

        plt.savefig("PLOTS/paper_"+prefix+"kPlanck.pdf",bbox_inches='tight')
    makeplot("std_")
    makeplot("x5.5_")
