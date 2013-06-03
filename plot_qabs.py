import numpy as np
import pylab as plt
from calcq import calcq
from dust_calculations import std_size_distr, calc_qabs, get_wavelen, get_densities, shock_size_distr

if __name__ == "__main__":
    RECALC_QABS = True
    MAKEPLOTS = False

    wavelen, logwavelen = get_wavelen()
    ac_rho, al2o3_rho, asi_rho, fe_rho, fes_rho, mg2sio4_rho, mg_rho, mgo_rho, sio2_rho, fe3o4_rho, mgsio3_rho, fe2sio4_rho = get_densities()

    if RECALC_QABS:
##         ac_a, ac_dnda = std_size_distr(ac_rho)
        ac_n, ac_k = np.load("opticalconsts/ac.npy")
        ac_opt = np.array([wavelen, ac_n, ac_k]).T
##         ac_qabs= calc_qabs(ac_opt,ac_a,ac_dnda,'std_ac')
        
##         al2o3_a, al2o3_dnda = std_size_distr(al2o3_rho)
        al2o3_n, al2o3_k = np.load("opticalconsts/al2o3.npy")
        al2o3_opt = np.array([wavelen, al2o3_n, al2o3_k]).T
##         al2o3_qabs= calc_qabs(al2o3_opt,al2o3_a,al2o3_dnda,'std_al2o3')
        
##         asi_a, asi_dnda = std_size_distr(asi_rho)
        asi_n, asi_k = np.load("opticalconsts/asi.npy")
        asi_opt = np.array([wavelen, asi_n, asi_k]).T
##         asi_qabs= calc_qabs(asi_opt,asi_a,asi_dnda,'std_asi')
        
##         fe_a, fe_dnda = std_size_distr(fe_rho)
        fe_n, fe_k = np.load("opticalconsts/fe.npy")
        fe_opt = np.array([wavelen, fe_n, fe_k]).T
##         fe_qabs= calc_qabs(fe_opt,fe_a,fe_dnda,'std_fe')
        
##         fes_a, fes_dnda = std_size_distr(fes_rho)
        fes_n, fes_k = np.load("opticalconsts/fes.npy")
        fes_opt = np.array([wavelen, fes_n, fes_k]).T
##         fes_qabs= calc_qabs(fes_opt,fes_a,fes_dnda,'std_fes')
        
##         mg2sio4_a, mg2sio4_dnda = std_size_distr(mg2sio4_rho)
        mg2sio4_n, mg2sio4_k = np.load("opticalconsts/mg2sio4.npy")
        mg2sio4_opt = np.array([wavelen, mg2sio4_n, mg2sio4_k]).T
##         mg2sio4_qabs= calc_qabs(mg2sio4_opt,mg2sio4_a,mg2sio4_dnda,'std_mg2sio4')
        
##         mg_a, mg_dnda = std_size_distr(mg_rho)
        mg_n, mg_k = np.load("opticalconsts/mg.npy")
        mg_opt = np.array([wavelen, mg_n, mg_k]).T
##         mg_qabs= calc_qabs(mg_opt,mg_a,mg_dnda,'std_mg')
        
##         mgo_a, mgo_dnda = std_size_distr(mgo_rho)
        mgo_n, mgo_k = np.load("opticalconsts/mgo.npy")
        mgo_opt = np.array([wavelen, mgo_n, mgo_k]).T
##         mgo_qabs= calc_qabs(mgo_opt,mgo_a,mgo_dnda,'std_mgo')
        
##         sio2_a, sio2_dnda = std_size_distr(sio2_rho)
        sio2_n, sio2_k = np.load("opticalconsts/sio2.npy")
        sio2_opt = np.array([wavelen, sio2_n, sio2_k]).T
##         sio2_qabs= calc_qabs(sio2_opt,sio2_a,sio2_dnda,'std_sio2')
        
##         fe3o4_a, fe3o4_dnda = std_size_distr(fe3o4_rho)
        fe3o4_n, fe3o4_k = np.load("opticalconsts/fe3o4.npy")
        fe3o4_opt = np.array([wavelen, fe3o4_n, fe3o4_k]).T
##         fe3o4_qabs= calc_qabs(fe3o4_opt,fe3o4_a,fe3o4_dnda,'std_fe3o4')
        
##         mgsio3_a, mgsio3_dnda = std_size_distr(mgsio3_rho)
        mgsio3_n, mgsio3_k = np.load("opticalconsts/mgsio3.npy")
        mgsio3_opt = np.array([wavelen, mgsio3_n, mgsio3_k]).T
##         mgsio3_qabs= calc_qabs(mgsio3_opt,mgsio3_a,mgsio3_dnda,'std_mgsio3')
        
##         fe2sio4_a, fe2sio4_dnda = std_size_distr(fe2sio4_rho)
        fe2sio4_n, fe2sio4_k = np.load("opticalconsts/fe2sio4.npy")
        fe2sio4_opt = np.array([wavelen, fe2sio4_n, fe2sio4_k]).T
##         fe2sio4_qabs= calc_qabs(fe2sio4_opt,fe2sio4_a,fe2sio4_dnda,'std_fe2sio4')

##         ac_a, ac_dnda = shock_size_distr(ac_rho,5.5)
##         ac_qabs_55 = calc_qabs(ac_opt,ac_a,ac_dnda,'x5.5_ac')
##         al2o3_a, al2o3_dnda = shock_size_distr(al2o3_rho,5.5)
##         al2o3_qabs_55 = calc_qabs(al2o3_opt,al2o3_a,al2o3_dnda,'x5.5_al2o3')
##         asi_a, asi_dnda = shock_size_distr(asi_rho,5.5)
##         asi_qabs_55 = calc_qabs(asi_opt,asi_a,asi_dnda,'x5.5_asi')
##         fe_a, fe_dnda = shock_size_distr(fe_rho,5.5)
##         fe_qabs_55 = calc_qabs(fe_opt,fe_a,fe_dnda,'x5.5_fe')
##         fes_a, fes_dnda = shock_size_distr(fes_rho,5.5)
##         fes_qabs_55 = calc_qabs(fes_opt,fes_a,fes_dnda,'x5.5_fes')
##         mg2sio4_a, mg2sio4_dnda = shock_size_distr(mg2sio4_rho,5.5)
##         mg2sio4_qabs_55 = calc_qabs(mg2sio4_opt,mg2sio4_a,mg2sio4_dnda,'x5.5_mg2sio4')
##         mg_a, mg_dnda = shock_size_distr(mg_rho,5.5)
##         mg_qabs_55 = calc_qabs(mg_opt,mg_a,mg_dnda,'x5.5_mg')
##         mgo_a, mgo_dnda = shock_size_distr(mgo_rho,5.5)
##         mgo_qabs_55 = calc_qabs(mgo_opt,mgo_a,mgo_dnda,'x5.5_mgo')
##         sio2_a, sio2_dnda = shock_size_distr(sio2_rho,5.5)
##         sio2_qabs_55 = calc_qabs(sio2_opt,sio2_a,sio2_dnda,'x5.5_sio2')
##         fe3o4_a, fe3o4_dnda = shock_size_distr(fe3o4_rho,5.5)
##         fe3o4_qabs_55 = calc_qabs(fe3o4_opt,fe3o4_a,fe3o4_dnda,'x5.5_fe3o4')
##         mgsio3_a, mgsio3_dnda = shock_size_distr(mgsio3_rho,5.5)
##         mgsio3_qabs_55 = calc_qabs(mgsio3_opt,mgsio3_a,mgsio3_dnda,'x5.5_mgsio3')
##         fe2sio4_a, fe2sio4_dnda = shock_size_distr(fe2sio4_rho,5.5)
##         fe2sio4_qabs_55 = calc_qabs(fe2sio4_opt,fe2sio4_a,fe2sio4_dnda,'x5.5_fe2sio4')

##         ac_a, ac_dnda = shock_size_distr(ac_rho,3.5)
##         ac_qabs_35 = calc_qabs(ac_opt,ac_a,ac_dnda,'x3.5_ac')
##         al2o3_a, al2o3_dnda = shock_size_distr(al2o3_rho,3.5)
##         al2o3_qabs_35 = calc_qabs(al2o3_opt,al2o3_a,al2o3_dnda,'x3.5_al2o3')
##         asi_a, asi_dnda = shock_size_distr(asi_rho,3.5)
##         asi_qabs_35 = calc_qabs(asi_opt,asi_a,asi_dnda,'x3.5_asi')
##         fe_a, fe_dnda = shock_size_distr(fe_rho,3.5)
##         fe_qabs_35 = calc_qabs(fe_opt,fe_a,fe_dnda,'x3.5_fe')
##         fes_a, fes_dnda = shock_size_distr(fes_rho,3.5)
##         fes_qabs_35 = calc_qabs(fes_opt,fes_a,fes_dnda,'x3.5_fes')
##         mg2sio4_a, mg2sio4_dnda = shock_size_distr(mg2sio4_rho,3.5)
##         mg2sio4_qabs_35 = calc_qabs(mg2sio4_opt,mg2sio4_a,mg2sio4_dnda,'x3.5_mg2sio4')
##         mg_a, mg_dnda = shock_size_distr(mg_rho,3.5)
##         mg_qabs_35 = calc_qabs(mg_opt,mg_a,mg_dnda,'x3.5_mg')
##         mgo_a, mgo_dnda = shock_size_distr(mgo_rho,3.5)
##         mgo_qabs_35 = calc_qabs(mgo_opt,mgo_a,mgo_dnda,'x3.5_mgo')
##         sio2_a, sio2_dnda = shock_size_distr(sio2_rho,3.5)
##         sio2_qabs_35 = calc_qabs(sio2_opt,sio2_a,sio2_dnda,'x3.5_sio2')
##         fe3o4_a, fe3o4_dnda = shock_size_distr(fe3o4_rho,3.5)
##         fe3o4_qabs_35 = calc_qabs(fe3o4_opt,fe3o4_a,fe3o4_dnda,'x3.5_fe3o4')
##         mgsio3_a, mgsio3_dnda = shock_size_distr(mgsio3_rho,3.5)
##         mgsio3_qabs_35 = calc_qabs(mgsio3_opt,mgsio3_a,mgsio3_dnda,'x3.5_mgsio3')
##         fe2sio4_a, fe2sio4_dnda = shock_size_distr(fe2sio4_rho,3.5)
##         fe2sio4_qabs_35 = calc_qabs(fe2sio4_opt,fe2sio4_a,fe2sio4_dnda,'x3.5_fe2sio4')

##         ac_a, ac_dnda = shock_size_distr(ac_rho,7.5)
##         ac_qabs_75 = calc_qabs(ac_opt,ac_a,ac_dnda,'x7.5_ac')
##         al2o3_a, al2o3_dnda = shock_size_distr(al2o3_rho,7.5)
##         al2o3_qabs_75 = calc_qabs(al2o3_opt,al2o3_a,al2o3_dnda,'x7.5_al2o3')
##         asi_a, asi_dnda = shock_size_distr(asi_rho,7.5)
##         asi_qabs_75 = calc_qabs(asi_opt,asi_a,asi_dnda,'x7.5_asi')
##         fe_a, fe_dnda = shock_size_distr(fe_rho,7.5)
##         fe_qabs_75 = calc_qabs(fe_opt,fe_a,fe_dnda,'x7.5_fe')
##         fes_a, fes_dnda = shock_size_distr(fes_rho,7.5)
##         fes_qabs_75 = calc_qabs(fes_opt,fes_a,fes_dnda,'x7.5_fes')
##         mg2sio4_a, mg2sio4_dnda = shock_size_distr(mg2sio4_rho,7.5)
##         mg2sio4_qabs_75 = calc_qabs(mg2sio4_opt,mg2sio4_a,mg2sio4_dnda,'x7.5_mg2sio4')
##         mg_a, mg_dnda = shock_size_distr(mg_rho,7.5)
##         mg_qabs_75 = calc_qabs(mg_opt,mg_a,mg_dnda,'x7.5_mg')
##         mgo_a, mgo_dnda = shock_size_distr(mgo_rho,7.5)
##         mgo_qabs_75 = calc_qabs(mgo_opt,mgo_a,mgo_dnda,'x7.5_mgo')
##         sio2_a, sio2_dnda = shock_size_distr(sio2_rho,7.5)
##         sio2_qabs_75 = calc_qabs(sio2_opt,sio2_a,sio2_dnda,'x7.5_sio2')
##         fe3o4_a, fe3o4_dnda = shock_size_distr(fe3o4_rho,7.5)
##         fe3o4_qabs_75 = calc_qabs(fe3o4_opt,fe3o4_a,fe3o4_dnda,'x7.5_fe3o4')
##         mgsio3_a, mgsio3_dnda = shock_size_distr(mgsio3_rho,7.5)
##         mgsio3_qabs_75 = calc_qabs(mgsio3_opt,mgsio3_a,mgsio3_dnda,'x7.5_mgsio3')
##         fe2sio4_a, fe2sio4_dnda = shock_size_distr(fe2sio4_rho,7.5)
##         fe2sio4_qabs_75 = calc_qabs(fe2sio4_opt,fe2sio4_a,fe2sio4_dnda,'x7.5_fe2sio4')

        ac_a, ac_dnda = shock_size_distr(ac_rho,4.5)
        ac_qabs_45 = calc_qabs(ac_opt,ac_a,ac_dnda,'x4.5_ac')
        al2o3_a, al2o3_dnda = shock_size_distr(al2o3_rho,4.5)
        al2o3_qabs_45 = calc_qabs(al2o3_opt,al2o3_a,al2o3_dnda,'x4.5_al2o3')
        asi_a, asi_dnda = shock_size_distr(asi_rho,4.5)
        asi_qabs_45 = calc_qabs(asi_opt,asi_a,asi_dnda,'x4.5_asi')
        fe_a, fe_dnda = shock_size_distr(fe_rho,4.5)
        fe_qabs_45 = calc_qabs(fe_opt,fe_a,fe_dnda,'x4.5_fe')
        fes_a, fes_dnda = shock_size_distr(fes_rho,4.5)
        fes_qabs_45 = calc_qabs(fes_opt,fes_a,fes_dnda,'x4.5_fes')
        mg2sio4_a, mg2sio4_dnda = shock_size_distr(mg2sio4_rho,4.5)
        mg2sio4_qabs_45 = calc_qabs(mg2sio4_opt,mg2sio4_a,mg2sio4_dnda,'x4.5_mg2sio4')
        mg_a, mg_dnda = shock_size_distr(mg_rho,4.5)
        mg_qabs_45 = calc_qabs(mg_opt,mg_a,mg_dnda,'x4.5_mg')
        mgo_a, mgo_dnda = shock_size_distr(mgo_rho,4.5)
        mgo_qabs_45 = calc_qabs(mgo_opt,mgo_a,mgo_dnda,'x4.5_mgo')
        sio2_a, sio2_dnda = shock_size_distr(sio2_rho,4.5)
        sio2_qabs_45 = calc_qabs(sio2_opt,sio2_a,sio2_dnda,'x4.5_sio2')
        fe3o4_a, fe3o4_dnda = shock_size_distr(fe3o4_rho,4.5)
        fe3o4_qabs_45 = calc_qabs(fe3o4_opt,fe3o4_a,fe3o4_dnda,'x4.5_fe3o4')
        mgsio3_a, mgsio3_dnda = shock_size_distr(mgsio3_rho,4.5)
        mgsio3_qabs_45 = calc_qabs(mgsio3_opt,mgsio3_a,mgsio3_dnda,'x4.5_mgsio3')
        fe2sio4_a, fe2sio4_dnda = shock_size_distr(fe2sio4_rho,4.5)
        fe2sio4_qabs_45 = calc_qabs(fe2sio4_opt,fe2sio4_a,fe2sio4_dnda,'x4.5_fe2sio4')
    else:
        ac_qabs = np.load("DATA/std_ac_qabs.npy")
        al2o3_qabs = np.load("DATA/std_al2o3_qabs.npy")
        asi_qabs = np.load("DATA/std_asi_qabs.npy")
        fe_qabs = np.load("DATA/std_fe_qabs.npy")
        fes_qabs = np.load("DATA/std_fes_qabs.npy")
        mg2sio4_qabs = np.load("DATA/std_mg2sio4_qabs.npy")
        mg_qabs = np.load("DATA/std_mg_qabs.npy")
        mgo_qabs = np.load("DATA/std_mgo_qabs.npy")
        sio2_qabs = np.load("DATA/std_sio2_qabs.npy")
        fe3o4_qabs = np.load("DATA/std_fe3o4_qabs.npy")
        mgsio3_qabs = np.load("DATA/std_mgsio3_qabs.npy")
        fe2sio4_qabs = np.load("DATA/std_fe2sio4_qabs.npy")

        ac_qabs_55 = np.load("DATA/x5.5_ac_qabs.npy")
        al2o3_qabs_55 = np.load("DATA/x5.5_al2o3_qabs.npy")
        asi_qabs_55 = np.load("DATA/x5.5_asi_qabs.npy")
        fe_qabs_55 = np.load("DATA/x5.5_fe_qabs.npy")
        fes_qabs_55 = np.load("DATA/x5.5_fes_qabs.npy")
        mg2sio4_qabs_55 = np.load("DATA/x5.5_mg2sio4_qabs.npy")
        mg_qabs_55 = np.load("DATA/x5.5_mg_qabs.npy")
        mgo_qabs_55 = np.load("DATA/x5.5_mgo_qabs.npy")
        sio2_qabs_55 = np.load("DATA/x5.5_sio2_qabs.npy")
        fe3o4_qabs_55 = np.load("DATA/x5.5_fe3o4_qabs.npy")
        mgsio3_qabs_55 = np.load("DATA/x5.5_mgsio3_qabs.npy")
        fe2sio4_qabs_55 = np.load("DATA/x5.5_fe2sio4_qabs.npy")
    
        ac_qabs_35 = np.load("DATA/x3.5_ac_qabs.npy")
        al2o3_qabs_35 = np.load("DATA/x3.5_al2o3_qabs.npy")
        asi_qabs_35 = np.load("DATA/x3.5_asi_qabs.npy")
        fe_qabs_35 = np.load("DATA/x3.5_fe_qabs.npy")
        fes_qabs_35 = np.load("DATA/x3.5_fes_qabs.npy")
        mg2sio4_qabs_35 = np.load("DATA/x3.5_mg2sio4_qabs.npy")
        mg_qabs_35 = np.load("DATA/x3.5_mg_qabs.npy")
        mgo_qabs_35 = np.load("DATA/x3.5_mgo_qabs.npy")
        sio2_qabs_35 = np.load("DATA/x3.5_sio2_qabs.npy")
        fe3o4_qabs_35 = np.load("DATA/x3.5_fe3o4_qabs.npy")
        mgsio3_qabs_35 = np.load("DATA/x3.5_mgsio3_qabs.npy")
        fe2sio4_qabs_35 = np.load("DATA/x3.5_fe2sio4_qabs.npy")

        ac_qabs_75 = np.load("DATA/x7.5_ac_qabs.npy")
        al2o3_qabs_75 = np.load("DATA/x7.5_al2o3_qabs.npy")
        asi_qabs_75 = np.load("DATA/x7.5_asi_qabs.npy")
        fe_qabs_75 = np.load("DATA/x7.5_fe_qabs.npy")
        fes_qabs_75 = np.load("DATA/x7.5_fes_qabs.npy")
        mg2sio4_qabs_75 = np.load("DATA/x7.5_mg2sio4_qabs.npy")
        mg_qabs_75 = np.load("DATA/x7.5_mg_qabs.npy")
        mgo_qabs_75 = np.load("DATA/x7.5_mgo_qabs.npy")
        sio2_qabs_75 = np.load("DATA/x7.5_sio2_qabs.npy")
        fe3o4_qabs_75 = np.load("DATA/x7.5_fe3o4_qabs.npy")
        mgsio3_qabs_75 = np.load("DATA/x7.5_mgsio3_qabs.npy")
        fe2sio4_qabs_75 = np.load("DATA/x7.5_fe2sio4_qabs.npy")

        ac_qabs_45 = np.load("DATA/x4.5_ac_qabs.npy")
        al2o3_qabs_45 = np.load("DATA/x4.5_al2o3_qabs.npy")
        asi_qabs_45 = np.load("DATA/x4.5_asi_qabs.npy")
        fe_qabs_45 = np.load("DATA/x4.5_fe_qabs.npy")
        fes_qabs_45 = np.load("DATA/x4.5_fes_qabs.npy")
        mg2sio4_qabs_45 = np.load("DATA/x4.5_mg2sio4_qabs.npy")
        mg_qabs_45 = np.load("DATA/x4.5_mg_qabs.npy")
        mgo_qabs_45 = np.load("DATA/x4.5_mgo_qabs.npy")
        sio2_qabs_45 = np.load("DATA/x4.5_sio2_qabs.npy")
        fe3o4_qabs_45 = np.load("DATA/x4.5_fe3o4_qabs.npy")
        mgsio3_qabs_45 = np.load("DATA/x4.5_mgsio3_qabs.npy")
        fe2sio4_qabs_45 = np.load("DATA/x4.5_fe2sio4_qabs.npy")

    ### Plots
    if MAKEPLOTS:
        which_a = [2,4,6,8]
    
        plt.subplot(6,2,1)
        for i in which_a:
            plt.plot(wavelen,ac_qabs[:,i])
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'amorphous C $\rho=$'+str(ac_rho))
    
        plt.subplot(6,2,2)
        for i in which_a:
            plt.plot(wavelen,al2o3_qabs[:,i])
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'$Al_2O_3$ $\rho=$'+str(al2o3_rho))
    
        plt.subplot(6,2,3)
        for i in which_a:
            plt.plot(wavelen,asi_qabs[:,i])
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'amorphous Si $\rho=$'+str(asi_rho))
    
        plt.subplot(6,2,4)
        for i in which_a:
            plt.plot(wavelen,fe_qabs[:,i])
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'Fe $\rho=$'+str(fe_rho))
    
        plt.subplot(6,2,5)
        for i in which_a:
            plt.plot(wavelen,fes_qabs[:,i])
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'FeS $\rho=$'+str(fes_rho))
    
        plt.subplot(6,2,6)
        for i in which_a:
            plt.plot(wavelen,mg2sio4_qabs[:,i])
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'$Mg_2SiO_4$ $\rho=$'+str(mg2sio4_rho))
    
        plt.subplot(6,2,7)
        for i in which_a:
            plt.plot(wavelen,mg_qabs[:,i])
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'Mg $\rho=$'+str(mg_rho))
    
        plt.subplot(6,2,8)
        for i in which_a:
            plt.plot(wavelen,mgo_qabs[:,i])
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'MgO $\rho=$'+str(mgo_rho))
    
        plt.subplot(6,2,9)
        for i in which_a:
            plt.plot(wavelen,sio2_qabs[:,i])
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'$SiO_2$ $\rho=$'+str(sio2_rho))
    
        plt.subplot(6,2,10)
        for i in which_a:
            plt.plot(wavelen,fe3o4_qabs[:,i])
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'$Fe_3O_4$ $\rho=$'+str(fe3o4_rho))
    
        plt.subplot(6,2,11)
        for i in which_a:
            plt.plot(wavelen,mgsio3_qabs[:,i])
        plt.xlabel(r'wavelength ($\mu$m)')
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'$MgSiO_3$ $\rho=$'+str(mgsio3_rho))
    
        plt.legend(('1e-2','1e-1','1e0','3e0'),loc="upper right")
    
        plt.subplot(6,2,12)
        for i in which_a:
            plt.plot(wavelen,fe2sio4_qabs[:,i])
        plt.xlabel(r'wavelength ($\mu$m)')
        plt.ylabel(r'$Q_{abs}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.ylim((10**-11, 10**1))
        plt.title(r'$Fe_2SiO_4$ $\rho=$'+str(fe2sio4_rho))
    
        plt.gcf().set_size_inches(12,28)
        plt.savefig("PLOTS/qabs_std.pdf")

        ######################
        def plot_qabs(wavelen,qabslist,filename):
            ac_rho, al2o3_rho, asi_rho, fe_rho, fes_rho, mg2sio4_rho, mg_rho, mgo_rho, sio2_rho, fe3o4_rho, mgsio3_rho, fe2sio4_rho = get_densities()
            which_a = [2,4,6,8]
            plt.clf()
            titles = [r'amorphous C $\rho=$'+str(ac_rho),
                      r'$Al_2O_3$ $\rho=$'+str(al2o3_rho),
                      r'amorphous Si $\rho=$'+str(asi_rho),
                      r'Fe $\rho=$'+str(fe_rho),
                      r'FeS $\rho=$'+str(fes_rho),
                      r'$Mg_2SiO_4$ $\rho=$'+str(mg2sio4_rho),
                      r'Mg $\rho=$'+str(mg_rho),
                      r'MgO $\rho=$'+str(mgo_rho),
                      r'$SiO_2$ $\rho=$'+str(sio2_rho),
                      r'$Fe_3O_4$ $\rho=$'+str(fe3o4_rho),
                      r'$MgSiO_3$ $\rho=$'+str(mgsio3_rho),
                      r'$Fe_2SiO_4$ $\rho=$'+str(fe2sio4_rho)]

            for j in xrange(12):
                plt.subplot(6,2,j+1)
                thisqabs = qabslist[j]
                for i in which_a:
                    plt.plot(wavelen,thisqabs[:,i])
                plt.ylabel(r'$Q_{abs}$')
                plt.gca().set_xscale('log')
                plt.gca().set_yscale('log')
                plt.ylim((10**-11, 10**1))
                plt.title(titles[j])
                if j+1==11:
                    plt.legend(('1e-2','1e-1','1e0','3e0'),loc="upper right")
    
            plt.gcf().set_size_inches(12,28)
            plt.savefig(filename)

        plot_qabs(wavelen,
                  [ac_qabs_35,al2o3_qabs_35,asi_qabs_35,fe_qabs_35,fes_qabs_35,
                   mg2sio4_qabs_35,mg_qabs_35,mgo_qabs_35,sio2_qabs_35,
                   fe3o4_qabs_35,mgsio3_qabs_35,fe2sio4_qabs_35],
                  "PLOTS/qabs_x3.5.pdf")

        plot_qabs(wavelen,
                  [ac_qabs_45,al2o3_qabs_45,asi_qabs_45,fe_qabs_45,fes_qabs_45,
                   mg2sio4_qabs_45,mg_qabs_45,mgo_qabs_45,sio2_qabs_45,
                   fe3o4_qabs_45,mgsio3_qabs_45,fe2sio4_qabs_45],
                  "PLOTS/qabs_x4.5.pdf")

        plot_qabs(wavelen,
                  [ac_qabs_55,al2o3_qabs_55,asi_qabs_55,fe_qabs_55,fes_qabs_55,
                   mg2sio4_qabs_55,mg_qabs_55,mgo_qabs_55,sio2_qabs_55,
                   fe3o4_qabs_55,mgsio3_qabs_55,fe2sio4_qabs_55],
                  "PLOTS/qabs_x5.5.pdf")

        plot_qabs(wavelen,
                  [ac_qabs_75,al2o3_qabs_75,asi_qabs_75,fe_qabs_75,fes_qabs_75,
                   mg2sio4_qabs_75,mg_qabs_75,mgo_qabs_75,sio2_qabs_75,
                   fe3o4_qabs_75,mgsio3_qabs_75,fe2sio4_qabs_75],
                  "PLOTS/qabs_x7.5.pdf")
