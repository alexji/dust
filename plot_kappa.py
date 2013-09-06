import numpy as np
import pylab as plt
from dust_calculations import calc_klambda, Blambda, std_size_distr, shock_size_distr
from dust_calculations import get_wavelen, get_densities, get_mass_fractions
from dust_calculations import kappaPlanck, Sgeometric

if __name__ == "__main__":
    RECALC_KAPPALAMBDA = False
    RECALC_DUSTMODELS = True
    MAKEPLOTS_KAPPALAMBDA = True
    MAKEPLOTS_KAPPAPLANCK = True
    # Pick ONE of these
    STANDARD = True
    SHOCK35  = False
    SHOCK55  = False
    SHOCK75  = False

    wavelen, logwavelen = get_wavelen()
    wavelen_cm = wavelen/10**4
    ac_rho, al2o3_rho, asi_rho, fe_rho, fes_rho, mg2sio4_rho, mg_rho, mgo_rho, sio2_rho, fe3o4_rho, mgsio3_rho, fe2sio4_rho = get_densities()

    if RECALC_KAPPALAMBDA:
        ## Dummy opt, since the only thing you need is the wavelength
        opt = np.array([wavelen, wavelen, wavelen]).T

        if STANDARD:
            prefix="std_"
            ac_a, ac_dnda = std_size_distr(ac_rho)
            ac_klambda= calc_klambda(opt,ac_a,ac_dnda,'std_ac')
            al2o3_a, al2o3_dnda = std_size_distr(al2o3_rho)
            al2o3_klambda= calc_klambda(opt,al2o3_a,al2o3_dnda,'std_al2o3')
            asi_a, asi_dnda = std_size_distr(asi_rho)
            asi_klambda= calc_klambda(opt,asi_a,asi_dnda,'std_asi')
            fe_a, fe_dnda = std_size_distr(fe_rho)
            fe_klambda= calc_klambda(opt,fe_a,fe_dnda,'std_fe')
            fes_a, fes_dnda = std_size_distr(fes_rho)
            fes_klambda= calc_klambda(opt,fes_a,fes_dnda,'std_fes')
            mg2sio4_a, mg2sio4_dnda = std_size_distr(mg2sio4_rho)
            mg2sio4_klambda= calc_klambda(opt,mg2sio4_a,mg2sio4_dnda,'std_mg2sio4')
            mg_a, mg_dnda = std_size_distr(mg_rho)
            mg_klambda= calc_klambda(opt,mg_a,mg_dnda,'std_mg')
            mgo_a, mgo_dnda = std_size_distr(mgo_rho)
            mgo_klambda= calc_klambda(opt,mgo_a,mgo_dnda,'std_mgo')
            sio2_a, sio2_dnda = std_size_distr(sio2_rho)
            sio2_klambda= calc_klambda(opt,sio2_a,sio2_dnda,'std_sio2')
            fe3o4_a, fe3o4_dnda = std_size_distr(fe3o4_rho)
            fe3o4_klambda= calc_klambda(opt,fe3o4_a,fe3o4_dnda,'std_fe3o4')
            mgsio3_a, mgsio3_dnda = std_size_distr(mgsio3_rho)
            mgsio3_klambda= calc_klambda(opt,mgsio3_a,mgsio3_dnda,'std_mgsio3')
            fe2sio4_a, fe2sio4_dnda = std_size_distr(fe2sio4_rho)
            fe2sio4_klambda= calc_klambda(opt,fe2sio4_a,fe2sio4_dnda,'std_fe2sio4')
        else:
            if SHOCK35:
                prefix="x3.5_"
                x = 3.5
            if SHOCK55:
                prefix="x5.5_"
                x = 5.5
            if SHOCK75:
                prefix="x7.5_"
                x = 7.5
            ac_a, ac_dnda = shock_size_distr(ac_rho,x)
            ac_klambda= calc_klambda(opt,ac_a,ac_dnda,prefix+'ac')
            al2o3_a, al2o3_dnda = shock_size_distr(al2o3_rho,x)
            al2o3_klambda= calc_klambda(opt,al2o3_a,al2o3_dnda,prefix+'al2o3')
            asi_a, asi_dnda = shock_size_distr(asi_rho,x)
            asi_klambda= calc_klambda(opt,asi_a,asi_dnda,prefix+'asi')
            fe_a, fe_dnda = shock_size_distr(fe_rho,x)
            fe_klambda= calc_klambda(opt,fe_a,fe_dnda,prefix+'fe')
            fes_a, fes_dnda = shock_size_distr(fes_rho,x)
            fes_klambda= calc_klambda(opt,fes_a,fes_dnda,prefix+'fes')
            mg2sio4_a, mg2sio4_dnda = shock_size_distr(mg2sio4_rho,x)
            mg2sio4_klambda= calc_klambda(opt,mg2sio4_a,mg2sio4_dnda,prefix+'mg2sio4')
            mg_a, mg_dnda = shock_size_distr(mg_rho,x)
            mg_klambda= calc_klambda(opt,mg_a,mg_dnda,prefix+'mg')
            mgo_a, mgo_dnda = shock_size_distr(mgo_rho,x)
            mgo_klambda= calc_klambda(opt,mgo_a,mgo_dnda,prefix+'mgo')
            sio2_a, sio2_dnda = shock_size_distr(sio2_rho,x)
            sio2_klambda= calc_klambda(opt,sio2_a,sio2_dnda,prefix+'sio2')
            fe3o4_a, fe3o4_dnda = shock_size_distr(fe3o4_rho,x)
            fe3o4_klambda= calc_klambda(opt,fe3o4_a,fe3o4_dnda,prefix+'fe3o4')
            mgsio3_a, mgsio3_dnda = shock_size_distr(mgsio3_rho,x)
            mgsio3_klambda= calc_klambda(opt,mgsio3_a,mgsio3_dnda,prefix+'mgsio3')
            fe2sio4_a, fe2sio4_dnda = shock_size_distr(fe2sio4_rho,x)
            fe2sio4_klambda= calc_klambda(opt,fe2sio4_a,fe2sio4_dnda,prefix+'fe2sio4')
    else:
        if STANDARD:
            prefix="std_"
            ac_a, ac_dnda = std_size_distr(ac_rho)
            al2o3_a, al2o3_dnda = std_size_distr(al2o3_rho)
            asi_a, asi_dnda = std_size_distr(asi_rho)
            fe_a, fe_dnda = std_size_distr(fe_rho)
            fes_a, fes_dnda = std_size_distr(fes_rho)
            mg2sio4_a, mg2sio4_dnda = std_size_distr(mg2sio4_rho)
            mg_a, mg_dnda = std_size_distr(mg_rho)
            mgo_a, mgo_dnda = std_size_distr(mgo_rho)
            sio2_a, sio2_dnda = std_size_distr(sio2_rho)
            fe3o4_a, fe3o4_dnda = std_size_distr(fe3o4_rho)
            mgsio3_a, mgsio3_dnda = std_size_distr(mgsio3_rho)
            fe2sio4_a, fe2sio4_dnda = std_size_distr(fe2sio4_rho)
        else:
            if SHOCK35:
                prefix="x3.5_"
                x=3.5
            elif SHOCK55:
                prefix="x5.5_"
                x=5.5
            elif SHOCK75:
                prefix="x7.5_"
                x=7.5
            ac_a, ac_dnda = shock_size_distr(ac_rho,x)
            al2o3_a, al2o3_dnda = shock_size_distr(al2o3_rho,x)
            asi_a, asi_dnda = shock_size_distr(asi_rho,x)
            fe_a, fe_dnda = shock_size_distr(fe_rho,x)
            fes_a, fes_dnda = shock_size_distr(fes_rho,x)
            mg2sio4_a, mg2sio4_dnda = shock_size_distr(mg2sio4_rho,x)
            mg_a, mg_dnda = shock_size_distr(mg_rho,x)
            mgo_a, mgo_dnda = shock_size_distr(mgo_rho,x)
            sio2_a, sio2_dnda = shock_size_distr(sio2_rho,x)
            fe3o4_a, fe3o4_dnda = shock_size_distr(fe3o4_rho,x)
            mgsio3_a, mgsio3_dnda = shock_size_distr(mgsio3_rho,x)
            fe2sio4_a, fe2sio4_dnda = shock_size_distr(fe2sio4_rho,x)
        ac_klambda = np.load("DATA/"+prefix+"ac_klambda.npy")
        al2o3_klambda= np.load("DATA/"+prefix+"al2o3_klambda.npy")
        asi_klambda= np.load("DATA/"+prefix+"asi_klambda.npy")
        fe_klambda= np.load("DATA/"+prefix+"fe_klambda.npy")
        fes_klambda= np.load("DATA/"+prefix+"fes_klambda.npy")
        mg2sio4_klambda= np.load("DATA/"+prefix+"mg2sio4_klambda.npy")
        mg_klambda= np.load("DATA/"+prefix+"mg_klambda.npy")
        mgo_klambda= np.load("DATA/"+prefix+"mgo_klambda.npy")
        sio2_klambda= np.load("DATA/"+prefix+"sio2_klambda.npy")
        fe3o4_klambda= np.load("DATA/"+prefix+"fe3o4_klambda.npy")
        mgsio3_klambda= np.load("DATA/"+prefix+"mgsio3_klambda.npy")
        fe2sio4_klambda= np.load("DATA/"+prefix+"fe2sio4_klambda.npy")

    if RECALC_DUSTMODELS:
        Sgeomarr = [Sgeometric(ac_a, ac_dnda),
                    Sgeometric(al2o3_a, al2o3_dnda),
                    Sgeometric(asi_a, asi_dnda),
                    Sgeometric(fe_a, fe_dnda),
                    Sgeometric(fes_a, fes_dnda),
                    Sgeometric(mg2sio4_a, mg2sio4_dnda),
                    Sgeometric(mg_a, mg_dnda),
                    Sgeometric(mgo_a, mgo_dnda),
                    Sgeometric(sio2_a, sio2_dnda),
                    Sgeometric(fe3o4_a, fe3o4_dnda),
                    Sgeometric(mgsio3_a, mgsio3_dnda),
                    Sgeometric(fe2sio4_a, fe2sio4_dnda)]
        klambdaarr = [ac_klambda,
                      al2o3_klambda,
                      asi_klambda,
                      fe_klambda,
                      fes_klambda,
                      mg2sio4_klambda,
                      mg_klambda,
                      mgo_klambda,
                      sio2_klambda,
                      fe3o4_klambda,
                      mgsio3_klambda,
                      fe2sio4_klambda]

        ### If made some AC dust from unmixed C: CD10 caveat for dust models 1 and 2
        farr = get_mass_fractions(mixed=False,depleted=True,mass=20,carbon=True)
        klambdaACUMD20 = np.zeros(ac_klambda.shape)
        SgeomACUMD20 = 0.0
        for i in xrange(len(farr)):
            klambdaACUMD20 = klambdaACUMD20 + farr[i] * klambdaarr[i]
            SgeomACUMD20 = SgeomACUMD20 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"AC-UM-D-20_klambda.npy",(klambdaACUMD20,SgeomACUMD20))
    
        farr = get_mass_fractions(mixed=False,depleted=False,mass=20,carbon=True)
        klambdaACUMND20 = np.zeros(ac_klambda.shape)
        SgeomACUMND20 = 0.0
        for i in xrange(len(farr)):
            klambdaACUMND20 = klambdaACUMND20 + farr[i] * klambdaarr[i]
            SgeomACUMND20 = SgeomACUMND20 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"AC-UM-ND-20_klambda.npy",(klambdaACUMND20,SgeomACUMND20))
        ### End AC dust

        farr = get_mass_fractions(mixed=False,depleted=True,mass=20)
        klambdaUMD20 = np.zeros(ac_klambda.shape)
        SgeomUMD20 = 0.0
        for i in xrange(len(farr)):
            klambdaUMD20 = klambdaUMD20 + farr[i] * klambdaarr[i]
            SgeomUMD20 = SgeomUMD20 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"UM-D-20_klambda.npy",(klambdaUMD20,SgeomUMD20))
    
        farr = get_mass_fractions(mixed=False,depleted=False,mass=20)
        klambdaUMND20 = np.zeros(ac_klambda.shape)
        SgeomUMND20 = 0.0
        for i in xrange(len(farr)):
            klambdaUMND20 = klambdaUMND20 + farr[i] * klambdaarr[i]
            SgeomUMND20 = SgeomUMND20 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"UM-ND-20_klambda.npy",(klambdaUMND20,SgeomUMND20))

        farr = get_mass_fractions(mixed=False,depleted=True,mass=170)
        klambdaUMD170 = np.zeros(ac_klambda.shape)
        SgeomUMD170 = 0.0
        for i in xrange(len(farr)):
            klambdaUMD170 = klambdaUMD170 + farr[i] * klambdaarr[i]
            SgeomUMD170 = SgeomUMD170 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"UM-D-170_klambda.npy",(klambdaUMD170,SgeomUMD170))

        farr = get_mass_fractions(mixed=False,depleted=False,mass=170)
        klambdaUMND170 = np.zeros(ac_klambda.shape)
        SgeomUMND170 = 0.0
        for i in xrange(len(farr)):
            klambdaUMND170 = klambdaUMND170 + farr[i] * klambdaarr[i]
            SgeomUMND170 = SgeomUMND170 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"UM-ND-170_klambda.npy",(klambdaUMND170,SgeomUMND170))
    
        farr = get_mass_fractions(mixed=True,depleted=True,mass=20)
        klambdaMD20 = np.zeros(ac_klambda.shape)
        SgeomMD20 = 0.0
        for i in xrange(len(farr)):
            klambdaMD20 = klambdaMD20 + farr[i] * klambdaarr[i]
            SgeomMD20 = SgeomMD20 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"M-D-20_klambda.npy",(klambdaMD20,SgeomMD20))
    
        farr = get_mass_fractions(mixed=True,depleted=False,mass=20)
        klambdaMND20 = np.zeros(ac_klambda.shape)
        SgeomMND20 = 0.0
        for i in xrange(len(farr)):
            klambdaMND20 = klambdaMND20 + farr[i] * klambdaarr[i]
            SgeomMND20 = SgeomMND20 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"M-ND-20_klambda.npy",(klambdaMND20,SgeomMND20))

        farr = get_mass_fractions(mixed=True,depleted=True,mass=170)
        klambdaMD170 = np.zeros(ac_klambda.shape)
        SgeomMD170 = 0.0
        for i in xrange(len(farr)):
            klambdaMD170 = klambdaMD170 + farr[i] * klambdaarr[i]
            SgeomMD170 = SgeomMD170 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"M-D-170_klambda.npy",(klambdaMD170,SgeomMD170))

        farr = get_mass_fractions(mixed=True,depleted=False,mass=170)
        klambdaMND170 = np.zeros(ac_klambda.shape)
        SgeomMND170 = 0.0
        for i in xrange(len(farr)):
            klambdaMND170 = klambdaMND170 + farr[i] * klambdaarr[i]
            SgeomMND170 = SgeomMND170 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"M-ND-170_klambda.npy",(klambdaMND170,SgeomMND170))
    
        farr = get_mass_fractions(caffau=True,mass=20,carbon=True)
        klambdaCaff20 = np.zeros(ac_klambda.shape)
        SgeomCaff20 = 0.0
        for i in xrange(len(farr)):
            klambdaCaff20 = klambdaCaff20 + farr[i] * klambdaarr[i]
            SgeomCaff20 = SgeomCaff20 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"Caff20_klambda.npy",(klambdaCaff20,SgeomCaff20))
    
        farr = get_mass_fractions(caffau=True,mass=35,carbon=True)
        klambdaCaff35 = np.zeros(ac_klambda.shape)
        SgeomCaff35 = 0.0
        for i in xrange(len(farr)):
            klambdaCaff35 = klambdaCaff35 + farr[i] * klambdaarr[i]
            SgeomCaff35 = SgeomCaff35 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"Caff35_klambda.npy",(klambdaCaff35,SgeomCaff35))
    
        farr = get_mass_fractions(caffau=True,mass=20,carbon=False)
        klambdaNoCCaff20 = np.zeros(ac_klambda.shape)
        SgeomNoCCaff20 = 0.0
        for i in xrange(len(farr)):
            klambdaNoCCaff20 = klambdaNoCCaff20 + farr[i] * klambdaarr[i]
            SgeomNoCCaff20 = SgeomNoCCaff20 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"NoC_Caff20_klambda.npy",(klambdaNoCCaff20,SgeomNoCCaff20))
    
        farr = get_mass_fractions(caffau=True,mass=35,carbon=False)
        klambdaNoCCaff35 = np.zeros(ac_klambda.shape)
        SgeomNoCCaff35 = 0.0
        for i in xrange(len(farr)):
            klambdaNoCCaff35 = klambdaNoCCaff35 + farr[i] * klambdaarr[i]
            SgeomNoCCaff35 = SgeomNoCCaff35 + farr[i] * Sgeomarr[i]
        np.save("DATA/"+prefix+"NoC_Caff35_klambda.npy",(klambdaNoCCaff35,SgeomNoCCaff35))
    
        klambdaSOIFPISN = .0641   * ac_klambda + \
                          .000582 * al2o3_klambda + \
                          .217    * mg2sio4_klambda + \
                          .705    * sio2_klambda + \
                          .0133   * fe3o4_klambda
        SgeomSOIFPISN = .0641   * Sgeometric(ac_a, ac_dnda) + \
                        .000582 * Sgeometric(al2o3_a, al2o3_dnda) + \
                        .217    * Sgeometric(mg2sio4_a, mg2sio4_dnda) + \
                        .705    * Sgeometric(sio2_a, sio2_dnda) + \
                        .0133   * Sgeometric(fe3o4_a, fe3o4_dnda)
        np.save("DATA/"+prefix+"SOIF06-PISN_klambda.npy",(klambdaSOIFPISN,SgeomSOIFPISN))

        klambdaSOIFCCSN = .0311   * ac_klambda + \
                          .000582 * al2o3_klambda + \
                          .0524   * mgsio3_klambda + \
                          .343    * mg2sio4_klambda + \
                          .000    * sio2_klambda + \
                          .0288   * fe3o4_klambda
        SgeomSOIFCCSN = .0311   * Sgeometric(ac_a, ac_dnda) + \
                        .000582 * Sgeometric(al2o3_a, al2o3_dnda) + \
                        .0524   * Sgeometric(mgsio3_a, mgsio3_dnda) + \
                        .343    * Sgeometric(mg2sio4_a, mg2sio4_dnda) + \
                        .000    * Sgeometric(sio2_a, sio2_dnda) + \
                        .0288   * Sgeometric(fe3o4_a, fe3o4_dnda)
        np.save("DATA/"+prefix+"SOIF06-CCSN_klambda.npy",(klambdaSOIFCCSN,SgeomSOIFCCSN))
    else:
        if STANDARD:
            prefix="std_"
        elif SHOCK35:
            prefix="x3.5_"
        elif SHOCK55:
            prefix="x5.5_"
        elif SHOCK75:
            prefix="x7.5_"
        klambdaACUMD20,SgeomACUMD20 = np.load("DATA/"+prefix+"AC-UM-D-20_klambda.npy")
        klambdaACUMND20,SgeomACUMND20= np.load("DATA/"+prefix+"AC-UM-ND-20_klambda.npy")
        klambdaUMD20,SgeomUMD20 = np.load("DATA/"+prefix+"UM-D-20_klambda.npy")
        klambdaUMND20,SgeomUMND20= np.load("DATA/"+prefix+"UM-ND-20_klambda.npy")
        klambdaUMD170,SgeomUMD170 = np.load("DATA/"+prefix+"UM-D-170_klambda.npy")
        klambdaUMND170,SgeomUMND170= np.load("DATA/"+prefix+"UM-ND-170_klambda.npy")
        klambdaMD20,SgeomMD20 = np.load("DATA/"+prefix+"M-D-20_klambda.npy")
        klambdaMND20,SgeomMND20= np.load("DATA/"+prefix+"M-ND-20_klambda.npy")
        klambdaMD170,SgeomMD170 = np.load("DATA/"+prefix+"M-D-170_klambda.npy")
        klambdaMND170,SgeomMND170= np.load("DATA/"+prefix+"M-ND-170_klambda.npy")
        klambdaSOIFPISN,SgeomSOIFPISN= np.load("DATA/"+prefix+"SOIF06-PISN_klambda.npy")
        klambdaSOIFCCSN,SgeomSOIFCCSN= np.load("DATA/"+prefix+"SOIF06-CCSN_klambda.npy")
        klambdaCaff20,SgeomCaff20 = np.load("DATA/"+prefix+"Caff20_klambda.npy")
        klambdaCaff35,SgeomCaff35 = np.load("DATA/"+prefix+"Caff35_klambda.npy")
        klambdaNoCCaff20,SgeomNoCCaff20 = np.load("DATA/"+prefix+"NoC_Caff20_klambda.npy")
        klambdaNoCCaff35,SgeomNoCCaff35 = np.load("DATA/"+prefix+"NoC_Caff35_klambda.npy")

    if MAKEPLOTS_KAPPALAMBDA:
        plt.subplot(7,2,1)
        plt.plot(wavelen,ac_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'amorphous C $\rho=$'+str(ac_rho))
        
        plt.subplot(7,2,2)
        plt.plot(wavelen,al2o3_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'$Al_2O_3$ $\rho=$'+str(al2o3_rho))
        
        plt.subplot(7,2,3)
        plt.plot(wavelen,asi_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'amorphous Si $\rho=$'+str(asi_rho))
        
        plt.subplot(7,2,4)
        plt.plot(wavelen,fe_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'Fe $\rho=$'+str(fe_rho))
        
        plt.subplot(7,2,5)
        plt.plot(wavelen,fes_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'FeS $\rho=$'+str(fes_rho))
        
        plt.subplot(7,2,6)
        plt.plot(wavelen,mg2sio4_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'$Mg_2SiO_4$ $\rho=$'+str(mg2sio4_rho))
        
        plt.subplot(7,2,7)
        plt.plot(wavelen,mg_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'Mg $\rho=$'+str(mg_rho))
        
        plt.subplot(7,2,8)
        plt.plot(wavelen,mgo_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'MgO $\rho=$'+str(mgo_rho))
        
        plt.subplot(7,2,9)
        plt.plot(wavelen,sio2_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'$SiO_2$ $\rho=$'+str(sio2_rho))
        
        plt.subplot(7,2,10)
        plt.plot(wavelen,fe3o4_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'$Fe_3O_4$ $\rho=$'+str(fe3o4_rho))
        
        plt.subplot(7,2,11)
        plt.plot(wavelen,mgsio3_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'$MgSiO_3$ $\rho=$'+str(mgsio3_rho))
        
        plt.subplot(7,2,12)
        plt.plot(wavelen,fe2sio4_klambda)
        plt.ylabel(r'$\kappa_{\lambda}$')
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r'$Fe_2SiO_4$ $\rho=$'+str(fe2sio4_rho))
        
        plt.subplot(7,2,13)
        plt.plot(wavelen,klambdaUMD20)
        plt.plot(wavelen,klambdaUMND20)
        plt.plot(wavelen,klambdaMD20)
        plt.plot(wavelen,klambdaMND20)
        plt.plot(wavelen,klambdaSOIFCCSN,'m:')
        plt.plot(wavelen,klambdaCaff20,'k-.')
        plt.plot(wavelen,klambdaCaff35,'k--')
        plt.plot(wavelen,klambdaACUMD20,'b:')
        plt.plot(wavelen,klambdaACUMND20,'g:')
        plt.gca().loglog()
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r"UM/M D/ND 20M$_\odot$ and S+06 CCSN and Caff20/35")

        plt.subplot(7,2,14)
        plt.plot(wavelen,klambdaUMD170)
        plt.plot(wavelen,klambdaUMND170)
        plt.plot(wavelen,klambdaMD170)
        plt.plot(wavelen,klambdaMND170)
        plt.plot(wavelen,klambdaSOIFPISN,'k:')
        plt.gca().loglog()
        plt.xlim((10**-1,10**3))
        plt.ylim((10**-3,10**6))
        plt.title(r"UM/M D/ND 170M$_\odot$ and S+06 PISN")

        plt.gcf().set_size_inches(12,32)
        plt.savefig("PLOTS/"+prefix+"klambda.pdf")
        
    if MAKEPLOTS_KAPPAPLANCK:
        plt.clf()
        Tarr = 10**np.linspace(1,4,150)
        kPlanckACUMD20_arr = np.zeros(len(Tarr))
        kPlanckACUMND20_arr= np.zeros(len(Tarr))
        kPlanckUMD20_arr = np.zeros(len(Tarr))
        kPlanckUMND20_arr= np.zeros(len(Tarr))
        kPlanckUMD170_arr = np.zeros(len(Tarr))
        kPlanckUMND170_arr= np.zeros(len(Tarr))
        kPlanckMD20_arr = np.zeros(len(Tarr))
        kPlanckMND20_arr= np.zeros(len(Tarr))
        kPlanckMD170_arr = np.zeros(len(Tarr))
        kPlanckMND170_arr= np.zeros(len(Tarr))
        kPlanckSOIFPISN_arr= np.zeros(len(Tarr))
        kPlanckSOIFCCSN_arr= np.zeros(len(Tarr))
        kPlanckCaff20_arr= np.zeros(len(Tarr))
        kPlanckCaff35_arr= np.zeros(len(Tarr))
        kPlanckNoCCaff20_arr= np.zeros(len(Tarr))
        kPlanckNoCCaff35_arr= np.zeros(len(Tarr))
        for i,T in enumerate(Tarr):
            kPlanckACUMD20_arr[i] = kappaPlanck(T, wavelen_cm, klambdaACUMD20)
            kPlanckACUMND20_arr[i]= kappaPlanck(T, wavelen_cm, klambdaACUMND20)
            kPlanckUMD20_arr[i] = kappaPlanck(T, wavelen_cm, klambdaUMD20)
            kPlanckUMND20_arr[i]= kappaPlanck(T, wavelen_cm, klambdaUMND20)
            kPlanckUMD170_arr[i] = kappaPlanck(T, wavelen_cm, klambdaUMD170)
            kPlanckUMND170_arr[i]= kappaPlanck(T, wavelen_cm, klambdaUMND170)
            kPlanckMD20_arr[i] = kappaPlanck(T, wavelen_cm, klambdaMD20)
            kPlanckMND20_arr[i]= kappaPlanck(T, wavelen_cm, klambdaMND20)
            kPlanckMD170_arr[i] = kappaPlanck(T, wavelen_cm, klambdaMD170)
            kPlanckMND170_arr[i]= kappaPlanck(T, wavelen_cm, klambdaMND170)
            kPlanckSOIFPISN_arr[i] = kappaPlanck(T, wavelen_cm, klambdaSOIFPISN)
            kPlanckSOIFCCSN_arr[i] = kappaPlanck(T, wavelen_cm, klambdaSOIFCCSN)
            kPlanckCaff20_arr[i] = kappaPlanck(T, wavelen_cm, klambdaCaff20)
            kPlanckCaff35_arr[i] = kappaPlanck(T, wavelen_cm, klambdaCaff35)
            kPlanckNoCCaff20_arr[i] = kappaPlanck(T, wavelen_cm, klambdaNoCCaff20)
            kPlanckNoCCaff35_arr[i] = kappaPlanck(T, wavelen_cm, klambdaNoCCaff35)

        plt.plot(Tarr,kPlanckUMD20_arr)
        plt.plot(Tarr,kPlanckUMND20_arr)
        plt.plot(Tarr,kPlanckMD20_arr)
        plt.plot(Tarr,kPlanckMND20_arr)
        plt.plot(Tarr,kPlanckSOIFCCSN_arr,'m:')
        plt.plot(Tarr,kPlanckCaff20_arr,'k--')
        plt.plot(Tarr,kPlanckCaff35_arr,'k--')
        plt.plot(Tarr,kPlanckNoCCaff20_arr,'k:')
        plt.plot(Tarr,kPlanckNoCCaff35_arr,'k:')
        plt.plot([1500,1500],[.1,10**4],'k--')
        plt.plot([2000,2000],[.1,10**4],'k--')
        plt.plot(Tarr,kPlanckACUMD20_arr,'b:')
        plt.plot(Tarr,kPlanckACUMND20_arr,'g:')
        plt.ylim((.1,10**4))
        plt.xlabel(r'$T_{dust}$')
        plt.ylabel(r'$\kappa_P$')
        plt.title(r'20M$_\odot$, CCSN')
        plt.legend([r'UM D 20M$_\odot$',r'UM ND 20M$_\odot$',r'M D 20M$_\odot$',r'M ND 20M$_\odot$',r'S+06 CCSN',r'Caff20',r'Caff35'],loc="lower right")
        plt.gca().loglog()
        plt.gcf().set_size_inches(6,6)
        plt.savefig("PLOTS/"+prefix+"kPlanck_20.pdf")

        plt.clf()
        plt.plot(Tarr,kPlanckUMD170_arr)
        plt.plot(Tarr,kPlanckUMND170_arr)
        plt.plot(Tarr,kPlanckMD170_arr)
        plt.plot(Tarr,kPlanckMND170_arr)
        plt.plot(Tarr,kPlanckSOIFPISN_arr,'k:')
        plt.plot([1500,1500],[.1,10**4],'k--')
        plt.plot([2000,2000],[.1,10**4],'k--')
        plt.ylim((.1,10**4))
        plt.xlabel(r'$T_{dust}$')
        plt.ylabel(r'$\kappa_P$')
        plt.title(r'170M$_\odot$, PISN')
        plt.legend([r'UM D 170M$_\odot$',r'UM ND 170M$_\odot$',r'M D 170M$_\odot$',r'M ND 170M$_\odot$',r'S+06 PISN'],loc="lower right")
        plt.gca().loglog()
        plt.gcf().set_size_inches(6,6)
        plt.savefig("PLOTS/"+prefix+"kPlanck_170.pdf")
