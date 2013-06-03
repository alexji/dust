import numpy as np
import pylab as plt
from scipy import interpolate

if __name__ == "__main__":
    ac      = np.loadtxt(open("opticalconsts/ac.csv",'r'),       delimiter=',',skiprows=1)
    al2o3   = np.loadtxt(open("opticalconsts/al2o3.csv",'r'),    delimiter=',',skiprows=1)
    asi     = np.loadtxt(open("opticalconsts/asi.csv",'r'),      delimiter=',',skiprows=1)
    fe      = np.loadtxt(open("opticalconsts/fe.csv",'r'),       delimiter=',',skiprows=1)
    fes     = np.loadtxt(open("opticalconsts/fes.csv",'r'),      delimiter=',',skiprows=1)
    mg2sio4 = np.loadtxt(open("opticalconsts/mg2sio4.csv",'r'),  delimiter=',',skiprows=1)
    mg      = np.loadtxt(open("opticalconsts/mg.csv",'r'),       delimiter=',',skiprows=1)
    mgo     = np.loadtxt(open("opticalconsts/mgo.csv",'r'),      delimiter=',',skiprows=1)
    sio2    = np.loadtxt(open("opticalconsts/sio2.csv",'r'),     delimiter=',',skiprows=1)
    fe3o4   = np.loadtxt(open("opticalconsts/fe3o4.csv",'r'),    delimiter=',',skiprows=1)
    mgsio3  = np.loadtxt(open("opticalconsts/mgsio3.csv",'r'),   delimiter=',',skiprows=1)
    fe2sio4 = np.loadtxt(open("opticalconsts/fe2sio4.csv",'r'),  delimiter=',',skiprows=1)

    ## Let this be the standardized wavelength array
    wavelen = 10**np.linspace(-2,5,1400)
    logwavelen = np.log10(wavelen)

    ####### A-Carbon
    plt.subplot(6,2,1)
    plt.plot(ac[:,0],ac[:,1],'bo-',markersize=2)
    plt.plot(ac[:,0],ac[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.ylabel('n, k')
    plt.title('amorphous C')
    
    int_ac_n_x     = ac[ac[:,1] > 0,0]
    int_ac_n_y     = ac[ac[:,1] > 0,1]
    int_ac_n       = interpolate.interp1d(np.log10(int_ac_n_x),np.log10(int_ac_n_y))
    spl_ac_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_ac_n_x[0:5]),np.log10(int_ac_n_y[0:5]),k=1)
    spl_ac_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_ac_n_x[-5:]),np.log10(int_ac_n_y[-5:]),k=1)
    def ac_opt_n(logwavelen):
        def ac_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_ac_n_x[0]):
                return 10**spl_ac_n_small(logwavelen)
            elif logwavelen > np.log10(int_ac_n_x[-1]):
                return 10**spl_ac_n_large(logwavelen)
            else:
                return 10**int_ac_n(logwavelen)
        return [ac_opt_n_helper(x) for x in logwavelen]
    ac_n = ac_opt_n(logwavelen)
    plt.plot(wavelen,ac_n,'r--')

    int_ac_k_x     = ac[ac[:,2] > 0,0]
    int_ac_k_y     = ac[ac[:,2] > 0,2]
    int_ac_k       = interpolate.interp1d(np.log10(int_ac_k_x),np.log10(int_ac_k_y))
    spl_ac_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_ac_k_x[0:5]),np.log10(int_ac_k_y[0:5]),k=1)
    spl_ac_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_ac_k_x[-5:]),np.log10(int_ac_k_y[-5:]),k=1)
    def ac_opt_k(logwavelen):
        def ac_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_ac_k_x[0]):
                return 10**spl_ac_k_small(logwavelen)
            elif logwavelen > np.log10(int_ac_k_x[-1]):
                return 10**spl_ac_k_large(logwavelen)
            else:
                return 10**int_ac_k(logwavelen)
        return [ac_opt_k_helper(x) for x in logwavelen]
    ac_k = ac_opt_k(logwavelen)
    plt.plot(wavelen,ac_k,'y--')
    plt.ylim(10**-4,10**2)
    np.save("opticalconsts/ac.npy",(ac_n,ac_k))

    ####### Al2O3
    plt.subplot(6,2,2)
    plt.plot(al2o3[:,0],al2o3[:,1],'bo-',markersize=2)
    plt.plot(al2o3[:,0],al2o3[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.ylabel('n, k')
    plt.title(r'$Al_2O_3$')
    
    int_al2o3_n_x     = al2o3[al2o3[:,1] > 0,0]
    int_al2o3_n_y     = al2o3[al2o3[:,1] > 0,1]
    int_al2o3_n       = interpolate.interp1d(np.log10(int_al2o3_n_x),np.log10(int_al2o3_n_y))
    spl_al2o3_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_al2o3_n_x[0:5]),np.log10(int_al2o3_n_y[0:5]),k=1)
    spl_al2o3_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_al2o3_n_x[-5:]),np.log10(int_al2o3_n_y[-5:]),k=1)
    def al2o3_opt_n(logwavelen):
        def al2o3_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_al2o3_n_x[0]):
                return 10**spl_al2o3_n_small(logwavelen)
            elif logwavelen > np.log10(int_al2o3_n_x[-1]):
                return 10**spl_al2o3_n_large(logwavelen)
            else:
                return 10**int_al2o3_n(logwavelen)
        return [al2o3_opt_n_helper(x) for x in logwavelen]
    al2o3_n = al2o3_opt_n(logwavelen)
    plt.plot(wavelen,al2o3_n,'r--')

    int_al2o3_k_x     = al2o3[al2o3[:,2] > 0,0]
    int_al2o3_k_y     = al2o3[al2o3[:,2] > 0,2]
    int_al2o3_k       = interpolate.interp1d(np.log10(int_al2o3_k_x),np.log10(int_al2o3_k_y))
    spl_al2o3_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_al2o3_k_x[0:5]),np.log10(int_al2o3_k_y[0:5]),k=1)
    spl_al2o3_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_al2o3_k_x[-5:]),np.log10(int_al2o3_k_y[-5:]),k=1)
    def al2o3_opt_k(logwavelen):
        def al2o3_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_al2o3_k_x[0]):
                return 10**spl_al2o3_k_small(logwavelen)
            elif logwavelen > np.log10(int_al2o3_k_x[-1]):
                return 10**spl_al2o3_k_large(logwavelen)
            else:
                return 10**int_al2o3_k(logwavelen)
        return [al2o3_opt_k_helper(x) for x in logwavelen]
    al2o3_k = al2o3_opt_k(logwavelen)
    plt.plot(wavelen,al2o3_k,'y--')
    np.save("opticalconsts/al2o3.npy",(al2o3_n,al2o3_k))

    ####### A-Silicon
    plt.subplot(6,2,3)
    plt.plot(asi[:,0],asi[:,1],'bo-',markersize=2)
    plt.plot(asi[:,0],asi[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.ylabel('n, k')
    plt.title('amorphous Si')
    
    int_asi_n_x     = asi[asi[:,1] > 0,0]
    int_asi_n_y     = asi[asi[:,1] > 0,1]
    int_asi_n       = interpolate.interp1d(np.log10(int_asi_n_x),np.log10(int_asi_n_y))
    spl_asi_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_asi_n_x[0:5]),np.log10(int_asi_n_y[0:5]),k=1)
    spl_asi_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_asi_n_x[-5:]),np.log10(int_asi_n_y[-5:]),k=1)
    def asi_opt_n(logwavelen):
        def asi_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_asi_n_x[0]):
                return 10**spl_asi_n_small(logwavelen)
            elif logwavelen > np.log10(int_asi_n_x[-1]):
                return 10**spl_asi_n_large(logwavelen)
            else:
                return 10**int_asi_n(logwavelen)
        return [asi_opt_n_helper(x) for x in logwavelen]
    asi_n = asi_opt_n(logwavelen)
    plt.plot(wavelen,asi_n,'r--')

    int_asi_k_x     = asi[asi[:,2] > 0,0]
    int_asi_k_y     = asi[asi[:,2] > 0,2]
    int_asi_k       = interpolate.interp1d(np.log10(int_asi_k_x),np.log10(int_asi_k_y))
    spl_asi_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_asi_k_x[0:5]),np.log10(int_asi_k_y[0:5]),k=1)
    spl_asi_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_asi_k_x[-10:-6]),np.log10(int_asi_k_y[-10:-6]),k=1)
    def asi_opt_k(logwavelen):
        def asi_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_asi_k_x[0]):
                return 10**spl_asi_k_small(logwavelen)
            elif logwavelen > np.log10(int_asi_k_x[-6]):
                return 10**spl_asi_k_large(logwavelen)
            else:
                return 10**int_asi_k(logwavelen)
        return [asi_opt_k_helper(x) for x in logwavelen]
    asi_k = asi_opt_k(logwavelen)
    plt.plot(wavelen,asi_k,'y--')
    plt.ylim(10**-4,10**1)
    np.save("opticalconsts/asi.npy",(asi_n,asi_k))

    ####### Fe
    plt.subplot(6,2,4)
    plt.plot(fe[:,0],fe[:,1],'bo-',markersize=2)
    plt.plot(fe[:,0],fe[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.ylabel('n, k')
    plt.title('Fe')
    
    int_fe_n_x     = fe[fe[:,1] > 0,0]
    int_fe_n_y     = fe[fe[:,1] > 0,1]
    int_fe_n       = interpolate.interp1d(np.log10(int_fe_n_x),np.log10(int_fe_n_y))
    spl_fe_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe_n_x[0:5]),np.log10(int_fe_n_y[0:5]),k=1)
    spl_fe_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe_n_x[-5:]),np.log10(int_fe_n_y[-5:]),k=1)
    def fe_opt_n(logwavelen):
        def fe_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_fe_n_x[0]):
                return 10**spl_fe_n_small(logwavelen)
            elif logwavelen > np.log10(int_fe_n_x[-1]):
                return 10**spl_fe_n_large(logwavelen)
            else:
                return 10**int_fe_n(logwavelen)
        return [fe_opt_n_helper(x) for x in logwavelen]
    fe_n = fe_opt_n(logwavelen)
    plt.plot(wavelen,fe_n,'r--')

    int_fe_k_x     = fe[fe[:,2] > 0,0]
    int_fe_k_y     = fe[fe[:,2] > 0,2]
    int_fe_k       = interpolate.interp1d(np.log10(int_fe_k_x),np.log10(int_fe_k_y))
    spl_fe_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe_k_x[0:5]),np.log10(int_fe_k_y[0:5]),k=1)
    spl_fe_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe_k_x[-5:]),np.log10(int_fe_k_y[-5:]),k=1)
    def fe_opt_k(logwavelen):
        def fe_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_fe_k_x[0]):
                return 10**spl_fe_k_small(logwavelen)
            elif logwavelen > np.log10(int_fe_k_x[-1]):
                return 10**spl_fe_k_large(logwavelen)
            else:
                return 10**int_fe_k(logwavelen)
        return [fe_opt_k_helper(x) for x in logwavelen]
    fe_k = fe_opt_k(logwavelen)
    plt.plot(wavelen,fe_k,'y--')
    np.save("opticalconsts/fe.npy",(fe_n,fe_k))

    ####### FeS
    plt.subplot(6,2,5)
    plt.plot(fes[:,0],fes[:,1],'bo-',markersize=2)
    plt.plot(fes[:,0],fes[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.ylabel('n, k')
    plt.title('FeS')
    
    int_fes_n_x     = fes[fes[:,1] > 0,0]
    int_fes_n_y     = fes[fes[:,1] > 0,1]
    int_fes_n       = interpolate.interp1d(np.log10(int_fes_n_x),np.log10(int_fes_n_y))
    spl_fes_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_fes_n_x[0:5]),np.log10(int_fes_n_y[0:5]),k=1)
    spl_fes_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_fes_n_x[-5:]),np.log10(int_fes_n_y[-5:]),k=1)
    def fes_opt_n(logwavelen):
        def fes_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_fes_n_x[0]):
                return 10**spl_fes_n_small(logwavelen)
            elif logwavelen > np.log10(int_fes_n_x[-1]):
                return 10**spl_fes_n_large(logwavelen)
            else:
                return 10**int_fes_n(logwavelen)
        return [fes_opt_n_helper(x) for x in logwavelen]
    fes_n = fes_opt_n(logwavelen)
    plt.plot(wavelen,fes_n,'r--')

    int_fes_k_x     = fes[fes[:,2] > 0,0]
    int_fes_k_y     = fes[fes[:,2] > 0,2]
    int_fes_k       = interpolate.interp1d(np.log10(int_fes_k_x),np.log10(int_fes_k_y))
    spl_fes_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_fes_k_x[0:5]),np.log10(int_fes_k_y[0:5]),k=1)
    spl_fes_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_fes_k_x[-5:]),np.log10(int_fes_k_y[-5:]),k=1)
    def fes_opt_k(logwavelen):
        def fes_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_fes_k_x[0]):
                return 10**spl_fes_k_small(logwavelen)
            elif logwavelen > np.log10(int_fes_k_x[-1]):
                return 10**spl_fes_k_large(logwavelen)
            else:
                return 10**int_fes_k(logwavelen)
        return [fes_opt_k_helper(x) for x in logwavelen]
    fes_k = fes_opt_k(logwavelen)
    plt.plot(wavelen,fes_k,'y--')
    np.save("opticalconsts/fes.npy",(fes_n,fes_k))

    ####### Mg2SiO4
    plt.subplot(6,2,6)
    plt.plot(mg2sio4[:,0],mg2sio4[:,1],'bo-',markersize=2)
    plt.plot(mg2sio4[:,0],mg2sio4[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.ylabel('n, k')
    plt.title(r'$Mg_2SiO_4$')
    
    int_mg2sio4_n_x     = mg2sio4[mg2sio4[:,1] > 0,0]
    int_mg2sio4_n_y     = mg2sio4[mg2sio4[:,1] > 0,1]
    int_mg2sio4_n       = interpolate.interp1d(np.log10(int_mg2sio4_n_x),np.log10(int_mg2sio4_n_y))
    spl_mg2sio4_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_mg2sio4_n_x[0:5]),np.log10(int_mg2sio4_n_y[0:5]),k=1)
    spl_mg2sio4_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_mg2sio4_n_x[-5:]),np.log10(int_mg2sio4_n_y[-5:]),k=1)
    def mg2sio4_opt_n(logwavelen):
        def mg2sio4_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_mg2sio4_n_x[0]):
                return 10**spl_mg2sio4_n_small(logwavelen)
            elif logwavelen > np.log10(int_mg2sio4_n_x[-1]):
                return 10**spl_mg2sio4_n_large(logwavelen)
            else:
                return 10**int_mg2sio4_n(logwavelen)
        return [mg2sio4_opt_n_helper(x) for x in logwavelen]
    mg2sio4_n = mg2sio4_opt_n(logwavelen)
    plt.plot(wavelen,mg2sio4_n,'r--')

    int_mg2sio4_k_x     = mg2sio4[mg2sio4[:,2] > 0,0]
    int_mg2sio4_k_y     = mg2sio4[mg2sio4[:,2] > 0,2]
    int_mg2sio4_k       = interpolate.interp1d(np.log10(int_mg2sio4_k_x),np.log10(int_mg2sio4_k_y))
    spl_mg2sio4_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_mg2sio4_k_x[0:5]),np.log10(int_mg2sio4_k_y[0:5]),k=1)
    spl_mg2sio4_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_mg2sio4_k_x[-5:]),np.log10(int_mg2sio4_k_y[-5:]),k=1)
    def mg2sio4_opt_k(logwavelen):
        def mg2sio4_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_mg2sio4_k_x[0]):
                return 10**spl_mg2sio4_k_small(logwavelen)
            elif logwavelen > np.log10(int_mg2sio4_k_x[-1]):
                return 10**spl_mg2sio4_k_large(logwavelen)
            else:
                return 10**int_mg2sio4_k(logwavelen)
        return [mg2sio4_opt_k_helper(x) for x in logwavelen]
    mg2sio4_k = mg2sio4_opt_k(logwavelen)
    plt.plot(wavelen,mg2sio4_k,'y--')
    np.save("opticalconsts/mg2sio4.npy",(mg2sio4_n,mg2sio4_k))

    ####### Mg
    plt.subplot(6,2,7)
    plt.plot(mg[:,0],mg[:,1],'bo-',markersize=2)
    plt.plot(mg[:,0],mg[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.ylabel('n, k')
    plt.title('Mg')
    
    int_mg_n_x     = mg[mg[:,1] > 0,0]
    int_mg_n_y     = mg[mg[:,1] > 0,1]
    int_mg_n       = interpolate.interp1d(np.log10(int_mg_n_x),np.log10(int_mg_n_y))
    spl_mg_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_mg_n_x[0:5]),np.log10(int_mg_n_y[0:5]),k=1)
    spl_mg_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_mg_n_x[-5:]),np.log10(int_mg_n_y[-5:]),k=1)
    def mg_opt_n(logwavelen):
        def mg_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_mg_n_x[0]):
                return 10**spl_mg_n_small(logwavelen)
            elif logwavelen > np.log10(int_mg_n_x[-1]):
                return 10**spl_mg_n_large(logwavelen)
            else:
                return 10**int_mg_n(logwavelen)
        return [mg_opt_n_helper(x) for x in logwavelen]
    mg_n = mg_opt_n(logwavelen)
    plt.plot(wavelen,mg_n,'r--')

    int_mg_k_x     = mg[mg[:,2] > 0,0]
    int_mg_k_y     = mg[mg[:,2] > 0,2]
    int_mg_k       = interpolate.interp1d(np.log10(int_mg_k_x),np.log10(int_mg_k_y))
    spl_mg_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_mg_k_x[0:5]),np.log10(int_mg_k_y[0:5]),k=1)
    spl_mg_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_mg_k_x[-5:]),np.log10(int_mg_k_y[-5:]),k=1)
    def mg_opt_k(logwavelen):
        def mg_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_mg_k_x[0]):
                return 10**spl_mg_k_small(logwavelen)
            elif logwavelen > np.log10(int_mg_k_x[-1]):
                return 10**spl_mg_k_large(logwavelen)
            else:
                return 10**int_mg_k(logwavelen)
        return [mg_opt_k_helper(x) for x in logwavelen]
    mg_k = mg_opt_k(logwavelen)
    plt.plot(wavelen,mg_k,'y--')
    np.save("opticalconsts/mg.npy",(mg_n,mg_k))

    ####### MgO
    plt.subplot(6,2,8)
    plt.plot(mgo[:,0],mgo[:,1],'bo-',markersize=2)
    plt.plot(mgo[:,0],mgo[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.ylabel('n, k')
    plt.title('MgO')
    
    int_mgo_n_x     = mgo[mgo[:,1] > 0,0]
    int_mgo_n_y     = mgo[mgo[:,1] > 0,1]
    int_mgo_n       = interpolate.interp1d(np.log10(int_mgo_n_x),np.log10(int_mgo_n_y))
    spl_mgo_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_mgo_n_x[0:5]),np.log10(int_mgo_n_y[0:5]),k=1)
    spl_mgo_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_mgo_n_x[-5:]),np.log10(int_mgo_n_y[-5:]),k=1)
    def mgo_opt_n(logwavelen):
        def mgo_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_mgo_n_x[0]):
                return 10**spl_mgo_n_small(logwavelen)
            elif logwavelen > np.log10(int_mgo_n_x[-1]):
                return 10**spl_mgo_n_large(logwavelen)
            else:
                return 10**int_mgo_n(logwavelen)
        return [mgo_opt_n_helper(x) for x in logwavelen]
    mgo_n = mgo_opt_n(logwavelen)
    plt.plot(wavelen,mgo_n,'r--')

    int_mgo_k_x     = mgo[mgo[:,2] > 0,0]
    int_mgo_k_y     = mgo[mgo[:,2] > 0,2]
    int_mgo_k       = interpolate.interp1d(np.log10(int_mgo_k_x),np.log10(int_mgo_k_y))
    spl_mgo_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_mgo_k_x[0:5]),np.log10(int_mgo_k_y[0:5]),k=1)
    spl_mgo_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_mgo_k_x[-5:]),np.log10(int_mgo_k_y[-5:]),k=1)
    def mgo_opt_k(logwavelen):
        def mgo_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_mgo_k_x[0]):
                return 10**spl_mgo_k_small(logwavelen)
            elif logwavelen > np.log10(int_mgo_k_x[-1]):
                return 10**spl_mgo_k_large(logwavelen)
            else:
                return 10**int_mgo_k(logwavelen)
        return [mgo_opt_k_helper(x) for x in logwavelen]
    mgo_k = mgo_opt_k(logwavelen)
    plt.plot(wavelen,mgo_k,'y--')
    np.save("opticalconsts/mgo.npy",(mgo_n,mgo_k))

    ####### SiO2
    plt.subplot(6,2,9)
    plt.plot(sio2[:,0],sio2[:,1],'bo-',markersize=2)
    plt.plot(sio2[:,0],sio2[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.ylabel('n, k')
    plt.title(r'$SiO_2$')

    int_sio2_n_x     = sio2[sio2[:,1] > 0,0]
    int_sio2_n_y     = sio2[sio2[:,1] > 0,1]
    int_sio2_n       = interpolate.interp1d(np.log10(int_sio2_n_x),np.log10(int_sio2_n_y))
    spl_sio2_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_sio2_n_x[0:5]),np.log10(int_sio2_n_y[0:5]),k=1)
    spl_sio2_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_sio2_n_x[-5:]),np.log10(int_sio2_n_y[-5:]),k=1)
    def sio2_opt_n(logwavelen):
        def sio2_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_sio2_n_x[0]):
                return 10**spl_sio2_n_small(logwavelen)
            elif logwavelen > np.log10(int_sio2_n_x[-1]):
                return 10**spl_sio2_n_large(logwavelen)
            else:
                return 10**int_sio2_n(logwavelen)
        return [sio2_opt_n_helper(x) for x in logwavelen]
    sio2_n = sio2_opt_n(logwavelen)
    plt.plot(wavelen,sio2_n,'r--')

    int_sio2_k_x     = sio2[sio2[:,2] > 0,0]
    int_sio2_k_y     = sio2[sio2[:,2] > 0,2]
    int_sio2_k       = interpolate.interp1d(np.log10(int_sio2_k_x),np.log10(int_sio2_k_y))
    spl_sio2_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_sio2_k_x[0:5]),np.log10(int_sio2_k_y[0:5]),k=1)
    spl_sio2_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_sio2_k_x[-5:]),np.log10(int_sio2_k_y[-5:]),k=1)
    def sio2_opt_k(logwavelen):
        def sio2_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_sio2_k_x[0]):
                return 10**spl_sio2_k_small(logwavelen)
            elif logwavelen > np.log10(int_sio2_k_x[-1]):
                return 10**spl_sio2_k_large(logwavelen)
            else:
                return 10**int_sio2_k(logwavelen)
        return [sio2_opt_k_helper(x) for x in logwavelen]
    sio2_k = sio2_opt_k(logwavelen)
    plt.plot(wavelen,sio2_k,'y--')
    np.save("opticalconsts/sio2.npy",(sio2_n,sio2_k))

    ####### Fe3O4
    plt.subplot(6,2,10)
    plt.plot(fe3o4[:,0],fe3o4[:,1],'bo-',markersize=2)
    plt.plot(fe3o4[:,0],fe3o4[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.ylabel('n, k')
    plt.title(r'$Fe_3O_4$')

    int_fe3o4_n_x     = fe3o4[fe3o4[:,1] > 0,0]
    int_fe3o4_n_y     = fe3o4[fe3o4[:,1] > 0,1]
    int_fe3o4_n       = interpolate.interp1d(np.log10(int_fe3o4_n_x),np.log10(int_fe3o4_n_y))
    spl_fe3o4_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe3o4_n_x[0:5]),np.log10(int_fe3o4_n_y[0:5]),k=1)
    spl_fe3o4_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe3o4_n_x[-5:]),np.log10(int_fe3o4_n_y[-5:]),k=1)
    def fe3o4_opt_n(logwavelen):
        def fe3o4_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_fe3o4_n_x[0]):
                return 10**spl_fe3o4_n_small(logwavelen)
            elif logwavelen > np.log10(int_fe3o4_n_x[-1]):
                return 10**spl_fe3o4_n_large(logwavelen)
            else:
                return 10**int_fe3o4_n(logwavelen)
        return [fe3o4_opt_n_helper(x) for x in logwavelen]
    fe3o4_n = fe3o4_opt_n(logwavelen)
    plt.plot(wavelen,fe3o4_n,'r--')

    int_fe3o4_k_x     = fe3o4[fe3o4[:,2] > 0,0]
    int_fe3o4_k_y     = fe3o4[fe3o4[:,2] > 0,2]
    int_fe3o4_k       = interpolate.interp1d(np.log10(int_fe3o4_k_x),np.log10(int_fe3o4_k_y))
    spl_fe3o4_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe3o4_k_x[0:5]),np.log10(int_fe3o4_k_y[0:5]),k=1)
    spl_fe3o4_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe3o4_k_x[-5:]),np.log10(int_fe3o4_k_y[-5:]),k=1)
    def fe3o4_opt_k(logwavelen):
        def fe3o4_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_fe3o4_k_x[0]):
                return 10**spl_fe3o4_k_small(logwavelen)
            elif logwavelen > np.log10(int_fe3o4_k_x[-1]):
                return 10**spl_fe3o4_k_large(logwavelen)
            else:
                return 10**int_fe3o4_k(logwavelen)
        return [fe3o4_opt_k_helper(x) for x in logwavelen]
    fe3o4_k = fe3o4_opt_k(logwavelen)
    plt.plot(wavelen,fe3o4_k,'y--')
    np.save("opticalconsts/fe3o4.npy",(fe3o4_n,fe3o4_k))

    ####### MgSiO3
    plt.subplot(6,2,11)
    plt.plot(mgsio3[:,0],mgsio3[:,1],'bo-',markersize=2)
    plt.plot(mgsio3[:,0],mgsio3[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.xlabel(r'Wavelength ($\mu$m)')
    plt.ylabel('n, k')
    plt.title(r'$MgSiO_3$')
    
    int_mgsio3_n_x     = mgsio3[mgsio3[:,1] > 0,0]
    int_mgsio3_n_y     = mgsio3[mgsio3[:,1] > 0,1]
    int_mgsio3_n       = interpolate.interp1d(np.log10(int_mgsio3_n_x),np.log10(int_mgsio3_n_y))
    spl_mgsio3_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_mgsio3_n_x[0:5]),np.log10(int_mgsio3_n_y[0:5]),k=1)
    spl_mgsio3_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_mgsio3_n_x[-5:]),np.log10(int_mgsio3_n_y[-5:]),k=1)
    def mgsio3_opt_n(logwavelen):
        def mgsio3_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_mgsio3_n_x[0]):
                return 10**spl_mgsio3_n_small(logwavelen)
            elif logwavelen > np.log10(int_mgsio3_n_x[-1]):
                return 10**spl_mgsio3_n_large(logwavelen)
            else:
                return 10**int_mgsio3_n(logwavelen)
        return [mgsio3_opt_n_helper(x) for x in logwavelen]
    mgsio3_n = mgsio3_opt_n(logwavelen)
    plt.plot(wavelen,mgsio3_n,'r--')

    int_mgsio3_k_x     = mgsio3[mgsio3[:,2] > 0,0]
    int_mgsio3_k_y     = mgsio3[mgsio3[:,2] > 0,2]
    int_mgsio3_k       = interpolate.interp1d(np.log10(int_mgsio3_k_x),np.log10(int_mgsio3_k_y))
    spl_mgsio3_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_mgsio3_k_x[0:5]),np.log10(int_mgsio3_k_y[0:5]),k=1)
    spl_mgsio3_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_mgsio3_k_x[-5:]),np.log10(int_mgsio3_k_y[-5:]),k=1)
    def mgsio3_opt_k(logwavelen):
        def mgsio3_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_mgsio3_k_x[0]):
                return 10**spl_mgsio3_k_small(logwavelen)
            elif logwavelen > np.log10(int_mgsio3_k_x[-1]):
                return 10**spl_mgsio3_k_large(logwavelen)
            else:
                return 10**int_mgsio3_k(logwavelen)
        return [mgsio3_opt_k_helper(x) for x in logwavelen]
    mgsio3_k = mgsio3_opt_k(logwavelen)
    plt.plot(wavelen,mgsio3_k,'y--')
    np.save("opticalconsts/mgsio3.npy",(mgsio3_n,mgsio3_k))

    ####### Fe2SiO4
    plt.subplot(6,2,12)
    plt.plot(fe2sio4[:,0],fe2sio4[:,1],'bo-',markersize=2)
    plt.plot(fe2sio4[:,0],fe2sio4[:,2],'go-',markersize=2)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.xlabel(r'Wavelength ($\mu$m)')
    plt.ylabel('n, k')
    plt.title(r'$Fe_2SiO_4$')
    
    int_fe2sio4_n_x     = fe2sio4[fe2sio4[:,1] > 0,0]
    int_fe2sio4_n_y     = fe2sio4[fe2sio4[:,1] > 0,1]
    int_fe2sio4_n       = interpolate.interp1d(np.log10(int_fe2sio4_n_x),np.log10(int_fe2sio4_n_y))
    spl_fe2sio4_n_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe2sio4_n_x[0:5]),np.log10(int_fe2sio4_n_y[0:5]),k=1)
    spl_fe2sio4_n_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe2sio4_n_x[-5:]),np.log10(int_fe2sio4_n_y[-5:]),k=1)
    def fe2sio4_opt_n(logwavelen):
        def fe2sio4_opt_n_helper(logwavelen):
            if logwavelen < np.log10(int_fe2sio4_n_x[0]):
                return 10**spl_fe2sio4_n_small(logwavelen)
            elif logwavelen > np.log10(int_fe2sio4_n_x[-1]):
                return 10**spl_fe2sio4_n_large(logwavelen)
            else:
                return 10**int_fe2sio4_n(logwavelen)
        return [fe2sio4_opt_n_helper(x) for x in logwavelen]
    fe2sio4_n = fe2sio4_opt_n(logwavelen)
    plt.plot(wavelen,fe2sio4_n,'r--')

    int_fe2sio4_k_x     = fe2sio4[fe2sio4[:,2] > 0,0]
    int_fe2sio4_k_y     = fe2sio4[fe2sio4[:,2] > 0,2]
    int_fe2sio4_k       = interpolate.interp1d(np.log10(int_fe2sio4_k_x),np.log10(int_fe2sio4_k_y))
    spl_fe2sio4_k_small = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe2sio4_k_x[48:50]),np.log10(int_fe2sio4_k_y[48:50]),k=1)
    spl_fe2sio4_k_large = interpolate.InterpolatedUnivariateSpline(np.log10(int_fe2sio4_k_x[-5:]),np.log10(int_fe2sio4_k_y[-5:]),k=1)
    def fe2sio4_opt_k(logwavelen):
        def fe2sio4_opt_k_helper(logwavelen):
            if logwavelen < np.log10(int_fe2sio4_k_x[48]):
                return 10**spl_fe2sio4_k_small(logwavelen)
            elif logwavelen > np.log10(int_fe2sio4_k_x[-1]):
                return 10**spl_fe2sio4_k_large(logwavelen)
            else:
                return 10**int_fe2sio4_k(logwavelen)
        return [fe2sio4_opt_k_helper(x) for x in logwavelen]
    fe2sio4_k = fe2sio4_opt_k(logwavelen)
    plt.plot(wavelen,fe2sio4_k,'y--')
    np.save("opticalconsts/fe2sio4.npy",(fe2sio4_n,fe2sio4_k))

    plt.gcf().set_size_inches(12,30)
    plt.savefig("PLOTS/optical_consts.pdf")
