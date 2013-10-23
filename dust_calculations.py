import numpy as np
import pylab as plt
from scipy.integrate import simps, trapz
from calcq import calcq

def get_wavelen():
    ## wavelen in um
    wavelen    = 10**np.linspace(-2,5,1400)
    logwavelen = np.log10(wavelen)
    return wavelen, logwavelen

def get_densities():
    ## densities in g/cm**3
    ac_rho     = 2.0 ##Estimated from Wikipedia #HOIC 2.25 graphite
    al2o3_rho  = 4.0 ##Estimated from Wikipedia #HOIC 3.965
    asi_rho    = 2.33##http://www.mit.edu/~6.777/matprops/asi.htm #HOIC 2.33-2.53
    fe_rho     = 7.87##Semenov et al
    fes_rho    = 4.83##Semenov et al
    mg2sio4_rho= 3.59##Semenov et al
    mg_rho     = 1.6 ##Estimated from Wikipedia #HOIC 1.57-1.74
    mgo_rho    = 3.58##Estimated from Wikipedia #HOIC 3.58
    sio2_rho   = 2.6 ##Estimated from Wikipedia #HOIC 2.53-2.65 quartz (amorphous slightly lower, but impure)
    fe3o4_rho  = 5.17##Estimated from Wikipedia #HOIC 5.17
    mgsio3_rho = 3.20##Semenov et al
    fe2sio4_rho= 4.39##http://webmineral.com/data/Fayalite.shtml and wikipedia
    return ac_rho, al2o3_rho, asi_rho, fe_rho, fes_rho, mg2sio4_rho, mg_rho, mgo_rho, sio2_rho, fe3o4_rho, mgsio3_rho, fe2sio4_rho

def get_mass_fractions(mixed=False,depleted=False,mass=20,caffau=False,carbon=False,carbfrac=None):
    if (caffau and mass==35):
        ## Caff35
        ## SDSS J102915 Schneider+12
        if carbon:
            ac_M     = 0.022
        else:
            ac_M     = 0
        al2o3_M  = 0.007
        asi_M    = 0.0
        fe_M     = 0.0
        fes_M    = 0.0
        mg2sio4_M= 0.12
        mg_M     = 0.0
        mgo_M    = 0.0
        sio2_M   = 0.0
        fe3o4_M  = 0.22
        mgsio3_M = 0.042
        fe2sio4_M= 0.0
    elif (caffau and mass==20):
        ## Caff20
        ## SDSS J102915 Schneider+12
        if carbon:
            ac_M     = 0.003
        else:
            ac_M     = 0
        al2o3_M  = 0.006
        asi_M    = 0.0
        fe_M     = 0.0
        fes_M    = 0.0
        mg2sio4_M= 0.07
        mg_M     = 0.0
        mgo_M    = 0.0
        sio2_M   = 0.0
        fe3o4_M  = 0.09
        mgsio3_M = 0.21
        fe2sio4_M= 0.0
    elif (not mixed and not depleted and mass==20):
        ## UM ND M20
        if carbon:
            ac_M     = .0145
        else:
            ac_M     = 1.3e-13
        al2o3_M  = 0.000086
        asi_M    = 0.030
        fe_M     = 0.000046
        fes_M    = 0.033
        mg2sio4_M= 0.0
        mg_M     = 0.00039
        mgo_M    = 0.00049
        sio2_M   = 0.039
        fe3o4_M  = 0.0
        mgsio3_M = 0.0
        fe2sio4_M= 0.0
    elif (not mixed and depleted and mass==20):
        ## UM D M20
        if carbon:
            ac_M     = .0145
        else:
            ac_M     = 1.3e-13
        al2o3_M  = 0.000086
        asi_M    = 0.030
        fe_M     = 0.000046
        fes_M    = 0.033
        mg2sio4_M= 0.089
        mg_M     = 0.00039
        mgo_M    = 0.00049
        sio2_M   = 0.0
        fe3o4_M  = 0.0
        mgsio3_M = 0.0
        fe2sio4_M= 0.0
    elif (not mixed and not depleted and mass==170):
        ## UM ND M170
        ac_M     = 1.2e-13
        al2o3_M  = 0.0296
        asi_M    = 1.963
        fe_M     = 0.000067
        fes_M    = 0.011
        mg2sio4_M= 2.474
        mg_M     = 0.00025
        mgo_M    = 0.0000025
        sio2_M   = 2.577
        fe3o4_M  = 0.0
        mgsio3_M = 0.0
        fe2sio4_M= 0.0
    elif (not mixed and depleted and mass==170):
        ## UM D M170
        ac_M     = 1.2e-13
        al2o3_M  = 0.0297
        asi_M    = 1.963
        fe_M     = 0.000067
        fes_M    = 0.011
        mg2sio4_M= 0.0
        mg_M     = 0.00802
        mgo_M    = 0.0000025
        sio2_M   = 3.638
        fe3o4_M  = 0.0
        mgsio3_M = 0.0
        fe2sio4_M= 0.0
    elif (mixed and not depleted and mass==20):
        ## M ND M20
        ac_M     = 0.0
        al2o3_M  = 0.00088
        asi_M    = 0.049
        fe_M     = 0.000432
        fes_M    = 0.0
        mg2sio4_M= 0.0
        mg_M     = 0.0014
        mgo_M    = 0.0
        sio2_M   = 0.105
        fe3o4_M  = 0.0
        mgsio3_M = 0.0
        fe2sio4_M= 0.0
    elif (mixed and depleted and mass==20):
        ## M D M20
        ac_M     = 0.0
        al2o3_M  = 0.00088
        asi_M    = 0.049
        fe_M     = 0.0
        fes_M    = 0.0
        mg2sio4_M= 0.160
        mg_M     = 0.0
        mgo_M    = 0.0
        sio2_M   = 0.0
        fe3o4_M  = 0.0
        mgsio3_M = 0.0
        fe2sio4_M= 0.125
    elif (mixed and not depleted and mass==170):
        ## M ND M170
        ac_M     = 0.0
        al2o3_M  = 0.003
        asi_M    = 8.1
        fe_M     = 0.004
        fes_M    = 0.0
        mg2sio4_M= 0.0
        mg_M     = 0.0
        mgo_M    = 0.0
        sio2_M   =17.3
        fe3o4_M  = 0.0
        mgsio3_M = 0.0
        fe2sio4_M= 0.0
    elif (mixed and depleted and mass==170):
        ## M D M170
        ac_M     = 0.0
        al2o3_M  = 0.003
        asi_M    = 8.1
        fe_M     = 0.004
        fes_M    = 0.0
        mg2sio4_M= 5.7
        mg_M     = 0.0
        mgo_M    = 0.0
        sio2_M   =12.9
        fe3o4_M  = 0.0
        mgsio3_M = 0.0
        fe2sio4_M= 6.6

    if (carbfrac==None):
        total_M = ac_M+al2o3_M+asi_M+fe_M+fes_M+mg2sio4_M+mg_M+mgo_M+sio2_M+fe3o4_M+mgsio3_M+fe2sio4_M
        return (ac_M/total_M,al2o3_M/total_M,asi_M/total_M,fe_M/total_M,fes_M/total_M,mg2sio4_M/total_M,mg_M/total_M,mgo_M/total_M,sio2_M/total_M,fe3o4_M/total_M,mgsio3_M/total_M,fe2sio4_M/total_M)
    else: # don't add carbon mass from model, instead put it in artificially to a fixed level
        total_M = al2o3_M+asi_M+fe_M+fes_M+mg2sio4_M+mg_M+mgo_M+sio2_M+fe3o4_M+mgsio3_M+fe2sio4_M
        # ratios sums to 1 without carbon
        ratios = (0,al2o3_M/total_M,asi_M/total_M,fe_M/total_M,fes_M/total_M,mg2sio4_M/total_M,mg_M/total_M,mgo_M/total_M,sio2_M/total_M,fe3o4_M/total_M,mgsio3_M/total_M,fe2sio4_M/total_M)
        ratios = np.array(ratios)*(1-carbfrac) # ratios sums to 1 - carbfrac
        ratios[0] = carbfrac #ratios sums to 1
        return tuple(ratios)

def Blambda(wavelen, T):
    """
    @param wavelen: wavelength (array) in cm
    @param T: temperature in Kelvin
    @return: Blackbody intensity in cgs units
    """
    return .00001191 * (1/wavelen**5) * 1/(np.exp(1.43878/(wavelen*T))-1)

def Sgeometric(a,dNda):
    """
    @param a: size array in cm
    @param dNda: number distribution per unit dust mass
    """
    return np.pi*trapz(a**2 * dNda,a)

def find_normalization(dnda_dat, rho):
    a    = dnda_dat[:,0]/10**4 #cm
    dnda = 10**dnda_dat[:,1]   #assume unitless (all units in normalization)
    I3 = trapz(a**3 * dnda, a) #units cm**4
    return .75/(np.pi*rho*I3)  #units 1/(g cm)

def std_size_distr(rho):
    ##a in um
    ##4/23/13 edit: added .001
    a = np.array([.001,.005,.010,.050,.100,.500,1.00,2.00,3.00,4.00,5.00])
    def nstd(a): ##returns log10(n)
        P0 = .005
        if a<P0:
            return 0
        elif a<1.0:
            return 3.5*np.log10(P0/a)
        elif a<=5.0:
            return 2*np.log10(1/.005) + 5.5*np.log10(P0/a)
        else:
            return -100 #very small number
    n = np.array([nstd(thisa) for thisa in a])
    norm = find_normalization(np.array([a,n]).T, rho)
    return a/10**4, norm*10**n

def shock_size_distr(rho,x):
    # n ~ 1, n ~ a^-x
    a = np.array([.001,.005,.010,.050,.100,.500,1.00,2.00,3.00,4.00,5.00])
    def nshock(a):
        P0=.005
        if a<P0:
            return 0
        else:
            return x * (np.log10(P0) - np.log10(a))
    n = np.array([nshock(thisa) for thisa in a])
    norm = find_normalization(np.array([a,n]).T, rho)
    return a/10**4, norm*10**n ##returns a in cm

def calc_qabs(opt,a_arr,dNda,fileprefix):
    """
    @param opt: optical constants
    @param a: size array in cm
    @param dNda: number distribution per unit dust mass
    @param fileprefix: prefix of filename that data will be saved to

    Calculates Q_abs(lambda, a), which is the normalized absorption cross section
    """
    print "Calculating Qabs for "+fileprefix
    wavelen = opt[:,0] #um
    qabs = np.zeros( (len(wavelen),len(a_arr)) )
    for i,a in enumerate(a_arr):
        qabs[:,i] = calcq(opt,a*10**4) #unitless; input here is in um, so convert a to um
    np.save('DATA/'+fileprefix+'_qabs.npy',qabs)
    return qabs

def calc_klambda(opt,a_arr,dNda,fileprefix,qabs=np.array([])):
    """
    @param opt: optical constants
    @param a_arr: size array in cm
    @param dNda: number distribution per unit dust mass (1/(g cm))
    @param fileprefix: prefix of filename that data will be saved to
    @param qabs: Qabs calculated from opt. If empty, loads qabs from fileprefix
    """
    if qabs.size == 0:
        print "Loading Qabs for "+fileprefix
        qabs = np.load('DATA/'+fileprefix+'_qabs.npy')

    print "Calculating kappa_lambda for "+fileprefix
    klambda = np.zeros(len(opt[:,0]))
    area = (np.pi * a_arr**2) * np.ones( qabs.shape )
    dNda = dNda * np.ones( qabs.shape )
    klambda_a = qabs * area * dNda

    klambda = trapz(klambda_a, a_arr)
    np.save('DATA/'+fileprefix+'_klambda.npy',klambda)
    return klambda

def Ikappa(T, wavelen, klambda):
    """
    Evaluates the integral in kappaPlanck (used separately)
    T in Kelvin
    wavelen in cm
    klambda in cm^2/g
    Outputs ergs/s/cm^2 * cm^2/g
    """
    integrand = Blambda(wavelen,T) * klambda
    return trapz(integrand,wavelen)

def kappaPlanck(T, wavelen, klambda):
    """
    Outputs the blackbody-averaged opacity
    T in Kelvin
    wavelen in cm
    """
    #if T>2000: ##sublimation temperature
    #    return 1e-2 #just a very tiny number
    sigma = 5.6704 * 10**(-5) ##in cgs: ergs/(cm^2 s K)
    I = Ikappa(T, wavelen, klambda)
    return I*np.pi / (sigma * T**4)

def dopckekappa(Td,k0=900):
    if Td < 200:
        return k0*(Td/200)**2
    elif Td < 1500:
        return k0
    else:
        return k0*(Td/1500)**(-12)

def powerkappa(Td,Tflat=200,Tsubl=1500,k0=400):
    if Td < Tflat:
        return k0*(Td/Tflat)**2
    elif Td < Tsubl:
        return k0
    else:
        return 0 #k0*(Td/Tsubl)**(-12)
    
if __name__ == "__main__":
    labels = ['ac','al2o3','asi','fe','fes','mg2sio4','mg','mgo','sio2','fe3o4','mgsio3','fe2sio4']
    densities = get_densities()

    for i in xrange(len(labels)):
        a,dnda = std_size_distr(densities[i])
        print labels[i]+':', "rho =",densities[i]
        print "S =",Sgeometric(a,dnda)
