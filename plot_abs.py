from calcq import calcq
import numpy as np
import pylab as plt
from scipy.integrate import simps
from scipy.integrate import trapz

def Blambda(wavelen, T):        ##all in cgs units, wavelen is in cm
    return .00001191 * (1/wavelen**5) * (1/(np.exp(1.43878/(wavelen*T))-1))
#return 1.19*0.00001 * (1/wavelen**5) * (1/(np.exp(1.439/(wavelen*T))))

def plot_optical(opt,fileprefix="sio2",title='SiO2'):
    plt.clf()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(opt[:,0],opt[:,1],'o--',markersize=2)
    ax.plot(opt[:,0],opt[:,2],'o--',markersize=2)
    ax.legend(('n','k'),loc="lower left")
    plt.xlabel('wavelength (um)')
    plt.ylabel('n,k')
    plt.title(title+' Optical')
    plt.savefig(fileprefix+"_optical.pdf",bbox_inches=0)
    
def plot_qabs(opt,fileprefix="sio2",title='SiO2'):
    plt.clf()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    a_arr = [.001,.01,0.1,1.0]
    for a in a_arr:
        qabs = calcq(opt,a)
        ax.plot(opt[:,0],qabs,'o--',markersize=2)
    ax.legend(a_arr,loc="lower left")
    plt.xlabel('wavelength (um)')
    plt.ylabel('Qabs')
    plt.title(title+' Qabs(a)')
    plt.savefig(fileprefix+"_qabs.pdf",bbox_inches=0)

def find_normalization(dnda_dat, rho):
    a    = dnda_dat[:,0]/10**4 #um -> cm
    dnda = 10**dnda_dat[:,1]
    #I2 = simps(a**2 * dnda, a)
    #I3 = simps(a**3 * dnda, a)
    I3 = trapz(a**3 * dnda, a)
    #print "find_normalization: S=", (0.75*I2)/(I3*rho)
    return .75/(np.pi*rho*I3) #cm^3/g

def plot_klambda(opt,dnda_dat,normalization,fileprefix="sio2",title='SiO2',loadfile=False,savefile=True):
    print "Plotting klambda for "+fileprefix
    wavelen   = opt[:,0] #um
    if loadfile == False:
        a_arr = dnda_dat[:,0] #um
        klambda   = np.zeros(len(wavelen))
        
        qabs = np.zeros( (len(wavelen),len(a_arr)) )
        for i,a in enumerate(a_arr):
            qabs[:,i] = calcq(opt,a) # dimensionless
        
        area = (np.pi * a_arr**2) * np.ones( qabs.shape ) # um^2
        dnda = normalization * (10**dnda_dat[:,1]) * np.ones( qabs.shape ) # number/cm^3/um -> number/g/um
        print normalization * (10**dnda_dat[:,1])
        print dnda[0,:]
        print area[0,:]
        
        klambda_a = qabs * area * dnda # um / g
        for i in xrange(len(klambda)):
            #klambda[i] = simps(klambda_a[i,:],a_arr) # um^2 / g
            klambda[i] = trapz(klambda_a[i,:],a_arr) # um^2 / g
        klambda = klambda * 10**(-8) #um^2/g --> cm^2/g
        if savefile:
            print "Saving qabs, klambda, klambda_a..."
            np.save('DATA/'+fileprefix+'_klambda_a.npy',klambda_a)
            np.save('DATA/'+fileprefix+'_klambda.npy',klambda)
            np.save('DATA/'+fileprefix+'_qabs.npy',qabs)
    else:
        print "Loading data from previously saved file..."
        klambda_a = np.load(fileprefix+'_klambda_a.npy')
        klambda   = np.load(fileprefix+'_klambda.npy')
        qabs      = np.load(fileprefix+'_qabs.npy')


    plt.clf()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(wavelen,klambda,'o--',markersize=2)
    plt.xlabel('wavelength (um)')
    plt.ylabel('kappa (cm^2/g)')
    plt.title(title+' Kappa')
    plt.savefig(fileprefix+'_klambda.pdf',bbox_inches=0)

    plt.clf()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    for i in xrange(qabs.shape[1]):
        ax.plot(wavelen,qabs[:,i],'o--',markersize=2)
    plt.xlabel('wavelength (um)')
    plt.ylabel('Qabs')
    plt.title(title+' Qabs(lambda; a)')
    plt.savefig(fileprefix+'_klambda_qabsa.pdf',bbox_inches=0)

def kappaPlanck(Tdust, opt, klambda):
    wavelen = opt[:,0] * 10**(-4) ##um to cm
    sigma = 5.6704 * 10**(-5) ##in cgs: ergs/(cm^2 s K)
    integrand = Blambda(wavelen,Tdust) * klambda
    return np.pi / (sigma * Tdust**4) * simps(integrand,wavelen)

def plot_kappaPlanck(opt,fileprefix="sio2",title='SiO2'):
    import os.path
    savefile = 'DATA/'+fileprefix+'_klambda.npy'
    if os.path.exists(savefile):
        klambda = np.load(savefile)
    else:
        print "File not found: "+savefile
        import sys; sys.exit()

    plt.clf()
    ax = plt.gca()
    Tarr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000]
    kP = np.zeros(len(Tarr))
    for i,Tdust in enumerate(Tarr):
        kP[i] = kappaPlanck(Tdust, opt, klambda)
    
    ax.plot(Tarr,kP,'o--',markersize=2)
    plt.xlabel('Dust Temperature (K)')
    plt.ylabel('kappa Planck')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.title(title+' kappa Planck')
    plt.savefig(fileprefix+"_kappaPlanck.pdf",bbox_inches=0)

def plot_bbody(opt,fileprefix="bbody.pdf"):
    wavelen = opt[:,0] * 0.0001 #um to cm
    plt.clf()
    ax = plt.gca()
    ax.set_xscale('log')
    temperature_arr = [50,100,150,200]
    for t in temperature_arr:
        y = Blambda(wavelen,t)
        y = y/np.max(y) #normalize
        ax.plot(opt[:,0],y)
    ax.legend(temperature_arr,loc="upper left")
    plt.xlabel('wavelength (um)')
    plt.ylabel('B (normalized)')
    plt.title('Blackbody Radiation')
    plt.savefig(fileprefix,bbox_inches=0)

if __name__ == "__main__":
    rho_sio2 = 2.6
    rho_ac   = 2.0

    M22     = (0.23+0.27+0.03+0.21+.005)*1.989e33
    M195    = (13+10+3+.8+.04)*1.989e33

    sio2_optical = np.loadtxt(open("sio2.csv",'r'),delimiter=',',skiprows=1)
    sio2_dnda195 = np.loadtxt(open("sio2_dnda_195msun.csv",'r'),delimiter=',',skiprows=1)
    norm_sio2195 = find_normalization(sio2_dnda195,rho_sio2)
##     #print norm_sio2195
    plot_klambda(sio2_optical,sio2_dnda195, norm_sio2195, fileprefix='sio2_195', title='SiO2 M195',loadfile=False)
    plot_kappaPlanck(sio2_optical,fileprefix="sio2_195",title="SiO2 M195")

    ac_optical = np.loadtxt(open("ac.csv",'r'),delimiter=',',skiprows=1)
    ac_dnda22  = np.loadtxt(open("ac_dnda_22msun.csv",'r'), delimiter=',',skiprows=1)
    ac_dnda195 = np.loadtxt(open("ac_dnda_195msun.csv",'r'),delimiter=',',skiprows=1)
    norm_ac195   = find_normalization(ac_dnda195,rho_ac)
    norm_ac22    = find_normalization(ac_dnda22,rho_ac)
    #print norm_ac195
    #print norm_ac22
    #print ac_dnda195[:,1]
    #print ac_dnda22[:,1]
    #print norm_ac195 * 10**ac_dnda195[:,1]
    #print norm_ac22  * 10**ac_dnda22[:,1]
    plot_klambda(ac_optical,ac_dnda22, norm_ac22,  fileprefix='ac_22', title='AC M22',loadfile=False)
    plot_klambda(ac_optical,ac_dnda195,norm_ac195, fileprefix='ac_195',title='AC M195',loadfile=False)
    plot_kappaPlanck(ac_optical,fileprefix="ac_22",title="AC M22")
    plot_kappaPlanck(ac_optical,fileprefix="ac_195",title="AC M195")

    ## With DataThief (really long computation)
##     sio2_dnda195 = np.loadtxt(open("sio2_195_dt_dnda.txt",'r'),delimiter=',',skiprows=1)
##     norm_sio2195 = find_normalization(sio2_dnda195,rho_sio2)
##     plot_klambda(sio2_optical,sio2_dnda195, norm_sio2195, fileprefix='sio2_dt_195', title='SiO2 M195',loadfile=False)
##     plot_kappaPlanck(sio2_optical,fileprefix="sio2_dt_195",title="SiO2 M195")

##     ac_dnda22  = np.loadtxt(open("ac_22_dt_dnda.txt",'r'), delimiter=',',skiprows=1)
##     ac_dnda195 = np.loadtxt(open("ac_195_dt_dnda.txt",'r'),delimiter=',',skiprows=1)
##     norm_ac195   = find_normalization(ac_dnda195,rho_ac)
##     norm_ac22    = find_normalization(ac_dnda22,rho_ac)
##     plot_klambda(ac_optical,ac_dnda22, norm_ac22,  fileprefix='ac_dt_22', title='AC M22',loadfile=False)
##     plot_klambda(ac_optical,ac_dnda195,norm_ac195, fileprefix='ac_dt_195',title='AC M195',loadfile=False)
##     plot_kappaPlanck(ac_optical,fileprefix="ac_dt_22",title="AC M22")
##     plot_kappaPlanck(ac_optical,fileprefix="ac_dt_195",title="AC M195")

    
    #plot_optical(sio2_optical)
    #plot_qabs(sio2_optical)
    #plot_optical(ac_optical,fileprefix='ac',title='AC')
    #plot_qabs(ac_optical,fileprefix='ac',title='AC')
    #plot_bbody(ac_optical)
