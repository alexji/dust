import numpy as np
import pylab as plt

from dust_calculations import get_wavelen
from cooling import get_kappafn, find_Td, find_Td_D

def Gamma_ad(n,T):
    kB = 1.38*10**-16
    G  = 6.67*10**-8
    mu = 4.0/3.0 #1 + 4*yHe
    mP = 1.67*10**-24
    tff= np.sqrt(3*np.pi/(32*G*mu*mP*n))
    return 1.5*n*kB*T/tff

def calc_model(heating,dustmodel,verbose=False):
    wavelen,lw = get_wavelen()
    klambda,S = np.load("DATA/"+dustmodel+"_klambda.npy")
    kappafn = get_kappafn(wavelen,klambda,Tsub=1500)
    
    Td,Hd,Ld = find_Td(n,T,S,kappafn)

    Dcrit = heating/Ld
    if verbose:
        print "---",dustmodel,"---"
        print " Td=",Td
        print " S=",S
        print " kappa(Td)=",kappafn(Td)
        print " Ld/D=",Ld
        print " Dcrit=",Dcrit
    return Td,S,kappafn(Td),Dcrit

if __name__ == "__main__":
    #wavelen, logwavelen = get_wavelen()
#    dustmodels = ['UM-D-20','UM-ND-20','UM-D-170','UM-ND-170',
#                  'M-D-20', 'M-ND-20', 'M-D-170', 'M-ND-170']
    dustmodels = ['UM-ND-20','UM-D-20','M-ND-20', 'M-D-20',
                  'UM-ND-170','UM-D-170','M-ND-170', 'M-D-170']
                  #'NoC_Caff20','NoC_Caff35',
                  #'Caff20','Caff35',
                  #'AC-UM-ND-20','AC-UM-D-20',
                  #'SOIF06-PISN','SOIF06-CCSN']
    stdmodels = ['std_'+x for x in dustmodels]
    x55models = ['x5.5_'+x for x in dustmodels]

    N = len(dustmodels)
    
    n = 10**12 #cm^-3
    T = 1000   #K
    print 'n=%e cm^-3, T=%d K' % (n, T)

    heating = Gamma_ad(n,T) #erg/sec/cm^3
    print "heating %e erg/cm^3" % heating

    std_Td = np.zeros(N)
    x55_Td = np.zeros(N)

    std_S  = np.zeros(N)
    x55_S  = np.zeros(N)

    std_kP = np.zeros(N)
    x55_kP = np.zeros(N)

    std_Dc = np.zeros(N)
    x55_Dc = np.zeros(N)

    for i in xrange(N):
        std_Td[i],std_S[i],std_kP[i],std_Dc[i] = calc_model(heating,stdmodels[i])
        x55_Td[i],x55_S[i],x55_kP[i],x55_Dc[i] = calc_model(heating,x55models[i])
        print dustmodels[i],std_Dc[i],x55_Dc[i]

    ind = np.arange(N)+1
    width = 0.35
    
    std_S  = np.log10(std_S)
    x55_S  = np.log10(x55_S)
    std_kP = np.log10(std_kP)
    x55_kP = np.log10(x55_kP)
    std_Dc = np.log10(std_Dc)
    x55_Dc = np.log10(x55_Dc)

    fig = plt.figure(1,figsize=(5,10))
    fig.subplots_adjust(bottom=.25,wspace=0,hspace=.12)

    ax  = fig.add_subplot(411)
    ax.plot(ind,x55_Td,'s-',color='r',alpha=0.6,mec='r')
    ax.plot(ind,std_Td,'D-',color='g',alpha=0.6,mec='g',mfc='None',mew=1.2)
    plt.ylabel(r'$T_d$ [K]',fontsize=14)
    ax.set_xticklabels(tuple(['' for x in np.arange(len(dustmodels))]))
    ax.set_xlim((ind[0]-.5,ind[N-1]+.5))
    plt.yticks((300,400,500,600,700))

    ax = fig.add_subplot(412)
    ax.plot(ind,x55_kP,'s-',color='r',alpha=0.6,mec='r')
    ax.plot(ind,std_kP,'D-',color='g',alpha=0.6,mec='g',mfc='None',mew=1.2)
    plt.ylabel(r'$\log_{10}{\kappa_P}$ [cm$^2$/g]',fontsize=14)
    ax.set_ylim((2.0,3.0))
    ax.set_xticklabels(tuple(['' for x in np.arange(len(dustmodels))]))
    ax.set_xlim((ind[0]-.5,ind[N-1]+.5))
    ax.legend(('shock','standard'),prop={'size':14},frameon=False,loc="lower right")

    ax = fig.add_subplot(413)
    ax.plot(ind,x55_S,'s-',color='r',alpha=0.6,mec='r')
    ax.plot(ind,std_S,'D-',color='g',alpha=0.6,mec='g',mfc='None',mew=1.2)
    plt.ylabel(r'$\log_{10}{S}$ [cm$^2$/g]',fontsize=14)
    ax.set_xticklabels(tuple(['' for x in np.arange(len(dustmodels))]))
    ax.set_xlim((ind[0]-.5,ind[N-1]+.5))
    plt.yticks((4.5,5.0,5.5,6.0))

    ax = fig.add_subplot(414)
    ax.plot(ind,x55_Dc,'s-',color='r',alpha=0.6,mec='r')
    ax.plot(ind,std_Dc,'D-',color='g',alpha=0.6,mec='g',mfc='None',mew=1.2)
    plt.plot([-1,N+1],[-8.357,-8.357],'k:')
    plt.ylabel(r'$\log_{10}{\mathcal{D}_{crit}}$',fontsize=14)
    #ax.set_ylim((-8.5,-7))
    plt.yticks((-7.0,-7.5,-8.0,-8.5))
    ax.set_xlim((ind[0]-.5,ind[N-1]+.5))
    #ax.invert_yaxis()
    plt.xlabel(r'Dust Model Number')

    for ax in fig.get_axes():
        ax.tick_params(axis='y',which='major',labelsize=12)

    plt.savefig("PLOTS/paper_compare_sizedistr.pdf",bbox_inches='tight')
