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
##     dustmodels = ['UM-D-20','UM-ND-20','UM-D-170','UM-ND-170',
##                   'M-D-20', 'M-ND-20', 'M-D-170', 'M-ND-170']
    dustmodels = ['UM-ND-20','UM-D-20','M-ND-20', 'M-D-20',
                  'UM-ND-170','UM-D-170','M-ND-170', 'M-D-170']
                  #'SOIF06-PISN','SOIF06-CCSN',
                  #'Caff20','Caff35']
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

    ind = np.arange(N)
    width = 0.35
    
    std_S  = np.log10(std_S)
    x55_S  = np.log10(x55_S)
    std_kP = np.log10(std_kP)
    x55_kP = np.log10(x55_kP)
    std_Dc = np.log10(std_Dc)
    x55_Dc = np.log10(x55_Dc)

    fig = plt.figure(1,figsize=(11,8.5))
    fig.subplots_adjust(bottom=.15,wspace=.3,hspace=.3)

    ax  = fig.add_subplot(221)
    rects1 = ax.bar(ind,      std_Td,width,color='g')
    rects2 = ax.bar(ind+width,x55_Td,width,color='r')
    plt.ylabel(r'$T_d$ [K]',fontsize=16)
    ax.set_xticks(ind+width)
    ax.set_xticklabels(tuple(dustmodels),rotation=90,size='x-small')

    ax = fig.add_subplot(222)
    rects1 = ax.bar(ind,      std_S,width,color='g')
    rects2 = ax.bar(ind+width,x55_S,width,color='r')
    plt.ylabel(r'$\log_{10}{S}$ [cm$^2$/g]',fontsize=16)
    ax.set_xticks(ind+width)
    ax.set_xticklabels(tuple(dustmodels),rotation=90,size='x-small')
    ax.set_ylim((3,6))

    ax = fig.add_subplot(223)
    rects1 = ax.bar(ind,      std_kP,width,color='g')
    rects2 = ax.bar(ind+width,x55_kP,width,color='r')
    plt.ylabel(r'$\log_{10}{\kappa_P}$ [cm$^2$/g]',fontsize=16)
    ax.set_xticks(ind+width)
    ax.set_xticklabels(tuple(dustmodels),rotation=90,size='x-small')
    ax.set_ylim((1.5,3.5))
    ax.legend((rects1[0],rects2[0]),('std','shock'),prop={'size':16},frameon=False)

    ax = fig.add_subplot(224)
    rects1 = ax.bar(ind,      std_Dc,width,color='g')
    rects2 = ax.bar(ind+width,x55_Dc,width,color='r')
    plt.plot([0,N],[-8.357,-8.357],'k:')
    plt.ylabel(r'$\log_{10}{\mathcal{D}_{crit}}$',fontsize=16)
    ax.set_xticks(ind+width)
    ax.set_xticklabels(tuple(dustmodels),rotation=90,size='x-small')
    ax.set_ylim((-9,-6))

    for ax in fig.get_axes():
        ax.tick_params(axis='y',which='major',labelsize=16)

    plt.savefig("PLOTS/paper_compare_sizedistr.pdf",bbox_inches='tight')
