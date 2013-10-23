import numpy as np
from dust_calculations import get_wavelen
from cooling import get_kappafn
import asciitable

import pylab as plt

if __name__=="__main__":
    dustmodels = ['UM-ND-20','UM-D-20','M-ND-20', 'M-D-20',
                  'UM-ND-170','UM-D-170','M-ND-170', 'M-D-170']
    stdmodels = ['std_'+x for x in dustmodels]
    x55models = ['x5.5_'+x for x in dustmodels]

    logTarr = np.linspace(0,3.5,351)
    Tarr = 10**logTarr
    wavelen,lw = get_wavelen()
    for model in (stdmodels+x55models):
        print "Computing for",model
        klambda,S = np.load("DATA/"+model+"_klambda.npy")
        kappafn = get_kappafn(wavelen,klambda,Tsub=3500)
        kappa = np.array([kappafn(T) for T in Tarr])
        asciitable.write({'logT':logTarr,'kappa':kappa},'DATA/kappa_'+model+'.dat',names=['logT','kappa'])

    plt.figure()
    plt.subplot(121)
    for model in (stdmodels):
        d = asciitable.read("DATA/kappa_"+model+".dat")
        plt.plot(d['logT'],np.log10(d['kappa']))
    plt.xlim((1,4))
    plt.ylim((-1,4))
    plt.subplot(122)
    for model in (x55models):
        d = asciitable.read("DATA/kappa_"+model+".dat")
        plt.plot(d['logT'],np.log10(d['kappa']))
    plt.xlim((1,4))
    plt.ylim((-1,4))
    plt.show()
