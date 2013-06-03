import numpy as np
import asciitable

## grid[row][col], grid[logT][logn]
_global_kappa_gas_grid = asciitable.read("DATA/md05_opacity.csv",delimiter=',')

def kappa_gas(n,T):
    logn = np.log10(n)
    logT = np.log10(T)

    mu = 4./3.          #1+4/12
    mP = 1.67 * 10**-24 #grams
    rho2n = -1*np.log10(mu * mP)
    logTarr = np.arange(1.8,3.9+.1,.1)
    lognarr = np.arange(-16+rho2n,-2+rho2n+1,1)

    grid = _global_kappa_gas_grid

    ixT = np.searchsorted(logTarr,logT)
    ixn = np.searchsorted(lognarr,logn)
    try:
        ixT = ixT[0]
    except TypeError:
        pass
    except IndexError:
        pass

    if ixn == 0:
        return 10.**(-33)
    if ixT == 0:
        # m(x-x0)+y0 + n-n0, linear interpolation
        m = (grid[1][0] - grid[0][0])/(logTarr[1]-logTarr[0])
        return 10.**(m*(logT - logTarr[0]) + grid[0][0] + logn - lognarr[0])

    if ixT >= len(logTarr) - 1:
        print "kappa_gas: error, T too large"
        return 10.**(-33.)
    if ixn >= len(lognarr) - 1:
        print "kappa_gas: error, n too large"
        return 10.**(-33.)

    # interpolate in T: smaller n
    logT1 = logTarr[ixT]
    logT2 = logTarr[ixT+1]
    #print ixT,ixn
    logk1 = grid[ixT][ixn]
    logk2 = grid[ixT+1][ixn]
    m = (logk2-logk1)/(logT2-logT1)
    logksmall = m*(logT - logT1) + logk1
    # interpolate in T: larger n
    logk1 = grid[ixT][ixn+1]
    logk2 = grid[ixT+1][ixn+1]
    m = (logk2-logk1)/(logT2-logT1)
    logklarge = m*(logT - logT1) + logk1

    # interpolate in n
    lognsmall = lognarr[ixn]
    lognlarge = lognarr[ixn+1]
    m = (logklarge-logksmall)/(lognlarge-lognsmall)
    return 10.**(m*(logn-lognsmall) + logksmall)
