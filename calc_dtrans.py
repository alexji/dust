import numpy as np
import asciitable
import pylab as plt

def findllum(data):
    """
    Luminosity estimate based on Teff and logg
    """
    llum = np.zeros(len(data))
    #llum[np.isfinite(data['llum'])] = data['llum'][np.isfinite(data['llum'])]
    #indices = np.where(np.isnan(data['llum']))[0]
    mask = np.where(np.isnan(data['teff']))[0]
    llum[mask] = data['llum'][mask]
    indices = np.where(np.isfinite(data['teff']))[0]
    for i in indices:
        teff = data['teff'][i]
        logg = data['logg'][i]
        lum = (4*np.pi*6.673e-11*1.989e30*5.67e-8*teff**4 * 100)/(10**logg * 3.845e26)
        if teff < 4600:
            lum = lum #* 0.8/0.8
        elif teff < 5800:
            lum = lum * 0.7/0.8
        else:
            lum = lum * 0.6/0.8
        llum[i] = np.log10(lum)
    return llum

def findcorr(llum):
    """
    Correction to [C/H] based on luminosity
    """
    if llum < 2.0:
        return 0.0
    elif llum < 3.1:
        return -0.6 * llum + 1.1
    else:
        return -0.76

def calc_dtrans(data):
    def dtransfn(ch,oh):
        return np.log10(10.**ch + 0.9 * 10**oh)
    dtrans = np.empty(len(data))   * np.nan
    dtranslo = np.empty(len(data)) * np.nan
    dtranshi = np.empty(len(data)) * np.nan
    dtransul= np.zeros(len(data))
    llum = findllum(data)
    corr = np.zeros(len(llum))
    for i,ll in enumerate(llum):
        corr[i] = findcorr(ll)
    for i,row in enumerate(data):
        ch = row['ch']-corr[i]
        oh = row['oh']
        ulo= row['ulo']; ulc = row['ulc']
        ## Check upper limits
        if ulo == 1 or ulc == 1:
            dtransul[i] = 1
        ## Both C and O
        if np.isfinite(ch) and np.isfinite(oh):
            dtrans[i] = dtransfn(ch,oh)
            ## dtranslo and dtranshi are nan
        ## Only C
        if np.isfinite(ch) and np.isnan(oh):
            ohmin = ch - 0.2
            ohmax = ch + 0.7
            dtrans[i] = dtransfn(ch,ohmin)
            dtranslo[i] = dtrans[i]
            dtranshi[i] = dtransfn(ch,ohmax)
        ## Only O
        if np.isnan(ch) and np.isfinite(oh):
            chmin = oh - 0.7
            chmax = oh + 0.2
            dtrans[i] = dtransfn(chmin,oh)
            dtranslo[i] = dtrans[i]
            dtranshi[i] = dtransfn(chmax,oh)
    return dtrans, dtranslo, dtranshi, dtransul


def plot_dtrans(feh,dtrans,dtranslo,dtranshi,dtransul):
    fig = plt.figure(1,figsize=(7.5,7.5))
    ax = fig.add_subplot(111)
    ## Both C and O
    mask = np.isnan(dtranslo)
    ax.scatter(feh[mask],dtrans[mask],marker='D',c='black',s=25)
    ## Just C or O
    mask = np.isfinite(dtranslo)
    ax.scatter(feh[mask],dtrans[mask],marker='o',facecolors='none',edgecolors='red',s=16)
    ax.vlines(feh[mask],dtranslo[mask],dtranshi[mask],color='r',lw=.2)
    ## Add upper limits
    indices = np.where(dtransul == 1)[0]
    arrowdx = 0.0; arrowdy = -0.06
    arrowhw = .03; arrowhl = .06
    for ix in indices:
        ax.arrow(feh[ix],dtrans[ix],arrowdx,arrowdy,fc='k',ec='k',head_width=arrowhw,head_length=arrowhl)
    
    plt.plot([-6.0,-3.0],[-3.3,-3.3],'k:')
    plt.plot([-6.0,-3.0],[-3.5,-3.5],'k--')
    plt.plot([-6.0,-3.0],[-3.7,-3.7],'k:')
    def dtranssolar(feh):
        return np.log10(10**feh + 0.9*10**feh)
    plt.plot([-4.5,-3.0],[dtranssolar(-4.5),dtranssolar(-3.0)],'k-')

    plt.xlim((-6.0,-3.4))
    plt.ylim((-4.0, 0.0))
    plt.xlabel(r'[Fe/H]')
    plt.ylabel(r'$D_{trans}$')
    plt.savefig("PLOTS/dtrans.pdf")

def plot_dtranssih(sih,ulsi,dtrans,dtranslo,dtranshi,dtransul):
    fig = plt.figure(1,figsize=(7.5,7.5))
    ax = fig.add_subplot(111)
    ## Cut by presence of Si
    mask = np.isfinite(sih)
    sih = sih[mask]; ulsi = ulsi[mask]
    dtrans = dtrans[mask]; dtransul = dtransul[mask]
    dtranslo = dtranslo[mask]; dtranshi = dtranshi[mask]
    ## Both C and O
    mask = np.isnan(dtranslo)
    ax.scatter(sih[mask],dtrans[mask],marker='D',c='black',s=25)
    ## Just C or O
    mask = np.isfinite(dtranslo)
    ax.scatter(sih[mask],dtrans[mask],marker='o',facecolors='none',edgecolors='red',s=16)
    ax.vlines(sih[mask],dtranslo[mask],dtranshi[mask],color='r',lw=.2)
    ## Add upper limits
    arrowdx = 0.0; arrowdy = -0.06
    arrowhw = .03; arrowhl = .06
    arrowcol = 'k'
    indices = np.where(dtransul == 1)[0]
    for ix in indices:
        ax.arrow(sih[ix],dtrans[ix],arrowdx,arrowdy,fc=arrowcol,ec=arrowcol,head_width=arrowhw,head_length=arrowhl)
    indices = np.where(ulsi == 1)[0]
    for ix in indices:
        ax.arrow(sih[ix],dtrans[ix],arrowdy,arrowdx,fc=arrowcol,ec=arrowcol,head_width=arrowhw,head_length=arrowhl)
    ## Add more silicon upper limits
    barwidth = .02
    # HE0557
    arrowcol = 'orange'
    he0557dtrans = -2.94
    he0557dtransdiff = -2.40+2.94
    ax.add_patch(plt.Rectangle((-5.4,he0557dtrans),-4.84+5.4,he0557dtransdiff,color='orange',alpha=.5))
    ax.add_patch(plt.Rectangle((-5.4,he0557dtrans),-4.84+5.4,barwidth,color='orange'))
    ax.arrow(-5.4,he0557dtrans,arrowdy,arrowdx,fc=arrowcol,ec=arrowcol,head_width=arrowhw,head_length=arrowhl)
    ax.arrow(-5.4,he0557dtrans,arrowdx,arrowdy,fc=arrowcol,ec=arrowcol,head_width=arrowhw,head_length=arrowhl)
    # HE1327
    he1327dtrans = -1.37
    ax.add_patch(plt.Rectangle((-5.4,he1327dtrans-barwidth/2),-5.2+5.4,barwidth,color='orange'))
    ax.arrow(-5.4,he1327dtrans,arrowdy,arrowdx,fc=arrowcol,ec=arrowcol,head_width=arrowhw,head_length=arrowhl)

    
    ## Dtrans limit
    ax.plot([-6.0,-1.0],[-3.3,-3.3],'k:')
    ax.plot([-6.0,-1.0],[-3.5,-3.5],'k--')
    ax.plot([-6.0,-1.0],[-3.7,-3.7],'k:')
    ## [Si/H] limit
    ax.plot([-5.53,-5.53],[-4.0,0.0],'b:')
    ax.plot([-5.33,-5.33],[-4.0,0.0],'b--')
    ax.plot([-5.13,-5.13],[-4.0,0.0],'b:')
    ## Solar abundance line
    def dtranssolar(feh):
        return np.log10(10**feh + 0.9*10**feh)
    ax.plot([-4.5,-1.0],[dtranssolar(-4.5),dtranssolar(-1.0)],'k-')

    plt.xlim((-6.0,-1.0))
    plt.ylim((-4.0, 0.0))
    plt.xlabel(r'[Si/H]')
    plt.ylabel(r'$D_{trans}$')
    plt.savefig("PLOTS/dtranssih.pdf")

def plot_sife(feh,sih,ulsi):
    sife = sih - feh
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(feh,sife)

    ## Add upper limits
    arrowdx = 0.0; arrowdy = -0.1
    arrowhw = .04; arrowhl = .08
    arrowcol = 'red'
    indices = np.where(ulsi == 1)[0]
    for ix in indices:
        ax.arrow(feh[ix],sife[ix],arrowdx,arrowdy,fc=arrowcol,ec=arrowcol,head_width=arrowhw,head_length=arrowhl)

    barwidth = .02
    ax.add_patch(plt.Rectangle((-4.81-barwidth/2.,-5.4+4.81),barwidth,-4.84+5.4,color='orange'))
    ax.arrow(-4.81,-5.4+4.81,arrowdx,arrowdy,fc=arrowcol,ec=arrowcol,head_width=arrowhw,head_length=arrowhl)
    ax.add_patch(plt.Rectangle((-5.76-barwidth/2.,-5.4+5.76),barwidth,-5.2+5.4,color='orange'))
    ax.arrow(-5.76,-5.4+5.76,arrowdx,arrowdy,fc=arrowcol,ec=arrowcol,head_width=arrowhw,head_length=arrowhl)

    plt.plot([-6.0,-3.0],[0.4,0.4],'k:')
    plt.xlim((-6.0,-3.0))
    plt.xlabel('[Fe/H]')
    plt.ylabel('[Si/Fe]')
    #plt.show()
    fig.savefig("PLOTS/sife.pdf")

if __name__ == "__main__":
    data = asciitable.read("DATA/lowfe_sample.csv",delimiter=',')
    data = data[data['feh'] <= -3.5]
    #plot_sife(data['feh'],data['sih'],data['ulsi'])
    #print data[data['sih']-data['feh'] > 3]
    #print data['star']
    dtrans, dtranslo, dtranshi, dtransul = calc_dtrans(data)
    print dtrans
    #print data[data['star']=='HE1424-0241']
    #print dtrans[data['star']=='HE1424-0241'],dtranshi[data['star']=='HE1424-0241']
    #print dtrans[data['star']=='HE1327-23261D']
    #print dtranslo[data['star']=='HE0557-4840'],dtranshi[data['star']=='HE0557-4840']

    #plot_dtrans(data['feh'], dtrans, dtranslo, dtranshi, dtransul)
    plt.clf()
    plot_dtranssih(data['sih'], data['ulsi'], dtrans, dtranslo, dtranshi, dtransul)
