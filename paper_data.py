import numpy as np
import pylab as plt
import asciitable

from calc_dtrans import calc_dtrans

def plot_data_sife():
    # load data
    #data = asciitable.read("sih_data.csv",delimiter=',')
    data = asciitable.read("DATA/lowfe_sample.csv",delimiter=',')
    data = data[data['feh'] <= -3.5]
    data['sih'][18] = -4.85
    data['sih'][22] = -5.2
    data['ulsi'][18] = 0
    data['ulsi'][22] = 1
    x = data['feh']
    y = data['sih']
    ulindex = np.where(data['ulsi'] == 1)[0]
    
    # create figure
    plt.clf()
    fig = plt.figure(figsize=(7.5,7.5))
    
    ax = plt.subplot(111)
    ax.scatter(x, y, c='k', s=15.0, facecolors='none')
    
    # plot arrows
    arrowdx = 0.0
    arrowdy = -0.1
    arrowhw = .04
    arrowhl = .08
    for ix in ulindex:
        ax.arrow(x[ix],y[ix],arrowdx,arrowdy,fc='k',ec='k',head_width=arrowhw,head_length=arrowhl)
    
    # plot upper limit bars for HE0557 and HE1327
    #barwidth = .02
    #ax.add_patch(plt.Rectangle((-4.81-barwidth/2.,-5.4),barwidth,-4.84+5.4,color='orange'))
    #ax.arrow(-4.81,-5.4,arrowdx,arrowdy,fc='r',ec='r',head_width=arrowhw,head_length=arrowhl)
    #ax.add_patch(plt.Rectangle((-5.76-barwidth/2.,-5.4),barwidth,-5.2+5.4,color='orange'))
    #ax.arrow(-5.76,-5.4,arrowdx,arrowdy,fc='r',ec='r',head_width=arrowhw,head_length=arrowhl)
    ax.set_xlabel(r'[Fe/H]',fontsize=16)
    ax.set_ylabel(r'[Si/H]',fontsize=16)
    return fig,ax

def plot_data_dtranssih():
    data = asciitable.read("DATA/lowfe_sample.csv",delimiter=',')
    data = data[data['feh'] <= -3.5]
    data['sih'][18] = -4.85
    data['sih'][22] = -5.2
    data['ulsi'][18] = 0
    data['ulsi'][22] = 1
    dtrans, dtranslo, dtranshi, dtransul = calc_dtrans(data)
    sih = data['sih']; ulsi = data['ulsi']
    mask = np.isfinite(sih)
    sih = sih[mask]; ulsi = ulsi[mask]
    dtrans = dtrans[mask]; dtransul = dtransul[mask]
    dtranslo = dtranslo[mask]; dtranshi = dtranshi[mask]
    
    fig = plt.figure(figsize=(7.5,7.5))
    ax = fig.add_subplot(111)

    ## Both C and O
    mask = np.isnan(dtranslo)
    ax.scatter(sih[mask],dtrans[mask],marker='D',c='black',s=25)
    ## Just C or O
    mask = np.isfinite(dtranslo)
    ax.scatter(sih[mask],dtrans[mask],marker='o',facecolors='none',edgecolors='blue',s=16)
    ax.vlines(sih[mask],dtranslo[mask],dtranshi[mask],color='blue',lw=.2)
    ## Add upper limits
    arrowdx = 0.0; arrowdy = -0.1
    arrowhw = .04; arrowhl = .08
    arrowcol = 'k'
    indices = np.where(dtransul == 1)[0]
    for ix in indices:
        ax.arrow(sih[ix],dtrans[ix],arrowdx,arrowdy,fc=arrowcol,ec=arrowcol,head_width=arrowhw,head_length=arrowhl)
    indices = np.where(ulsi == 1)[0]
    for ix in indices:
        ax.arrow(sih[ix],dtrans[ix],arrowdy,arrowdx,fc=arrowcol,ec=arrowcol,head_width=arrowhw,head_length=arrowhl)

    ax.set_xlabel(r'[Si/H]',fontsize=16)
    ax.set_ylabel(r'$D_{trans}$',fontsize=16)
    return fig,ax

plt.clf()
fig,ax = plot_data_sife()
# plot [Si/H]crit lines
stdDcrit = [4.52e-8,5.28e-8,3.87e-8,5.49e-8,4.33e-8,3.88e-8,3.87e-8,4.46e-8]
x55Dcrit = [7.86e-9,1.09e-8,6.02e-9,9.84e-9,7.08e-9,6.08e-9,6.01e-9,7.17e-9]
MsiMdust = [.469,.312,.625,.293,.648,.519,.637,.486]
stdsihcrit = np.log10(stdDcrit)+np.log10(4./(3*28.1))+np.log10(MsiMdust)+12-7.51
x55sihcrit = np.log10(x55Dcrit)+np.log10(4./(3*28.1))+np.log10(MsiMdust)+12-7.51
for stdsih in stdsihcrit:
    ax.plot([-6.0,-3.0],[stdsih,stdsih],'g--')
for x55sih in x55sihcrit:
    ax.plot([-6.0,-3.0],[x55sih,x55sih],'r--')
ax.add_patch(plt.Polygon([[-6,-6],[-3,-3],[-3,-3+.4],[-6,-6+.4]],closed=True,alpha=.2,color='cyan'))

ax.set_xlim(-6.0,-3.0)
ax.set_ylim(-6.0,-2.0)
plt.savefig("PLOTS/paper_sih.pdf",bbox_inches='tight')

plt.clf()
fig,ax = plot_data_dtranssih()
#for stdsih in stdsihcrit:
#    ax.plot([stdsih,stdsih],[-4.0,0.0],'g--')
#for x55sih in x55sihcrit:
#    ax.plot([x55sih,x55sih],[-4.0,0.0],'r--')
## Dtrans limit
ax.plot([-6.0,-1.0],[-3.3,-3.3],'k:')
ax.plot([-6.0,-1.0],[-3.5,-3.5],'k--')
ax.plot([-6.0,-1.0],[-3.7,-3.7],'k:')

ax.add_patch(plt.Polygon([[min(stdsihcrit),0.0],[max(stdsihcrit),0.0],[max(stdsihcrit),-4.0],[min(stdsihcrit),-4.0]],closed=True,alpha=.6,color='green'))
ax.add_patch(plt.Polygon([[min(x55sihcrit),0.0],[max(x55sihcrit),0.0],[max(x55sihcrit),-4.0],[min(x55sihcrit),-4.0]],closed=True,alpha=.6,color='red'))
#ax.add_patch(plt.Polygon([[-6.0,-3.3],[-6.0,-3.7],[-2.0,-3.7],[-2.0,-3.3]],closed=True,alpha=.4,color='black'))

ax.text(-5.8,-2.0,"shock dust cooling fails",color='red',fontsize=16,
        horizontalalignment='center', verticalalignment='center',
        rotation='vertical')

ax.text(-4.7,-2.0,"standard dust cooling fails",color='green',fontsize=16,
        horizontalalignment='center', verticalalignment='center',
        rotation='vertical')

ax.text(-4.0,-3.8,"fine structure cooling fails",color='black',fontsize=16,
        horizontalalignment='left', verticalalignment='center')

ax.set_xlim(-6.0,-2.0)
ax.set_ylim(-4.0,0.0)
plt.savefig("PLOTS/paper_dtranssih.pdf",bbox_inches='tight')
