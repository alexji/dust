import numpy as np
import pylab as plt
from scipy.integrate import trapz
from scipy.integrate import simps

#columns are grain radius in um, log(n) in cm^-3
ac195   = np.loadtxt(open("ac_195_dt_dnda.txt"),delimiter=',',skiprows=1)
ac22    = np.loadtxt(open("ac_22_dt_dnda.txt"),delimiter=',',skiprows=1)
sio2195 = np.loadtxt(open("sio2_195_dt_dnda.txt"),delimiter=',',skiprows=1)

#convert to cm, n
ac195[:,0]   = ac195[:,0]/10**4
ac195[:,1]   = 10**ac195[:,1]
ac22[:,0]    = ac22[:,0]/10**4
ac22[:,1]    = 10**ac22[:,1]
sio2195[:,0] = sio2195[:,0]/10**4
sio2195[:,1] = 10**sio2195[:,1]

#pick densities in g cm^-3
rho_ac = 2
rho_si = 2.6

I2ac195 = trapz(ac195[:,0]**2 * ac195[:,1], ac195[:,0])
I3ac195 = trapz(ac195[:,0]**3 * ac195[:,1], ac195[:,0])
I2si195 = trapz(sio2195[:,0]**2 * sio2195[:,1], sio2195[:,0])
I3si195 = trapz(sio2195[:,0]**3 * sio2195[:,1], sio2195[:,0])
I2ac22  = trapz(ac22[:,0]**2 * ac22[:,1], ac22[:,0])
I3ac22  = trapz(ac22[:,0]**3 * ac22[:,1], ac22[:,0])

S_ac195 = .75*I2ac195/(rho_ac*I3ac195)
S_si195 = .75*I2si195/(rho_si*I3si195)
S_ac22  = .75*I2ac22 /(rho_ac*I3ac22)

print "S for AC 195:  ", S_ac195
print "S for SiO2 195:", S_si195
print "S for AC 22:   ", S_ac22


## plt.clf()
## plt.subplot(311)
## plt.plot(ac195[:,0],ac195[:,1],'o')
## plt.ylabel("n (cm^-3 cm^-1)")
## plt.xlim((0,.0001))
## plt.legend(["AC195"])
## plt.subplot(312)
## plt.plot(ac22[:,0],ac22[:,1],'o')
## plt.ylabel("n (cm^-3 cm^-1)")
## plt.xlim((0,.00002))
## plt.legend(["AC22"])
## plt.subplot(313)
## plt.plot(sio2195[:,0],sio2195[:,1],'o')
## plt.xlabel("a (cm)")
## plt.ylabel("n (cm^-3 cm^-1)")
## plt.legend(["SiO2195"])
## plt.xlim((0,.00002))
## plt.savefig("dnda.pdf")

## #individual dust mass in grams
## #Mac22   = 0.23*1.989e33
## #Mac195  = 3.00*1.989e33
## #Msio2195= 13.0*1.989e33

## #total dust mass in grams
## M22     = (0.23+0.27+0.03+0.21+.005)*1.989e33
## M195    = (13+10+3+.8+.04)*1.989e33

## #ejecta volume when dust models have stopped running
## V22     = 8e48
## V195    = 4.2e48

## print "--------------------------------------"
## print "22: ",V22/M22
## print "195: ",V195/M195
## print "--------------------------------------"
## print "--------------------------------------"
## print "Full trapz"
## print "--------------------------------------"
## I0ac195 = trapz(ac195[:,1],ac195[:,0])
## I2ac195 = np.pi * trapz(ac195[:,0]**2 * ac195[:,1], ac195[:,0])
## I3ac195 = np.pi * trapz(ac195[:,0]**3 * ac195[:,1], ac195[:,0]) * 4./3.
## I2si195 = np.pi * trapz(sio2195[:,0]**2 * sio2195[:,1], sio2195[:,0])
## I3si195 = np.pi * trapz(sio2195[:,0]**3 * sio2195[:,1], sio2195[:,0]) * 4./3.
## I2ac22  = np.pi * trapz(ac22[:,0]**2 * ac22[:,1], ac22[:,0])
## I3ac22  = np.pi * trapz(ac22[:,0]**3 * ac22[:,1], ac22[:,0]) * 4./3.
## print "AC, M195 (SFS04)"
## print "I0, I2, I3: ", I0ac195, I2ac195, I3ac195
## print "My S=I2*V/M: ", I2ac195 * V195/M195
## print "Actual S =8.06e4 / My S ", 8.06e4 / (I2ac195 * V195/M195)
## print "Their grain density (g/cm^3): ", I2ac195/I3ac195 / 8.06e4
## #print "NEW: I0*V^2/M*I2 = ", I0ac195*V195**2/M195 * I2ac195
## print ""
## print "SiO2, M195 (SFS04)"
## print "I2, I3: ", I2si195, I3si195
## print "My S=I2*V/M: ", I2si195 * V195/M195
## print "Actual S =5.58e4 / My S ", 5.58e4 / (I2si195 * V195/M195)
## print "Their grain density (g/cm^3): ", I2si195/I3si195 / 5.58e4
## print ""
## print "AC, M22 (TF01)"
## print "I2, I3: ", I2ac22, I3ac22
## print "My S=I2*V/M: ", I2ac22 * V22/M22
## print "Actual S =8.50e4 / My S ", 8.50e4 / (I2ac22 * V22/M22)
## print "Their grain density (g/cm^3): ", I2ac22/I3ac22 / 8.50e4


## ## print "--------------------------------------"
## ## print "Cut off trapz"
## ## print "--------------------------------------"
## ## I2ac195 = np.pi * trapz(ac195[ac195[:,1] > .0001,0]**2 * ac195[ac195[:,1] > .0001,1], ac195[ac195[:,1] > .0001,0])
## ## I3ac195 = np.pi * trapz(ac195[ac195[:,1] > .0001,0]**3 * ac195[ac195[:,1] > .0001,1], ac195[ac195[:,1] > .0001,0]) * 4./3.
## ## I2si195 = np.pi * trapz(sio2195[sio2195[:,1] > .0001,0]**2 * sio2195[sio2195[:,1] > .0001,1], sio2195[sio2195[:,1] > .0001,0])
## ## I3si195 = np.pi * trapz(sio2195[sio2195[:,1] > .0001,0]**3 * sio2195[sio2195[:,1] > .0001,1], sio2195[sio2195[:,1] > .0001,0]) * 4./3.
## ## I2ac22  = np.pi * trapz(ac22[ac22[:,1] > .0001,0]**2 * ac22[ac22[:,1] > .0001,1], ac22[ac22[:,1] > .0001,0])
## ## I3ac22  = np.pi * trapz(ac22[ac22[:,1] > .0001,0]**3 * ac22[ac22[:,1] > .0001,1], ac22[ac22[:,1] > .0001,0]) * 4./3.
## ## print "AC, M195 (SFS04)"
## ## print "I2, I3: ", I2ac195, I3ac195
## ## print "My S=I2*V/M: ", I2ac195 * V195/M195
## ## print "Actual S =8.06e4 / My S ", 8.06e4 / (I2ac195 * V195/M195)
## ## print "Their grain density (g/cm^3): ", I2ac195/I3ac195 / 8.06e4
## ## print ""
## ## print "SiO2, M195 (SFS04)"
## ## print "I2, I3: ", I2si195, I3si195
## ## print "My S=I2*V/M: ", I2si195 * V195/M195
## ## print "Actual S =5.58e4 / My S ", 5.58e4 / (I2si195 * V195/M195)
## ## print "Their grain density (g/cm^3): ", I2si195/I3si195 / 5.58e4
## ## print ""
## ## print "AC, M22 (TF01)"
## ## print "I2, I3: ", I2ac22, I3ac22
## ## print "My S=I2*V/M: ", I2ac22 * V22/M22
## ## print "Actual S =8.50e4 / My S ", 8.50e4 / (I2ac22 * V22/M22)
## ## print "Their grain density (g/cm^3): ", I2ac22/I3ac22 / 8.50e4


## ## print "--------------------------------------"
## ## print "Cut off simps"
## ## print "--------------------------------------"
## ## I2ac195 = np.pi * simps(ac195[ac195[:,1] > .0001,0]**2 * ac195[ac195[:,1] > .0001,1], ac195[ac195[:,1] > .0001,0])
## ## I3ac195 = np.pi * simps(ac195[ac195[:,1] > .0001,0]**3 * ac195[ac195[:,1] > .0001,1], ac195[ac195[:,1] > .0001,0]) * 4./3.
## ## ## I2si195 = np.pi * simps(sio2195[sio2195[:,1] > .0001,0]**2 * sio2195[sio2195[:,1] > .0001,1], sio2195[sio2195[:,1] > .0001,0])
## ## ## I3si195 = np.pi * simps(sio2195[sio2195[:,1] > .0001,0]**3 * sio2195[sio2195[:,1] > .0001,1], sio2195[sio2195[:,1] > .0001,0]) * 4./3.
## ## ## I2ac22  = np.pi * simps(ac22[ac22[:,1] > .0001,0]**2 * ac22[ac22[:,1] > .0001,1], ac22[ac22[:,1] > .0001,0])
## ## ## I3ac22  = np.pi * simps(ac22[ac22[:,1] > .0001,0]**3 * ac22[ac22[:,1] > .0001,1], ac22[ac22[:,1] > .0001,0]) * 4./3.
## ## print "AC, M195 (SFS04)"
## ## print "I2, I3: ", I2ac195, I3ac195
## ## print "My S=I2*V/M: ", I2ac195 * V195/M195
## ## print "Actual S =8.06e4 / My S ", 8.06e4 / (I2ac195 * V195/M195)
## ## print "Their grain density (g/cm^3): ", I2ac195/I3ac195 / 8.06e4
## ## ## print ""
## ## ## print "SiO2, M195 (SFS04)"
## ## ## print "I2, I3: ", I2si195, I3si195
## ## ## print "My S=I2*V/M: ", I2si195 * V195/M195
## ## ## print "Actual S =5.58e4 / My S ", 5.58e4 / (I2si195 * V195/M195)
## ## ## print "Their grain density (g/cm^3): ", I2si195/I3si195 / 5.58e4
## ## ## print ""
## ## ## print "AC, M22 (TF01)"
## ## ## print "I2, I3: ", I2ac22, I3ac22
## ## ## print "My S=I2*V/M: ", I2ac22 * V22/M22
## ## ## print "Actual S =8.50e4 / My S ", 8.50e4 / (I2ac22 * V22/M22)
## ## ## print "Their grain density (g/cm^3): ", I2ac22/I3ac22 / 8.50e4



## ############ Ordering of the constants makes no difference, good...
## ## I2ac195 = trapz(np.pi * ac195[:,0]**2 * ac195[:,1], ac195[:,0])
## ## I3ac195 = trapz(np.pi * 4./3. * ac195[:,0]**3 * ac195[:,1], ac195[:,0])
## ## I2si195 = trapz(np.pi * sio2195[:,0]**2 * sio2195[:,1], sio2195[:,0])
## ## I3si195 = trapz(np.pi * 4./3. * sio2195[:,0]**3 * sio2195[:,1], sio2195[:,0])
## ## I2ac22  = trapz(np.pi * ac22[:,0]**2 * ac22[:,1], ac22[:,0])
## ## I3ac22  = trapz(np.pi * 4./3. * ac22[:,0]**3 * ac22[:,1], ac22[:,0])

## ## print ""
## ## print "--------------------------------------"
## ## print "AC, M195 (SFS04)"
## ## print "I2, I3: ", I2ac195, I3ac195
## ## print "My S=I2*V/M: ", I2ac195 * V195/M195
## ## print "Actual S =8.06e4 / My S ", 8.06e4 / (I2ac195 * V195/M195)
## ## print "Their grain density (g/cm^3): ", I2ac195/I3ac195 / 8.06e4
## ## print ""
## ## print "SiO2, M195 (SFS04)"
## ## print "I2, I3: ", I2si195, I3si195
## ## print "My S=I2*V/M: ", I2si195 * V195/M195
## ## print "Actual S =5.58e4 / My S ", 5.58e4 / (I2si195 * V195/M195)
## ## print "Their grain density (g/cm^3): ", I2si195/I3si195 / 5.58e4
## ## print ""
## ## print "AC, M22 (TF01)"
## ## print "I2, I3: ", I2ac22, I3ac22
## ## print "My S=I2*V/M: ", I2ac22 * V22/M22
## ## print "Actual S =8.50e4 / My S ", 8.50e4 / (I2ac22 * V22/M22)
## ## print "Their grain density (g/cm^3): ", I2ac22/I3ac22 / 8.50e4

## ## get_ipython().system(u'ls -F --color ')
## ## numpy
## ## loadtxt
## ## dnda = loadtxt(open("ac_195_dt.txt"),delimiter=',',skiprows=1)
## ## dnda
## ## plot(dnda[:,0],10**dnda[:,1])
## ## a = dnda[:,0]
## ## n = dnda[:,1]
## ## a
## ## n
## ## n = 10**n
## ## n
## ## a = a/10**4
## ## a
## ## from scipy.integrate import trapz
## ## from scipy.integrate import simps
## ## np.pi * trapz(a**2 *n)
## ## trapz(np.pi * a**2 *n)
## ## trapz(np.pi * a**2 *n, a)
## ## np.pi * trapz(a**2 *n, a)
## ## 4/3 * np.pi * trapz(a**3 *n, a)
## ## np.pi * trapz(a**2 *n, a) / 8.06 e 4
## ## np.pi * trapz(a**2 *n, a) / 8.06e4
## ## 3e55e-21
## ## 3.55e-21
## ## 3.55e-21 / (4/3 * np.pi * trapz(a**3 *n, a))
## ## I2 = np.pi * trapz(a**2 * n,a)
## ## I3 = 4*np.pi/3 * trapz(a**3 * n,a)
## ## I2
## ## I3
## ## 4/3 * np.pi * trapz(a**3 * n,a)
## ## ((4*np.pi)/3) * trapz(a**3*n,a)
## ## 4*np.pi/3
## ## 4/3*np.pi
## ## I2/I3
## ## I2/I3 * 1/8.06e4
## ## 6.0e33/4.2e48 /I3
## ## si = loadtxt(open("sio2_195_dt.txt"),delimiter=',',skiprows=1)
## ## si
## ## si[:,0] = si[:,0]/10**4
## ## si[:,1] = 10**si[:,1]
## ## I2SiO2 = np.pi * trapz(si[:,0]**2 * si[:,1],si[:,0])
## ## I3SiO2 = 4.0*np.pi/3.0 * trapz(si[:,0]**3 * si[:,1],si[:,0])
## ## I2SiO2
## ## I3SiO2
## ## 5.58e4/I2SiO2
## ## I2SiO4 / I3SiO4 / 5.58e4
## ## I2SiO2 / I3SiO2 / 5.58e4
## ## 8.15/2.82
## ## 13/3
## ## 13/3.
