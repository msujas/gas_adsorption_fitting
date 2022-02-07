import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt

file = 'Ar_occupancy_1bar.txt'
T,ads1,ads2 = np.loadtxt(file,unpack = True, skiprows=1)
print(min(ads1),min(ads2))
min_value = min([min(ads1),min(ads2)])
def twosite_nonequiv_intra(Ha,Hb,Sa,Sb,Jab,minads,maxads,T=T):
    R = 8.314
    dA = -1*(Ha-T*Sa)
    dB = -1*(Hb-T*Sb)
    S = dA + dB
    D = dA - dB

    siga = (2*np.exp(-Jab/(R*T))*np.sinh(S/(R*T)) + 2*np.exp(Jab/(R*T))*np.sinh(D/(R*T)))/ \
    (2*np.exp(-Jab/(R*T))*np.cosh(S/(R*T))+2*np.exp(Jab/(R*T))*np.cosh(D/(R*T)))
    sigb = (2*np.exp(-Jab/(R*T))*np.sinh(S/(R*T)) - 2*np.exp(Jab/(R*T))*np.sinh(D/(R*T)))/ \
    (2*np.exp(-Jab/(R*T))*np.cosh(S/(R*T))+2*np.exp(Jab/(R*T))*np.cosh(D/(R*T)))

    gamma_a = (siga+1)/2
    gamma_b = (sigb+1)/2
    gamma_tot = (gamma_a + gamma_b)/2
    adsfit1 = gamma_a*(maxads-minads)+minads
    adsfit2 = gamma_b*(maxads-minads)+minads
    adsav = (adsfit1+adsfit2)/2
    sigtot = (siga+sigb)/2
    return adsfit1, adsfit2

def twosite_nonequiv_intra_optimise(x):
    global T
    global ads1
    global ads2
    Ha = x[0]
    Hb = x[1]
    Sa = x[2]
    Sb = x[3]
    Jab = x[4]
    minads = x[5]
    maxads = x[6]

    adsfit1,adsfit2 = twosite_nonequiv_intra(Ha,Hb,Sa,Sb,Jab,minads,maxads)


    yresid  = np.append(np.array([ads1-adsfit1]),np.array([ads2-adsfit2]))
    return yresid

def twosite_nonequiv_intra_optimise_entropy_bounds(x):
    global T
    global ads1
    global ads2
    Ha = x[0]
    Hb = x[1]
    Sa = x[2]
    Sb = Sa*x[3]
    Jab = x[4]
    minads = x[5]
    maxads = x[6]

    adsfit1,adsfit2 = twosite_nonequiv_intra(Ha,Hb,Sa,Sb,Jab,minads,maxads)


    yresid  = np.append(np.array([ads1-adsfit1]),np.array([ads2-adsfit2]))
    return yresid

initial = np.array([-5000,-5000,-50,-50,1500,0,16.7])
bounds = (np.array([-50000,-50000,-1000,-1000,-100,-1,14]),
          np.array([0,0,1000,1000,20000,0,20]))
yresid = optimize.least_squares(twosite_nonequiv_intra_optimise,initial,
                               bounds= bounds)
print(yresid)
Ha = yresid['x'][0]
Hb = yresid['x'][1]
Sa = yresid['x'][2]
Sb = yresid['x'][3]
Jab = yresid['x'][4]
minads = yresid['x'][5]
maxads = yresid['x'][6]

print('residuals =',np.sum(yresid['fun']))

print('dHa =',Ha)
print('dHb =',Hb)
print('dSa =', Sa)
print('dSb =', Sb)
print('Jab =', Jab)
print('min. ads. =', minads)
print('max. ads. =', maxads)

print(yresid['x'])
adsfit1,adsfit2 = twosite_nonequiv_intra(Ha,Hb,Sa,Sb,Jab,minads,maxads)

adsav = (ads1+ads2)/2
adsfitav = (adsfit1+adsfit2)/2



'''
entropy bound fit
'''

initial2 = np.array([-5000,-5000,-50, 1,5000,0,16.7])
bounds2 = (np.array([-50000,-50000,-1000,0.67,-100,-1,14]),
          np.array([0,0,1000,1.5,20000,0,20]))
yresid2 = optimize.least_squares(twosite_nonequiv_intra_optimise_entropy_bounds,
                                initial2, bounds = bounds2)

print(yresid2)
Ha2 = yresid2['x'][0]
Hb2 = yresid2['x'][1]
Sa2 = yresid2['x'][2]
Sb2 = yresid2['x'][3]*Sa
Jab2 = yresid2['x'][4]
minads2 = yresid2['x'][5]
maxads2 = yresid2['x'][6]

print('dHa =',Ha2)
print('dHb =',Hb2)
print('dSa =', Sa2)
print('dSb =', Sb2)
print('Jab =', Jab2)
print('min. ads. =', minads2)
print('max. ads. =', maxads2)

adsfit1_2,adsfit2_2 = twosite_nonequiv_intra(Ha2,Hb2,Sa2,Sb2,Jab2,minads2,maxads2)
adsfitav2 = (adsfit1_2+adsfit2_2)/2
'''
plotting
'''

fontsize = 16

fig,(ax1,ax2) = plt.subplots(2,1)

ax1.set_title('no entropy or enthalpy bounds')
ax1.plot(T,ads1,'o',label='Ar1')
ax1.plot(T,ads2,'o',label = 'Ar2')
ax1.plot(T,adsav,'o',label = 'Ar av.')
ax1.plot(T,adsfit1, label = 'Ar1 fit',linewidth = 3)
ax1.plot(T,adsfit2,label = 'Ar2 fit', linewidth = 3)
ax1.plot(T,adsfitav, label = 'fit av.', linewidth = 3)
ax1.set_ylabel('adsorption',fontsize = fontsize)
ax1.set_xlabel('temperature (K)',fontsize = fontsize)


ax1.legend(fontsize = fontsize)

ax2.set_title('Entropy bound 0.67x - 1.5x')
ax2.plot(T,ads1,'o',label='Ar1')
ax2.plot(T,ads2,'o',label = 'Ar2')
ax2.plot(T,adsav,'o',label = 'Ar av.')
ax2.plot(T,adsfit1_2, label = 'Ar1 fit',linewidth = 3)
ax2.plot(T,adsfit2_2,label = 'Ar2 fit', linewidth = 3)
ax2.plot(T,adsfitav2, label = 'fit av.', linewidth = 3)
ax2.set_ylabel('adsorption',fontsize = fontsize)
ax2.set_xlabel('temperature (K)',fontsize = fontsize)

ax2.legend(fontsize = fontsize)
plt.show()
