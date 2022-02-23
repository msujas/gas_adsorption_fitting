import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize


R = 8.31446
temperatures = np.arange(10,301)
T = np.arange(10,301)
Ha = -8.15532658e+03
Hb = -4.54085862e+03
Sa = -7.01540102e+01
Sb = -3.59307747e+01
Jab = 2000
J1 =  600

def two_site_inter_iterative(x,Ha,Hb,Sa,Sb,Jab,J1):
    Ga = Ha-T*Sa
    Gb = Hb-T*Sb
    lenx_half = int(len(x)/2)
    sigmaa_guess = x[:lenx_half]
    sigmab_guess = x[lenx_half:]
    sigmaa = np.tanh((-Ga-2*Jab*sigmab_guess-2*J1*sigmaa_guess)/(R*T))
    sigmab = np.tanh((-Gb-2*Jab*sigmaa_guess-2*J1*sigmab_guess)/(R*T))

    return np.append(sigmaa-sigmaa_guess,sigmab-sigmab_guess)

def two_site_inter_iterative_fit(sigmaa_guess,sigmab_guess):
    Ga = Ha-T*Sa
    Gb = Hb-T*Sb

    sigmaa = np.tanh((-Ga-2*Jab*sigmab_guess-2*J1*sigmaa_guess)/(R*T))
    sigmab = np.tanh((-Gb-2*Jab*sigmaa_guess-2*J1*sigmab_guess)/(R*T))

    return sigmaa,sigmab

def twosite_nonequiv(T,Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi,minadsi2,maxadsi2):
    R = 8.31446
    dA = -1*(Hai-T*Sai)
    dB = -1*(Hbi-T*Sbi)
    S = dA + dB
    D = dA - dB

    siga = (2*np.exp(-Jabi/(R*T))*np.sinh(S/(R*T)) + \
    2*np.exp(Jabi/(R*T))*np.sinh(D/(R*T)))/ \
    (2*np.exp(-Jabi/(R*T))*np.cosh(S/(R*T))+\
    2*np.exp(Jabi/(R*T))*np.cosh(D/(R*T)))


    sigb = (2*np.exp(-Jabi/(R*T))*np.sinh(S/(R*T)) - \
    2*np.exp(Jabi/(R*T))*np.sinh(D/(R*T)))/ \
    (2*np.exp(-Jabi/(R*T))*np.cosh(S/(R*T))+\
    2*np.exp(Jabi/(R*T))*np.cosh(D/(R*T)))

    gamma_a = (siga+1)/2
    gamma_b = (sigb+1)/2
    gamma_tot = (gamma_a + gamma_b)/2
    adsfit1 = gamma_a*(maxadsi-minadsi)+minadsi
    adsfit2 = gamma_b*(maxadsi2-minadsi2)+minadsi2
    #adsav = (adsfit1+adsfit2)/2
    #sigtot = (siga+sigb)/2
    return siga, sigb
sigmaa_guess,sigmab_guess = twosite_nonequiv(T,Ha,Hb,Sa,Sb,Jab,0,1,0,1)

x = np.append(sigmaa_guess,sigmab_guess)
yopt = optimize.least_squares(two_site_inter_iterative,x0=x,args = (Ha,Hb,Sa,Sb,Jab,J1))
print(yopt['fun'])
lenopt_half = int(len(yopt['x'])/2)
sigmaa_fit = yopt['x'][:lenopt_half]
sigmab_fit = yopt['x'][lenopt_half:]
sigmaa,sigmab = two_site_inter_iterative_fit(sigmaa_fit,sigmab_fit)



plt.plot(temperatures,sigmaa,'o')
plt.plot(temperatures,sigmab,'o')
plt.plot(temperatures,sigmaa_fit)
plt.plot(temperatures,sigmab_fit)
plt.show()
