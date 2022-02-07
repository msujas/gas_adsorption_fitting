import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt

file = 'Ar_occupancy_1bar.txt'



class Adsorption_profile:
    '''
    built in values: temperature (T), adsorption 1 (ads1), adsorption 2 (ads2),
    dH1 (Ha), dH2 (Hb), dS1 (Sa), dS2 (Sb), Jab, minimum adsorption (minads),
    maximum adsorption (maxads), array of thermodynamic values (values),
    another array where Hb and Sb are replaced proportions of Ha and Sa (value_proportion),
    fitted data 1 (fit1), fitted data 2 (fit2)

    example usage:
    adsorption_data = Adsorption_profile()
    adsorption_data.read_file('<text file name>')
    adsorption_data.run_optimise(bounds = <optional tuple of bounds>)
    adsorption_data.adsorption_plot()
    '''
    def __init__(self):
        self.T = np.array([])
        self.ads1 = np.array([])
        self.ads2 = np.array([])
        self.Ha = -5000
        self.Hb = -5000
        self.Sa = -50
        self.Sb = -50
        self.Jab = 0
        self.minads = 0
        self.maxads = 1
        self.values = np.array([self.Ha,self.Hb,self.Sa,self.Sb,self.Jab,self.minads,self.maxads])
        self.value_proportion = np.array([self.Ha,1,self.Sa,1,self.Jab,self.minads,self.maxads])

        self.fit1 = np.array([])
        self.fit2 = np.array([])
    def update_values(self):
        self.values = np.array([self.Ha,self.Hb,self.Sa,self.Sb,self.Jab,self.minads,self.maxads])
    def update_proportion_values(self):
        self.value_proportion = np.array([self.Ha,self.Hb/self.Ha,self.Sa,self.Sb/self.Sa,self.Jab,self.minads,self.maxads])
    def read_file(self,file):
        self.T,self.ads1,self.ads2 = np.loadtxt(file,unpack = True, skiprows = 1)
        self.maxads = max([max(self.ads1),max(self.ads2)])*1.1
        self.update_values()
        self.update_proportion_values()
    def twosite_nonequiv(self,Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi):
        R = 8.314
        dA = -1*(Hai-self.T*Sai)
        dB = -1*(Hbi-self.T*Sbi)
        S = dA + dB
        D = dA - dB

        siga = (2*np.exp(-Jabi/(R*self.T))*np.sinh(S/(R*self.T)) + \
        2*np.exp(Jabi/(R*self.T))*np.sinh(D/(R*self.T)))/ \
        (2*np.exp(-Jabi/(R*self.T))*np.cosh(S/(R*self.T))+\
        2*np.exp(Jabi/(R*self.T))*np.cosh(D/(R*self.T)))


        sigb = (2*np.exp(-Jabi/(R*self.T))*np.sinh(S/(R*self.T)) - \
        2*np.exp(Jabi/(R*self.T))*np.sinh(D/(R*self.T)))/ \
        (2*np.exp(-Jabi/(R*self.T))*np.cosh(S/(R*self.T))+\
        2*np.exp(Jabi/(R*self.T))*np.cosh(D/(R*self.T)))

        gamma_a = (siga+1)/2
        gamma_b = (sigb+1)/2
        gamma_tot = (gamma_a + gamma_b)/2
        adsfit1 = gamma_a*(maxadsi-minadsi)+minadsi
        adsfit2 = gamma_b*(maxadsi-minadsi)+minadsi
        adsav = (adsfit1+adsfit2)/2
        sigtot = (siga+sigb)/2
        return adsfit1, adsfit2

    def twosite_nonequiv_optimise(self, x=np.array([None])):
        if x.any() == None:
            x = self.values

        Hai = x[0]
        Hbi = x[1]
        Sai = x[2]
        Sbi = x[3]
        Jabi = x[4]
        minadsi = x[5]
        maxadsi = x[6]

        adsfit1,adsfit2 = self.twosite_nonequiv(Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi)


        yresid  = np.append(np.array([self.ads1-adsfit1]),np.array([self.ads2-adsfit2]))
        return yresid
    def run_optimise(self,initial=np.array([None]),bounds = ([-1*np.inf]*7,[np.inf]*7)):
        if initial.any() == None:
            initial = self.values
        yopt = optimize.least_squares(self.twosite_nonequiv_optimise,initial,
                                       bounds = bounds)
        self.Ha = yopt['x'][0]
        self.Hb = yopt['x'][1]
        self.Sa = yopt['x'][2]
        self.Sb = yopt['x'][3]
        self.Jab = yopt['x'][4]
        self.minads = yopt['x'][5]
        self.maxads = yopt['x'][6]
        #print('residuals =',np.sum(yopt['fun']))
        self.update_values()
        print('dHa =',       self.Ha)
        print('dHb =',       self.Hb)
        print('dSa =',       self.Sa)
        print('dSb =',       self.Sb)
        print('Jab =',       self.Jab)
        print('min. ads. =', self.minads)
        print('max. ads. =', self.maxads)
        self.fit1, self.fit2 = self.twosite_nonequiv(*self.values)
        print('residual =',np.sum(yopt['fun']))
        return yopt

    def twosite_nonequiv_intra_optimise_proportion_bounds(self,x=np.array([None])):
        '''
        x is list of arguements: Ha, Hb/Ha, Sa, Sb/Sa, Jab, minads, maxads2
        Same as twosite_nonequiv_intra_optimise, but 2nd and 4th arguements are
        Hb/Ha and Sb/Sa respectively as it makes it more convenient to bound Hb
        and Sb to propotions of Ha and Sa.
        '''
        if x.any() == None:
            x = self.value_proportion
        Hai = x[0]
        Hbi = Ha*x[1]
        Sai = x[2]
        Sbi = Sa*x[3]
        Jabi = x[4]
        minadsi = x[5]
        maxadsi = x[6]

        adsfit1,adsfit2 = self.twosite_nonequiv(Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi)


        yresid  = np.append(np.array([self.ads1-adsfit1]),np.array([self.ads2-adsfit2]))
        return yresid
    def run_optimise_proportion(self,initial=np.array([None]),bounds=([-1*np.inf]*7,[np.inf]*7)):
        if initial.any() == None:
            initial = self.value_proportion
        yopt = optimize.least_squares(self.twosite_nonequiv_intra_optimise_proportion_bounds,
                                        initial, bounds = bounds)
        print(yopt)
        self.Ha = yopt['x'][0]
        self.Hb = yopt['x'][1]
        self.Sa = yopt['x'][2]
        self.Sb = yopt['x'][3]*self.Sa
        self.Jab = yopt['x'][4]
        self.minads = yopt['x'][5]
        self.maxads = yopt['x'][6]
        self.update_values()
        self.update_proportion_values()
        print('dHa =',       self.Ha)
        print('dHb =',       self.Hb)
        print('dSa =',       self.Sa)
        print('dSb =',       self.Sb)
        print('Jab =',       self.Jab)
        print('min. ads. =', self.minads)
        print('max. ads. =', self.maxads)
        self.fit1, self.fit2 = self.twosite_nonequiv(*self.values)
        print('residual =',np.sum(yopt['fun']))
        return yopt
    def adsorption_plot(self):
        adsav = (self.ads1 + self.ads2)/2
        fitav = (self.fit1 + self.fit2)/2
        plt.plot(self.T,self.ads1,'o',label = 'adsorption 1')
        plt.plot(self.T,self.ads2,'o', label = 'adsorption 2')
        plt.plot(self.T,adsav,'o', label = 'adsorption av.')
        plt.plot(self.T,self.fit1, label = 'fit 1')
        plt.plot(self.T,self.fit2, label = 'fit 2')
        plt.plot(self.T,fitav, label = 'fit av.')
        plt.xlabel('Temperature (K)')
        plt.ylabel('adsorption')
        plt.xlim(min(self.T),max(self.T))
        plt.show()



bounds = (np.array([-50000,-50000,-1000,-1000,-100,-1,14]),
          np.array([0,0,1000,1000,20000,0,20]))

Ar_ads_data = Adsorption_profile()

Ar_ads_data.read_file('Ar_occupancy_1bar.txt')
print(Ar_ads_data.values)
yopt = Ar_ads_data.run_optimise(bounds = bounds)

print(Ar_ads_data.values)

Ar_ads_data.adsorption_plot()
