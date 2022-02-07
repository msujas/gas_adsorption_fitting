import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt

file = 'Ar_occupancy_1bar.txt'



class Adsorption_profile:
    '''
    Attributes: temperature (T), adsorption 1 (ads1), adsorption 2 (ads2),
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
        self.value_proportion = np.array([self.Ha,self.Hb/self.Ha,self.Sa,self.Sb/self.Sa,self.Jab,self.minads,self.maxads])
    def read_file(self,file,delimiter=None,skiprows = 0,usecols = None):
        '''
        uses numpy.loadtxt. Positional arguments: file. Keyword arguments:
        delimiter=None, skiprows=0, usecols=None.
        '''
        self.T,self.ads1,self.ads2 = np.loadtxt(file,unpack = True, skiprows = skiprows,
        delimiter = delimiter,usecols=usecols)
        self.maxads = max([max(self.ads1),max(self.ads2)])*1.1
        self.fit1 = np.zeros(len(self.T))
        self.fit2 = np.zeros(len(self.T))
        self.update_values()

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
        #adsav = (adsfit1+adsfit2)/2
        #sigtot = (siga+sigb)/2
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

        yresid  = np.append(self.ads1-adsfit1,self.ads2-adsfit2)
        return yresid

    def twosite_nonequiv_minmax2(self,Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi,minadsi2,maxadsi2):
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
        #gamma_tot = (gamma_a + gamma_b)/2
        adsfit1 = gamma_a*(maxadsi-minadsi)+minadsi
        adsfit2 = gamma_b*(maxadsi2-minadsi2)+minadsi2
        #adsav = (adsfit1+adsfit2)/2
        #sigtot = (siga+sigb)/2
        return adsfit1, adsfit2

    def twosite_nonequiv_minmax2_optimise(self, x=np.array([None])):
        if x.any() == None:
            x = self.values

        Hai = x[0]
        Hbi = x[1]
        Sai = x[2]
        Sbi = x[3]
        Jabi = x[4]
        minadsi = x[5]
        maxadsi = x[6]
        minadsi2 = x[7]
        maxadsi2 = x[8]

        adsfit1,adsfit2 = self.twosite_nonequiv_minmax2(Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi,minadsi2,maxadsi2)

        yresid  = np.append(self.ads1-adsfit1,self.ads2-adsfit2)
        return yresid

    def run_optimise(self,initial=np.array([None]),bounds = ([-1*np.inf]*7,[np.inf]*7)):
        if initial.any() == None:
            initial = self.values
        if len(initial) == 7:
            yopt = optimize.least_squares(self.twosite_nonequiv_optimise,initial,
                                       bounds = bounds)
        elif len(initial) == 9:
            yopt = optimize.least_squares(self.twosite_nonequiv_minmax2_optimise,initial,
                       bounds = bounds)
        self.Ha = yopt['x'][0]
        self.Hb = yopt['x'][1]
        self.Sa = yopt['x'][2]
        self.Sb = yopt['x'][3]
        self.Jab = yopt['x'][4]
        self.minads = yopt['x'][5]
        self.maxads = yopt['x'][6]
        if len(initial) == 9:
            self.minads2 = yopt['x'][7]
            self.maxads2 = yopt['x'][8]

        self.update_values()
        rw = np.sum(yopt['fun'])/np.sum(np.append(self.ads1,self.ads2))**2
        string = f'''dHa = {self.Ha}
dHb = {self.Hb}
dSa = {self.Sa}
dSb = {self.Sb}
Jab = {self.Jab}
min. ads. = {self.minads}
max. ads. = {self.maxads}
rw = {rw}
'''
        if len(initial) == 9:
            string += f'min. ads. 2 = {self.minads2}\nmax. ads. 2 = {self.maxads2}'

        print(string)

        f = open('fitted_parameters.txt','w')
        f.write(string)
        f.close()
        if len(initial == 7):
            self.fit1, self.fit2 = self.twosite_nonequiv(*self.values)
        if len(initial == 9):
            self.values = np.append(self.values,self.minads2)
            self.values = np.append(self.values,self.maxads2)
            self.fit1, self.fit2 = self.twosite_nonequiv_minmax2(*self.values)

        np.savetxt('adsorption_fit.txt',np.array([self.T,self.ads1,self.fit1,self.ads2,self.fit2]).transpose(),
        header = 'T(K) adsorption1 fit1 adsorption2 fit2')


        return yopt

    def twosite_nonequiv_intra_optimise_proportion_bounds(self,x=np.array([None])):
        '''
        x is list of arguments: Ha, Hb/Ha, Sa, Sb/Sa, Jab, minads, maxads2
        Same as twosite_nonequiv_intra_optimise, but 2nd and 4th arguments are
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


        yresid  = np.append(self.ads1-adsfit1,self.ads2-adsfit2)
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

        rw = np.sum(yopt['fun'])/np.sum(np.append(self.ads1,self.ads2))**2
        string = f'''dHa = {self.Ha}
dHb = {self.Hb}
dSa = {self.Sa}
dSb = {self.Sb}
Jab = {self.Jab}
min. ads. = {self.minads}
max. ads. = {self.maxads}
rw = {rw}
'''
        print(string)
        f = open('fitted_parameters.txt','w')
        f.write(string)
        f.close()
        self.fit1, self.fit2 = self.twosite_nonequiv(*self.values)
        np.savetxt('adsorption_fit.txt',np.array([self.T,self.ads1,self.fit1,
        self.ads2,self.fit2]).transpose(),
        header = 'T(K) adsorption1 fit1 adsorption2 fit2')

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
        plt.legend()
        plt.show()


initial2 = np.array([-5000,-5000,-50,-50,0,0,15,0,15])
bounds = (np.array([-50000,-50000,-1000,-1000,-100,-1,14]),
          np.array([0,0,1000,1000,20000,0,20]))
bounds2 = (np.array([-50000,-50000,-1000,-1000,-100,-1,14,-1,14]),
          np.array([0,0,1000,1000,20000,0,30,0,17]))
Ar_ads_data = Adsorption_profile()

Ar_ads_data.read_file('Ar_occupancy_1bar.txt',skiprows = 1)
#print(Ar_ads_data.values)
yopt = Ar_ads_data.run_optimise(initial = initial2,bounds = bounds2)

#print(Ar_ads_data.values)

Ar_ads_data.adsorption_plot()
