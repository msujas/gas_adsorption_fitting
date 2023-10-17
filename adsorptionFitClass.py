import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import os
import matplotlib
matplotlib.rcParams.update({'font.size': 14})


class One_site_adsorption_profile:
    def __init__(self):
        self.T = np.array([])
        self.ads = np.array([])
        self.H = -5000
        self.S = -50
        self.J = 0
        self.p = 1
        self.minads = 0
        self.maxads = 1
        self.values = np.array([self.H,self.S,self.J,self.minads,self.maxads])
        self.fit = np.array([])
        self.minmax_refine = False
    def update_values(self):
        self.values = np.array([self.H,self.S,self.J,self.minads,self.maxads])
        self.gamma = (self.ads - self.minads)/(self.maxads-self.minads)
    def read_file(self,file,delimiter=None,usecols = (0,1),skiprows = None,unit = 'K'):
        '''
        uses numpy.loadtxt. Positional arguments: file. Keyword arguments:
        delimiter=None, skiprows=0, usecols=None.
        '''

        self.filename = file
        self.directory = os.path.split(file)[0]+'/'

        if self.directory == '/':
            self.directory = os.path.curdir +'/'
        R = 8.31446
        if file.endswith('.csv'):
            delimiter = ','
        if skiprows == None:
            for i in range(20):
                try:
                    self.T,self.ads = np.loadtxt(file,unpack = True, skiprows = i,
                    delimiter = delimiter,usecols=usecols)
                    break
                except:
                    if i == 19:
                        raise RuntimeError
                    else:
                        continue

        else:

            self.T,self.ads = np.loadtxt(file,unpack = True, skiprows = skiprows,
            delimiter = delimiter,usecols=usecols)

        self.maxads = max(self.ads)*1.1
        if unit == 'C':
            self.T += 273.15

        self.update_values()
        self.fit = self.one_site_simple(self.H,self.S,self.minads,self.maxads)

        self.Keq = self.gamma/(1-self.gamma)
        self.Keq_fit = np.exp(-self.H/(R*self.T) + self.S/R)

    def one_site_simple(self,Hi,Si,minadsi,maxadsi):
        deltaG = Hi - self.T*Si
        R = 8.31446
        gamma = np.exp(-deltaG/(R*self.T))*self.p/(1+np.exp(-deltaG/(R*self.T))*self.p)
        adsorption = gamma*(maxadsi-minadsi) + minadsi
        return adsorption
    def one_site_simple_optimise(self,x=np.array([None])):
        if x.any() == None:
            x = [self.H,self.S,self.minads,self.maxads]
        if len(x) == 4:
            adsfit = self.one_site_simple(*x)
        elif len(x) == 2:
            adsfit = self.one_site_simple(x[0],x[1],self.minads,self.maxads)
        return self.ads - adsfit
    def run_optimise(self,x = np.array([None]),bounds = ([np.inf]*4,[np.inf]*4)):
        if x.any() == None and len(x) == 4:
            x = np.array([self.H,self.S,self.minads,self.maxads])
        elif x.any() == None and len(x) == 2:
            x = np.array([self.H,self.S])
            bounds = (bounds[0][:2],bounds[1][:2])

        yopt = least_squares(self.one_site_simple_optimise,x,bounds=bounds)
        self.H = yopt['x'][0]
        self.S = yopt['x'][1]

        if len(yopt['x']) == 4:
            self.minads = yopt['x'][2]
            self.maxads = yopt['x'][3]
        R = 8.31446
        self.fit = self.one_site_simple(self.H,self.S,self.minads,self.maxads)
        rw = np.sum(self.ads-self.fit)**2/np.sum(self.ads**2)
        string = f'''dH = {self.H}
dS = {self.S}
min. ads. = {self.minads}
max. ads = {self.maxads}
rw = {rw}
'''
        print(string)
        output_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_fit_pars.txt'
        f = open(self.directory+output_filename,'w')
        f.write(string)
        f.close()
        print('fitted parameters saved to',output_filename)

        self.update_values()
        self.Keq = self.gamma/(1-self.gamma)
        self.Keq_fit = np.exp(-self.H/(R*self.T) + self.S/R)
        fit_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_ads_fit.txt'
        np.savetxt(fit_filename,np.array([self.T,self.ads,self.fit]).transpose(),
        header = 'T(K) adsorption fit')
        print('fitted data saved to',fit_filename)

        return yopt

    def one_site_coop(self,Hi,Si,Ji,minadsi,maxadsi):

        gamma = (self.ads-minadsi)/(maxadsi-minadsi)
        R = 8.31446

        #gamma = np.arange(0.001,0.999,0.001)
        Tfit = (Hi-2*Ji*(2*gamma-1))/(Si-R*np.log(gamma/(self.p*(1-gamma))))
        return Tfit

    def one_site_coop_optimise(self,x):
        Hi = x[0]
        Si = x[1]
        Ji = x[2]
        minadsi = self.minads
        maxadsi = self.maxads
        if len(x) == 5:
            minadsi = x[3]
            maxadsi = x[4]

        Tfit = self.one_site_coop(Hi,Si,Ji,minadsi,maxadsi)
        return self.T - Tfit
    def run_coop_optimise(self,x,bounds):
        yopt = least_squares(self.one_site_coop_optimise,x,bounds = bounds)
        params = yopt['x']
        self.H = params[0]
        self.S = params[1]
        self.J = params[2]
        if len(x) == 5:
            self.minads = params[3]
            self.maxads = params[4]
        self.update_values()
        self.fit = self.one_site_coop(self.H,self.S,self.J,self.minads,self.maxads)
        rw = np.sum(self.T-self.fit)**2/np.sum(self.T**2)
        output = f'''deltaH = {self.H}
deltaS = {self.S}
J = {self.J}
min. ads. = {self.minads}
max. ads. = {self.maxads}
rw = {rw}'''
        print(output)
        output_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_1sitecoop_fit_pars.txt'
        f = open(self.directory+output_filename,'w')
        f.write(output)
        f.close()
        print('fitted parameters saved to',output_filename)
        fit_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_1sitecoop_ads_fit.txt'
        np.savetxt(fit_filename,np.array([self.T,self.ads,self.fit]).transpose(),
        header = 'T(K) adsorption fit')
        print('fitted data saved to',fit_filename)
        return yopt

    def twosite_nonequiv(self,Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi):
        R = 8.31446
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
        #sigtot = (siga+sigb)/2
        return adsav
    def twosite_optimise(self,x):
        Hai = x[0]
        Hbi = x[1]
        Sai = x[2]
        Sbi = x[3]
        Jabi = x[4]
        minadsi = self.minads
        maxadsi = self.maxads
        if len(x) == 7:
            minadsi = x[5]
            maxadsi = x[6]
        adsavfit = self.twosite_nonequiv(Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi)
        return self.ads - adsavfit

    def run_twosite_optimise(self,x,bounds):
        yopt = least_squares(self.twosite_optimise,x,bounds = bounds)
        params = yopt['x']
        self.H = params[0]
        self.Hb = params[1]
        self.S = params[2]
        self.Sb = params[3]
        self.J = params[4]
        if len(x) == 7:
            self.minads = params[5]
            self.maxads = params[6]
        self.update_values()
        self.fit = self.twosite_nonequiv(self.H,self.Hb,self.S,self.Sb,self.J,self.minads,self.maxads)
        rw = np.sum(self.ads - self.fit)**2/np.sum(self.ads**2)
        string = (f"deltaHa = {self.H}\n"
        f"deltaHb = {self.Hb}\n"
        f"deltaSa = {self.S}\n"
        f"deltaSb {self.Sb}\n"
        f"Jab = {self.J}\n"
        f"min. ads. = {self.minads}\n"
        f"max. ads. = {self.maxads}\n"
        f"rw = {rw}")
        print(string)
        output_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_fit_pars.txt'
        f = open(self.directory+output_filename,'w')
        f.write(string)
        f.close()
        print('fitted parameters saved to',output_filename)
        fit_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_ads_fit.txt'
        np.savetxt(fit_filename,np.array([self.T,self.ads,self.fit]).transpose(),
        header = 'T(K) adsorption fit')
        print('fitted data saved to',fit_filename)
        return yopt
    def adsorption_plot(self):
        plt.plot(self.T,self.ads,'o',label = 'adsorption')
        plt.plot(self.T,self.fit, label = 'fit')
        plt.xlim(min(self.T),max(self.T))
        plt.xlabel('Temperature (K)')
        plt.ylabel('Adsorption')
        plt.legend()
        plt.show()
    def vant_hoff_plot(self):
        plt.plot(1/self.T,np.log(self.Keq),'o',label = 'measured (min. max. normalised)')
        plt.plot(1/self.T,np.log(self.Keq_fit),label = 'fit')
        plt.xlim(min(1/self.T),max(1/self.T))
        plt.xlabel('1/T (K$^{-1}$)')
        plt.ylabel('ln(K)')
        plt.legend()
        plt.show()
    def ads_vant_hoff_plot(self):
        fig, (ax1,ax2) = plt.subplots(2,1)
        ax1.plot(self.T,self.ads,'o',label = 'adsorption')
        ax1.plot(self.T,self.fit, label = 'fit')
        ax1.set_xlim(min(self.T),max(self.T))
        ax1.set_xlabel('Temperature (K)')
        ax1.set_ylabel('Adsorption')
        ax1.legend()
        ax2.plot(1/self.T,np.log(self.Keq),'o',label = 'measured (min. max. normalised)')
        ax2.plot(1/self.T,np.log(self.Keq_fit),label = 'fit')
        ax2.set_xlim(min(1/self.T),max(1/self.T))
        ax2.set_xlabel('1/T (K$^{-1}$)')
        ax2.set_ylabel('ln(K)')
        ax2.legend()
        plt.show()

    def one_site_coop_plot(self):
        plt.plot(self.T,self.ads,'o',label = 'adsorption')
        plt.plot(self.fit,self.ads,label = 'fit')
        plt.xlabel('Temperature (K)')
        plt.ylabel('adsorption')
        plt.xlim(min(self.T),max(self.T))
        plt.legend()
        plt.show()








class Two_site_adsorption_profile:
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
        self.minads2 = 0
        self.maxads2 = 1
        self.values = np.array([self.Ha,self.Hb,self.Sa,self.Sb,self.Jab,self.minads,self.maxads,self.minads2,self.maxads2])
        self.value_proportion = np.array([self.Ha,1,self.Sa,1,self.Jab,self.minads,self.maxads])
        self.fit1 = np.array([])
        self.fit2 = np.array([])

    def update_values(self):
        self.values = np.array([self.Ha,self.Hb,self.Sa,self.Sb,self.Jab,self.minads,self.maxads,self.minads2,self.maxads2])
        self.value_proportion = np.array([self.Ha,self.Hb/self.Ha,self.Sa,self.Sb/self.Sa,self.Jab,self.minads,self.maxads])
    def read_file(self,file,delimiter=None,skiprows = 0,usecols = (0,1,2), unit = 'K'):
        '''
        uses numpy.loadtxt. Positional arguments: file. Keyword arguments:
        delimiter=None, skiprows=0, usecols=None.
        '''
        self.filename = file
        self.directory = os.path.split(file)[0]+'/'
        if self.directory == '/':
            self.directory = os.path.curdir +'/'
        if file.endswith('.csv'):
            delimiter = ','
        for i in range(20):
            try:
                self.T,self.ads1,self.ads2 = np.loadtxt(file,unpack = True, skiprows = i,
                delimiter = delimiter,usecols=usecols)
                break
            except:
                continue
        if unit == 'C':
            self.T += 273.15
        self.maxads = max([max(self.ads1),max(self.ads2)])*1.1
        self.update_values()
        self.fit1,self.fit2 = self.twosite_nonequiv(*self.values)
    def twosite_nonequiv(self,Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi,minadsi2,maxadsi2):
        R = 8.31446
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
        adsfit2 = gamma_b*(maxadsi2-minadsi2)+minadsi2
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
        minadsi = self.minads
        maxadsi = self.maxads
        minadsi2 = self.minads
        maxadsi2 = self.maxads
        if len(x) >= 7:
            minadsi = x[5]
            maxadsi = x[6]
            minadsi2 = x[5]
            maxadsi2 = x[6]
        if len(x) == 9:
            minadsi2 = x[7]
            maxadsi2 = x[8]


        adsfit1,adsfit2 = self.twosite_nonequiv(Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi,minadsi2,maxadsi2)

        yresid  = np.append(self.ads1-adsfit1,self.ads2-adsfit2)
        return yresid



    def run_optimise(self,x=np.array([None]),bounds = ([None],[None])):
        if x.any() == None:
            x = self.values[:7]
        #if any(bounds[0]) == False:
#
        #    bounds = ([-1*np.inf]*len(x),[np.inf]*len(x))

        yopt = least_squares(self.twosite_nonequiv_optimise,x,
                                   bounds = bounds)
        self.Ha = yopt['x'][0]
        self.Hb = yopt['x'][1]
        self.Sa = yopt['x'][2]
        self.Sb = yopt['x'][3]
        self.Jab = yopt['x'][4]
        if len(x) >= 7:
            self.minads = yopt['x'][5]
            self.maxads = yopt['x'][6]
            self.minads2 = yopt['x'][5]
            self.maxads2 = yopt['x'][6]
        if len(x) == 9:
            yopt = least_squares(self.twosite_nonequiv_minmax2_optimise,initial,
                       bounds = bounds)
            self.minads2 = yopt['x'][7]
            self.maxads2 = yopt['x'][8]




        self.update_values()

        rw = np.sum(yopt['fun'])**2/np.sum(np.append(self.ads1,self.ads2)**2)
        string = f'''dHa = {self.Ha}
dHb = {self.Hb}
dSa = {self.Sa}
dSb = {self.Sb}
Jab = {self.Jab}
min. ads. = {self.minads}
max. ads. = {self.maxads}
rw = {rw}
'''





        self.fit1, self.fit2 = self.twosite_nonequiv(*self.values)
        if len(x) == 9:

            string += f'min. ads. 2 = {self.minads2}\nmax. ads. 2 = {self.maxads2}'
        output_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_fit_pars.txt'
        print(string)
        f = open(self.directory+output_filename,'w')
        f.write(string)
        f.close()
        print('fitted parameters saved to',output_filename)

        fit_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_intrapore_ads_fit.txt'
        np.savetxt(fit_filename,np.array([self.T,self.ads1,self.fit1,self.ads2,self.fit2]).transpose(),
        header = 'T(K) adsorption1 fit1 adsorption2 fit2')
        print('fitted data saved to',fit_filename)

        return yopt



    def two_site_inter_sigma_solver(self,x,Hai,Hbi,Sai,Sbi,Jabi,J1i,minadsi,maxadsi):
        R = 8.31446
        Ga = Hai-self.T*Sai
        Gb = Hbi-self.T*Sbi
        lenx_half = int(len(x)/2)
        adsa_guess = x[:lenx_half]
        adsb_guess = x[lenx_half:]
        sigmaa_guess = 2*((adsa_guess-minadsi)/(maxadsi-minadsi))-1
        sigmab_guess = 2*((adsb_guess-minadsi)/(maxadsi-minadsi))-1
        sigmaa = np.tanh((-Ga-2*Jabi*sigmab_guess-2*J1i*sigmaa_guess)/(R*self.T))
        sigmab = np.tanh((-Gb-2*Jabi*sigmaa_guess-2*J1i*sigmab_guess)/(R*self.T))

        return np.append(sigmaa-sigmaa_guess,sigmab-sigmab_guess)

    def two_site_inter_fit(self,adsa_guess,adsb_guess,Hai,Hbi,Sai,Sbi,Jabi,J1i,minadsi,maxadsi):
        R = 8.31446
        Ga = Hai-self.T*Sai
        Gb = Hbi-self.T*Sbi
        sigmaa_guess = 2*((adsa_guess-minadsi)/(maxadsi-minadsi))-1
        sigmab_guess = 2*((adsb_guess-minadsi)/(maxadsi-minadsi))-1
        sigmaa = np.tanh((-Ga-2*Jabi*sigmab_guess-2*J1i*sigmaa_guess)/(R*self.T))
        sigmab = np.tanh((-Gb-2*Jabi*sigmaa_guess-2*J1i*sigmab_guess)/(R*self.T))
        gammaa = (sigmaa+1)/2
        gammab = (sigmab+1)/2
        adsa = gammaa*(maxadsi-minadsi) + minadsi
        adsb = gammab*(maxadsi-minadsi) + minadsi
        return adsa,adsb

    def two_site_inter_optimise(self,x):
        Hai = x[0]
        Hbi = x[1]
        Sai = x[2]
        Sbi = x[3]
        Jabi = x[4]
        J1i = x[5]
        minadsi = self.minads
        maxadsi = self.maxads
        if len(x) == 8:
            minadsi = x[6]
            maxadsi = x[7]
        adsa_guess,adsb_guess = self.twosite_nonequiv(Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi,minadsi,maxadsi)
        initial = np.append(adsa_guess,adsb_guess)
        yopt_sigma = least_squares(self.two_site_inter_sigma_solver,initial,args = (Hai,Hbi,Sai,Sbi,Jabi,J1i,minadsi,maxadsi))
        lenopt_half = int(len(yopt_sigma['x'])/2)
        adsa_fit = yopt_sigma['x'][:lenopt_half]
        adsb_fit = yopt_sigma['x'][lenopt_half:]
        self.adsa_guess = adsa_fit
        self.adsb_guess = adsb_fit
        return np.append(self.ads1 - adsa_fit,self.ads2 - adsb_fit)

    def run_twosite_inter_optimise(self,x,bounds):

        yopt = least_squares(self.two_site_inter_optimise,x,bounds = bounds)


        self.Ha = yopt['x'][0]
        self.Hb = yopt['x'][1]
        self.Sa = yopt['x'][2]
        self.Sb = yopt['x'][3]
        self.Jab = yopt['x'][4]
        self.J1 = yopt['x'][5]
        if len(x) == 8:
            self.minads = yopt['x'][6]
            self.maxads = yopt['x'][7]
        adsa_guess,adsb_guess = self.twosite_nonequiv(self.Ha,self.Hb,self.Sa,self.Sb,self.Jab,self.minads,self.maxads,self.minads,self.maxads)
        initial = np.append(adsa_guess,adsb_guess)
        yopt_sigma = least_squares(self.two_site_inter_sigma_solver,initial,args = (self.Ha,self.Hb,self.Sa,self.Sb,self.Jab,self.J1,self.minads,self.maxads))
        lenopt_half = int(len(yopt_sigma['x'])/2)
        self.adsa_guess = yopt_sigma['x'][:lenopt_half]
        self.adsb_guess = yopt_sigma['x'][lenopt_half:]
        #print(yopt_sigma['fun'])

        self.fit1,self.fit2 = self.two_site_inter_fit(self.adsa_guess,self.adsb_guess,self.Ha,self.Hb,self.Sa,self.Sb,self.Jab,self.J1,self.minads,self.maxads)
        rw = np.sum(yopt['fun'])**2/np.sum(np.append(self.ads1,self.ads2)**2)
        string = f'''dHa = {self.Ha}
dHb = {self.Hb}
dSa = {self.Sa}
dSb = {self.Sb}
Jab = {self.Jab}
J1 = {self.J1}
min. ads. = {self.minads}
max. ads. = {self.maxads}
rw = {rw}
'''
        print(string)
        output_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_interpore_fit_pars.txt'
        f = open(output_filename,'w')
        f.write(string)
        f.close()

        fit_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_interpore_ads_fit.txt'
        np.savetxt(fit_filename,np.array([self.T,self.ads1,self.fit1,self.ads2,self.fit2]).transpose(),
        header = 'T(K) adsorption1 fit1 adsorption2 fit2')
        print('fitted data saved to',fit_filename)
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

    def adsorption_plot_inter(self):
        fig,(ax1,ax2) = plt.subplots(2,1)
        adsav = (self.ads1 + self.ads2)/2
        fitav = (self.fit1 + self.fit2)/2
        ax1.plot(self.T,self.fit1,'o',label = 'fit 1')
        ax1.plot(self.T,self.fit2,'o', label = 'fit 2')
        ax1.plot(self.T,self.adsa_guess, label = 'guess 1')
        ax1.plot(self.T,self.adsb_guess, label = 'guess 2')
        ax1.set_xlabel('Temperature (K)')
        ax1.set_ylabel('adsorption')
        ax1.set_xlim(min(self.T),max(self.T))
        ax1.legend()

        ax2.plot(self.T,self.ads1,'o',label = 'adsorption 1')
        ax2.plot(self.T,self.ads2,'o', label = 'adsorption 2')
        ax2.plot(self.T,adsav,'o', label = 'adsorption av.')
        ax2.plot(self.T,self.fit1, label = 'fit 1')
        ax2.plot(self.T,self.fit2, label = 'fit 2')
        ax2.plot(self.T,fitav, label = 'fit av.')
        ax2.set_xlabel('Temperature (K)')
        ax2.set_ylabel('adsorption')
        ax2.set_xlim(min(self.T),max(self.T))
        ax2.legend()

        plt.show()


class AdsorptionProfile:
    def __init__(self):
        self.T = np.array([])
        self.H = -5000
        self.S = -50
        self.J = 0
        self.p = 1
        self.minads = 0
        self.maxads = 1
        self.values = np.array([self.H,self.S,self.J,self.minads,self.maxads])
        self.fit = np.array([])
        self.minmax_refine = False
    def update_values(self):
        self.values = np.array([self.H,self.S,self.J,self.minads,self.maxads])
        self.gamma = (self.ads[0] - self.minads)/(self.maxads-self.minads)
    def read_file(self,file,delimiter=None,skiprows = None,unit = 'K'):
        '''
        uses numpy.loadtxt. Positional arguments: file. Keyword arguments:
        delimiter=None, skiprows=0, usecols=None.
        '''

        self.filename = file
        self.directory = os.path.split(file)[0]+'/'

        if self.directory == '/':
            self.directory = os.path.curdir +'/'
        R = 8.31446
        if file.endswith('.csv'):
            delimiter = ','
        if skiprows == None:
            for i in range(30):
                try:
                    profileArray = np.loadtxt(file,unpack = True, skiprows = i,
                    delimiter = delimiter,usecols=usecols)
                    break
                except:
                    if i == 29:
                        raise RuntimeError
                    else:
                        continue

        else:

            profileArray = np.loadtxt(file,unpack = True, skiprows = skiprows,
            delimiter = delimiter,usecols=usecols)
        self.T = profileArray[0]
        self.ads = {}
        for i,p in enumerate(profileArray[1:]):
            self.ads[i] = p
        self.maxads = max(self.ads)*1.1
        if unit == 'C':
            self.T += 273.15

        self.update_values()
        self.fit = self.one_site_simple(self.H,self.S,self.minads,self.maxads)

        self.Keq = self.gamma/(1-self.gamma)
        self.Keq_fit = np.exp(-self.H/(R*self.T) + self.S/R)

    def one_site_simple(self,Hi,Si,minadsi,maxadsi):
        deltaG = Hi - self.T*Si
        R = 8.31446
        gamma = np.exp(-deltaG/(R*self.T))*self.p/(1+np.exp(-deltaG/(R*self.T))*self.p)
        adsorption = gamma*(maxadsi-minadsi) + minadsi
        return adsorption
    def one_site_simple_optimise(self,x=np.array([None])):
        if x.any() == None:
            x = [self.H,self.S,self.minads,self.maxads]
        if len(x) == 4:
            adsfit = self.one_site_simple(*x)
        elif len(x) == 2:
            adsfit = self.one_site_simple(x[0],x[1],self.minads,self.maxads)
        return self.ads[1] - adsfit
    def run_optimise(self,x = np.array([None]),bounds = ([np.inf]*4,[np.inf]*4)):
        if x.any() == None and len(x) == 4:
            x = np.array([self.H,self.S,self.minads,self.maxads])
        elif x.any() == None and len(x) == 2:
            x = np.array([self.H,self.S])
            bounds = (bounds[0][:2],bounds[1][:2])

        yopt = least_squares(self.one_site_simple_optimise,x,bounds=bounds)
        self.H = yopt['x'][0]
        self.S = yopt['x'][1]

        if len(yopt['x']) == 4:
            self.minads = yopt['x'][2]
            self.maxads = yopt['x'][3]
        R = 8.31446
        self.fit = self.one_site_simple(self.H,self.S,self.minads,self.maxads)
        rw = np.sum(self.ads[1]-self.fit)**2/np.sum(self.ads[1]**2)
        string = f'''dH = {self.H}
dS = {self.S}
min. ads. = {self.minads}
max. ads = {self.maxads}
rw = {rw}
'''
        print(string)
        output_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_fit_pars.txt'
        f = open(self.directory+output_filename,'w')
        f.write(string)
        f.close()
        print('fitted parameters saved to',output_filename)

        self.update_values()
        self.Keq = self.gamma/(1-self.gamma)
        self.Keq_fit = np.exp(-self.H/(R*self.T) + self.S/R)
        fit_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_ads_fit.txt'
        np.savetxt(fit_filename,np.array([self.T,self.ads[1],self.fit]).transpose(),
        header = 'T(K) adsorption fit')
        print('fitted data saved to',fit_filename)

        return yopt

    def one_site_coop(self,Hi,Si,Ji,minadsi,maxadsi):

        gamma = (self.ads[1]-minadsi)/(maxadsi-minadsi)
        R = 8.31446

        #gamma = np.arange(0.001,0.999,0.001)
        Tfit = (Hi-2*Ji*(2*gamma-1))/(Si-R*np.log(gamma/(self.p*(1-gamma))))
        return Tfit

    def one_site_coop_optimise(self,x):
        Hi = x[0]
        Si = x[1]
        Ji = x[2]
        minadsi = self.minads
        maxadsi = self.maxads
        if len(x) == 5:
            minadsi = x[3]
            maxadsi = x[4]

        Tfit = self.one_site_coop(Hi,Si,Ji,minadsi,maxadsi)
        return self.T - Tfit
    def run_coop_optimise(self,x,bounds):
        yopt = least_squares(self.one_site_coop_optimise,x,bounds = bounds)
        params = yopt['x']
        self.H = params[0]
        self.S = params[1]
        self.J = params[2]
        if len(x) == 5:
            self.minads = params[3]
            self.maxads = params[4]
        self.update_values()
        self.fit = self.one_site_coop(self.H,self.S,self.J,self.minads,self.maxads)
        rw = np.sum(self.T-self.fit)**2/np.sum(self.T**2)
        output = f'''deltaH = {self.H}
deltaS = {self.S}
J = {self.J}
min. ads. = {self.minads}
max. ads. = {self.maxads}
rw = {rw}'''
        print(output)
        output_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_1sitecoop_fit_pars.txt'
        f = open(self.directory+output_filename,'w')
        f.write(output)
        f.close()
        print('fitted parameters saved to',output_filename)
        fit_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_1sitecoop_ads_fit.txt'
        np.savetxt(fit_filename,np.array([self.T,self.ads[1],self.fit]).transpose(),
        header = 'T(K) adsorption fit')
        print('fitted data saved to',fit_filename)
        return yopt

    def twosite_1profile(self,Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi):
        R = 8.31446
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
        #sigtot = (siga+sigb)/2
        return adsav
    def twosite_1profile_optimise(self,x):
        Hai = x[0]
        Hbi = x[1]
        Sai = x[2]
        Sbi = x[3]
        Jabi = x[4]
        minadsi = self.minads
        maxadsi = self.maxads
        if len(x) == 7:
            minadsi = x[5]
            maxadsi = x[6]
        adsavfit = self.twosite_1profile(Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi)
        return self.ads[1] - adsavfit

    def run_twosite_optimise(self,x,bounds):
        yopt = least_squares(self.twosite_1profile_optimise,x,bounds = bounds)
        params = yopt['x']
        self.H = params[0]
        self.Hb = params[1]
        self.S = params[2]
        self.Sb = params[3]
        self.J = params[4]
        if len(x) == 7:
            self.minads = params[5]
            self.maxads = params[6]
        self.update_values()
        self.fit = self.twosite_nonequiv(self.H,self.Hb,self.S,self.Sb,self.J,self.minads,self.maxads)
        rw = np.sum(self.ads[1] - self.fit)**2/np.sum(self.ads[1]**2)
        string = (f"deltaHa = {self.H}\n"
        f"deltaHb = {self.Hb}\n"
        f"deltaSa = {self.S}\n"
        f"deltaSb {self.Sb}\n"
        f"Jab = {self.J}\n"
        f"min. ads. = {self.minads}\n"
        f"max. ads. = {self.maxads}\n"
        f"rw = {rw}")
        print(string)
        output_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_fit_pars.txt'
        f = open(self.directory+output_filename,'w')
        f.write(string)
        f.close()
        print('fitted parameters saved to',output_filename)
        fit_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_ads_fit.txt'
        np.savetxt(fit_filename,np.array([self.T,self.ads[1],self.fit]).transpose(),
        header = 'T(K) adsorption fit')
        print('fitted data saved to',fit_filename)
        return yopt
    def adsorption_plot(self):
        plt.plot(self.T,self.ads[1],'o',label = 'adsorption')
        plt.plot(self.T,self.fit, label = 'fit')
        plt.xlim(min(self.T),max(self.T))
        plt.xlabel('Temperature (K)')
        plt.ylabel('Adsorption')
        plt.legend()
        plt.show()
    def vant_hoff_plot(self):
        plt.plot(1/self.T,np.log(self.Keq),'o',label = 'measured (min. max. normalised)')
        plt.plot(1/self.T,np.log(self.Keq_fit),label = 'fit')
        plt.xlim(min(1/self.T),max(1/self.T))
        plt.xlabel('1/T (K$^{-1}$)')
        plt.ylabel('ln(K)')
        plt.legend()
        plt.show()
    def ads_vant_hoff_plot(self):
        fig, (ax1,ax2) = plt.subplots(2,1)
        ax1.plot(self.T,self.ads[1],'o',label = 'adsorption')
        ax1.plot(self.T,self.fit, label = 'fit')
        ax1.set_xlim(min(self.T),max(self.T))
        ax1.set_xlabel('Temperature (K)')
        ax1.set_ylabel('Adsorption')
        ax1.legend()
        ax2.plot(1/self.T,np.log(self.Keq),'o',label = 'measured (min. max. normalised)')
        ax2.plot(1/self.T,np.log(self.Keq_fit),label = 'fit')
        ax2.set_xlim(min(1/self.T),max(1/self.T))
        ax2.set_xlabel('1/T (K$^{-1}$)')
        ax2.set_ylabel('ln(K)')
        ax2.legend()
        plt.show()

    def one_site_coop_plot(self):
        plt.plot(self.T,self.ads[1],'o',label = 'adsorption')
        plt.plot(self.fit,self.ads[1],label = 'fit')
        plt.xlabel('Temperature (K)')
        plt.ylabel('adsorption')
        plt.xlim(min(self.T),max(self.T))
        plt.legend()
        plt.show()

    #two site models


    def twosite_intra(self,Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi,minadsi2,maxadsi2):
        R = 8.31446
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
        adsfit2 = gamma_b*(maxadsi2-minadsi2)+minadsi2
        #adsav = (adsfit1+adsfit2)/2
        #sigtot = (siga+sigb)/2
        return adsfit1, adsfit2


    def twosite_intra_optimise(self, x=np.array([None])):
        if x.any() == None:
            x = self.values

        Hai = x[0]
        Hbi = x[1]
        Sai = x[2]
        Sbi = x[3]
        Jabi = x[4]
        minadsi = self.minads
        maxadsi = self.maxads
        minadsi2 = self.minads
        maxadsi2 = self.maxads
        if len(x) >= 7:
            minadsi = x[5]
            maxadsi = x[6]
            minadsi2 = x[5]
            maxadsi2 = x[6]
        if len(x) == 9:
            minadsi2 = x[7]
            maxadsi2 = x[8]


        adsfit1,adsfit2 = self.twosite_intra(Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi,minadsi2,maxadsi2)

        yresid  = np.append(self.ads[1]-adsfit1,self.ads[2]-adsfit2)
        return yresid



    def run_optimise(self,x=np.array([None]),bounds = ([None],[None])):
        if x.any() == None:
            x = self.values[:7]
        #if any(bounds[0]) == False:
#
        #    bounds = ([-1*np.inf]*len(x),[np.inf]*len(x))

        yopt = least_squares(self.twosite_intra_optimise,x,
                                   bounds = bounds)
        self.Ha = yopt['x'][0]
        self.Hb = yopt['x'][1]
        self.Sa = yopt['x'][2]
        self.Sb = yopt['x'][3]
        self.Jab = yopt['x'][4]
        if len(x) >= 7:
            self.minads = yopt['x'][5]
            self.maxads = yopt['x'][6]
            self.minads2 = yopt['x'][5]
            self.maxads2 = yopt['x'][6]
        if len(x) == 9:
            yopt = least_squares(self.twosite_nonequiv_minmax2_optimise,initial,
                       bounds = bounds)
            self.minads2 = yopt['x'][7]
            self.maxads2 = yopt['x'][8]




        self.update_values()

        rw = np.sum(yopt['fun'])**2/np.sum(np.append(self.ads[1],self.ads[2])**2)
        string = f'''dHa = {self.Ha}
dHb = {self.Hb}
dSa = {self.Sa}
dSb = {self.Sb}
Jab = {self.Jab}
min. ads. = {self.minads}
max. ads. = {self.maxads}
rw = {rw}
'''





        self.fit1, self.fit2 = self.twosite_nonequiv(*self.values)
        if len(x) == 9:

            string += f'min. ads. 2 = {self.minads2}\nmax. ads. 2 = {self.maxads2}'
        output_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_fit_pars.txt'
        print(string)
        f = open(self.directory+output_filename,'w')
        f.write(string)
        f.close()
        print('fitted parameters saved to',output_filename)

        fit_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_intrapore_ads_fit.txt'
        np.savetxt(fit_filename,np.array([self.T,self.ads[1],self.fit1,self.ads[2],self.fit2]).transpose(),
        header = 'T(K) adsorption1 fit1 adsorption2 fit2')
        print('fitted data saved to',fit_filename)

        return yopt



    def two_site_inter_sigma_solver(self,x,Hai,Hbi,Sai,Sbi,Jabi,J1i,minadsi,maxadsi):
        R = 8.31446
        Ga = Hai-self.T*Sai
        Gb = Hbi-self.T*Sbi
        lenx_half = int(len(x)/2)
        adsa_guess = x[:lenx_half]
        adsb_guess = x[lenx_half:]
        sigmaa_guess = 2*((adsa_guess-minadsi)/(maxadsi-minadsi))-1
        sigmab_guess = 2*((adsb_guess-minadsi)/(maxadsi-minadsi))-1
        sigmaa = np.tanh((-Ga-2*Jabi*sigmab_guess-2*J1i*sigmaa_guess)/(R*self.T))
        sigmab = np.tanh((-Gb-2*Jabi*sigmaa_guess-2*J1i*sigmab_guess)/(R*self.T))

        return np.append(sigmaa-sigmaa_guess,sigmab-sigmab_guess)

    def two_site_inter_fit(self,adsa_guess,adsb_guess,Hai,Hbi,Sai,Sbi,Jabi,J1i,minadsi,maxadsi):
        R = 8.31446
        Ga = Hai-self.T*Sai
        Gb = Hbi-self.T*Sbi
        sigmaa_guess = 2*((adsa_guess-minadsi)/(maxadsi-minadsi))-1
        sigmab_guess = 2*((adsb_guess-minadsi)/(maxadsi-minadsi))-1
        sigmaa = np.tanh((-Ga-2*Jabi*sigmab_guess-2*J1i*sigmaa_guess)/(R*self.T))
        sigmab = np.tanh((-Gb-2*Jabi*sigmaa_guess-2*J1i*sigmab_guess)/(R*self.T))
        gammaa = (sigmaa+1)/2
        gammab = (sigmab+1)/2
        adsa = gammaa*(maxadsi-minadsi) + minadsi
        adsb = gammab*(maxadsi-minadsi) + minadsi
        return adsa,adsb

    def two_site_inter_optimise(self,x):
        Hai = x[0]
        Hbi = x[1]
        Sai = x[2]
        Sbi = x[3]
        Jabi = x[4]
        J1i = x[5]
        minadsi = self.minads
        maxadsi = self.maxads
        if len(x) == 8:
            minadsi = x[6]
            maxadsi = x[7]
        adsa_guess,adsb_guess = self.twosite_nonequiv(Hai,Hbi,Sai,Sbi,Jabi,minadsi,maxadsi,minadsi,maxadsi)
        initial = np.append(adsa_guess,adsb_guess)
        yopt_sigma = least_squares(self.two_site_inter_sigma_solver,initial,args = (Hai,Hbi,Sai,Sbi,Jabi,J1i,minadsi,maxadsi))
        lenopt_half = int(len(yopt_sigma['x'])/2)
        adsa_fit = yopt_sigma['x'][:lenopt_half]
        adsb_fit = yopt_sigma['x'][lenopt_half:]
        self.adsa_guess = adsa_fit
        self.adsb_guess = adsb_fit
        return np.append(self.ads[1] - adsa_fit,self.ads[2] - adsb_fit)

    def run_twosite_inter_optimise(self,x,bounds):

        yopt = least_squares(self.two_site_inter_optimise,x,bounds = bounds)


        self.Ha = yopt['x'][0]
        self.Hb = yopt['x'][1]
        self.Sa = yopt['x'][2]
        self.Sb = yopt['x'][3]
        self.Jab = yopt['x'][4]
        self.J1 = yopt['x'][5]
        if len(x) == 8:
            self.minads = yopt['x'][6]
            self.maxads = yopt['x'][7]
        adsa_guess,adsb_guess = self.twosite_nonequiv(self.Ha,self.Hb,self.Sa,self.Sb,self.Jab,self.minads,self.maxads,self.minads,self.maxads)
        initial = np.append(adsa_guess,adsb_guess)
        yopt_sigma = least_squares(self.two_site_inter_sigma_solver,initial,args = (self.Ha,self.Hb,self.Sa,self.Sb,self.Jab,self.J1,self.minads,self.maxads))
        lenopt_half = int(len(yopt_sigma['x'])/2)
        self.adsa_guess = yopt_sigma['x'][:lenopt_half]
        self.adsb_guess = yopt_sigma['x'][lenopt_half:]
        #print(yopt_sigma['fun'])

        self.fit1,self.fit2 = self.two_site_inter_fit(self.adsa_guess,self.adsb_guess,self.Ha,self.Hb,self.Sa,self.Sb,self.Jab,self.J1,self.minads,self.maxads)
        rw = np.sum(yopt['fun'])**2/np.sum(np.append(self.ads[1],self.ads[2])**2)
        string = f'''dHa = {self.Ha}
dHb = {self.Hb}
dSa = {self.Sa}
dSb = {self.Sb}
Jab = {self.Jab}
J1 = {self.J1}
min. ads. = {self.minads}
max. ads. = {self.maxads}
rw = {rw}
'''
        print(string)
        output_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_interpore_fit_pars.txt'
        f = open(output_filename,'w')
        f.write(string)
        f.close()

        fit_filename = os.path.splitext(os.path.split(self.filename)[-1])[0] + '_interpore_ads_fit.txt'
        np.savetxt(fit_filename,np.array([self.T,self.ads[1],self.fit1,self.ads[2],self.fit2]).transpose(),
        header = 'T(K) adsorption1 fit1 adsorption2 fit2')
        print('fitted data saved to',fit_filename)
        return yopt

    def adsorption_plot(self):
        adsav = (self.ads[1] + self.ads[2])/2
        fitav = (self.fit1 + self.fit2)/2
        plt.plot(self.T,self.ads[1],'o',label = 'adsorption 1')
        plt.plot(self.T,self.ads[2],'o', label = 'adsorption 2')
        plt.plot(self.T,adsav,'o', label = 'adsorption av.')
        plt.plot(self.T,self.fit1, label = 'fit 1')
        plt.plot(self.T,self.fit2, label = 'fit 2')
        plt.plot(self.T,fitav, label = 'fit av.')
        plt.xlabel('Temperature (K)')
        plt.ylabel('adsorption')
        plt.xlim(min(self.T),max(self.T))
        plt.legend()
        plt.show()

    def adsorption_plot_inter(self):
        fig,(ax1,ax2) = plt.subplots(2,1)
        adsav = (self.ads[1] + self.ads[2])/2
        fitav = (self.fit1 + self.fit2)/2
        ax1.plot(self.T,self.fit1,'o',label = 'fit 1')
        ax1.plot(self.T,self.fit2,'o', label = 'fit 2')
        ax1.plot(self.T,self.adsa_guess, label = 'guess 1')
        ax1.plot(self.T,self.adsb_guess, label = 'guess 2')
        ax1.set_xlabel('Temperature (K)')
        ax1.set_ylabel('adsorption')
        ax1.set_xlim(min(self.T),max(self.T))
        ax1.legend()

        ax2.plot(self.T,self.ads[1],'o',label = 'adsorption 1')
        ax2.plot(self.T,self.ads[2],'o', label = 'adsorption 2')
        ax2.plot(self.T,adsav,'o', label = 'adsorption av.')
        ax2.plot(self.T,self.fit1, label = 'fit 1')
        ax2.plot(self.T,self.fit2, label = 'fit 2')
        ax2.plot(self.T,fitav, label = 'fit av.')
        ax2.set_xlabel('Temperature (K)')
        ax2.set_ylabel('adsorption')
        ax2.set_xlim(min(self.T),max(self.T))
        ax2.legend()

        plt.show()
