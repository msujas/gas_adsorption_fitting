
import adsorptionFitClass
import numpy as np

initial2 = np.array([-5000,-5000,-50,-50,0,0,15,0,15])
bounds = (np.array([-50000,-50000,-1000,-1000,-100,-1,14]),
          np.array([0,0,1000,1000,20000,0,20]))
bounds2 = (np.array([-50000,-50000,-1000,-1000,-100,-1,14,-1,14]),
          np.array([0,0,1000,1000,20000,0,30,0,17]))

Ar_ads_data = adsorptionFitClass.One_site_adsorption_profile()
direc2 = r'C:\Users\kenneth1a\Documents\gas_absorption_lectures/'
Ar_ads_data.read_file(direc2+'generation of plastic crystal methane rotator.csv',unit = 'C')
#print(Ar_ads_data.values)

bounds1site = (np.array([-50000,-1000,-1,Ar_ads_data.ads.max()]),
          np.array([0,1000,0,20]))
yopt = Ar_ads_data.run_optimise(x = np.array([-5000,-50,0,2.2]),bounds = bounds1site)

#print(Ar_ads_data.values)

#Ar_ads_data.adsorption_plot()
Ar_ads_data.ads_vant_hoff_plot()
