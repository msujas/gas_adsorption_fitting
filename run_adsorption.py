
import Ar_adsorption_fit_class

initial2 = np.array([-5000,-5000,-50,-50,0,0,15,0,15])
bounds = (np.array([-50000,-50000,-1000,-1000,-100,-1,14]),
          np.array([0,0,1000,1000,20000,0,20]))
bounds2 = (np.array([-50000,-50000,-1000,-1000,-100,-1,14,-1,14]),
          np.array([0,0,1000,1000,20000,0,30,0,17]))
Ar_ads_data = Two_site_adsorption_profile()

Ar_ads_data.read_file('Ar_occupancy_1bar.txt',skiprows = 1)
#print(Ar_ads_data.values)
yopt = Ar_ads_data.run_optimise(bounds = bounds)

#print(Ar_ads_data.values)

Ar_ads_data.adsorption_plot()
