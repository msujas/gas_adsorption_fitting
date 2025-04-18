# gas_adsorption_fitting
scripts for fitting gas adsorption profiles
Includes 4 models. To install do 'pip install -e .' in the repository. Then to run do "adsFitGui" in a terminal, select a file containing temperature vs adsorption data 
(some example files given, can be with one or two site data, i.e. x,y or x,y1,y2), select a model, starting values and click Run Fit. Will plot a fit 
and print thermodynamic parameters on side. The Update initial values button will put the refined values into the initial values.

Uses SciPy to run fits, GUI made with PyQt5.

The GUI

![image](https://github.com/msujas/gas_adsorption_fitting/assets/79653376/f0efb3e6-f14b-475a-a478-2ae43956dca0)

Single-site, simple fit

![image](https://github.com/msujas/gas_adsorption_fitting/assets/79653376/9a123e8a-9913-48f1-becc-c4b0a29c3e4f)


