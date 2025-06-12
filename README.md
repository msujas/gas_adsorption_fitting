# gas_adsorption_fitting
scripts for fitting gas adsorption profiles
Includes 4 models. To install do 'pip install -e .' in the repository. Then to run do "adsFitGui" in a terminal, select a file containing temperature vs adsorption data 
(some example files given, can be with one or two site data, i.e. x,y or x,y1,y2), select a model, starting values and click Run Fit. Will plot a fit 
and print thermodynamic parameters on side. The Update initial values button will put the refined values into the initial values.

Uses SciPy to run fits, GUI made with PyQt5.

The GUI

![image](https://github.com/user-attachments/assets/e193e838-4004-4f00-8cf9-22423dd69952)

Single-site, simple fit

![image](https://github.com/user-attachments/assets/9dd8d55f-87de-46b4-9651-2d1676477ca5)



