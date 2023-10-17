# gas_adsorption_fitting
scripts for fitting gas adsorption profiles
Includes 4 models. To run do "python adsorptionFitGui.py" in a terminal, select a file containing temperature vs adsorption data 
(some example files given, can be with one or two site data, i.e. x,y or x,y1,y2), select a model, starting values and click Run Fit. Will plot a fit 
and print thermodynamic parameters on side. The Update initial values button will put the refined values into the initial values.

Uses SciPy to run fits, GUI made with PyQt5.
