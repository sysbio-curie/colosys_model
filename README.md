# README for 'colosys_model' folder
### July/2020
### Mihaly Koltai

## Folders

## data



## exastolog

This folder contains the files to run [ExaStoLog](https://github.com/sysbio-curie/exact-stoch-log-mod/) that requires MATLAB (>=2016a).
It is the file *colosys_model.m* that contains the scripts specific to COlOSYS and the subfolders *colosys_plots* and *model_files*.
It is mainly the model '*krasreduced_dnarep_simplif.bnet*' that I worked with.
In *colosys_model.m* the initial conditions are set up for this model version.
First an individual simulation and a 1-dimensional parameter scan is done. This is followed by script for parameter fitting for the transition rates that were found to have a significant (>0.1 or >0.2 or some other threshold) effect on the stationary solutions (of at least some node).
The parameter fitting is done by lsqnonlin (trust region least squares) which shows the fastest convergence.
Since this is a local search method I performed multistart (different initial values for parameters) searches and plotted the results by boxplots and histograms. This is to see which parameters converge consistently to a certain value, ie. are identifiable.

Then there are scripts to pull in the data from the *'data/'* folder (*'INTEGRATE cell line data'*): viability, WES, Cytof. At the *viability* data there are also some lines of code where I tried to perform nonlinear regression (with sigmoidals in the transition rates) on the stationary solutions, but have not pursued this further.

I eventually switched to pymaboss since ExaStoLog was slowed down by the presence of numerous large cycle of the STG of the models that are due to negative feedbacks.

## network_paths


## pymaboss
