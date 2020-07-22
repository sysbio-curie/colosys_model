# README for 'colosys_model' folder
### July/2020
### Mihaly Koltai

## Folders

The R files in the root folder are from earlier work when I was running MaBoSS from R on the command line with a smaller KRAS model (15 nodes). I changed the approach and methods in 2019, and switched to ExaStoLog and PyMaboss. I leave the R files here as an archive.

*Gene list.bublitz.xlsx* contains a list of genes

## data
The subfolders
*crispri/*  
*cytof/*  
*viability/*  
*wes_facs/*

contain the different data types. In *plots/* the cytof, wes and viability (at 48h) data are summarized on the 3 pdf files. These are the ones used for fitting.
The tidy format data tables per subfolder are:
- cytof: *'All means_DMSO_alllog2.csv'*, *'All means_DMSO_allraw.csv'*
- viability: *'viab_df_tidy_all.csv'*
- wes: *'wes_facs_data_tidy.csv'*

The crispri I haven't yet (20/july/2020) worked with.

## exastolog

This folder contains the files to run [ExaStoLog](https://github.com/sysbio-curie/exact-stoch-log-mod/) that requires MATLAB (>=2016a).
It is the file *colosys_model.m* that contains the scripts specific to COlOSYS and the subfolders *colosys_plots* and *model_files*.
It is mainly the model '*krasreduced_dnarep_simplif.bnet*' that I worked with.
In *colosys_model.m* the initial conditions are set up for this model version.
First an individual simulation and a 1-dimensional parameter scan is done. This is followed by script for parameter fitting for the transition rates that were found to have a significant (>0.1 or >0.2 or some other threshold) effect on the stationary solutions (of at least some node).
The parameter fitting is done by lsqnonlin (trust region least squares) which shows the fastest convergence.
Since this is a local search method I performed multistart (different initial values for parameters) searches and plotted the results by boxplots and histograms. This is to see which parameters converge consistently to a certain value, ie. are identifiable.

Then there are scripts to pull in the data from the *'data/'* folder (*'INTEGRATE cell line data'*): viability, WES, Cytof. At the *viability* data there are also some lines of code where I tried to perform nonlinear regression (with sigmoidals in the transition rates) on the stationary solutions, but have not pursued this further.

I eventually switched to pymaboss since ExaStoLog was slowed down by the presence of numerous large cycles in the STG of the models that are due to the negative feedbacks.

## network_paths

*'slides_scCRC1.pdf'* is the network drawing from Charite that was the basis of our models. There is a readme within this folder that explains files.  

## pymaboss

*'COLOSYS.ipynb'* is the main file for running PyMaboss simulations. I do ensemble modeling in the script, ie. there are a few nodes for which there are a number of possible logical rules, then I generate all model permutations from this. These different model versions are simulated and scored by how many 'behavioral criteria' they satisfy. These are for instance how proliferation/cell death change due to the presence of the mk2i/chk1i inhibitions.
*model_ens* contains the different model versions.
*figs* contains the figures generated from simulations.
*outputs* contains csv files that are model outputs from parameter scans.
*'own_functions.py'* contains function for simulations.
