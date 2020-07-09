# ExaStoLog

A MATLAB toolbox for the exact calculation of stationary states + parameter sensitivity analysis & fitting of stochastic logical models.  
Author: Mihály Koltai, [Computational Systems Biology of Cancer group](https://github.com/sysbio-curie) at Institut Curie.

## Read the [tutorial here](https://github.com/mbkoltai/exact-stoch-log-mod/tree/master/doc).

### Requirements

#### - MATLAB version 2015b or later

#### - clone the repository

#### - unzip the file 'toolboxes.zip' for the external MATLAB libraries used:

- [export_fig](https://mathworks.com/matlabcentral/fileexchange/23629-export_fig) (export figures as EPS/PDF/PNG as they appear on screen)  
- [Customizable heatmaps](https://mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps) (for heatmaps)  
- [Redblue colormap](https://mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap) (for heatmaps)  
- [distinguishable_colors](https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors) (for lineplots of one dimensional parameter scans. Requires Image Processing toolbox!)  
- [tight subplots](https://mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w) (for customizable gaps between subplots)  
- [Simulated annealing](https://mathworks.com/matlabcentral/fileexchange/10548-general-simulated-annealing-algorithm) (parameter fitting by simulated annealing)  

These libraries are for visualizations and parameter fitting, not for the calculations themselves, so they are optional, but we recommend using them for all the features of the toolbox.

#### - enter the directory and add the folders to the path by running 'add_toolboxes_paths'

For the detailed examples of calculations read the [tutorial here](https://github.com/mbkoltai/exact-stoch-log-mod/tree/master/doc).
