# Clasts-Pattern-Mars
The codes are to reproduce the figures and data analysis within our manuscript titled "Hyperuniformity on Mars: Pebbles scattered on sand" by Zheng Zhu, Bernard Hallet, Andras Sipos, Gabor Domokos, Quan-Xing Liu*, submitted to 2024. arXiv [Preprint] (2023). https://doi.org/10.48550/arXiv.2312.13818

The codes are classified by their function on model simulation, statistical analysis, and spatial hyperuniformity. 


# Code:

The fold within "sk and df" were modified from our previous pulished paper,

    sk and df # this fold icnludes all the statistical analysis on hyperuniform states. 
    Huang MJ et al. 2021.Circular swimming motility and disordered hyperuniform state in an algae system, PNAS 118: e2100493118 

The fold "Model simulation"  implement the clast reptation on sand ripples.

    Model simulation # this file contains all simulation functions 
    Werner BT. 1995. Eolian dunes: computer simulation and attractor interpretation. Geology 23: 1107–1110.

It is sorted by the figures in the article, so there may be duplicates in different folders.
The data and code list reads as follows: codes, FIG_xx_yy, MSL image raw data, and Gobi image raw data, where xx and yy represent the figure number and data source, respectively. 

`density_fluctuation.m` :Two-dimensional density_fluctuation of a set of points.

`struct_factor_sq.m` :Two-dimensional static structure factor of a set of points.

`clastsrun.m` : individual-based model of the clasts movements on aeolian sands written in MATLAB 2019a.

`spectral2D.m` :Two-dimensional spectral analysis of an image file, from http://www.nioo.knaw.nl/homepages/koppel

`main.m` :Computation of obeserved data and possion point process in Fig.2.

`rippleformation.m` :sand dune model of Werner 1995.

`sinuosity.m` :Mean sinuosity of the sand ripples from a binary image.

`plot*.m` :Plot of the figure in each folder.

`powerspectrum2d.m` :Two-dimensional powerspectrum of an image file.

`DJI_*.JPG`  :Photos taken by DJI drones at the Gobi desert.

# Dataset

`MSL image Raw data` :Mars images in Gale Crater downloaded from https://mars.nasa.gov/msl/home/.

`*.csv` and `*.mat`:files that save the results of the analyses in each folder for further ploting used in main figures and SI figures.

`Gobi image Raw data` :Images of Gobi clasts were collected from Alxa Left Banner. 

`Movies ripple migration` :Images and movies show the spatial scale and migration speed of ripples  alongside Gobi clasts in Alxa Left Banner. 
