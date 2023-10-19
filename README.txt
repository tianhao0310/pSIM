--- Codes for the simulation studies of the piecewise SIM

--- Authors: Tianhao Wang & Yingcun Xia

--- DESCRIPTION: 

This folder contains the simulation files and the auxiliary functions which are used to estimate the piecewise SIM

--- FUNCTION LIST:
------ cvadap.m: a function to infer the smoothing bandwith though the cross-validation method
------ knn.m: a function for k-nearest neighbor classification
------ ksLLadap.m: local linear kernel smoothing function
------ m2.m:  a function to calculate the distances of the colomns between two matrixs
------ pSIM.m: a function to cluster the pointwise local gradients into a given number of groups
------ pSIMfit.m: fit a piecewise SIM for a given data
------ rMAVE.m: the rMAVE method of Xia et al (2002)
------ rMAVEfit.m: find the e.d.r. directions with rMAVE
------ rOPGadap.m: function to find the pointwise local gradients
------ rOPGfit.m: find the e.d.r. directions with rOPGadap

--- DATASET LIST:
------ cars.mat: a dataset that studies the fuel efficiency for automobiles; it is from the ASA Data Exposition dataset (1983) collected by Professor Ernesto Ramos and Professor David Donoho.
------ LAozone.mat: a dataset that studies the atmospheric ozone concentration in Los Angeles basin (Breiman and Friedman 1985).
------ salaryRemo3.mat: a dataset for Hitters' salary data with the outliers removed, which was originally given at 1988 ASA Graphics Poster Session (Chaudhuri et al., 1994).
------ Irnk258.mat: standard random numbers ranging from 1 to 258 for out-of-sample prediction exercise on dataset "salaryRemo3"
------ Irnk330.mat: standard random numbers ranging from 1 to 330 for out-of-sample prediction exercise on dataset "LAozone"
------ Irnk392.mat: standard random numbers ranging from 1 to 392 for out-of-sample prediction exercise on dataset "cars"

--- SIMULATION FILES LIST
------ simuCarsFit.m: fit the piecewise single index model for the dataset "cars"
------ simuCarsPF.m: Random out-of-sample prediction exercise for the dataset "cars"
------ simuLAozoneFit.m: fit the piecewise single index model for the dataset "LAozone"
------ simuLAozonePF.m: Random out-of-sample prediction exercise for the dataset "LAozone"
------ simuLinPF.m: the simulation file for the Example 1 of the paper
------ simuNlinPF.m: the simulation file for the Example 2 of the paper
------ simuSalaryFit.m: fit the piecewise single index model for the dataset "salaryRemo3"
------ simuSalaryPF.m: Random out-of-sample prediction exercise for the dataset "salaryRemo3"

--- USAGE
------ Add the path of this folder in the "addpath(...)" command at the top of the simulation files, then run the simulation files. 
------ The results of "simuLAozoneFit.m" and "simuSalaryFit.m" are plotted directly or displayed in the "Command window"
------ The results of  "simuLinPF.m", "simuNlinPF.m", "simuLAozonePF.m" and "simuSalaryPF.m" are saved with the names defined at the end of each simualtion file.

--- REMARKS
------ All codes in this folder are tested on Matlab R2010a; some commands within the functions may not be supported by earlier versions
------ The files "simuLinPF.m" and "simuNlinPF.m" use the parallel toolbox in Matlab; it may have some problem to run them on a PC.

--- REFERENCE
------ Breiman, L. and Friedman, J. H. (1985) Estimating optimal transformations for multiple regression and correlation. Journal of the American Statistical Association 80, 580-597.
------ Chaudhuri, P., Huang, M. C., Loh, W. Y., and Yao, R. (1994) Piecewise polynomial regression trees. Statistica Sinica 4, 143-167.

--- Last update: 19/06/2013