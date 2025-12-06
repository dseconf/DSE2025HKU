Instructions for Estimation Code
********************************


General Information
*******************

All code is written in R and C++.  The C++ code is called from the Rcpp package, so you do not need to install anything separately.

Every estimation run or counterfactual has a starting file which calls subfiles and produces output (estimation_master*.r, fit_data.r, fit_data_clust*.r).  Before running a starting file, you need to set the following path variables:

datapath is the directory where the estimation data is
fpath is the file path which points to where the code is being run


Generating the main model estimates
***********************************

The code to generate the main set of estimates is in the folder called "main estimates".  You need to edit the file estimation_master2.r, fill in the path variables, and then run this file.  Since it take a long time to estimate the model, our suggestion is to run it in batch mode.  The main output file will be called est2_mcmcout.RData.  Note that the main file loads in some starting information which is in four RData files we have supplied.  To generate these files from scratch, follow the steps listed under the heading "Generating initial files" below.

The main code will first call estimation_setup2.r, which will initialize parameters needed for running the MCMC loop.  It then will run run_initial_ml_data2.r, which will calculate the inverse information matrix at the starting parameters, which is used as a proposal for the population-fixed parameters.  The main MCMC loop is executed in run_ijc_data2.r.  Note that the code to actually run the MCMC loop is written in C++, and is in the file ijc_Rcpp_par.cpp.  The C++ function MCMCLoops will run the MCMC loop and return the results in an R list object.  This object is what gets saved in est2_mcmcout.RData.  To see an example of how to work with the data, in the file make_tables.r under the "estimation tables" directory we provide sample code to produce the main tables of estimates and discount factor heterogeneity plot. To run this you need to copy over a file named mcmc_starting_params2.RData, which is produced by running the main estimation file, to the same directory as make_tables.r.

Note that the 2 file suffix indexes a particular model specification, which is defined by the spec.index variable in the starting file.  Index 1 is an initial run that is used to generate starting information (see below).

Generating initial files
*************************

1. In the folder step 1:
a) run fit_data.r.  This will run an initial maximum likelihood estimation assuming individuals are myopic.
b) this will produce a file called ml_fit.RData.  Copy this to the step 2 folder.

2. In the folder step 2:
a) run estimation_master1.r.  This runs a simplified model without parameter heterogeneity, which is used in the full model to generate the variance matrix for drawing out the fixed parameters.
b) the main output file you need is called bdrawstart.RData.  This will be accessed by all the estimation files.  The file will also produce a larger file called est1_mcmcout.RData which is used to construct bdrawstart.RData.

3. In the folder step 3:
a) run the files fit_data_clust1_1.r, fit_data_clust1_2.r and fit_data_clust1_3.r.  These will generate the files out_ml_A_fix1_3_*.RData which are also used to start the main MCMC algorithm.


Generating the additional model estimates
*****************************************

In the paper we also present estimation results from a model specification that assumes the discount factor is calibrated to the interest rate, and a myopic specification.  The code to run these specifications can be found in the folders "beta 0.9995 estimates" and "myopic estimates", respectively.  They are indexed with specification numbers 3 and 4.  The results can be produced by running the estimation_master files in the same way as the main files.


Generating the counterfactuals
******************************

The code to generate the counterfactuals is located under the directory counterfactuals.  This directory has two subdirectories named "main" and "fixed discount factor", which run the counterfactual calculations for each scenario.  Each directory has five subdirectories, which correspond to different counterfactual simulations.  The directory named "base" runs the counterfactual simulation at the estimated parameters and observed data.  To run this, you should copy the results files (est2_mcmcout.RData and estml2.RData) from the main specification over to this folder, edit counterfactual_master1.r, and insert the path variables.  Then run the R file.  It will produce a results file called cfout1.RData.  Note that there is a variable called spec.index in the counterfactual master file that tells the code which estimation input file to use.  There is also a variable called cf.index which is a number for the particular counterfactual.  The other counterfactuals (increased promotional frequency and depth for long-term and short-term) are in folders with the appropriate names.  Counterfactuals for the fixed discount factor specification are under the folder "fixed discount factor" and have the same subfolder names.  They load in estimation results from specification 3.

To see how to work with the counterfactual output, in the "estimation tables" directory we include a file called counterfactual_tables.r that will create counterfactual tables.  This file has a variable called cfpath that should point the program to the location of the counterfactual output files, which are named cfout*.RData.  The * will refer to one of the indexes of a counterfactual, which corresponds to the file name.  The correspondence is also shown in the table file.  After running the counterfactual code, we would suggest one to move all the cfout files to a single directory.






