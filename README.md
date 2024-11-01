# Prescreening-Code-and-Examples
Example code for simulations and creating plots as seen in the paper


The file "Prescreening_functions.R" contains functions that will be used for the simulations.

The file "Bivariate_simulations.R" contains code to run prescreening simulations for the bivariate model. There is examples for the prosecution case with common evidence and the defense cases moving the trace towards the mean of the source level mixture. There is also an example for using the rare prosecution case when the control and trace are moved away from the other subpopulations, and the intersection point is the .99 probability contour of the component of the subpopulation of interest. The direction vector and contour can be changed to explore more cases.

The file "Making_plots.R" contains code to make the plots to visualize the simulation results using the pointwise average and standard deviation of the normal- and density-based LRs as in the paper. For convenience we have included the results from the bivariate simulations when the prosecution is true and the control and trace are rare for the subpopulation and moved away from the other subpopulation. These results are in the "Bivariate_Sim_Results_mean_away_K3" folder.

