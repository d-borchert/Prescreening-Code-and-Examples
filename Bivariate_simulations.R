# simulations for bivariate model
# Author: Dylan Borchert
# updated 11/1/2024


# set working directory to save files in correct directory
# Gerneally setting working directory to source file location works best
setwd()


# load the functions needed for the simulations
source("Prescreening_functions.R")


## generate the mxiure with parameters used in paper
set.seed(23232)
mix <- MixSim(BarOmega = 0.0001, K = 3, p = 2, hom=TRUE)



# Simulation parameters

# this controls the number of replications, B in the paper
n_iter <- 3

# this controls the number of prescreening steps, L in the paper
n_la <- 10

# This controls the number of background objects
n_a <- 50

# this controls the number of observations within each background object
n_w <- 5

# this controls the subpopulation which acts as the subpopulation of interest
subpop_id <- 1


# simulations for three clusters
# K = 3

# parameters of the model
biv_mu <- mix$Mu
biv_S_b <- mix$S
biv_S_w <- .1*mix$S[,,1]
biv_pi <- mix$Pi

# first and second mean
biv_mu1 <- mix$Mu[1,]
biv_mu2 <- mix$Mu[2,]

# this is the means of the trace for the 
# common evidence prosecution and 4 defense cases
moving_means <- get_means(means = biv_mu, subpop_id = 1, covariance = biv_S_b[,,1], pi = biv_pi, seq_number = 5)


for (i in 1:5){
  
  # for monitoring how far we are in the simualtion
  print(paste0('starting sim ', i))
  
  # generate the trace
  e_u <- e_u_dat(10, moving_means[i,], biv_S_w)
  # and control
  e_s <- e_s_dat(10, biv_mu1, biv_S_w)
  
  # specify where to save the results
  foldername_biv <- 'Bivariate_Sim_Results_K3'
  # specify the folder when there is no subpopulation structure
  simname_biv1 <- paste0('biv_sim_single_', i)
  
  # prescreening simulation
  biv_sim1 <- thresholding_simulation(e_u = e_u, 
                                      e_s = e_s,
                                      n_iter = n_iter,
                                      n_la = n_la,
                                      n_a = n_a,
                                      n_w = n_w,
                                      pi = biv_pi,
                                      mu = biv_mu,
                                      sig_b = biv_S_b,
                                      sig_w = biv_S_w,
                                      is_supopulations = FALSE,
                                      subpopid = subpop_id,
                                      foldername = foldername_biv,
                                      simname = simname_biv1)
  
  # specify the folder when there is a subpopulation structure
  simname_biv2 <- paste0('biv_sim_all_', i)
  
  # prescreening simulation
  biv_sim2 <- thresholding_simulation(e_u = e_u, 
                                      e_s = e_s,
                                      n_iter = n_iter,
                                      n_la = n_la,
                                      n_a = n_a,
                                      n_w = n_w,
                                      pi = biv_pi,
                                      mu = biv_mu,
                                      sig_b = biv_S_b,
                                      sig_w = biv_S_w,
                                      is_supopulations = TRUE,
                                      subpopid = subpop_id,
                                      foldername = foldername_biv,
                                      simname = simname_biv2)
  
  
}




# Example simulation moving the control and trace away from the mean


# get mean between the two subpoulations, on the 99% contour of the second subpop.
# this the the direction from mu1 away from mu2 where it intersects the 99% contour
source_mean_away <- move_mean((biv_mu1 - biv_mu2), biv_mu1, biv_S_b[,,1], .99)

# generate trace
e_u <- e_u_dat(10, source_mean_away, biv_S_w)
# and control
e_s <- e_s_dat(10, source_mean_away, biv_S_w)

# specify folder to save results
foldername_biv <- 'Bivariate_Sim_Results_mean_away_K3'

# simulation with no subpopulations
simname_biv3 <- 'biv_sim_single'
biv_sim1 <- thresholding_simulation(e_u = e_u, 
                                    e_s = e_s,
                                    n_iter = n_iter,
                                    n_la = n_la,
                                    n_a = n_a,
                                    n_w = n_w,
                                    pi = biv_pi,
                                    mu = biv_mu,
                                    sig_b = biv_S_b,
                                    sig_w = biv_S_w,
                                    is_supopulations = FALSE,
                                    subpopid = subpop_id,
                                    foldername = foldername_biv,
                                    simname = simname_biv3)
# simulation with subpopulations
simname_biv4 <- 'biv_sim_all'
biv_sim2 <- thresholding_simulation(e_u = e_u, 
                                    e_s = e_s,
                                    n_iter = n_iter,
                                    n_la = n_la,
                                    n_a = n_a,
                                    n_w = n_w,
                                    pi = biv_pi,
                                    mu = biv_mu,
                                    sig_b = biv_S_b,
                                    sig_w = biv_S_w,
                                    is_supopulations = TRUE,
                                    subpopid = subpop_id,
                                    foldername = foldername_biv,
                                    simname = simname_biv4)


