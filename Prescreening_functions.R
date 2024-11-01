##########################################################################
# Author: Dylan Borchert
# This file contains functions used in the prescreening simulations and
# exploratory examples
# updated 11/1/2024
##########################################################################

# required packaged
library(MASS)
library(mvtnorm)
library(comparison)
library(MixSim)



############################################################################
# Functions for moving e_u
############################################################################

# function to move the mean of the trace for simualtion
# means: matrix of means
# subpop_id: which subpopulation is the subpopulation of interest
# covariance: coviarance matrix
# pi: mixing proportion of the mixture model
# seq_number: how many means you want in the sequence 
get_means <- function(means, subpop_id, covariance, pi, seq_number){
  
  # dimension
  p = ncol(means)
  
  # center of the means
  center = pi %*% means
  
  # direction from current mean to center
  direction_vec = center - means[subpop_id,]
  
  # norm squred of the direction vector
  dir_vec_norm_sq = sum(direction_vec^2)
  
  # chi squared quantile for level
  c_sq = qchisq(.33, p)
  
  # numerator and denominator for scale length
  num = c_sq * dir_vec_norm_sq
  denom = direction_vec %*% solve(covariance) %*% t(direction_vec)
  
  # the scale factor
  max_scale = sqrt(num/denom)
  
  scale_seq = seq(from = 0, to = as.numeric(max_scale), length.out = seq_number)
  
  points_mat = matrix(NA, nrow = seq_number, ncol = p)
  # point in the direction of direction vector
  # scale away from the mean
  for (i in 1:seq_number){
    points_mat[i,] <- means[subpop_id,] + (as.numeric(scale_seq[i]/sqrt(dir_vec_norm_sq))*direction_vec)
  }
  return(points_mat)
}





# function to find location of a mean using a direction vector and a covarince matrix
# for the contour ellipse
# parameters:
# direction_vector: vector that has the direction from the population mean
# mean: location to move from
# covariance: covariance matrix
# level: contour level the point should lie on.
move_mean <- function(direction_vec, mean, covariance, level){
  
  # dimension
  p <- length(direction_vec)
  
  # norm squred of the direction vector
  dir_vec_norm_sq <- sum(direction_vec^2)
  
  # chi squared quantile for level
  c_sq <- qchisq(level, p)
  
  # numerator and denominator for scale length
  num <- c_sq * dir_vec_norm_sq
  denom <- t(direction_vec) %*% solve(covariance) %*% direction_vec
  
  # the scale factor
  scale <- sqrt(num/denom)
  
  # point in the direction of direction vector
  # scale away from the mean
  point <- mean + (as.numeric(scale/sqrt(dir_vec_norm_sq))*direction_vec)
  return(point)
}




############################################################################
# Functions for use in simulations
############################################################################

############################################################################
# Function to simulate a background population of hierarchical normal sample
# parameters:
# n_a: number of sources in the background population
# n_w: number of samples withing each source
# theta: mean of the source distribution
# sig_b: between source covariance
# sig_w : within source cavariance
############################################################################
# output a data frame where the first column is titled 'obj' and the next 
# p columns are the measurements in p dimensions
############################################################################
sim_dat_fun <- function(n_a, n_w, mu, sig_b, sig_w){
  # initialize storer
  dat_tmp = NULL
  # n_a times
  for(obj in 1:n_a){
    # sample a source mean
    obj_mu = rmvnorm(1, mu, sig_b)
    # draw samples from the source
    obj_samps = rmvnorm(n_w, obj_mu, sig_w)
    # store in a data frame
    dat_tmp = data.frame(rbind(dat_tmp, cbind(obj, obj_samps)))
  }
  return(dat_tmp)
}






################################################################################
# Function to simulate a background population with a mixture of normals
# between sources and normal within source
# parameters:
# n_a: number of sources in the background population
# n_w: number of samples withing each source
# pi: vector of mixing proportions
# mu: matrix of means for between source mxiture model with
# rows corresponding to mean of a cluster
# S: Array of covariance matrices for between source
# sig_w : within source cavariance
################################################################################
# data frame where the first column is titled 'obj'
# the next p columns are measurements
# the last column is 'id' which is the subpopulation id
################################################################################
sim_dat_fun2 <- function(n_a, n_w, pi, mu, S, sig_w){
  # sample n_a sources from the between source distribution
  mixture = simdataset(n_a, pi, mu, S)
  means = mixture$X
  ids = mixture$id
  dat_tmp <- NULL
  for(obj in 1:n_a){
    # draw samples from each source
    obj_samps = rmvnorm(n_w, means[obj,], sig_w)
    dat_tmp = data.frame(rbind(dat_tmp, cbind(obj, obj_samps, id = ids[obj])))
  }
  return(dat_tmp)
}


##############################################
# Function to generate the trace samples
# parameters:
# n_u: number of samples
# mu_u: mean 
# sig_u: covariance
##############################################
# data frame where the first column is obj
# the next p columns are measurements
##############################################
e_u_dat <- function(n_u, mu_u, sig_u){
  dat_u <- rmvnorm(n_u, mu_u, sig_u)
  return(data.frame(obj='E_u', dat_u))
}


##############################################
# Function to generate the control samples
# parameters:
# n_s: number of samples
# mu_s: mean 
# sig_s: covariance
##############################################
# data frame where the first column is obj
# the next p columns are measurements
##############################################
e_s_dat <- function(n_s, mu_s, sig_s){
  dat_s <- rmvnorm(n_s, mu_s, sig_s)
  return(data.frame(obj='E_s', dat_s))
}





# Function to do a chi_squared test of equality of means in normal samples
# parameters
# x,y: samples to compare mean
# S: common covariance matrix
chi_squared_stat <- function(x, y, S){
  # number of samples and dimensions of x and y
  nx = nrow(x)
  ny = nrow(y)
  px = ncol(x)
  py = ncol(y)
  
  # check the dimension of x and y
  if (px != py)
    stop("Both samples must have the same number of variables (columns)")
  
  # set dimension
  p = px
  
  # calculate the sample mean vectors
  mx = as.matrix(apply(x, 2, mean), ncol = 1)
  my = as.matrix(apply(y, 2, mean), ncol = 1)
  
  # calculate the test stat
  test_stat = ((nx*ny)/(nx+ny))*t(mx - my)%*%solve(S)%*%(mx - my)
  
  # under the null (equality of means) the test stat 
  # follows a Chi Squared distribution with p df.
  p_val = pchisq(test_stat, p, lower.tail = FALSE)
  
  #return the p value and test_stat
  return(list(pval = p_val, teststat = test_stat))
}



# Function to get the within covariance matrix from a collection
# of samples from different sources
# parameters:
# e_a: data set (need a variable titled `obj`)
# dat_cols: vector containing the columns with the measurements
get_Sw <- function(e_a, dat_cols){
  uni_obj = unique(e_a$obj)
  S_w_sum = matrix(0, nrow=length(dat_cols), ncol=length(dat_cols))
  for (i in uni_obj){
    tmp_dat = e_a[e_a$obj==i,]
    tmp_cov = cov(tmp_dat[,dat_cols])
    S_w_sum = S_w_sum + tmp_cov
  }
  
  S_w = S_w_sum / length(uni_obj)
  
  return(S_w)
}





# Function to return test stat comparing each source against e_u
# parameters:
# e_a: data (need variable titles `obj`)
# e_u: sample to test against
# dat_cols: columns with measurements (should be same for e_a and e_u) 
condit_test <- function(e_a, e_u, dat_cols){
  
  # get the pooled covariance from the given background
  # population
  tmp_Sw <- get_Sw(e_a, dat_cols)
  
  test_stat_stor <- NULL
  uni_obj <- unique(e_a$obj)
  for (i in uni_obj){
    # current object
    tmp_dat <- e_a[e_a$obj==i,]
    # test stat from chi squared test multiplied by negative one
    test_stat <- chi_squared_stat(e_u[,dat_cols], tmp_dat[,dat_cols], tmp_Sw)$teststat
    test_stat_stor <- c(test_stat_stor, test_stat)
  }
  return(cbind(uni_obj, test_stat_stor))
}





#####################################################################
# two.level.normal.LR() from comparison package but edited to return
# the numerator and denominator of the LR
####################################################################

two.level.normal.LR.edited <- function (control, recovered, background) 
{
  U = background$v.within
  C = background$v.between
  mu = background$overall.means
  cont.means = control$item.means
  n.cont = control$n.replicates
  rec.means = recovered$item.means
  n.rec = recovered$n.replicates
  y.star = ((n.cont * cont.means) + (n.rec * rec.means))/(n.cont + 
                                                            n.rec)
  diff.cont.rec = cont.means - rec.means
  diff.cont.mu = cont.means - mu
  diff.rec.mu = rec.means - mu
  diff.y.star.mu = y.star - mu
  nom1 = exp(-1/2 * t(diff.cont.rec) %*% solve(U/n.cont + U/n.rec) %*% 
               (diff.cont.rec)) * (det(U/n.cont + U/n.rec))^(-1/2)
  nom2 = exp(-1/2 * t(diff.y.star.mu) %*% solve(U/(n.cont + 
                                                     n.rec) + C) %*% (diff.y.star.mu)) * (det(U/(n.cont + 
                                                                                                   n.rec) + C))^(-1/2)
  nom = nom1 * nom2
  denom1 = exp(-1/2 * t(diff.cont.mu) %*% solve(U/n.cont + 
                                                  C) %*% (diff.cont.mu)) * (det(U/n.cont + C))^(-1/2)
  denom2 = exp(-1/2 * t(diff.rec.mu) %*% solve(U/n.rec + C) %*% 
                 (diff.rec.mu)) * (det(U/n.rec + C))^(-1/2)
  denom = denom1 * denom2
  LR = as.numeric(nom/denom)
  return(list(LR = LR, nom = nom, denom = denom))
}






#function for thresholding simulation with hierachical normal structure
# parameters:
# e_u: data set of recovered evidence
# e_S: data set of control evidence
# n_iter: number of simulations
# n_la: numer of cutoffs to use for prescreening
# n_a: number of sources in the background population
# n_w: number of samples within each source
# pi: mixing proportions
# mu: matrix where ith row is the mean of the ith subpopulation
# sig_b: between source covariance matrix array
# sig_w: within source covariance matrix
# is_subpopulations: True or False indicating if the background population should
# be generated from multiple subpopulations or not
# subpopid: integer indicating which subpopulation the evidence was generated from
# foldername: name of folder
# simname: name of folder to store sim results in
thresholding_simulation <- function(e_u, e_s, n_iter, n_la, n_a, n_w, pi, mu, sig_b, 
                                    sig_w, is_supopulations = FALSE, subpopid, foldername, simname){
  
  if(dir.exists(foldername) == FALSE){
    dir.create(foldername)
  }
  
  filename = paste0(foldername, '/', simname, '/')
  
  
  #create a folder
  if(dir.exists(paste0(foldername, '/', simname)) == FALSE){
    dir.create(paste0(foldername, '/', simname))
  }
  # create files to store results
  file.create(paste0(filename, 'NumberCut.txt'))
  file.create(paste0(filename, 'props.txt'))
  file.create(paste0(filename, 'NormalLRs.txt'))
  file.create(paste0(filename, 'DensityLRs.txt'))
  file.create(paste0(filename, 'NormalNumerators.txt'))
  file.create(paste0(filename, 'NormalDenominators.txt'))
  file.create(paste0(filename, 'TotalWithin.txt'))
  file.create(paste0(filename, 'TotalBetween.txt'))
  
  write.table(rbind(e_u, e_s), file = paste0(filename, 'Evidence.txt'))
  
  param_df <- data.frame(n_iter = n_iter,
                         n_cuts = n_la,
                         source_num = n_a,
                         within_num = n_w)
  
  write.table(param_df, file = paste0(filename, 'SimParameters.txt'))
  
  # dimensionality
  p = nrow(sig_w)
  dat_cols = 2:(p+1)
  
  NLR_mat = NULL
  DLR_mat = NULL
  num_mat = NULL
  denom_mat = NULL
  number_cut_mat = NULL
  total_within_mat = NULL
  total_between_mat = NULL
  prop_mat = NULL
  # loop for simulation
  for (i in 1:n_iter){
    
    # generate the background population
    if (is_supopulations == FALSE){
      
      # generate e_a from just the subpopulation of interest
      e_a = sim_dat_fun(n_a, n_w, mu[subpopid,], sig_b[,,subpopid], sig_w)
      
    } else{
      
      # generate e_a from all subpopulations
      e_a = sim_dat_fun2(n_a, n_w, pi, mu, sig_b, sig_w)
      
    }
    
    write.table(e_a, file = paste0(filename, 'Background.txt'))
    
    # get stat comparing e_s to e_u
    S_w_hat = get_Sw(e_a, dat_cols)
    given_evidence_stat = chi_squared_stat(e_u[,dat_cols], e_s[,dat_cols], S_w_hat)
    
    # get test stat for all background objects against e_u
    results = condit_test(e_a, e_u, dat_cols)
    
    # the objects and their test stat values
    e_a_obj = results[,1]
    e_a_stats = results[,2]
    
    # the stats that are greater than the given evidence stat
    e_a_stats_considered = e_a_stats[e_a_stats > given_evidence_stat]
    e_a_obj_considered = e_a_obj[e_a_stats > given_evidence_stat]
    
    # order e_a based on the stat values
    order_of_stats = order(e_a_stats_considered, decreasing = TRUE)
    order_objs = e_a_obj_considered[order_of_stats]
    order_stats = e_a_stats_considered[order_of_stats]
    
    # sequence of number of objects to exclude from e_a
    num_reject_sec = ceiling(seq(from = 1, to = length(order_objs), length.out = n_la))
    
    
    #storers
    NLR_stor = rep(NA, n_la)
    DLR_stor = rep(NA, n_la)
    num_stor = rep(NA, n_la)
    denom_stor = rep(NA, n_la)
    number_cut_stor = rep(NA, n_la)
    total_between_stor = rep(NA, n_la)
    total_within_stor = rep(NA, n_la)
    prop_stor = rep(NA, n_la)
    # loop through this sequence
    for (j in num_reject_sec){
      
      ind_j = which(num_reject_sec == j)
      
      # remove objects from the background population
      objs_to_reject = order_objs[1:j]
      e_a_subset = e_a[!e_a$obj %in% objs_to_reject,]
      
      prop_sub = sum(e_a_subset$id == subpopid) / nrow(e_a_subset)
      
      if (n_a - j < p + 1){
        break
      }
      
      # calcualte desntiy based and normal-based LR
      e_u_format = two.level.comparison.items(e_u, dat_cols)
      e_s_format = two.level.comparison.items(e_s, dat_cols)
      Z = two.level.components(e_a_subset, dat_cols, 1)
      
      if (det(Z$v.within/(2*n_w) + Z$v.between) < 1e-40| 
          is.nan(det(Z$v.between)) == TRUE |
          prod(eigen(Z$v.between)$values) < 0){
        break
      } 
      
      # calculate LRs
      NLR_output = two.level.normal.LR.edited(e_s_format, e_u_format, Z)
      NLR = NLR_output$LR
      num = NLR_output$nom
      denom = NLR_output$denom
      
      DLR = two.level.density.LR(e_s_format, e_u_format, Z)
      
      NLR_stor[ind_j] = NLR
      DLR_stor[ind_j] = DLR
      num_stor[ind_j] = num
      denom_stor[ind_j] = denom
      number_cut_stor[ind_j] = j
      total_between_stor[ind_j] = det(Z$v.between)
      total_within_stor[ind_j] = det(Z$v.within)
      prop_stor[ind_j] = prop_sub
      
      
    }
    
    print(i)
    
    NLR_mat = rbind(NLR_mat, NLR_stor)
    DLR_mat = rbind(DLR_mat, DLR_stor)
    num_mat = rbind(num_mat, num_stor)
    denom_mat = rbind(denom_mat, denom_stor)
    number_cut_mat = rbind(number_cut_mat, number_cut_stor)
    total_within_mat = rbind(total_within_mat, total_within_stor)
    total_between_mat = rbind(total_between_mat, total_between_stor)
    prop_mat = rbind(prop_mat, prop_stor)
    
    
    #write to files
    ncols <- length(prop_stor)
    write(number_cut_stor, file = paste0(filename, 'NumberCut.txt'), append = TRUE, ncolumns = ncols)
    write(prop_stor, file = paste0(filename, 'props.txt'), append = TRUE, ncolumns = ncols)
    write(NLR_stor, file = paste0(filename, 'NormalLRs.txt'), append = TRUE, ncolumns = ncols)
    write(DLR_stor, file = paste0(filename, 'DensityLRs.txt'), append = TRUE, ncolumns = ncols)
    write(num_stor, file = paste0(filename, 'NormalNumerators.txt'), append = TRUE, ncolumns = ncols)
    write(denom_stor, file = paste0(filename, 'NormalDenominators.txt'), append = TRUE, ncolumns = ncols)
    write(total_within_stor, file = paste0(filename, 'TotalWithin.txt'), append = TRUE, ncolumns = ncols)
    write(total_between_stor, file = paste0(filename, 'TotalBetween.txt'), append = TRUE, ncolumns = ncols)
    
  }
  return(list(NLR = NLR_mat,
              DLR = DLR_mat,
              num = num_mat,
              denom = denom_mat,
              n_c = number_cut_mat,
              tw = total_within_mat,
              tb = total_between_mat,
              prop = prop_mat))
}






