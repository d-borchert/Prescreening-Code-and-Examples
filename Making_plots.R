##########################################################################
# Author: Dylan Borchert
# updated: 11/1/2024
# This file contains code to create plots of final results
##########################################################################

library(latex2exp)
library(mvtnorm)

# set wd to read files
setwd()



###################
# Functions
###################

# calcualte the true LR in the subpopulation of interest
True_subpop_LR <- function(ybar_u, ybar_s, mu, S_a, S_e, n){
  
  # dimensionality
  p = length(mu)
  
  # number of trace and control samples
  n_u = n
  n_s = n
  
  # block matrix for prosecution model
  big_S_p = rbind(cbind(S_a + S_e/n_u, S_a), cbind(S_a, S_a + S_e/n_s))
  
  # block matrix for defense model
  big_S_d = rbind(cbind(S_a + S_e/n_u, matrix(0, nrow = p, ncol = p)),
                  cbind(matrix(0, nrow = p, ncol = p), S_a + S_e/n_s))
  
  # the numerator and denominator of LR
  num = mvtnorm::dmvnorm(c(ybar_u, ybar_s), c(mu, mu), big_S_p)
  denom = mvtnorm::dmvnorm(c(ybar_u, ybar_s), c(mu, mu), big_S_d)
  
  # calculate LR
  LR = num/denom

  return(list(LR = LR, num = num, denom = denom))
}


# plotting LR means and shading for error
prescreening_plot <- function(folder, simname, legend_location, params, title, skew = FALSE){
  
  # plotting nargins
  par(mar=c(4,4,1,0.1))
  
  # what was the number of background samples used
  n_a <- 100
  
  # read files for normal- and density-based LRs
  NLR <- log10(as.matrix(read.delim(paste0(folder, '/', simname, '/NormalLRs.txt'), header = FALSE, sep = ' ')))
  DLR <- log10(as.matrix(read.delim(paste0(folder, '/', simname, '/DensityLRs.txt'), header = FALSE, sep = ' ')))
  
  # read file for number of cuts at each position
  cuts <- as.matrix(read.delim(paste0(folder, '/', simname, '/NumberCut.txt'), header = FALSE, sep = ' '))
  
  # find number of non NA columns (NAs in the files are due to not being able to calculate LR due to too
  # few backgorund observations at that step)
  # basically we are getting rid of runs where we were not complete
  DLR_inds <- colSums(is.na(DLR)) == 0
  NLR <- NLR[1:200, colSums(is.na(NLR)) == 0]
  DLR <- DLR[1:200, DLR_inds]
  cuts <- cuts[1:200, DLR_inds]
  
  # how many points to calculate the pointwise means and sds
  num_cuts <- ncol(cuts)
  
  # initialize storers
  mean_NLR <- NULL
  mean_DLR <- NULL
  sd_NLR <- NULL
  sd_DLR <- NULL
  
  # loop through points to calcualte pointwise values
  for(j in 1:num_cuts){
    
    # Current values
    tmp_NLR <- NLR[,j]
    tmp_DLR <- DLR[,j]
    
    # pointwise means
    mean_NLR <- c(mean_NLR, mean(tmp_NLR))
    mean_DLR <- c(mean_DLR, mean(tmp_DLR))
    
    # pointwise sds
    sd_NLR <- c(sd_NLR, sd(tmp_NLR))
    sd_DLR <- c(sd_DLR, sd(tmp_DLR))
  }
  
  # do we have any inifity values
  num_finite <- min(sum(is.infinite(mean_NLR) == 0), sum(is.infinite(mean_DLR) == 0))
  
  # get rid of them if we do
  mean_NLR <- mean_NLR[1:num_finite]
  mean_DLR <- mean_DLR[1:num_finite]
  sd_NLR <- sd_NLR[1:num_finite]
  sd_DLR <- sd_DLR[1:num_finite]
  cut_seq <- cuts[1,1:num_finite]
  
  # convert raw counts to propotion
  cut_seq <- (n_a - cut_seq)/n_a
  
  # make empty plotting region
  plot(cut_seq, mean_NLR + 3*sd_NLR, type = 'n', xlim = c(1,0), xlab = TeX('$p$'), ylab = 'LLR',
       ylim = c(min(c(mean_NLR - 3*sd_NLR, rev(mean_NLR + 3*sd_NLR), mean_DLR - 3*sd_DLR, rev(mean_DLR + 3*sd_DLR))), 
                max(c(mean_NLR - 3*sd_NLR, rev(mean_NLR + 3*sd_NLR), mean_DLR - 3*sd_DLR, rev(mean_DLR + 3*sd_DLR)))))
  
  # plot the title if given
  mtext(title, 3)
  
  # use this if you want to shade
  # polygon(c(cut_seq, rev(cut_seq)), 
  #         c(mean_NLR - 3*sd_NLR, rev(mean_NLR + 3*sd_NLR)), 
  #         col = rgb(0,0,1,0.2), border = FALSE)
  
  # use this if you want to use dashed lines over shading
  lines(cut_seq, mean_NLR - 2*sd_NLR, col = 'blue', lty = 2, lwd = 2)
  lines(cut_seq, mean_NLR + 2*sd_NLR, col = 'blue', lty = 2, lwd = 2)
  
  # plot the pointwise average mean line or normal LR
  lines(cut_seq, mean_NLR, type='l', col ='blue', lwd = 2)
  
  
  # point wise average mean for density LR
  lines(cut_seq, mean_DLR, lty = 1, col = 'red', lwd = 2)
  
  # use this if you want to shade
  # polygon(c(cut_seq, rev(cut_seq)), 
  #         c(mean_DLR - 3*sd_DLR, rev(mean_DLR + 3*sd_DLR)), 
  #         col = rgb(1,0,0,0.2), border = FALSE)
  
  # use this if you want dashed lines
  lines(cut_seq, mean_DLR - 2*sd_DLR, col = 'red', lty = 2, lwd = 2)
  lines(cut_seq, mean_DLR + 2*sd_DLR, col = 'red', lty = 2, lwd = 2)
  
  
  
  # read in the files with trace and control samples
  evidence <- read.delim(paste0(folder, '/', simname, '/Evidence.txt'), sep = '', header=TRUE)
  # numbers of samples in trace and control samples
  n <- sum(evidence$obj == 'E_u')
  
  # sample mean
  ybar1 <- colMeans(evidence[evidence$obj == 'E_u', -1 ])
  ybar2 <- colMeans(evidence[evidence$obj == 'E_s', -1 ])
  
  # calculate the target LR
  target <- True_subpop_LR(ybar1, ybar2, params$mu, params$Sigma_a, params$Sigma_epsilon, n)
  # add in the corresponding target LR
  abline(h = log10(target$LR))

  # make the legend
  legend(legend_location, legend = c('MVN', 'KDE', 'Target'), 
         col = c('blue', 'red', 'black'), lty = c(1,1), lwd = 2)
  
}




#########################################
# Making plots
#########################################

# bivariate K = 3

# generate mixture used in the paper 
set.seed(23232)
mix <- MixSim(BarOmega = 0.0001, K = 3, p = 2, hom=TRUE)


# parameters
biv_mu <- mix$Mu
biv_S_b <- mix$S
biv_S_w <- .1*mix$S[,,1]
biv_pi <- mix$Pi

# put into the list needed
params <- list(mu = biv_mu[1,], Sigma_a = biv_S_b[,,1], Sigma_epsilon = biv_S_w)




# plot for subpopulation case
prescreening_plot('Bivariate_Sim_Results_mean_away_K3', 
                  'biv_sim_all',
                  'topleft',
                  params,
                  '')

# plot for no subpopulation case
prescreening_plot('Bivariate_Sim_Results_mean_away_K3', 
                  'biv_sim_single',
                  'topleft',
                  params,
                  '')







