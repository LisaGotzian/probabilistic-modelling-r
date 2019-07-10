#----------------- 0) Simulating the data ----------------
# Purpose: simulate the data and parameters in the model by
# Agarwal et al. 2015 on organic search results.

# Erik Schau, Lisa Gotzian, August 2018

library(readxl)
library(tidyverse)
library(mvtnorm)
library(truncnorm)
library(MCMCpack)

set.seed(66)

## instantiate the values to estimate
num_keywords <- 36
num_days <- 40
num_observations <- 1440

## generate estimates for the variables used
# poisson impressions
impressions <- matrix(round(rpois(n=num_observations, lambda = rep(72.8, num_observations)), 2),
                       nrow=num_keywords, ncol=num_days)
impressions[impressions > 1666] <- 1666
impressions[impressions < 1] <- 1

# poisson clicks
clicks <- matrix(round(rpois(n=num_observations, lambda = rep(1.1, num_observations)),2),
                      nrow=num_keywords, ncol=num_days)
clicks[clicks > 24] <- 24
clicks[clicks < 0] <- 0

# poission orders
orders <- matrix(round(rpois(n=num_observations, lambda = rep(0.03, num_observations)),2),
                      nrow=num_keywords, ncol=num_days)
orders[orders > 3] <- 3
orders[orders < 0] <- 0

adpos <- matrix(round(rtruncnorm(n=num_observations, a=1, b=9.78, mean=3.47, sd=1.7),2),
                nrow=num_keywords, ncol=num_days)
organic <- matrix(rbinom(n=num_observations, size = 1, p=0.15),
                nrow=num_keywords, ncol=num_days)
organic_comp <- matrix(round(rtruncnorm(n=num_observations, a=0, b=3.08, mean=0.78, sd=0.54),2),
                                nrow=num_keywords, ncol=num_days)
iv_organic <- matrix(round(rtruncnorm(n=num_observations, a=0.1, b=2.34, mean=0.42, sd=0.23),2),
                nrow=num_keywords, ncol=num_days)
sponsored_comp <- matrix(round(rtruncnorm(n=num_observations, a=0, b=5.59, mean=1.79, sd=0.89),2),
                nrow=num_keywords, ncol=num_days)
iv_sponsored <- matrix(round(rtruncnorm(n=num_observations, a=0.13, b=5.19, mean=1.49, sd=1),2),
                nrow=num_keywords, ncol=num_days)
lqscore <- matrix(round(rtruncnorm(n=num_observations, a=6, b=10, mean=8, sd=1.5),0),
                nrow=num_keywords, ncol=num_days)

brand <- matrix(rbinom(n=num_keywords, size = 1, p=0.6),
                nrow=num_keywords, ncol=1)

# turned this to a binomial distribution
specificity <- matrix(rbinom(n=num_keywords, size = 1, p=0.4),
                nrow=num_keywords, ncol=1)
bid <- matrix(round(rtruncnorm(n=num_observations, a=0.08, b=2, mean=0.5, sd=0.3),2),
                nrow=num_keywords, ncol=num_days)

# Calculating clickthrough rate based on the data
CTR <- clicks/impressions

# Calculating conversion rate based on the data
CONV <- orders/impressions

#------------calculate the error terms-------------------
r1 <- c(0.370,	-0.036,	0.000,	-0.005)
r2 <- c(-0.036,	0.232,	0.007,	0.02)
r3 <- c(0,	0.007,	0.064,	-0.005)
r4 <- c(-0.005,	0.02,	-0.005,	0.087)

o_mat_cols <- c('CONV','CTR','Pos','Organic_comp')
omega_matrix <- matrix(c(r1, r2, r3, r4), nrow=4, ncol=4, byrow=TRUE,
                       dimnames = list(o_mat_cols))
colnames(omega_matrix) <- o_mat_cols

error_terms <- rmvnorm(n = num_observations, sigma=omega_matrix)
errors <- vector(mode="list", length=4)
names(errors) <- c('theta_error','beta_error','gamma_error','alpha_error')
errors[[1]] <- matrix(error_terms[,1], nrow=num_keywords, ncol=num_days)
errors[[2]] <- matrix(error_terms[,2], nrow=num_keywords, ncol=num_days)
errors[[3]] <- matrix(error_terms[,3], nrow=num_keywords, ncol=num_days)
errors[[4]] <- matrix(error_terms[,4], nrow=num_keywords, ncol=num_days)

#------------/calculate the error terms------------------#
#----fct to calculate heterogeneous variables------------

compare_coeff_paper <- function(keywords, coefficient, v) {
  l <- coefficient + rnorm(keywords, 0, v)
  return(l)
}

calc_coeff_matr <- function(keywords, delta_matr, z, v_matr) {
  l = matrix(nrow= 4, ncol=keywords)
  for (k in 1:keywords){
    u <- matrix(rmvnorm(n=1, mean=rep(0,nrow(delta_matr)), v_matr),
                nrow=nrow(delta_matr))
    d <- delta_matr %*% z[,k] + u
    l[,k] <-  d
  }
  return(l)
}

#---/fct to calculate heterogeneous variables------------#
#-----instantiate variables to use in the equations------

# create k z's
z <- matrix(nrow=2, ncol=num_keywords)
for (k in 1:num_keywords){
  z[,k] <- c(brand[k], specificity[k])
}

# The covariance matrices: basically the error "within" a set of parameters
# 4 thetas depend on k, so VTheta is 4x4. All in all, there are 4 V matrices,
# one for each parameter set, two of size 4x4, two of size 2x2.
# As there are no estimates given in the paper, we sample from the IW distribution,
# as stated in the appendix, with degrees of freedom equal to the columns of
# the matrix and the identity matrix.

VTheta <- riwish(4, diag(1,nrow=4))
VBeta <- riwish(4, diag(1,nrow=4))
VGamma <- riwish(2, diag(1,nrow=2))
VAlpha <- riwish(2, diag(1,nrow=2))


# The theta_delta matrix is given incompletely, all missing values are treated as 0,
# see table 4 in the main paper. It is all parameters depending on k in interaction
# with a keyword's brand and specificity, giving a 4x2 (or 2x2) matrix.
theta_delta <- matrix(c(c(0.25, -0.51), # const x brand and spec
                        c(0,0), # adpos x brand and spec
                        c(-0.09, 0.64), # organic_comp x brand and spec
                        c(-0.01,0.55)), # sponsored_comp x brand and spec
                      nrow=4, ncol=2, byrow=TRUE)

beta_delta <- matrix(c(c(-0.77, -0.25), # const x brand and spec
                        c(0,0), # adpos x brand and spec
                        c(1.27, 0.52), # organic_comp x brand and spec
                        c(-0.31,0.58)), # sponsored_comp x brand and spec
                      nrow=4, ncol=2, byrow=TRUE)

# The estimates for the gamma or alpha matrix are not given at all. As they
# are normally distributed (see appendix), we draw them with a mean of 0
# and a standard sd of 0.5^2.
gamma_delta <- matrix(rnorm(4, 0, 0.25), # const and bid x brand and spec
                     nrow=2, ncol=2, byrow=TRUE)

alpha_delta <- matrix(rnorm(4, 0, 0.25), # const and iv_organic x brand and spec
                      nrow=2, ncol=2, byrow=TRUE)


# The thetas for CTR
random_thetas <- calc_coeff_matr(num_keywords, theta_delta, z, VTheta)

random_thetas_paper <- matrix(c(compare_coeff_paper(num_keywords, -4.22, 0.55), #const
                         compare_coeff_paper(num_keywords, -1.42, 0.18), # adpos
                         compare_coeff_paper(num_keywords, -0.41, 0.15), #organiccomp
                         compare_coeff_paper(num_keywords, -0.17, 0.12)), #sponsoredcomp
                         nrow = 4, ncol = 36, byrow = TRUE)

# Took the mean values for the thetas that do not depend on k.
theta4 <- -0.2; theta5 <- 0.3; theta6 <- -0.01;
theta_error <- errors$theta_error


# The betas for CONV
random_betas <- calc_coeff_matr(num_keywords, beta_delta, z, VBeta)

random_betas_paper <- matrix(c(compare_coeff_paper(num_keywords, -2.75, 0.68), #const
                                compare_coeff_paper(num_keywords, 0.81, 0.28), # adpos
                                compare_coeff_paper(num_keywords, 0.7, 0.18), #organiccomp
                                compare_coeff_paper(num_keywords, -0.27, 0.17)), #sponsoredcomp
                              nrow = 4, ncol = 36, byrow = TRUE)

beta4 <- 0.17; beta5 <- -0.01; beta6 <- -0.01;
beta_error <- errors$beta_error


# The gammas for AdPos
random_gammas <- calc_coeff_matr(num_keywords, gamma_delta, z, VGamma)

random_gammas_paper <- matrix(c(compare_coeff_paper(num_keywords, 1.03, 0.16), #const
                               compare_coeff_paper(num_keywords, -0.45, 0.08)), # log(bid)
                             nrow = 2, ncol = 36, byrow = TRUE)
gamma2 <- -0.049; gamma3 <- -0.002;
gamma_error <- errors$gamma_error

# The alphas for organic_comp
random_alphas <- calc_coeff_matr(num_keywords, alpha_delta, z, VAlpha)

random_gammas_paper <- matrix(c(compare_coeff_paper(num_keywords, 0.07, 0.17), #const
                                compare_coeff_paper(num_keywords, 1.86, 0.21)), # log(bid)
                              nrow = 2, ncol = 36, byrow = TRUE)
alpha2 <- -0.004
alpha_error <- errors$alpha_error

#-----\instantiate variables to use in the equations------#
#-----instantiate matrices to store results--------------
u_ctr <- matrix(nrow=num_keywords, ncol=num_days)
CTR_estimate <- matrix(nrow=num_keywords, ncol=num_days)
u_conv <- matrix(nrow=num_keywords, ncol=num_days)
CONV_estimate <- matrix(nrow=num_keywords, ncol=num_days)
adpos_eq <- matrix(nrow=num_keywords, ncol=num_days)
#lnadpos <- matrix(nrow=num_keywords, ncol=num_days)
orgcomp_eq <- matrix(nrow=num_keywords, ncol=num_days)
#-----\instantiate matrices to store results--------------#
#---------------RUN EQUATIONS!---------------------------
for (t in 1:num_days) {
  for (k in 1:num_keywords){ # Attention, due to the indexing, random_gamma[[k]][1]
    # corresponds to gamma0
    adpos_eq[k,t] <- exp(random_gammas[1,k] + random_gammas[2,k]*log(bid[k,t]) + gamma2*lqscore[k,t] +
      gamma3*t + gamma_error[k,t])
    #lnadpos[k,t] <- log(adpos_eq[k,t])

    orgcomp_eq[k,t] <- random_alphas[1,k] + random_alphas[2,k]*iv_organic[k,t] +
      alpha2*t + alpha_error[k,t]

    u_ctr[k,t] <- random_thetas[1,k] + random_thetas[2,k]*adpos_eq[k,t] +
      random_thetas[3,k]*orgcomp_eq[k,t] + random_thetas[4,k]*sponsored_comp[k,t] +
      theta4*organic[k,t] + theta5*lqscore[k,t] +
      theta6*t + theta_error[k,t]
    CTR_estimate[k,t] <- exp(u_ctr[k,t]) / (1 + exp(u_ctr[k,t]))

    u_conv[k,t] <- random_betas[1,k] + random_betas[2,k]*adpos_eq[k,t] +
      random_betas[3,k]*orgcomp_eq[k,t] + random_betas[4,k]*sponsored_comp[k,t] +
      beta4*organic[k,t] + beta5*lqscore[k,t] +
      beta6*t + beta_error[k,t]
    CONV_estimate[k,t] <- exp(u_conv[k,t]) / (1 + exp(u_conv[k,t]))
  }
}
#--------------/RUN EQUATIONS!--------------------------#
#---------------------- Plots -----------------------------------
# Compare with values from the paper
# CTR_estimate and CONV_estimate # vs. CTR and CONV based on the data
# theta etc. with what we have
plot(CTR)
plot(CTR_estimate)
plot(CONV)
plot(CONV_estimate)

