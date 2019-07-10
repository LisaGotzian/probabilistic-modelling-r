## TO FIGURE OUT:
# - how to calculate z's
# - how to calcualte delta matrices
# - how to calculate V matrices

setwd('C:/Users/5SQQJ72/Documents/leuphana/semester_2/probabilistic_modeling/group_project')

library(readxl)
library(tidyverse)
library(mvtnorm)
library(truncnorm)

## instantiate the values to estimate
num_keywords <- 36
num_days <- 40
num_observations <- 1440

## generate estimates for the variables used
# poisson impressions
impressions <- matrix(round(rpois(n=num_observations, lambda = rep(72.8, num_observations)), 0),
                      nrow=num_keywords, ncol=num_days)
impressions[impressions > 1666] <- 1666
impressions[impressions < 1] <- 1

# turned this to poisson
#impressions <- matrix(round(rtruncnorm(n=num_observations, a=1, b=16666, mean=72.8, sd=159),0),
#                nrow=num_keywords, ncol=num_days)

# poisson clicks
clicks <- matrix(round(rpois(n=num_observations, lambda = rep(1.1, num_observations)),0),
                 nrow=num_keywords, ncol=num_days)
clicks[clicks > 24] <- 24
clicks[clicks < 0] <- 0

#turned this to poisson
#clicks <- matrix(round(rtruncnorm(n=num_observations, a=0, b=24, mean=1.1, sd=2.2),0),
#                nrow=num_keywords, ncol=num_days)

# poission orders
orders <- matrix(round(rpois(n=num_observations, lambda = rep(0.03, num_observations)),2),
                 nrow=num_keywords, ncol=num_days)
orders[orders > 3] <- 3
orders[orders < 0] <- 0

#turned this to poisson
#orders <- matrix(round(rtruncnorm(n=num_observations, a=0, b=3, mean=0.03, sd=0.2),0),
#                nrow=num_keywords, ncol=num_days)

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
# also make this normal?
brand <- matrix(rbinom(n=num_keywords, size = 1, p=0.6),
                nrow=num_keywords, ncol=1)
# normal distribution specificity
# specificity <- matrix(rtruncnorm(n=num_keywords, a=0, b=1, mean=0.4, sd=0.7),
#                 nrow=num_keywords, ncol=1)
# binomial distribution specificity
specificity <- matrix(rbinom(n=num_keywords, size=1, p=0.4),
                      nrow=num_keywords, ncol=1)
bid <- matrix(round(rtruncnorm(n=num_observations, a=0.08, b=2, mean=0.5, sd=0.3),2),
              nrow=num_keywords, ncol=num_days)
#------------calculate the error terms-------------------#
r1 <- c(0.370,	-0.036,	0.000,	-0.005)
r2 <- c(-0.036,	0.232,	0.007,	0.02)
r3 <- c(0,	0.007,	0.064,	-0.005)
r4 <- c(-0.005,	0.02,	-0.005,	0.087)

o_mat_cols <- c('CONV','CTR','Pos','Organic_comp')
omega_matrix <- matrix(c(r1, r2, r3, r4), nrow=4, ncol=4, byrow=TRUE,
                       dimnames = list(o_mat_cols))
colnames(omega_matrix) <- o_mat_cols

## I think this needs to be changed to generate 1 n x m matrix for each error term
error_terms <- rmvnorm(n = num_observations, sigma=omega_matrix)
errors <- vector(mode="list", length=4)
names(errors) <- c('theta_error','beta_error','gamma_error','alpha_error')
errors[[1]] <- matrix(error_terms[,1], nrow=num_keywords, ncol=num_days)
errors[[2]] <- matrix(error_terms[,2], nrow=num_keywords, ncol=num_days)
errors[[3]] <- matrix(error_terms[,3], nrow=num_keywords, ncol=num_days)
errors[[4]] <- matrix(error_terms[,4], nrow=num_keywords, ncol=num_days)

#------------/calculate the error terms------------------#
#----fct to calculate heterogeneous variables------------#
#is this correct?

calc_coeff <- function(keywords, coefficient, v) {
  l <- coefficient + rnorm(keywords, 0, v)
  return(l)
}

calc_coeff_matr <- function(keywords, delta_matr, z, v_matr) {
  l = list(n=keywords)
  for (k in 1:keywords){
    u <- matrix(rmvnorm(n=1, mean=c(0,0), v_matr),nrow=2)
    d <- delta_matr %*% z[[k]] + u
    l[k] <-  list(d)
  }
  return(l)
}
#---/fct to calculate heterogeneous variables------------#
#-----instantiate variables to use in the equations------#
# this we'll have to do: compare our thetas from calc_weights to the paper (sd given)
# create k z's
z <- list()
for (k in 1:num_keywords){
  z[k] <- list(c(brand[k], specificity[k]))
}

# create V covarince matrix
V <- matrix(c(c(-0.41, 0.4),c(0.4, -0.17)), nrow=2)
# generate variables for each equation
# - CTR
theta_delta <- matrix(c(c(-0.09, 0.64),c(-0.01,0.55)), nrow=2, ncol=2, byrow=TRUE)
random_thetas <- calc_coeff_matr(num_keywords, theta_delta, z, V)
theta0 <- calc_coeff(num_keywords, -4.22, 0.55)
theta1 <- calc_coeff(num_keywords, -1.42, 0.18)
theta2 <- list()
theta3 <- list()
for (k in 1:num_keywords) {
  theta2[k] <- list(random_thetas[[k]][1,])
  theta3[k] <- list(random_thetas[[k]][2,])
}
theta4 <- -0.2; theta5 <- 0.3; theta6 <- -0.01;
theta_brand <- 0.25; theta_spec <- -0.51; theta_error <- errors$theta_error

# - CONV
beta_delta <- matrix(c(c(1.27, 0.52), c(-0.31,0.58)), nrow=2, ncol=2, byrow=TRUE)
random_betas <- calc_coeff_matr(num_keywords, beta_delta, z, V)
beta0 <- calc_coeff(num_keywords, -2.75, 0.68)
beta1 <- calc_coeff(num_keywords, 0.81, 0.28)
beta2 <- list()
beta3 <- list()
for (k in 1:num_keywords) {
  beta2[k] <- list(random_betas[[k]][1,])
  beta3[k] <- list(random_betas[[k]][2,])
}
beta4 <- 0.17; beta5 <- -0.01; beta6 <- -0.01;
beta_brand <- -0.77; beta_spec <- -0.25;
beta_error <- errors$beta_error

# - AdPos
gamma0 <- calc_coeff(num_keywords, 1.03, 0.16)
gamma1 <- calc_coeff(num_keywords, -0.45, 0.08)
gamma2 <- -0.049; gamma3 <- -0.002;
gamma_error <- errors$gamma_error

# - Organic Comp
alpha0 <- calc_coeff(num_keywords, 0.07, 0.17)
alpha1 <- calc_coeff(num_keywords, 1.86, 0.21)
alpha2 <- -0.004
alpha_error <- errors$alpha_error
#----/instantiate variables to use in the equations------#
#-----instantiate matrices to store results--------------#
u_ctr <- matrix(nrow=num_keywords, ncol=num_days)
ctr <- matrix(nrow=num_keywords, ncol=num_days)
u_conv <- matrix(nrow=num_keywords, ncol=num_days)
conv <- matrix(nrow=num_keywords, ncol=num_days)
adpos_eq <- matrix(nrow=num_keywords, ncol=num_days)
lnadpos <- matrix(nrow=num_keywords, ncol=num_days)
orgcomp_eq <- matrix(nrow=num_keywords, ncol=num_days)
#----/instantiate matrices to store results--------------#
#---------------RUN EQUATIONS!---------------------------#
for (t in 1:num_days) {
  for (k in 1:num_keywords){
    adpos_eq[k,t] <- exp(gamma0[k] + gamma1[k]*log(bid[k,t]) + gamma2*lqscore[k,t] +
                           gamma3*t + gamma_error[k,t])
    
    orgcomp_eq[k,t] <- alpha0[k] + alpha1[k]*iv_organic[k,t] +
      alpha2*t + alpha_error[k,t]
    
    u_ctr[k,t] <- theta0[k] + theta1[k]*adpos_eq[k,t] +
      theta2[[k]]*orgcomp_eq[k,t] + theta3[[k]]*sponsored_comp[k,t] +
      theta4*organic[k,t] + theta5*lqscore[k,t] +
      theta6*t + theta_error[k,t]
    ctr[k,t] <- exp(u_ctr[k,t]) / (1 + exp(u_ctr[k,t]))
    
    u_conv[k,t] <- beta0[k] + beta1[k]*adpos_eq[k,t] +
      beta2[[k]]*orgcomp_eq[k,t] + beta3[[k]]*sponsored_comp[k,t] +
      beta4*organic[k,t] + beta5*lqscore[k,t] +
      beta6*t + beta_error[k,t]
    conv[k,t] <- exp(u_conv[k,t]) / (1 + exp(u_conv[k,t]))
  }
}
#--------------/RUN EQUATIONS!--------------------------#
#--------estimate clicks and orders---------------------#
click_est <- round(impressions * ctr, 0)
#click_est[click_est > 24] <- 24

orders_est <- round(click_est * conv, 0)
#orders_est[orders_est > 3] <- 3
#-------/estimate clicks and orders---------------------#
mean(click_est)
mean(clicks)

library('rjags') 
library('runjags') # MODEL DOES NOT WORK! How to make it work? 
# do we have to source one of the files by Kruschke? this way, we could use his functions.

delta = 0.5
data_list = list(
  z = z,
  N = num_days,
  M = num_keywords,
  adpos = adpos,
  organic = organic,
  orgcomp = organic_comp,
  sponsored_comp = sponsored_comp,
  lqscore = lqscore,
  theta_error = theta_error,
  V = V,
  delta = delta
)

modelString = "
model {
for ( i in 1:NumKeywords ){
theta0[i] <- delta*z[i] + mu[i]
theta1[i] <- delta*z[i] + mu[i]
theta2[i] <- delta*z[i] + mu[i]
theta3[i] <- delta*z[i] + mu[i]

mu[i] ~ rnorm(0, V)
}

for ( t in 1:N ) {
for ( k in 1:M ){
ctr[k, t] ~ dbern( p[k, t] )
p[k, t] <- exp( z[k, t] ) / ( 1 + exp( z[k, t] ) )
z[k, t] <- theta0[k] + theta1[k]*adpos[k, t] + theta2[k]*orgcomp[k, t] +
theta3[k]*sponsored_comp[k, t] + theta4*organic[k] +
theta5*lqscore[k, t] + theta6*t + theta_error[k, t]
}
}
theta4 <- -0.2
theta5 <- 0.3
theta6 <- -0.01
}
"
writeLines( modelString, con="TEMPmodel.txt")

jagsModel = rjags::jags.model( file = "TEMPmodel.txt",
                       data = data_list,
                       #inits = initsList,
                       n.chains = 3,
                       n.adapt = 500 )

update( jagsModel, n.iter = 500 )
#
# codaSamples = coda.samples( jagsModel,
#                             variable.names = c("ctr","conv"),
#                             n.iter = 4000 )

nChains = 3
nAdaptSteps = 500
nBurninSteps = 1000

runJagsOut <- runjags::run.jags( method = 'parallel',
                                 model = 'TEMPmodel.txt',
                                 monitor = c('theta0','theta1','theta2','theta3'),
                                 data = data_list,
                                 #inits = initsList,
                                 n.chains = nChains,
                                 adapt = nAdaptSteps,
                                 burnin = nBurninSteps,
                                 summarise = FALSE,
                                 plots = FALSE )

codaSamples = as.mcmc.list( runJagsOut )
