#-------------------- 1) JAGS attempt -------------------------
# Purpose:
# 1) simulate the data and parameters in the model by
#   Agarwal et al. 2015 on organic search results.
# 2) Run the model in JAGS. As it is recursive, JAGS doesn't work.
#   Read more in probabilistic_modelling_agarwal.pdf.

# Lisa Gotzian, Erik Schau, August 2018

#----------------- Simulating the data ----------------

rm(list=ls())  # Careful! This clears all of R's memory!

set.seed(66)
library(R2OpenBUGS) #for BUGS
library('runjags') #for JAGS

source("DBDA2E-utilities.R") #for JAGS
source("data_simulation.R")  # the simulated data
num_days = 40
num_keywords = 36

# yields a 2x2, so only brand and spec. Apparently, it's the same for every
# beta0, beta1 etc.
deltaq <- solve(solve(z %*% t(z)) + diag(0.1, nrow = 2, ncol = 2))

data_list = list(
  z = z,
  N = num_days,
  M = num_keywords,
  #AdPos = AdPos,
  organic = organic,
  #organicComp = organicComp,
  sponsored_comp = sponsored_comp,
  lqscore = lqscore,
  bid = bid,
  iv_organic = iv_organic,
  CTR = CTR,
  CONV = CONV,
  diag4 = diag(1,nrow=4), # as JAGS doesn't know the function diag()
  diag2 = diag(1,nrow=2),
  zero4 = rep(0,4),
  zero3 = rep(0,3),
  zero2 = rep(0,2),
  threesigma = 100 * diag(1,nrow=3),  
  twosigma = 100 * diag(1,nrow=2),
  deltaq = deltaq)




modelString <-"model{
  for ( t in 1:N ) {
    for ( k in 1:M ){ #likelihood
      CTR[k,t] ~ dbern(pCTR[k,t]) 
      CONV[k,t] ~ dbern(pCONV[k,t])
    }}
  for ( t in 1:N ) {
      for ( k in 1:M ){

      pCTR[k,t] <- exp(UCtr[k,t])/(1+exp(UCtr[k,t]))
      pCONV[k,t] <- exp(UConv[k,t])/(1+exp(UConv[k,t]))

      UCtr[k,t] <- thetak[1,k] + thetak[2,k] * AdPos[k,t] + thetak[3,k] * organicComp[k,t] +
                  thetak[4,k] * sponsored_comp[k,t] + theta[1] * organic[k,t] +
                  theta[2] * lqscore[k,t] + theta[3] * t + epsilon[2,k,t]

      UConv[k,t] <- betak[1,k] + betak[2,k] * AdPos[k,t] + betak[3,k] * organicComp[k,t] +
                  betak[4,k] * sponsored_comp[k,t] + beta[1] * organic[k,t] +
                  beta[2] * lqscore[k,t] + beta[3] * t + epsilon[1,k,t]

      log(AdPos[k,t]) <- gammak[1,k] + gammak[2,k] * log(bid[k,t]) +
                        gamma[1] * lqscore[k,t] + gamma[2] * t + epsilon[3,k,t]

      log(organicComp[k,t]) <- alphak[1,k] + alphak[2,k] * iv_organic[k,t] +
                              alpha[1] * t + epsilon[4,k,t]

      
      epsilon[1:4,k,t] ~ dmnorm(zero4[], invOmega[,])
    }
  }

  for ( k in 1:M){
      thetak[1:4,k] <- deltaTheta %*% z[1:2,k] + uTheta[1:4,k]
      betak[1:4,k] <- deltaBeta %*% z[1:2,k] + uBeta[1:4,k]
      gammak[1:2,k] <- deltaGamma %*% z[1:2,k] + uGamma[1:2,k]
      alphak[1:2,k] <- deltaAlpha %*% z[1:2,k] + uAlpha[1:2,k]

      uTheta[1:4,k] ~ dmnorm(zero4[], invVtheta[,])
      uBeta[1:4,k] ~ dmnorm(zero4[], invVbeta[,])
      uGamma[1:2,k] ~ dmnorm(zero2[], invVgamma[,])
      uAlpha[1:2,k] ~ dmnorm(zero2[], invValpha[,])
  }


    # Disregard deltaBar as it is 0 (or is that just the Initial value?)
    deltaTildeBeta <- deltaq %*% (z %*% t(betak))
    deltaTildeTheta <- deltaq %*% (z %*% t(thetak))
    deltaTildeGamma <- deltaq %*% (z %*% t(gammak))
    deltaTildeAlpha <- deltaq %*% (z %*% t(alphak))

    for (i in 1:4){ #not done
      deltaBeta[i,1:2] ~ dmnorm(deltaTildeBeta[1:2,i],  deltaq[,])
      deltaTheta[i,1:2] ~ dmnorm(deltaTildeTheta[1:2,i],  deltaq[,])
    }
    for (i in 1:2){
      deltaGamma[i,1:2] ~ dmnorm(deltaTildeGamma[1:2,i],  deltaq[,])
      deltaAlpha[i,1:2] ~ dmnorm(deltaTildeAlpha[1:2,i],  deltaq[,])
    }

    theta[1:3] ~ dmnorm(zero3[], threesigma[,])
    beta[1:3] ~ dmnorm(zero3[], threesigma[,])
    gamma[1:2] ~ dmnorm(zero2[], twosigma[,])
    alpha ~ dnorm(0, 100)

    # For the covariance matrices, I used a weakly informed prior
    # with degrees of freedom correpsonding to the dimension of the 
    # outcome instead of recreating the appendix which is prone to errors
    invVtheta[1:4, 1:4] ~ dwish(diag4[,], 4)
    invVbeta[1:4, 1:4] ~ dwish(diag4[,], 4)
    invVgamma[1:2, 1:2] ~ dwish(diag2[,], 2)
    invValpha[1:2, 1:2] ~ dwish(diag2[,], 2)

    invOmega[1:4, 1:4] ~ dwish(diag4[,], 4)
}
"

writeLines( modelString, con="TEMPmodel.txt")


#---------------------- For JAGS ---------------------------
nChains = 3
nAdaptSteps = 1000
nBurninSteps = 500
nUseSteps=10000  # total number of used steps
nThinSteps=2

runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=c("CTR", "CONV", "AdPos", "organicComp",
                                  "thetak", "betak", "gammak", "alphak",
                                  "theta", "beta", "gamma", "alpha",
                                  "deltaTheta", "deltaBeta", "deltaGamma", "deltaAlpha") ,
                        data=data_list ,
                        #inits=initsList ,
                        n.chains=nChains ,
                        adapt=nAdaptSteps ,
                        burnin=nBurninSteps ,
                        sample=ceiling(nUseSteps/nChains) ,
                        thin=nThinSteps ,
                        summarise=FALSE ,
                        plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )
