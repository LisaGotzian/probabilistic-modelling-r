#-------------------- 2) BUGS attempt -------------------------
# Purpose:
# 1) simulate the data and parameters in the model by
#   Agarwal et al. 2015 on organic search results.
# 2) Run the model in JAGS. As it is recursive, JAGS doesn't work.
#   Therefore, BUGS is tried, but also shows its fallacies. Read more in
#   probabilistic_modelling_agarwal.pdf.

# Lisa Gotzian, Erik Schau, August 2018

rm(list=ls())  # Careful! This clears all of R's memory!

set.seed(66)
library(R2OpenBUGS) #for BUGS
library('runjags') #for JAGS

source("data_simulation.R")  # the simulated data
source("DBDA2E-utilities.R") #for JAGS
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




modelString <-"model{ # Careful - reversed many indices b/c BUGS requires it
for ( t in 1:N ) {
for ( k in 1:M ){
CTR[k,t] ~ dnorm(0.5,pCTR[k,t]) # I give up on making this any better
CONV[k,t] ~ dnorm(0.5,pCONV[k,t])

}
}
for ( t in 1:N ) {
for ( k in 1:M ){
pCTR[k,t] <- exp(UCtr[k,t])/(1+exp(UCtr[k,t]))
pCONV[k,t] <- exp(UConv[k,t])/(1+exp(UConv[k,t]))

UCtr[k,t] <- thetak[k,1] + thetak[k,2] * AdPos[k,t] + thetak[k,3] * organicComp[k,t] +
thetak[k,4] * sponsored_comp[k,t] + theta[1] * organic[k,t] +
theta[2] * lqscore[k,t] + theta[3] * t + epsilon[k,t,2]

UConv[k,t] <- betak[k,1] + betak[k,2] * AdPos[k,t] + betak[k,3] * organicComp[k,t] +
betak[k,4] * sponsored_comp[k,t] + beta[1] * organic[k,t] +
beta[2] * lqscore[k,t] + beta[3] * t + epsilon[k,t,1]

log(AdPos[k,t]) <- gammak[k,1] + gammak[k,2] * log(bid[k,t]) +
gamma[1] * lqscore[k,t] + gamma[2] * t + epsilon[k,t,3]

log(organicComp[k,t]) <- alphak[k,1] + alphak[k,2] * iv_organic[k,t] +
alpha * t + epsilon[k,t,4]


epsilon[k,t,1:4] ~ dmnorm(zero4[], invOmega[,])
}
}

for ( k in 1:M){
for (i in 1:4){ # BUGS cannot do matrix multiplication :(
# also, the delta matrices like deltaTheta are the only variables
# that have stayed with the same dimensions.
deltaThetaMult[k,i] <- deltaTheta[i,1] * z[1,k] + deltaTheta[i,2] * z[2,k]
deltaBetaMult[k,i] <- deltaBeta[i,1] * z[1,k] + deltaBeta[i,2] * z[2,k]
}
for (i in 1:2){
deltaGammaMult[k,i] <- deltaGamma[i,1] * z[1,k] + deltaGamma[i,2] * z[2,k]
deltaAlphaMult[k,i] <- deltaAlpha[i,1] * z[1,k] + deltaAlpha[i,2] * z[2,k]
}

for (i in 1:4){ # as BUGS doesn't let me define vectors in this way:
# thetak[k, 1:4] <- deltaThetaMult[k,1:4] + uTheta[k,1:4]
# for distributions, this apparently is ok. It behaves like a little kid.
# Have I ever cussed mmore often in code? I think not.
thetak[k,i] <- deltaThetaMult[k,i] + uTheta[k,i]
betak[k, i] <- deltaBetaMult[k,i] + uBeta[k,i]
}
for (i in 1:2){
gammak[k, i] <- deltaGammaMult[k,i] + uGamma[k,i]
alphak[k, i] <- deltaAlphaMult[k,i] + uAlpha[k,i]
}


uTheta[k, 1:4] ~ dmnorm(zero4[], invVtheta[,])
uBeta[k, 1:4] ~ dmnorm(zero4[], invVbeta[,])
uGamma[k, 1:2] ~ dmnorm(zero2[], invVgamma[,])
uAlpha[k, 1:2] ~ dmnorm(zero2[], invValpha[,])
}

##### This entire part does only a matrix multiplication by hand. And yes, I was thinking about doing
# it in an array, but it seemed to me this is even more hairy.
for (k in 1:M){
for (i in 1:4){ #both z by hand, the 4 theta in this loop
# for the motherland, for beta! Matrix multiplication, what a nightmare you are.
# betak is actually transposed, that's why its index looks off. bugs doesn't know t()
newzbeta[k,i,1] <- z[1,k] * betak[k,i]
newzbeta[k,i,2] <- z[2,k] * betak[k,i]
zbetak[1,i,k] <- step(1.5-k)*newzbeta[k,i,1] + step(-1.5+k) * zbetak[1,i,k] + newzbeta[k,i,1]
zbetak[2,i,k] <- step(1.5-k)*newzbeta[k,i,2] + step(-1.5+k) * zbetak[2,i,k] + newzbeta[k,i,2]

# coming to the theta...
newztheta[k,i,1] <- z[1,k] * thetak[k,i]
newztheta[k,i,2] <- z[2,k] * thetak[k,i]
zthetak[1,i,k] <- step(1.5-k)*newztheta[k,i,1] + step(-1.5+k) * zthetak[1,i,k] + newztheta[k,i,1]
zthetak[2,i,k] <- step(1.5-k)*newztheta[k,i,2] + step(-1.5+k) * zthetak[2,i,k] + newztheta[k,i,2]
}

for (i in 1:2){
# looking at you, gamma
newzgamma[k,i,1] <- z[1,k] * gammak[k,i]
newzgamma[k,i,2] <- z[2,k] * gammak[k,i]

# the step function implements a conditional if statement
zgammak[1,i,k] <- step(1.5-k)*newzgamma[k,i,1] + step(-1.5+k) * zgammak[1,i,k] + newzgamma[k,i,1]
zgammak[2,i,k] <- step(1.5-k)*newzgamma[k,i,2] + step(-1.5+k) * zgammak[2,i,k] + newzgamma[k,i,2]


# all the hate now targeted towards alpha
newzalpha[k,i,1] <- z[1,k] * alphak[k,i]
newzalpha[k,i,2] <- z[2,k] * alphak[k,i]
zalphak[1,i,k] <- step(1.5-k)*newzalpha[k,i,1] + step(-1.5+k) * zalphak[1,i,k] + newzalpha[k,i,1]
zalphak[2,i,k] <- step(1.5-k)*newzalpha[k,i,2] + step(-1.5+k) * zalphak[2,i,k] + newzalpha[k,i,2]
}
}


for (j in 1:2){
for (i in 1:4){ #the last index on zbetak is a dummy index
deltaTildeBeta[j,i] <- deltaq[j,1] * zbetak[1,i,1] + deltaq[j,2] * zbetak[2,i,1]
deltaTildeTheta[j,i] <- deltaq[j,1] * zthetak[1,i,1] + deltaq[j,2] * zthetak[2,i,1]
}
for (i in 1:2){
deltaTildeGamma[j,i] <- deltaq[j,1] * zgammak[1,i,1] + deltaq[j,2] * zgammak[2,i,1]
deltaTildeAlpha[j,i] <- deltaq[j,1] * zalphak[1,i,1] + deltaq[j,2] * zalphak[2,i,1]
}
}

# and all the hustle to do this!!
# deltaTildeBeta <- deltaq %*% (z %*% t(betak))
# deltaTildeTheta <- deltaq %*% (z %*% t(thetak))
# deltaTildeGamma <- deltaq %*% (z %*% t(gammak))
# deltaTildeAlpha <- deltaq %*% (z %*% t(alphak))
##### Matrix multiplication end

for (i in 1:4){
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
# Uncomment later
nChains = 3
nAdaptSteps = 1000
nBurninSteps = 500
nUseSteps=10000  # total number of used steps
nThinSteps=2

# runJagsOut <- run.jags( method="parallel" ,
#                         model="TEMPmodel.txt" ,
#                         monitor=c("CTR", "CONV", "AdPos", "organicComp",
#                                   "thetak", "betak", "gammak", "alphak",
#                                   "theta", "beta", "gamma", "alpha",
#                                   "deltaTheta", "deltaBeta", "deltaGamma", "deltaAlpha") ,
#                         data=data_list ,
#                         #inits=initsList ,
#                         n.chains=nChains ,
#                         adapt=nAdaptSteps ,
#                         burnin=nBurninSteps ,
#                         sample=ceiling(nUseSteps/nChains) ,
#                         thin=nThinSteps ,
#                         summarise=FALSE ,
#                         plots=FALSE )
# codaSamples = as.mcmc.list( runJagsOut )

#---------------------- For BUGS ---------------------------
parameters = c("CTR", "CONV", "AdPos", "organicComp",
               "thetak", "betak", "gammak", "alphak",
               "theta", "beta", "gamma", "alpha",
               "deltaTheta", "deltaBeta", "deltaGamma", "deltaAlpha",
               "invOmega", "invVtheta", "invVbeta", "invVgamma", "invValpha",
               "deltaTheta", "deltaBeta", "deltaGamma", "deltaAlpha",
               "uTheta", "uBeta", "uGamma", "uAlpha")
in1 = list(beta = as.vector(matrix(1, nrow = 3, ncol = 1)),
           theta = as.vector(matrix(1, nrow = 3, ncol = 1)),
           gamma = as.vector(matrix(1, nrow = 2, ncol = 1)),
           alpha = 1,
           
           deltaBeta = matrix(0, nrow = 4, ncol = 2),
           deltaTheta = matrix(0, nrow = 4, ncol = 2),
           deltaGamma = matrix(0, nrow = 2, ncol = 2),
           deltaAlpha = matrix(0, nrow = 2, ncol = 2),
           
           invVbeta = diag(x = 1, nrow = 4, ncol = 4),
           invVtheta = diag(x = 1, nrow = 4, ncol = 4),
           invVgamma = diag(x = 1, nrow = 2, ncol = 2),
           invValpha = diag(x = 1, nrow = 2, ncol = 2),
           
           invOmega = diag(x = 1, nrow = 4, ncol = 4),
           
           # careful: these are weirdly defined b/c BUGS needs the same k uThetas
           uTheta = matrix(0, nrow = k, ncol = 4),
           uBeta = matrix(0, nrow = k, ncol = 4),
           uGamma = matrix(0, nrow = k, ncol = 2),
           uAlpha = matrix(0, nrow = k, ncol = 2),
           
           epsilon = array(0, dim=c(k,t,4)))
in2 = in1
in3 = in1
inits = list(in1, in2, in3)

bugs.inits(inits, nChains, digits = 4)

bugsdata <- bugs.data(data_list, data.file = "bugsdata.txt")


# set the WINE working directory and the directory to OpenBUGS - change the OpenBUGS.exe location as necessary
WINE="/usr/local/bin/wine"
#WINEPATH="/usr/local/bin/winepath"
OpenBUGS.pgm="/Users/lisagotzian/.wine/drive_c/Program Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe"

bugstest <- bugs(data = "bugsdata.txt", inits = inits, parameters = parameters, model.file = "TEMPmodel.txt",
                 n.chains = nChains, n.iter = nUseSteps, debug = T, OpenBUGS.pgm=OpenBUGS.pgm, WINE=WINE,useWINE=T)
print(bugstest)


