#-------------------- 3) The raw Metropolis Algorithm -------------------------
# Purpose:
# 1) simulate the data and parameters in the model by
#   Agarwal et al. 2015 on organic search results.
# 2) Run the model in JAGS. As it is recursive, JAGS doesn't work.
#   Therefore, BUGS is tried, but also shows its fallacies. That's why the raw
#   Metropolis is tried out in this script. Read more in
#   probabilistic_modelling_agarwal.pdf.
#
# Exact steps:
# - Following Agarwal et al.'s approach, a random walk Metropolis-Hastings
#   algorithm is built. Acknowledgements to Rossi et al. 2005.
# - priors, likelihood and the steps are listed in the author's appendix.
# - It is possible to source this script.

# Lisa Gotzian, August 2018

#------------------------ Preliminaries --------------------------
rm(list=ls())  # Careful! This clears all of R's memory!

library(pacman)
p_load(SciViews, # for the convenience of ln()
       blockmatrix, # for dealing with a matrix of matrices
       MCMCpack, #for calling the inverse wishart distribution
       tidyverse,
       mvtnorm,
       truncnorm, # to use the truncated norm
       matrixcalc, # to use is.positive.definite
       beepr) # rather daring, but for the supermario sound when it's done
              
## Values from the paper
set.seed(1)
k = 36 # number of keywords
t = 40 # number of days
num_observations = k*t # number of observations
nIter = 100
burnIn = 40000


#-------------------------- The data -----------------------------
# poisson impressions
impressions <- matrix(round(rpois(n=num_observations, lambda = rep(72.8, num_observations)), 0),
                      nrow=k, ncol=t)
impressions[impressions > 1666] <- 1666
impressions[impressions < 1] <- 1


# poisson clicks
clicks <- matrix(round(rpois(n=num_observations, lambda = rep(1.1, num_observations)),0),
                 nrow=k, ncol=t)
clicks[clicks > 24] <- 24
clicks[clicks < 0] <- 0

# poission orders
orders <- matrix(round(rpois(n=num_observations, lambda = rep(0.03, num_observations)),2),
                 nrow=k, ncol=t)
orders[orders > 3] <- 3
orders[orders < 0] <- 0


# Maybe take out AdPos/organic_comp as they are endogeneous
AdPos <- matrix(round(rtruncnorm(n=num_observations, a=1, b=9.78, mean=3.47, sd=1.7),2),
                nrow=k, ncol=t)
organic <- matrix(as.numeric(rbinom(n=num_observations, size = 1, p=0.15)),
                  nrow=k, ncol=t)
organic_comp <- matrix(round(rtruncnorm(n=num_observations, a=0, b=3.08, mean=0.78, sd=0.54),2),
                      nrow=k, ncol=t)
iv_organic <- matrix(round(rtruncnorm(n=num_observations, a=0.1, b=2.34, mean=0.42, sd=0.23),2),
                     nrow=k, ncol=t)
sponsored_comp <- matrix(round(rtruncnorm(n=num_observations, a=0, b=5.59, mean=1.79, sd=0.89),2),
                         nrow=k, ncol=t)
iv_sponsored <- matrix(round(rtruncnorm(n=num_observations, a=0.13, b=5.19, mean=1.49, sd=1),2),
                       nrow=k, ncol=t)
lqscore <- matrix(round(rtruncnorm(n=num_observations, a=6, b=10, mean=8, sd=1.5),0),
                  nrow=k, ncol=t)
brand <- matrix(rbinom(n=k, size = 1, p=0.6),
                nrow=k, ncol=1)

# binomial distribution specificity
specificity <- matrix(rbinom(n=k, size=1, p=0.4),
                      nrow=k, ncol=1)
bid <- matrix(round(rtruncnorm(n=num_observations, a=0.08, b=2, mean=0.5, sd=0.3),2),
              nrow=k, ncol=t)

# We understand time as time t/40 as big values yield too big values for the omega matrices.


#-----instantiate variables to use in the equations------#
# this we'll have to do: compare our thetas from calc_weights to the paper (sd given)
# create k z's
z <- matrix(nrow = 2, ncol = k)
for (i in 1:k){
  z[,i] <- as.numeric(c(brand[i], specificity[i]))
}

# Calculating clickthrough rate based on the data
CTR <- clicks/impressions

# Calculating conversion rate based on the data
CONV <- orders/impressions

#--------------------------- The MCMC ------------------------

##--------------------- Initialize the chains ----------------
# Initialize the latent utilities
# If it has to be initialized, it's a mistake in the paper.
uCtrOld <- matrix(0, nrow = k, ncol = t) 
uConvOld <- matrix(0, nrow = k, ncol = t)

uCtrNew <- uConvNew <- 
LambdaCtrOld <- LambdaConvOld <- 
LambdaCtrNew <- LambdaConvNew <- 
lUOldmat <- lUNewmat <- 
Xthetak <- Xtheta <- Xctr <- 
Xbetak <- Xbeta <- Xconv <- 
Xgammak <- Xgamma <- XadPos <- eAdPos <- 
Xalphak <- Xalpha <- Xorganic_comp <- eorganic_comp <- matrix(nrow = k, ncol = t)

epsilon <- array(dim = c(k, t, 4))

# betak is the matrix for the betas depending on k, therefore with k columns,
# beta is the vector uniform for all k's.
betak <- matrix(0, nrow = 4, ncol = k)
beta <- as.vector(matrix(0, nrow = 3, ncol = 1))
thetak <- matrix(0, nrow = 4, ncol = k)
theta <- as.vector(matrix(0, nrow = 3, ncol = 1))
gammak <- matrix(0, nrow = 2, ncol = k)
gamma <- as.vector(matrix(0, nrow = 2, ncol = 1))
alphak <- matrix(1, nrow = 2, ncol = k) # to avoid ln(orgcomp) being negative
alpha <- 1

# The delta matrices have the same amount of rows as there are coefficients
# depending on k, e.g. for beta, that's 4. The 2 columns denote the effect
# on brand and specificity.
deltaBeta <- matrix(0, nrow = 4, ncol = 2)
deltaTheta <- matrix(0, nrow = 4, ncol = 2)
deltaGamma <- matrix(0, nrow = 2, ncol = 2)
deltaAlpha <- matrix(0, nrow = 2, ncol = 2)

# mistake: wasn't given in the paper
Omega <- diag(x = 1, nrow = 4, ncol = 4)


## Values for step 2, the keyword specific thetas
# These four xk arrays are written in one matrix with values only on the diagonal, probably to draw the inverse
# of each of them. They are seperated here for easier calculation.
# the matrices are transposed. The first column is the constant bias like beta0.
ones <- matrix(1, nrow = k, ncol = t)

step2xkbetat <- array(dim = c(t, 4, k))
step2xkthetat <- array(dim = c(t, 4, k))
step2xkgammat <- array(dim = c(t, 2, k))
step2xkalphat <- array(dim = c(t, 2, k))

# in the second step, the transposed matrices like step2xkbetat are transposed again
step2xkbeta <- array(dim = c(4, t, k))
step2xktheta <- array(dim = c(4, t, k))
step2xkgamma <- array(dim = c(2, t, k))
step2xkalpha <- array(dim = c(2, t, k))

# the only relevant data for keyword-specific thetas
for (i in 1:k){
  step2xkbetat[,, i] <- array(c(ones[i,], AdPos[i,], organic_comp[i,], sponsored_comp[i,]), dim = c(t, 4))
  step2xkthetat[,, i] <- array(c(ones[i,], AdPos[i,], organic_comp[i,], sponsored_comp[i,]), dim = c(t, 4))
  step2xkgammat[,, i] <- array(c(ones[i,], bid[i,]), dim = c(t, 2))
  step2xkalphat[,, i] <- array(c(ones[i,], iv_organic[i,]), dim = c(t, 2))
  
  step2xkbeta[,, i] <- t(step2xkbetat[,, i])
  step2xktheta[,, i] <- t(step2xkthetat[,, i])
  step2xkgamma[,, i] <- t(step2xkgammat[,, i])
  step2xkalpha[,, i] <- t(step2xkalphat[,, i])
}

# Values for V
Vbeta <- diag(x = 1, nrow = 4, ncol = 4)
Vtheta <- diag(x = 1, nrow = 4, ncol = 4)
Vgamma <- diag(x = 1, nrow = 2, ncol = 2)
Valpha <- diag(x = 1, nrow = 2, ncol = 2)

## Values for step 3, the non-keyword specific thetas
step3xbetat <- array(dim = c(t, 3, k))
step3xthetat <- array(dim = c(t, 3, k))
step3xgammat <- array(dim = c(t, 2, k))
step3xalphat <- array(dim = c(t, 1, k))

# in the second step, the transposed matrices like step2xkbetat are transposed again
step3xbeta <- array(dim = c(3, t, k))
step3xtheta <- array(dim = c(3, t, k))
step3xgamma <- array(dim = c(2, t, k))
step3xalpha <- array(dim = c(1, t, k))

#uniform thetas (different from step 2)
for (i in 1:k){
  step3xbetat[,, i] <- array(c(organic[i,], lqscore[i,], (1:t)/40), dim = c(t, 3))
  step3xthetat[,, i] <- array(c(organic[i,], lqscore[i,], (1:t)/40), dim = c(t, 3))
  step3xgammat[,, i] <- array(c(lqscore[i,], (1:t)/40), dim = c(t, 2))
  step3xalphat[,, i] <- array(c((1:t)/40), dim = c(t, 1))
  
  step3xbeta[,, i] <- t(step3xbetat[,, i])
  step3xtheta[,, i] <- t(step3xthetat[,, i])
  step3xgamma[,, i] <- t(step3xgammat[,, i])
  step3xalpha[,, i] <- t(step3xalphat[,, i])
}


## Values for step 4
vOmega = 10

# For monitoring
Uarray <- array(dim = c(k,t,2,nIter))

print("Chain initialized")

# start nIter iterations and plan to have 40,000 as burn-in, 80,000 total.
#-------------------- Step 1: Draw the main part ----------------------
for (g in 1:nIter){
  for ( j in 1:t ) {
   for ( i in 1:k ){
      
      deltaCtr <-  rnorm(1,0,0.2) # draw 1 with mean 0 and sd 0.2
      deltaConv <-  rnorm(1,0,0.2)
      
      ## Aus alt mach neu
      uCtrNew[i,j] = uCtrOld[i,j] + deltaCtr
      uConvNew[i,j] = uConvOld[i,j] + deltaConv
      
      ## The logit model for the probability of clicking
      LambdaCtrOld[i,j] = exp(uCtrOld[i,j])/(1 + exp(uCtrOld[i,j]))
      LambdaConvOld[i,j] = exp(uConvOld[i,j])/(1 + exp(uConvOld[i,j]))
      
      LambdaCtrNew[i,j] = exp(uCtrNew[i,j])/(1 + exp(uCtrNew[i,j]))
      LambdaConvNew[i,j] = exp(uConvNew[i,j])/(1 + exp(uConvNew[i,j]))
      
      ## The likelihood for orders and clicks
      # It's the product (colloseum). It's later handled as one value by prod()
      lUOldmat[i,j] = ( (LambdaConvOld[i,j] * LambdaCtrOld[i,j])^orders[i,j] *
                   ((1-LambdaConvOld[i,j])*LambdaCtrOld[i,j])^(clicks[i,j]-orders[i,j]) * (1-LambdaCtrOld[i,j])^(impressions[i,j]-clicks[i,j]))
      lUNewmat[i,j] = ( (LambdaConvNew[i,j] * LambdaCtrNew[i,j])^orders[i,j] *
                       ((1-LambdaConvNew[i,j])*LambdaCtrNew[i,j])^(clicks[i,j]-orders[i,j]) * (1-LambdaCtrNew[i,j])^(impressions[i,j]-clicks[i,j]))
  
      ## The linear model formulation for the latent utility
      # Attention: due to indexing, thetak[1,] is actually theta0
      
      Xthetak[i,j] = thetak[1, i] + thetak[2, i] * AdPos[i,j] +
        thetak[3, i] * organic_comp[i,j] +
        thetak[4, i] * sponsored_comp[i,j]
      Xtheta[i,j] = theta[1] * organic[i,j] +
        theta[2] * lqscore[i,j] +
        theta[3] * j/t
      Xctr[i,j] =  Xthetak[i,j] + Xtheta[i,j]
      
      Xbetak[i,j] = betak[1, i] + betak[2, i] * AdPos[i,j] +
        betak[3, i] * organic_comp[i,j] +
        betak[4, i] * sponsored_comp[i,j]
      Xbeta[i,j] = beta[1] * organic[i,j] +
        beta[2] * lqscore[i,j] +
        beta[3] * j/t
      Xconv[i,j] = Xbetak[i,j] + Xbeta[i,j]
        
      
      ## The instrumental variables
      # Draw the epsilons
      epsilon[i,j, ] <- rmvnorm(n = 1, mean = rep(0, 4), sigma = solve(Omega))
      
      # As from this distribution, the epsilon can sometimes become too big (even with the identity
      # matrix, this "cheat" was introduced to ensure positive values of AdPos and organic_comp):
      # replace negative values by 0.01
      
      Xgammak[i,j] = gammak[1, i] + gammak[2, i] * ln ( bid[i,j] )
      Xgamma[i,j] = gamma[1] * lqscore[i,j] + gamma[2] * j/t
      XadPos[i,j] =  Xgammak[i,j] + Xgamma[i,j]
      AdPos[i,j] = exp ( XadPos[i,j] + epsilon[i,j,3] )
      if (any(AdPos[i,j] <= 0, is.infinite(AdPos[i,j])) ){
        AdPos[i,j] = 0.01
      }
      eAdPos[i,j] = ln ( AdPos[i,j] ) - XadPos[i,j]
      
      Xalphak[i,j] = alphak[1, i] + alphak[2, i] * iv_organic[i,j]
      Xalpha[i,j] =  alpha * j/t
      Xorganic_comp[i,j] =  Xalphak[i,j] + Xalpha[i,j]
      organic_comp[i,j] = Xorganic_comp[i,j] + epsilon[i,j,4]
      if (any(organic_comp[i,j] <= 0, is.infinite(organic_comp[i,j]))){
        organic_comp[i,j] = 0.01
      }
      eorganic_comp[i,j] = ln ( organic_comp[i,j] ) - Xorganic_comp[i,j]
    }
  }
  
  ## Preparation of other variables
  uNew <- array(c(uCtrNew, uConvNew), dim = c(k, t, 2))
  uOld <- array(c(uCtrOld, uConvOld), dim = c(k, t, 2))
  lUOld <- prod(lUOldmat)
  lUNew <- prod(lUNewmat)
  X <- array(c(Xctr, Xconv), dim = c(k, t, 2))
  e <- array(c(eAdPos, eorganic_comp), dim = c(k, t, 2))
  W11 <- Omega[1:2, 1:2]
  W22 <- Omega[3:4, 3:4]
  W12 <- Omega[1:2, 3:4]
  
  # matrix times array by hand
  matrixTimesArray <- function(matrix, array){ #the easy case for a 2x2 * kxtx2
    product <- array(c(matrix[1,1] * array[,,1] + matrix[1,2] * array[,,2],
                       matrix[2,1] * array[,,1] + matrix[2,2] * array[,,2]),
                     dim = c(k, t, 2))
    return(product)
  }
  
  Wprep <- W12 %*% solve(W22)
  E <- matrixTimesArray(matrix = Wprep, array = e)
  InvA <- W11 - W12 %*% solve(W22) %*% W12
  A <- solve(InvA)
  
  ## Acceptance of the jump
  # The proposed algorithm accepts the entire draw or rejects it completely.
  # Matrix multiplication by hand. UXE is a step in between.
  UXEnew = uNew - X - E
  UXEold = uOld - X - E
  UXEnewarray = matrixTimesArray(matrix = A, array = UXEnew)
  UXEoldarray = matrixTimesArray(matrix = A, array = UXEold)
  
  acceptance <- array(dim = c(k,t))
  for (j in 1:t){
    for (i in 1:k){
      acceptance[i,j] = ( exp( -0.5* (UXEnew[i,j,1] * UXEnewarray[i,j,1] + 
                               UXEnew[i,j,2] * UXEnewarray[i,j,2]) * lUNew))/
        ( exp( -0.5 * (UXEold[i,j,1] * UXEoldarray[i,j,2] +
                     UXEold[i,j,2] * UXEoldarray[i,j,2]) * lUOld))
     
       if (acceptance[i,j] < 1){ 
        uCtrOld[i,j] = uCtrNew[i,j]
        uConvOld[i,j] = uConvNew[i,j]
      }
    }
  }
  
  uOld <- array(c(uCtrOld, uConvOld), dim = c(k,t,2))
  
  
  #------- Step 2: Draw b[k] = [ beta[k], theta[k], gamma[k], alpha[k] ] ------
  ## Preliminaries
  yk <- array(c(
    uConvOld - Xbeta,
    uCtrOld - Xtheta,
    ln ( AdPos ) - Xgamma,
    ln ( organic_comp ) - Xalpha), dim = c(k,t,4))
  
  #make it 36 matrices to be multiplied with the rest
  yk <- aperm(yk)
  
  # one entry of omega times one entry of the entire array step2xkbetat
  Omegaxkbeta <- Omega[1,1] * step2xkbetat 
  Omegaxktheta <- Omega[2,2] * step2xkthetat
  Omegaxkgamma <- Omega[3,3] * step2xkgammat
  Omegaxkalpha <- Omega[4,4] * step2xkalphat
  
  
  ## Calculating the covariance matrix Qk
  # allocate space
  Qkbeta <- array(dim = c(4, 4, k))
  Qktheta <- array(dim = c(4, 4, k))
  Qkgamma <- array(dim = c(2, 2, k))
  Qkalpha <- array(dim = c(2, 2, k))
  
  # Inverse of V-1 is just an inverse of each of the V matrices as the inverse of a diagonal matrix is
  # the inverses of each of the diagonals
  for (i in 1:k) { # handling everything on the keyword level makes things more comprehensive
    Qkbeta[,, i] <- solve(solve(step2xkbeta[,, i] %*% Omegaxkbeta[,, i]) + # returns a 4x4 matrix
                            solve(Vbeta))
    Qktheta[,, i] <- solve(solve(step2xktheta[,, i] %*% Omegaxktheta[,, i]) + 
                            solve(Vtheta))
    Qkgamma[,, i] <- solve(solve(step2xkgamma[,, i] %*% Omegaxkgamma[,, i]) + 
                            solve(Vgamma))
    Qkalpha[,, i] <- solve(solve(step2xkalpha[,, i] %*% Omegaxkalpha[,, i]) + 
                            solve(Valpha))
  }
  
  ## Calculate the mean bTildek
  InvOmega <- solve(Omega)
  # Create the values within the bBark vector
  deltaBetaz <- deltaBeta %*% z # 4x2 times 2x36 gives a 4x36 matrix
  deltaThetaz <- deltaTheta %*% z
  deltaGammaz <- deltaGamma %*% z
  deltaAlphaz <- deltaAlpha %*% z
  
  # The last part of bTildek: 2 4x36 and 2 2x36 matrices
  # This cannot be in one array because gamma and alpha have
  # different dimensions from beta and theta.
  step2addition2Beta <- solve(Vbeta) %*% deltaBetaz
  step2addition2Theta <- solve(Vtheta) %*% deltaThetaz 
  step2addition2Gamma <- solve(Vgamma) %*% deltaGammaz
  step2addition2Alpha <- solve(Valpha) %*% deltaAlphaz
  
  # to have things sorted by k again, this loop is necessary...
  step2InverseProd <- array(dim = c(4,t,k))
  bTildekBeta <- array(dim = c(4,k))
  bTildekTheta <- array(dim = c(4,k))
  bTildekGamma <- array(dim = c(2,k))
  bTildekAlpha <- array(dim = c(2,k))
  
  for (i in 1:k){
    step2InverseProd[,,i] <- InvOmega %*% yk[,,i] # 4x4 times a 4x40 matrix
    # from now on, the first row belongs to beta, the second to theta etc.
    
    bTildekBeta[,i] <- Qkbeta[,,i] %*% (step2xkbeta[,,i] %*% step2InverseProd[1,,i]) +
      step2addition2Beta[,i]
    bTildekTheta[,i] <- Qktheta[,,i] %*% (step2xktheta[,,i] %*% step2InverseProd[2,,i]) +
      step2addition2Theta[,i]
    bTildekGamma[,i] <- Qkgamma[,,i] %*% (step2xkgamma[,,i] %*% step2InverseProd[3,,i]) +
      step2addition2Gamma[,i]
    bTildekAlpha[,i] <- Qkalpha[,,i] %*% (step2xkalpha[,,i] %*% step2InverseProd[4,,i]) +
      step2addition2Alpha[,i]
  }
  
  
  ## Final step: draw all betak, thetak, gammak and alphak
  for (i in 1:k){
    betak[,i] = rmvnorm(n = 1, mean = bTildekBeta[,i], sigma = Qkbeta[,,i])
    thetak[,i] = rmvnorm(n = 1, mean = bTildekTheta[,i], sigma = Qktheta[,,i])
    gammak[,i] = rmvnorm(n = 1, mean = bTildekGamma[,i], sigma = Qkgamma[,,i])
    alphak[,i] = rmvnorm(n = 1, mean = bTildekAlpha[,i], sigma = Qkalpha[,,i])
  }
  
  #------------- Step 3: Draw b = [ beta, theta, gamma, alpha ] --------------
  # This section was basically copied from above and only altered regarding k 
  # and a few minor changes such as Vbar, bBar and the y vector (including x
  # which was adjusted before the chains).
  ## Preliminaries
  y <- array(c(
    uConvOld - Xbetak,
    uCtrOld - Xthetak,
    ln ( AdPos ) - Xgammak,
    ln ( organic_comp ) - Xalphak), dim = c(k,t,4))
  
  y <- aperm(y)
  
  # one entry of omega times one entry of the entire array step3xbetat
  Omegaxbeta <- Omega[1,1] * step3xbetat 
  Omegaxtheta <- Omega[2,2] * step3xthetat
  Omegaxgamma <- Omega[3,3] * step3xgammat
  Omegaxalpha <- Omega[4,4] * step3xalphat
  
  
  ## Calculating the covariance matrix Q
  # allocate space - adjusted the dimensions in comparison to step 2
  QbetaArray <- array(dim = c(3, 3, k))
  QthetaArray <- array(dim = c(3, 3, k))
  QgammaArray <- array(dim = c(2, 2, k))
  QalphaArray <- array(dim = c(1, 1, k))
  
  # Inverse of V-1 is just an inverse of each of the V matrices as the inverse of a diagonal matrix is
  # the inverses of each of the diagonals - just did it as diagonals here
  for (i in 1:k) { # handling everything on the keyword level makes things more comprehensive
    QbetaArray[,, i] <- solve(solve(step3xbeta[,, i] %*% Omegaxbeta[,, i]) + # returns a 3x3 matrix
                            solve(diag(x = 100, nrow = 3, ncol = 3))) 
    QthetaArray[,, i] <- solve(solve(step3xtheta[,, i] %*% Omegaxtheta[,, i]) + 
                            solve(diag(x = 100, nrow = 3, ncol = 3)))
    QgammaArray[,, i] <- solve(solve(step3xgamma[,, i] %*% Omegaxgamma[,, i]) + 
                            solve(diag(x = 100, nrow = 2, ncol = 2)))
    QalphaArray[,, i] <- solve(solve(step3xalpha[,, i] %*% Omegaxalpha[,, i]) + 
                            solve(diag(x = 100, nrow = 1, ncol = 1)))
  }
  
  ## Calculate the mean bTilde
  # The bBar vector is just 0's, I suspect the part is just included
  # to show the similarity to step 2. It's left out here.
  
  # to have things sorted by k again, this loop is necessary...
  # changed the dimensions as compared to step 2
  step3InverseProd <- array(dim = c(4,t,k))
  bTildeBeta <- array(dim = c(3,k))
  bTildeTheta <- array(dim = c(3,k))
  bTildeGamma <- array(dim = c(2,k))
  bTildeAlpha <- array(dim = c(1,k))
  
  for (i in 1:k){
    step3InverseProd[,,i] <- InvOmega %*% y[,,i] # 4x4 times a 4x40 matrix
    # from now on, the first row belongs to beta, the second to theta etc.
    
    # disregarded solve(vBar)*bBar as this is 0
    bTildeBeta[,i] <- QbetaArray[,,i] %*% (step3xbeta[,,i] %*% step3InverseProd[1,,i])
    bTildeTheta[,i] <- QthetaArray[,,i] %*% (step3xtheta[,,i] %*% step3InverseProd[2,,i])
    bTildeGamma[,i] <- QgammaArray[,,i] %*% (step3xgamma[,,i] %*% step3InverseProd[3,,i])
    bTildeAlpha[,i] <- QalphaArray[,,i] %*% (step3xalpha[,,i] %*% step3InverseProd[4,,i])
  }
  
  
  ## Final step: draw all beta, theta, gamma and alpha
  # this is cheated: unconvential step in between: I took the mean out of all k's
  # because I just have data for each k.
  Qbeta <- rowMeans(QbetaArray, dims = 2)
  Qtheta <- rowMeans(QthetaArray, dims = 2)
  Qgamma <- rowMeans(QgammaArray, dims = 2)
  Qalpha <- rowMeans(QalphaArray, dims = 2)
  
  bTildeBeta <- rowMeans(bTildeBeta)
  bTildeTheta <- rowMeans(bTildeTheta)
  bTildeGamma <- rowMeans(bTildeGamma)
  bTildeAlpha <- rowMeans(bTildeAlpha)
  
  beta = rmvnorm(n = 1, mean = bTildeBeta, sigma = Qbeta)
  theta = rmvnorm(n = 1, mean = bTildeTheta, sigma = Qtheta)
  gamma = rmvnorm(n = 1, mean = bTildeGamma, sigma = Qgamma)
  alpha = rnorm(n = 1, mean = bTildeAlpha, sd = Qalpha)
  
  
  #------------------------ Step 4: Draw Omega ------------------------------
  ## Preliminaries
  ykt <- array(c(
    uConvOld - Xbetak - Xbeta,
    uCtrOld - Xthetak - Xtheta,
    ln ( AdPos ) - Xgammak - Xgamma,
    ln ( organic_comp ) - Xalphak - Xalpha
  ), dim = c(k,t,4))
  
  # Calculating the S matrix
  # The transpose in the formula corresponds to the dot product which has been
  # done by hand for the arrays in this case. 
  # Unconventionally, the transpose was used to enable the matrix
  # multiplication. After three ys, I believe the ys on this page
  # have been incorrectly as column vectors of matrices instead of
  # row vectors.
  Somega <- matrix(
    c(sum(ykt[,,1] %*% t(ykt[,,1])), sum(ykt[,,1] %*% t(ykt[,,2])), sum(ykt[,,1] %*% t(ykt[,,3])), sum(ykt[,,1] %*% t(ykt[,,4])),
      sum(ykt[,,2] %*% t(ykt[,,1])), sum(ykt[,,2] %*% t(ykt[,,2])), sum(ykt[,,2] %*% t(ykt[,,3])), sum(ykt[,,2] %*% t(ykt[,,4])),
      sum(ykt[,,3] %*% t(ykt[,,1])), sum(ykt[,,3] %*% t(ykt[,,2])), sum(ykt[,,3] %*% t(ykt[,,3])), sum(ykt[,,3] %*% t(ykt[,,4])),
      sum(ykt[,,4] %*% t(ykt[,,1])), sum(ykt[,,4] %*% t(ykt[,,2])), sum(ykt[,,4] %*% t(ykt[,,3])), sum(ykt[,,4] %*% t(ykt[,,4]))),
      nrow = 4, ncol = 4) +
    diag(10, nrow = 4, ncol = 4) # finish calculating Somega with this identity matrix
  
  # Make sure Somega is positive definite (rounding errors!)
  Somega <- round(Somega, 2)


  # Do I have to go for the density here?
  Omega <- riwish(v = vOmega + num_observations, S = Somega)

  #--------------------------- Step 7: Draw V ----------------------------------
  # Generate all the matrices for the scale matrix for the Inverse Wishart
  # distribution.
  # ? that's not correct, dimension will have to be 4x4. I am deliberately
  # changing the transpose again to get a 4x4 matrix. How *** did they write down
  # their vectors and matrices?
  
  ## Initialize here because dimensions are changing
  IWbeta <- array(dim = c(4,4,k))
  IWtheta <- array(dim = c(4,4,k))
  IWgamma <- array(dim = c(2,2,k))
  IWalpha <- array(dim = c(2,2,k))
  
  for (i in 1:k){
    IWbeta[,,i] <- (betak[,i] - deltaBetaz[,i]) %*% t(betak[,i] - deltaBetaz[,i]) + diag(10, nrow = 4, ncol = 4)
    IWtheta[,,i] <- (thetak[,i] - deltaThetaz[,i]) %*% t(thetak[,i] - deltaThetaz[,i]) + diag(10, nrow = 4, ncol = 4)
    IWgamma[,,i] <- (gammak[,i] - deltaGammaz[,i]) %*% t(gammak[,i] - deltaGammaz[,i]) + diag(10, nrow = 2, ncol = 2)
    IWalpha[,,i] <- (alphak[,i] - deltaAlphaz[,i]) %*% t(alphak[,i] - deltaAlphaz[,i]) + diag(10, nrow = 2, ncol = 2)
  }
  
  # then collapse by k
  IWbeta <- rowSums(IWbeta, dims = 2)
  IWtheta <- rowSums(IWtheta, dims = 2)
  IWgamma <- rowSums(IWgamma, dims = 2)
  IWalpha <- rowSums(IWalpha, dims = 2)
  
  Vbeta = riwish(v = vOmega + k, IWbeta)
  Vtheta = riwish(v = vOmega + k, IWtheta)
  Vgamma = riwish(v = vOmega + k, IWgamma)
  Valpha = riwish(v = vOmega + k, IWalpha)
  
  #-------------------------- Step 8: Draw Delta -------------------------------
  # Does the matrix for a mvt normal distribution has to be a certain way?
  # 2x2 is probably not correct
  # Though given differently in the paper, given from the notation, we only need
  # one deltaq.
  deltaq <- solve(solve(z %*% t(z)) + diag(0.1, nrow = 2, ncol = 2))
  # yields a 2x2, so only brand and spec. Apparently, it's the same for every
  # beta0, beta1 etc.
  
  # Disregard deltaBar as it is 0 (or is that just the Initial value?)
  deltaTildeBeta <- deltaq %*% (z %*% t(betak))
  deltaTildeTheta <- deltaq %*% (z %*% t(thetak))
  deltaTildeGamma <- deltaq %*% (z %*% t(gammak))
  deltaTildeAlpha <- deltaq %*% (z %*% t(alphak))
  
  # deltaBeta is a 4x2 matrix, has been initialized in the beginning.
  # Deltaq accounts on the brand/spec level (2x2) and is therefore not indexed.
  for (i in 1:4){
    deltaBeta[i,] <- rmvnorm(1, mean = deltaTildeBeta[,i], sigma = deltaq)
    deltaTheta[i,] <- rmvnorm(1, mean = deltaTildeTheta[,i], sigma = deltaq)
  }
  for (i in 1:2){
    deltaGamma[i,] <- rmvnorm(1, mean = deltaTildeGamma[,i], sigma = deltaq)
    deltaAlpha[i,] <- rmvnorm(1, mean = deltaTildeAlpha[,i], sigma = deltaq)
  }
  
  
  ### Monitor all the jumps for later. If this is running, it'll need the parameters etc. as well.
  Uarray[,,,i] <- uOld
  
  if (g %% 10 == 0){paste0(c("Iteration", g, "done!"), collapse = " ")}
}

beep(3)
#/////////////////////////////////////////////////////////////////////////////
#//////////////////////////// Done with the chain ////////////////////////////
#/////////////////////////////////////////////////////////////////////////////

#------------------------------- Diagnostics ----------------------------------
