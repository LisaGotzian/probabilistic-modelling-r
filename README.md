# Probabilistic modelling using JAGS, BUGS and the Metropolis algorithm in R
This project follows the bayesian hierarchical model by [Agarwal et al., 2015: Do organic results help or hurt sponsored search performance?](http://dx.doi.org/10.1287/isre.2015.0593), including its supplementary material, in R using JAGS, BUGS and a bottom-up approach.
`probabilistic_modelling_agarwal.pdf` is a complete documentation and description of the project including inline code. (I urge you, read it. It'll be worth the read!)

## Data simulation
0. As a first step, we simulate the data and parameters in the model by Agarwal et al. 2015 in `0_data_simulation.R`.

## MCMC using JAGS, BUGS and the raw Metropolis-Hastings algorithm
1. `JAGS_attempt.R` is a formulation of the model in JAGS. As it is recursive, JAGS doesn't work.
2. `BUGS_attempt.R` is a formulation of the same model in BUGS. BUGS allows for recursive models, but can lead to stackoverflows, as happened here.
3. `Metropolist_agarwal.R` is a (loose) interpretation of the authors' appendix as this is written vaguely. It is the raw Metropolos-Hastings algorithm without any libraries used.
4. `JAGS_running.R` is a running model in JAGS which has been simplified to not be recursive anymore.

We highly value clarification by the authors on their notation and model formulation.

## Built with
* MCMCpack
* JAGS from [here](https://sites.google.com/site/doingbayesiandataanalysis/software-installation), including Kruschke's DBDA2E-utilities
* R2OpenBUGS
* rjags
