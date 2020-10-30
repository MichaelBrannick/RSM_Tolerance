
##########################################################################################
# This R program runs a Monte Carlo study of meta-analysis methods.
#
#   We compare credibility and related interval estimates by  
#    (a) Schmidt-Hunter credibility interval without k correction (HS),
#    (b) Prediction interval(base on SH quantities) (PI)
#    (c) analytic tolerance interval from engineering (ATI)
#    (d) eliminated Schmidt-Hunter bootstrapped tolerance interval (BTI) for this run
#    (e) added Mee tolerance interval
#
#  The bootstrap fails on rare occasions.  When it does, all the remaining 
#   estimates are equal to -9, so be sure to check each run and rerun
#   any conditions during and after the failure.
#
##########################################################################################
#  The paramater set is...rho, SD.rho, number of studies(k), 
#   sample sizes (average N).
#  We follow the simulation design used by 
#  Sanchez-Meca & Marin-Martinez (2008).  Confidence intervals for the overall effect size
#   in random-effects meta-analysis.  Psychological Methods, 13, 31-48.
##########################################################################################

############################################################################
# If not installed, install packages for the libraries below
#options(scipen=999) # this suppresses scientific notation (exponential notation)
library(psych)      # for descriptive stats
library(MASS)       # for sampling correlations
library(beepr)      # sound at end to know MC is done
library(boot)       # for bootstrap TI
#############################Clear Console#################################
rm(list=ls()) # cleans (empties) global environment
dev.off()     # cleans plots/graphics
cat("\014")   # cleans out the console
###########################################################################
# set parameters for the simulation
###########################################################################
# number of meta-analysis (trials) in the simulation run 
# numbermetas = counter for the main simulation loop for Monte Carlo
numbermetas <- 2000             # use 2000 for actual study
# rho values at .1, .25 and .4
rho <- .4                      # underlying mean of infinite-sample effect sizes
# tau = .05, .10, .15, .20 (tau is SD.rho)05
SD.rho <-.20                    # underlying SD of infinite-sample effect sizes
BootIter <- 500                # Number of iterations for the bootstrap (use 500 for simulation)
mu <- rep(0,2)                  # multivariate means of zero
rhobox <- matrix(-9,2,2)+10*diag(2) # placeholder for correlation matrix for sampling
# number of independent samples in one meta-analysis (value of k)
# values are set to one of these:   10, to 100 by 10
CR80L.pop <- rho-1.282*SD.rho  # population 80% credibility interval lower bound
CR80U.pop <- rho+1.282*SD.rho  # upper bound
# Set up sample sizes.
ssizes <- 1:60
ssizes <- matrix(ssizes,12,5)
ssizes[1, ] <-c(12,16,18,20,84)    # average N = 30
ssizes[2, ] <- ssizes[1, ]+10      # average N = 40
ssizes[3, ] <-c(32,36,38,40,104)   # ave N = 50
ssizes[4, ] <-ssizes[3, ]+10       # ave N = 60
ssizes[5, ] <-ssizes[4, ]+10       # ave N = 70
ssizes[6, ] <-c(62,66,68,70,134)   # ave N = 80
ssizes[7, ] <-ssizes[6, ]+10       # ave N = 90
ssizes[8, ] <-c(82,86,88,90,154)   # ave N = 100
ssizes[9, ] <- ssizes[8, ]+100      # ave N = 200
ssizes[10, ] <- ssizes[8, ]+200     # ave N = 300
ssizes[11, ] <- ssizes[8, ]+300     # ave N = 400
ssizes[12, ] <- ssizes[8, ]+400     # ave N = 500
# ssizes
pop.r <- rnorm(1000000,rho,SD.rho)              # sample 1 million correlations from population
pop.rs <- sample(pop.r[pop.r > -.99 & pop.r <.99], size = 1000000, replace = T) # resample to be in bounds
#
ctrl.parms <- cbind(numbermetas, rho, SD.rho, BootIter)
Sys.time() # to print the starting time
#
########################################################################
# Tolerance function (engineering)
P <- .8
alpha <- .05
kfactor <-function(alpha = .05, P = .8, n = 100) # define function and arguments
{
  u <- qnorm((1+P)/2)*sqrt(1+1/n) # Howe 1969 calc u
  chival <- qchisq(alpha, (n-1))  # Howe 1969 needed chi-square value
  v <- sqrt((n-1)/chival)         # Howe 1969 calc v
  w <- sqrt(1+(n-3-chival)/(2*(n+1)^2)) # Howe 1969 calc w
  k.factor <- u*v*w               # Howe 1969 k factor
}
##########################################################################
iter <-1 # start
Output <- matrix(-9, 120, 17)
# number of samples loop
# ksam 1 to 10
# nsam 1 to 9
for (ksam in 1:10){
  numbersamples <- ksam*10
# sample sizes loop
  for (nsam in 1:12) {
    ssizes1 <- ssizes[nsam, ]
    avesamsize <- mean(ssizes1) # set Nbar, average sample size per study
#
Nreps <- numbersamples/5           # Nreps is the number of times each set of sample sizes (ssizes) is repeated
#
# set placeholders for the simuation to hold the results

Results <- matrix(-9,numbermetas,20) # placeholder for results - 
                                    # replications by the number of desired items 
############################################################################
# MAIN Loop - run as many times as you want replications
# set/reset placeholders for one meta

for (j in 1:numbermetas){
out1 <- 1:numbersamples
ns <- out1

# Compute 1 meta-analysis Loop - compute study data and meta-analyze
sam=0
for (rep in 1:Nreps) {
for (i in 1:5){  
# set the local study value of rho, sampled from the parameter set
rhobox[2, 1] <-sample(pop.rs, size = 1, replace = TRUE)
rhobox[1, 2] <- rhobox[2, 1]
sam = sam+1
sample.rs <-mvrnorm(ssizes1[i],mu=mu,Sigma=rhobox) # sample ni data for obs r from rho
cor1 <-cor(sample.rs)                              # compute unrestricted correlation from sampled data
                                
out1[sam] <- cor1[2,1]          # output r and sample size for each group
ns[sam] <- ssizes1[i]
}                               # end meta-analysis loop (inner loop)
}                               # end reps loop
                                # we have sampled data for 'numbersamples' or k studies
meta <- cbind(out1,ns)          # meta has r and n - what we need for a meta-analysis
k <- length(out1)               # number of effect sizes (studies)
mod.data <- data.frame(out1, ns)  # data for bootstrap
################################################Hunter & Schmidt####################################
# compute the H&S 'bare bones' meta-analysis in r 
#
ri <- meta[, 1]
ri <- ifelse(ri == 0, 0.0000001, ri)        # avoid dividing by zero later on
ni <- meta[, 2]
nR <- ri*ni
Nsum <- sum(ni)
rbar <- sum(nR)/Nsum
ei <- ((1-rbar^2)^2)/(ni-1)
Dev <- ri-rbar
Wt.Dev <- ni*(Dev*Dev)
V.obs <- sum(Wt.Dev)/Nsum
V.m <- V.obs/k
V.obs.k <- k/(k-1)*V.obs
V.err <- sum(ni*ei)/(Nsum)
V.rho <- V.obs-V.err
if(V.rho < 0) {V.rho = .0001}
Rhat <- V.err/V.rho
if(Rhat >5) {Rhat <- 5}
HStau <- sqrt(V.rho)
V.rho.k <- V.obs.k-V.err
if(V.rho.k < 0) {V.rho.k = 0}
HStauk <- sqrt(V.rho.k)
# compute credibility interval
CRl80 <- rbar - 1.282 * HStau  #HS lower bound
CRu80 <- rbar + 1.282 * HStau  #HS upper bound
pctCR <- pnorm(CRu80, rho, SD.rho)-pnorm(CRl80, rho, SD.rho)
#
#HSl80k <- rbar - 1.282 * HStauk #HS-k lower bound
#HSu80k <- rbar + 1.282 * HStauk #HS-k upper bound
#pctHSk <- pnorm(HSu80k, rho, SD.rho)-pnorm(HSl80k, rho, SD.rho)
###################################################
# Bootstrap Tolerance Interval function for Schmidt-Hunter
#  deleted and replaced by Mee estimator
################################################
  # find prediction interval PI
  res.s <- sqrt(V.m+V.rho)                # PI SD
  
  PILB <- rbar + qt(.1, df=(k-2))*res.s # 80 percent prediction interval
  PIUB <- rbar - qt(.1, df=(k-2))*res.s
  pctPI <- pnorm(PIUB, rho, SD.rho)-pnorm(PILB, rho, SD.rho) # area captured
 # find analytic tolerance interval 
  k.factor <- kfactor(n=k, alpha=.05)               # Howe 1969 k factor
  Atoll80 <- rbar-k.factor*HStau                     # engineering tolerance interval
  Atolu80 <- rbar+k.factor*HStau
  pctATI <- pnorm(Atolu80, rho, SD.rho)-pnorm(Atoll80, rho, SD.rho)
  # find Mee tolerance interval
  nstar <- k/(1+Rhat)^.5       # find adjusted n
  s.sqr <- V.obs*(k/(k-1)) # adjust for S-H omission of n-1 in the denominator
  lambda <- kfactor(n=nstar, alpha = .05)
  lamstar <- lambda/(1+Rhat)^.5
  MeeL <- rbar-lamstar*sqrt(s.sqr)
  MeeU <- rbar+lamstar*sqrt(s.sqr)
  pctMee <- pnorm(MeeU, rho, SD.rho)-pnorm(MeeL, rho, SD.rho)
  
  # the areas captured by pctHigt pctTI are just for reference; they are not the quantity that is of 
  #  interest for the prediction or tolerance interval.  That is computed later in the program.
  
  # Compute coverage for all intervals
newstudy <- rnorm(1, rho, SD.rho)

CR.cov <-0
if(pctCR >= .8) {CR.cov<-1}
#
#HSk.cov <-0
#if(pctHSk >= .8) {HSk.cov<-1}
Meetol.cov <-0
if(pctMee >=.8) {Meetol.cov <- 1}
#
Atol.cov <-0
if(pctATI >= .8) {Atol.cov<-1}
#
PI.cov <-0
if(PILB<= newstudy & newstudy <= PIUB) {PI.cov<-1}

##########################################
# Output results
##########################################
Results[j, ] <- cbind(rho, SD.rho, CR80L.pop, CR80U.pop, # 1, 2, 3, 4
                      CRl80, CRu80, MeeL, MeeU,       # 5, 6, 7, 8
                      Atoll80, Atolu80, PILB, PIUB,      # 9, 10, 11, 12, 
                      pctCR, pctMee, pctATI, pctPI,   # 13, 14, 15, 16
                      CR.cov, Meetol.cov, Atol.cov, PI.cov)  # 17, 18, 19, 20
} # END outer loop.

# Rename the columns of the Results
resultsnames = c("rho", "SD.rho", "true lower bound", "true upper bound", # 1, 2, 3, 4                        #1, 2, 3
                 "CRL80", "CRU80", "MeeLB", "MeeUB",          # 5, 6, 7, 8
                 "Atoll80", "Atolu80", "PILB", "PIUB",      # 9, 10, 11, 12,
                 "pctCR", "pctMee", "pctATI", "pctPI",   # 13, 14, 15, 16
                 "CR.cov", "Mee.cov", "Atol.cov", "PI.cov")  # 17, 18, 19, 20
colnames(Results) = resultsnames
# Results # Show data from MC
#
# Calculate the biases for the 3 estimators.
CR.Mean.pct <- sum(Results[,13])/numbermetas
Mee.Mean.pct <-sum(Results[,14])/numbermetas
Atol.Mean.pct <- sum(Results[,15])/numbermetas
PI.Mean.pct <- sum(Results[,16])/numbermetas
# Calculate RMSE for each estimator of the lower bound
CR.RMSE <-   sqrt(sum((Results[,13]-.8)^2)/numbermetas)
Mee.RMSE <- sqrt(sum((Results[,14]-.8)^2)/numbermetas)
Atol.RMSE <- sqrt(sum((Results[,15]-.8)^2)/numbermetas)
PI.RMSE <- sqrt(sum((Results[,16]-.8)^2)/numbermetas)
# Determine whether estimate is within 10 percent of desired .80
CR.cov <- sum(Results[, 17])/numbermetas
Mee.cov <- sum(Results[, 18])/numbermetas
Atol.cov <- sum(Results[, 19])/numbermetas
PI.cov <- sum(Results[, 20])/numbermetas

# Create database for output to collect across runs
Output[iter, ] <- c(rho, SD.rho, avesamsize, numbersamples, numbermetas,     # 1-5 
                     CR.Mean.pct, Mee.Mean.pct, Atol.Mean.pct, PI.Mean.pct, # 6 - 9
                     CR.RMSE, Mee.RMSE, Atol.RMSE, PI.RMSE,     # 10 - 13
                     CR.cov, Mee.cov, Atol.cov, PI.cov)           # 14 - 17
iter <- iter+1                     
} # end sample sizes loop
} # end number samples loop
outnames = c("rho", "SD.rho", "avesamsize", "numbersamples", "numbermetas", # 1 - 5
             "CR.Mean.pct", "Mee.Mean.pct", "Atol.Mean.pct", "PI.Mean.pct", # 6 - 9
             "CR.RMSE", "Mee.RMSE", "Atol.RMSE", "PI.RMSE",     # 10 - 13
             "CR.cov", "Mee.cov", "Atol.cov", "PI.cov")           # 14 - 17
colnames(Output) = outnames
describe(Output)
beep()

#Quick Plot to Check if Any Crashes Occured, Observations of -9 indicate crash
boxplot(Output[, 14], Output[, 15], Output[, 16], Output[, 17], 
 ylim = c(0, 1.1), ylab='Percent', main='Coverage')
abline(h=.8)
abline(h=.95, lty=2)
text(1,1.05, "CR")
text(2,1.05, "Mee")
text(3,1.05, "Atol")
text(4,1.05, "PI")
mtext(side=3, paste("rho=", rho, "SDrho=", SD.rho))

#



# THIS WILL EXPORT OUTPUT TO SHARED O DRIVE
# Make sure to set col.names to TRUE for one condition so headers are added just once
# write.table(Output, file = "out25_05.csv", row.names=F, append=T, col.names=T, sep=",")
write.table(Output, file = "output1.csv", row.names=F, append=T, col.names=T, sep=",")
# THIS WILL EXPORT RESULTS TO SHARED O DRIVE 
# NEED TO CHANGE RESULTS TO CONDITION NUMBER; i.e condition saved as results1.csv 
# write.table(Results, file = "res25_05.csv", row.names=F, append=T, col.names=T, sep=",")