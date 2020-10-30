
##########################################################################################
# This R program runs a Monte Carlo study of meta-analysis methods.
#
#   We compare the conventional Schmidt-Hunter credibility interval, 
#    prediction interval (t distribution),
#    bootstrap tolerance interval is omitted in this version, 
#    and analytic tolerance interval (Howe 1969 k factor),
#    but we added the Mee (1984) tolerance interval
#    for coverage of the 80% interval at alpha = .05
#    Looking to see coverage values at 95 percent for the
#    80 percent interval.  All based on S&H (2015) bare bones, ES = d.
#
#   We also compute captured proportions of the underlying population 
#    (analogous to bias and efficiency).
#
#   The paramater set is...delta, SD.delta, number of studies(k), 
#    and sample sizes (average N). For one run, set delta & SD delta.  Rest are all included.
#    We follow the simulation design used by
#    Sanchez-Meca & Marin-Martinez (2008).  Confidence intervals for the overall effect size
#    in random-effects meta-analysis.  Psychological Methods, 13, 31-48.
##########################################################################################
#############################Clear Console##################################
rm(list=ls()) # cleans (empties)global environment
dev.off()     # cleans plots/graphics
cat("\014")   # cleans out the console
############################################################################
# If not installed, install packages for the libraries below
#options(scipen=999) # this suppresses scientific notation (exponential notation)
# install.packages('psych')
# install.packages('metafor')
# install.packages('beepr')
# install.packages('boot')
library(psych)      # for descriptive stats
# library(metafor)    # for meta-analysis with REML or D
library(beepr)      # sound at end to know MC is done
# library(boot)       # for bootstrap tolerance interval



###########################################################################
# set parameters for the simulation
###########################################################################
# number of meta-analysis (trials) in the simulation run 
# numbermetas = counter for the main simulation loop for Monte Carlo
numbermetas <- 2000 # for actual Monte Carlo, set to 2000 (two thousand)
# delta values at .2, .5 and .8 (pick 1)
delta <- .8           # underlying mean of infinite-sample effect sizes
# tau-squared = .04, .08, .16, .32 (use values of tau in the simulation)
# tau = .2, .283, .4, .566 (tau is SD.delta; pick 1)
SD.delta <-.566               # underlying SD of infinite-sample effect sizes
 BootIter <- 1000             # number of iteractions for the bootstrap
# The program will cycle through all the combinations of k and Nbar in one run,
#  so you don't have to pick those.
# Compute actual lower and upper bounds of credibility interval
CR80L.pop <- delta-1.282*SD.delta  # population 80% credibility interval lower bound
CR80U.pop <- delta+1.282*SD.delta  # upper bound
## sample sizes
ssizes <- 1:60   #  average sample size: 30, to 100 by 10; 100 to 500 by 100
ssizes <- matrix(ssizes,12,5)      # shape the matrix
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
ctrl.parms <- cbind(numbermetas, delta, SD.delta, BootIter)
Sys.time() # to print the starting time
################################################################################
kfactor <-function(alpha = .05, P = .8, n = 100) # define function and arguments
{
  u <- qnorm((1+P)/2)*sqrt(1+1/n) # Howe 1969 calc u
  chival <- qchisq(alpha, (n-1))  # Howe 1969 needed chi-square value
  v <- sqrt((n-1)/chival)         # Howe 1969 calc v
  w <- sqrt(1+(n-3-chival)/(2*(n+1)^2)) # Howe 1969 calc w
  k.factor <- u*v*w               # Howe 1969 k factor
}
################################################################################
# Tolerance items
P <- .8
alpha <- .05
###############################################################################
iter <-1 # start counter for output matrix
Output <- matrix(-9, 120, 17) # placeholder for results set to -9
Output2 <- matrix(-9, 120, 8)
#########################################################################
##        begin the simulation
# number of samples loop
for (ksam in 1:10){
  numbersamples <- ksam*10  # set k, the number of effect sizes
# sample sizes loop
  for (nsam in 1:12) {
    ssizes1 <- ssizes[nsam, ]
    avesamsize <- mean(ssizes1) # set Nbar, average sample size per study
#
Nreps <- numbersamples/5        # Nreps is the number of times each set of sample sizes (ssizes) is repeated
#
# set placeholders for the simuation to hold the results

Results <- matrix(0,numbermetas,21) #placeholder for results - 
#  replications by the number of desired items 


############################################################################
# MAIN Loop - run as many times as you want replications
# set/reset placeholders for one meta
for (j in 1:numbermetas){
out1 <- 1:numbersamples
n1all<- 1:numbersamples
n2all<- 1:numbersamples
# Compute 1 meta-analysis Loop - compute study data and meta-analyze
sam=0
for (rep in 1:Nreps) {
for (i in 1:5){  
# set the local study value of delta, sampled from the parameter set
delta.i <-rnorm(1,delta,SD.delta)
sam = sam+1
n1 <-ssizes1[i]/2
n2 <-n1
sample1 <- rnorm(n1,delta.i,1) #experimental group data
sample2 <- rnorm(n2,0,1)       #control group data
df1 <- n1-1
df2 <- n2-1
Ms1 = mean(sample1)
Vs1 =var(sample1)
Ms2 = mean(sample2)
Vs2 = var(sample2)
d = (Ms1-Ms2)/sqrt((Vs1*df1+Vs2*df2)/(df1+df2)) # compute local value of d
#output d and sample size for each group
out1[sam] <- d
n1all[sam] <- n1
n2all[sam] <- n2
} # end meta-analysis loop (inner loop)
} # end reps loop
  # we have sampled data for 'numbersamples' or k studies
meta <- cbind(out1,n1all,n2all) # meta has d, n1 and n2 - what we need for a meta-analysis
k <- length(out1)               # number of effect sizes (studies)
#
Vdi <- ((meta[,2] + meta[,3]) / (meta[,2] * meta[,3])) + (meta[,1]^2 / (2 * (meta[,2] + meta[,3])))
gi <- (1 - (3 / (4 * (meta[,2]+meta[,3] - 3)-1))) * meta[,1]
Vgi <- (1 - (3 / (4 * (meta[,2]+meta[,3] - 3)-1)))^2 * Vdi

mod.data <- data.frame(gi, Vgi, n1all, n2all)

################################################Hunter & Schmidt####################################
# compute the H&S meta-analysis in d (use g, unbiasd estimator of standarized mean difference)
#
HSgi <- gi   # unbiased estimates of standardized mean difference computed above
HSVgi <- Vgi # unbiased variance computed above
HSwi <- n1all+n2all # total N
HSwigi <- HSwi*HSgi # weighted ES
HSsumwigi <- sum(HSwigi) # sum weighted ES
HSsumwi <- sum(HSwi)     # sum weights
HSmeang <- HSsumwigi/HSsumwi # mean d; dbar
Nbar <- sum(HSwi)/k          # mean N; nbar
HSdevg <- HSgi-HSmeang       # deviations from the mean
HSdevgsq <- HSdevg*HSdevg    # squared deviations
HSNdevgsq <- (HSwi)*HSdevgsq # weighted, squared deviations
HSsumNdevgsq <- sum(HSNdevgsq) # sum the weighted, squared deviations
HSVg <-HSsumNdevgsq/HSsumwi  # variance (SS/Sw)
HSVgk <- (k)/(k-1)*HSVg      # k corrected variance
HSVser <-((Nbar-1)/(Nbar-3))*((4/Nbar)*(1+HSmeang/8)) # sampling error
HStausq <- HSVg-HSVser       # REVC
HStausqk <- HSVgk-HSVser     # k corrected REVC
if(HStausq <= 0) {HStausq <- .0000001} # keep from zero
if(HStausqk <= 0) {HStausqk <- .0000001}
HStau <- HStausq^.5
HStauk <- HStausqk^.5
HSVg <- HSVg/numbersamples  # needed variance of the mean calculation
HSVgksqrt <- HSVgk^.5       # standard error of the mean
# end of one meta-analysis in d
#########################Compute HS Intervals##########################
#  Conventional credibility interval
HSl80 <- HSmeang - 1.282 * HStau  #HS lower bound
HSu80 <- HSmeang + 1.282 * HStau  #HS upper bound
# Percentage of population captured by conventional credibility interval
pctHS <- pnorm(HSu80, delta, SD.delta)-pnorm(HSl80, delta, SD.delta)
# k corrected version - omit during ORM revision
#HSl80k <- HSmeang - 1.282 * HStauk #HS-k lower bound
#HSu80k <- HSmeang + 1.282 * HStauk #HS-k upper bound
#pctHSk <- pnorm(HSu80k, delta, SD.delta)-pnorm(HSl80k, delta, SD.delta)
###################################################
#  Prediction Interval (Higgins interval with S&H values)
res.s <- sqrt(HStausq+HSVg)
HSPILB <- HSmeang + qt(.1, df=(k-2))*res.s
HSPIUB <- HSmeang - qt(.1, df=(k-2))*res.s
pctHSPI <- pnorm(HSPIUB, delta, SD.delta)-pnorm(HSPILB, delta, SD.delta)
#######################################################
# Howe (1969) tolerance interval
k.factor <- kfactor(n=k, alpha=.05)               # Howe 1969 k factor
TILB <- HSmeang-k.factor*HStau
TIUB <- HSmeang+k.factor*HStau
pctATI <- pnorm(TIUB, delta, SD.delta)-pnorm(TILB, delta, SD.delta)
#######################################################
# Mee(1984) tolerance interval
Rhat <- HSVser/HStausq # sampling variance/REVC
if(Rhat > 5) {Rhat <- 5} # maximum value of R hat if REVC approaches zero
nstar <- k/(1+Rhat)^.5       # find adjusted n
s.sqr <- HSVgk # adjust for S-H omission of n-1 in the denominator
lambda <- kfactor(n=nstar, alpha = .05)
lamstar <- lambda/(1+Rhat)^.5
MeeL <- HSmeang-lamstar*sqrt(s.sqr)
MeeU <- HSmeang+lamstar*sqrt(s.sqr)
pctMee <- pnorm(MeeU, delta, SD.delta)-pnorm(MeeL, delta, SD.delta)
####################################################################
# Compute whether interval contains parameters of interest for one meta
#  (coverage instance) for all intervals
newstudy <- rnorm(1, delta, SD.delta)
# conventional credibility interval
CR.cov <-0
if(pctHS >= .8) {CR.cov<-1}
 #HSk.cov <-0
 #if(pctHSk >= .8) {HSk.cov<-1}
# Mee tolerance interval
Meetol.cov <-0
if(pctMee >=.8) {Meetol.cov <- 1}
# analytic tolerance interval
ATI.cov <-0
if(pctATI >= .8) {ATI.cov<-1}
# prediction interval
PI.cov <-0
if(HSPILB<= newstudy & newstudy <= HSPIUB) {PI.cov<-1}

##########################################
# Output results
##########################################
Results[j, ] <- cbind(delta, SD.delta, CR80L.pop, CR80U.pop, # 1, 2, 3, 4
                      HSl80, HSu80, MeeL, MeeU, # 5, 6, 7, 8
                      HSPILB, HSPIUB, TILB, TIUB,     # 9, 10, 11, 12
                      pctHS, pctMee, pctHSPI, pctATI, # 13, 14, 15, 
                      CR.cov, Meetol.cov, ATI.cov, PI.cov, Rhat)  # 17, 18, 19, 20
} # END outer loop.

# Rename the columns of the Results
resultsnames = c("delta", "SD.delta", "true lower bound", "true upper bound", # 1, 2, 3, 4                        #1, 2, 3
                 "CRLB", "CRUB",  "MeeLB", "MeeUB",   # 5, 6, 7, 8
                 "HSPILB", "HSPIUB", "AtolLB", "AtolUB",         # 9, 10, 11, 12
                 "pctCR", "PctBtol", "pctPI", "pctAtol",       # 13, 14, 15, 16
                 "CR.cov", "MeeTol.cov", "ATI.cov", "PI.cov", "Rhat")  # 17, 18, 19, 20
colnames(Results) = resultsnames
# Results # Show data from MC
#
# Calculate the biases for the 3 estimators.
CR.Mean.pct <- sum(Results[,13])/numbermetas
#HS.K.Mean.pct <-sum(Results[,18])/numbermetas
Meetol.Mean.pct <- sum(Results[,14])/numbermetas
Atol.Mean.pct <-sum(Results[,15])/numbermetas
PI.Mean.pct <- sum(Results[,16])/numbermetas

# Calculate RMSE for each estimator of the lower bound
CR.RMSE <-     sqrt(sum((Results[,13]-.8)^2)/numbermetas)
#HS.K.RMSE <-   sqrt(sum((Results[,18]-.8)^2)/numbermetas)
Meetol.RMSE <- sqrt(sum((Results[,14]-.8)^2)/numbermetas)
Atol.RMSE <-     sqrt(sum((Results[,15]-.8)^2)/numbermetas)
PI.RMSE <-   sqrt(sum((Results[,16]-.8)^2)/numbermetas)

# Mean coverage
CR.covm <-    sum(Results[, 17])/numbermetas
#HSk.covm <-   sum(Results[, 21])/numbermetas
Meetol.covm <- sum(Results[, 18])/numbermetas
Atol.covm <-    sum(Results[, 19])/numbermetas
PI.covm <-  sum(Results[, 20])/numbermetas

# calculate the mean upper and lower bound estimates for each estimator
HS.lb.m <- sum(Results[,5])/numbermetas
HS.ub.m <- sum(Results[,6])/numbermetas
#HS.K.lb.m <-sum(Results[,11])/numbermetas
#HS.tol.ub.m <- sum(Results[,12])/numbermetas
Meetol.lb.m <-sum(Results[,7])/numbermetas
Meetol.ub.m <- sum(Results[,8])/numbermetas
Atol.lb.m <-sum(Results[,11])/numbermetas
Atol.ub.m <- sum(Results[,12])/numbermetas
PI.lb.m <-sum(Results[,9])/numbermetas
PI.ub.m <- sum(Results[,10])/numbermetas

# Create database for output to collect across runs
Output[iter, ] <- c(delta, SD.delta, avesamsize, numbersamples, numbermetas,   #1 - 5
                    CR.Mean.pct, Meetol.Mean.pct, Atol.Mean.pct, PI.Mean.pct,  #6 - 9
                     CR.RMSE, Meetol.RMSE, Atol.RMSE, PI.RMSE,               #10-13
                     CR.covm, Meetol.covm, Atol.covm, PI.covm)                  #14-17
Output2[iter, ] <- c(HS.lb.m,HS.ub.m,
                   Meetol.lb.m, Meetol.ub.m, 
                   Atol.lb.m, Atol.ub.m,
                   PI.lb.m, PI.ub.m)
#
iter <- iter+1                     
} # end sample sizes loop
} # end number samples loop
######################################################################


outnames = c("delta", "SD.delta", "avesamsize", "numbersamples", "numbermetas",     # 1- 5
             "CR.Mean.pct", "Meetol.Mean.pct", "Atol.Mean.pct", "PI.Mean.pct",        # 6 - 9
             "CR.RMSE", "Meetol.RMSE", "Atol.RMSE", "PI.RMSE",                        # 10 - 13
             "CR.covm", "Meetol.covm", "Atol.covm", "PI.covm")                        # 14 - 17
colnames(Output) = outnames
out2names = c("CRlow", "CRhi", "Meelow", "Meehi", "Atollow",
             "Atolhi", "PIlow", "PIhi")
colnames(Output2) = out2names
Output
describe(Output)
Output2
describe(Output2)
#Quick Plot to Check if Any Crashes Occured, Observations of -9 indicate crash
# mean captuer
#boxplot(Output[, 6], Output[, 7], Output[, 8], Output[, 9],
#        ylim = c(0, 1.1), ylab='Percent', main='Mean Capture Pct')
#abline(h=.8)
#abline(h=.9, lty=2)
#abline(h=.7, lty=2)
#text(1,1.05, "HS")
#text(2,1.05, "HS.k")
#text(3,1.05, "Tol")
#text(4,1.05, "Hig.t")
#text(5,1.05, "Hig.t")
#mtext(side=3, paste("delta=", delta, "SDdelta=", SD.delta))
#
# Coverage
boxplot(Output[, 14], Output[, 15], Output[, 16], Output[, 17], 
        ylim = c(0, 1.1), ylab='Percent', main='Coverage')
abline(h=.95)
abline(h=.9)
abline(h=.8)
abline(h=.5)
text(1,1.05, "CR")
text(2,1.05, "MeeTol")
text(3,1.05, "AnTol")
text(4,1.05, "PredInt")
mtext(side=3, paste("delta=", delta, "SDdelta=", SD.delta))
#
describe(Output)
ctrl.parms
Sys.time() # to print the ending time - to see how long it took to run the simulation

#outall <- data.frame(Output, Output2)
#describeBy(outall$HStol.covm, group = outall$avesamsize)
beep(3)
# Coverage
# THIS WILL EXPORT OUTPUT TO SHARED O DRIVE
# Make sure to set col.names to TRUE for one condition so headers are added just once
  write.table(Output, file = "Outd1.csv", row.names=F, append=T, col.names=T, sep=",")
#write.table(Output, file = "O:/CAS-PSY/output.csv", row.names=F, append=T, col.names=F, sep=",")

# THIS WILL EXPORT RESULTS TO SHARED O DRIVE 
# NEED TO CHANGE RESULTS TO CONDITION NUMBER; i.e condition saved as results1.csv 
# write.table(Results, file = "results5_20.csv", row.names=F, append=T, col.names=T, sep=",")
  