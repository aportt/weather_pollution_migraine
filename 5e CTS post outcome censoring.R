
#Title: Post-outcome censoring

################################################################################
# Updated version of the R code for the simulation study in:
#
#  "The case time series analysis"
#  Antonio Gasparrini
#  http://www.ag-myresearch.com/2021_gasparrini_epidemiol.html
#
# * The case time series analysis code is available at:
#   https://github.com/gasparrini/CaseTimeSeries
################################################################################

################################################################################
# SIMULATED SCENARIO: POST-OUTCOME CENSORED PERIOD
################################################################################
# Authors: Antonio Gasparrini, Andrea Portt
################################################################################

#This code has been revised from the 2021 publication to account for outcome-dependent follow-up to 
#simulate a post-outcome censored period. 
#It is based on Scenario 12 from the Case time series 2021 publication. 

#This is done by setting the outcome variable to NA in the period 1-10 days after an event.

###################################

# PACKAGES
library(splines)
library(gnm); library(survival)
library(data.table) ; library(scales)
library(foreach) ; library(doParallel)

# SET SEED
set.seed(13041975)

# SET SIMULATION SETTINGS
nsim <- 500

# SET PARAMETERS: NUMBER OF SUBJECTS, FOLLOW-UP, TIME VARIABLES
n <- 500
dstart <- as.Date("2019-01-01")
dend <- as.Date("2019-12-31")
date <- seq(dstart, dend, by=1)

################################################################################
# PREPARE THE DATA

# GENERATE DATASET
data <- data.table(id = factor(rep(seq(n), each=length(date))), date = date)

# TRUE EFFECT OF EXPOSURE EPISODE
# AP: Strongest point estimate is around 1.09 for ozone, so this is a healthy 
# over-estimation of risk, which will give a generous/cautious estimate 
#of bias
beta <- log(1.15)

################################################################################
# EXAMPLE

# SIMULATE EXPOSURE EPISODES (10% OF TIMES)
nexp <- round(length(date)*0.1)
data[, x := seq(length(date)) %in% sample(length(date), nexp) + 0, by=id]

# DEFINE EXPOSURE HISTORY AND THE CUMULATED EXPOSURE OVER LAG TIMES 0-3
exphist <- as.matrix(data[, shift(x, 0:3, fill=0), by=id][,2:5])
data$cumexp <- rowSums(exphist)

# COMPUTE RISK ASSOCIATED WITH EXPOSURE
data$rrexp <- exp(data$cumexp*beta)

# SIMULATE THE OUTCOME (5-20 EVENTS)
# per individual
data[, y := rmultinom(1, sample(5:20, 1), prob=rrexp), by=id]

# SIMULATE CENSORED PERIOD AFTER OUTCOME
#change so that there are only 2 days of censoring
ind <- rowSums(data[, shift(y, 1:10, fill=0), by=id][,2:3]) > 0
data$y[ind] <- NA

# PLOT
plot(date, date, ylim=c(0.5,10+0.5), yaxt="n", ylab="", xlab="Follow-up",
  frame.plot=F, main="Post-outcome censoring")
axis(1, at=range(date)+c(-500,500), labels=F, lwd.ticks=0)
axis(2, at=rev(seq(10)), labels=paste("Sub",seq(10)), cex.axis=0.8, lwd=0, las=1)
for(i in seq(10)) {
  sub <- subset(data, id==rev(seq(10))[i])
  lines(date, rep(i, length(date)), col=grey(0.9))
  lines(date, ifelse(is.na(sub$y), NA, i))
  for(t in seq(length(sub$date))[!is.na(sub$y) & sub$y>0]) {
    points(rep(sub$date[t], sub$y[t]), i+(seq(sub$y[t])-mean(seq(sub$y[t])))/4,
      pch=1, cex=1.5, bg=2)
  }
}

# RUN THE MODEL
mod <- gnm(y ~ cumexp, eliminate=id, data=data, family="poisson")

# ESTIMATED EFFECT
beta ; coef(mod)[[1]]

################################################################################
# SIMULATIONS

# SET THE PARALLELIZATION
ncores <- detectCores()
cl <- makeCluster(max(1,ncores-2))
registerDoParallel(cl)

# ITERATIONS USING FOREACH
pkg <- c("splines","gnm","survival","data.table")
eff <- foreach(i=seq(nsim), .packages=pkg, .combine=rbind) %dopar% {
  
  # SET SEED
  set.seed(13041975+i)
  
  # SIMULATE EXPOSURE, RISKS, OUTCOME (AS ABOVE)
  data[, x := seq(length(date)) %in% sample(length(date), nexp) + 0, by=id]
  exphist <- as.matrix(data[, shift(x, 0:3, fill=0), by=id][,2:5])
  data$cumexp <- rowSums(exphist)
  data$rrexp <- exp(data$cumexp*beta)
  data[, y := rmultinom(1, sample(5:20, 1), prob=rrexp), by=id]
  
  # SIMULATE CENSORED PERIOD AFTER OUTCOME
  ind <- rowSums(data[, shift(y, 1:10, fill=0), by=id][,2:3]) > 0
  data$y[ind] <- NA
  
  # RUN THE MODEL AND EXTRACT COEF-VCOV-CI
  mod <- gnm(y ~ cumexp, eliminate=id, data=data, family="poisson")
  coef <- coef(mod)[[1]]
  vcov <- as.matrix(vcov(mod))[1,1]
  cilow <- (coef-qnorm(0.975)*sqrt(vcov)) 
  cihigh <- (coef+qnorm(0.975)*sqrt(vcov))
  
  # STORE THE RESULTS
  cbind(est = coef, cov = cilow < beta & cihigh > beta)
}

# STOP CLUSTER
stopCluster(cl)

################################################################################
# RESULTS

# RELATIVE BIAS (%)
(bias <- (mean(eff[,1])-beta)/beta * 100)
#0.05624837


# COVERAGE
# ratio of times the CI overlaps the true value
(cov <- mean(eff[,2]))
#0.938

# RELATIVE RMSE (%) 
#root mean square error - standard deviation of the residuals (prediction errors). 
# How concentrated the data points are around the regression line
(rmse <- sqrt(mean((eff[,1]-beta)^2))/beta * 100)
# 16.13292
