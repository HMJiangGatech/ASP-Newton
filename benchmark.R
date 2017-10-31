## Notice that the precision have to be fine tunned for getting comparable results
setwd("~/Desktop/2018 MPB Prox-Newton/2018 mpb prox-newton/testcode")
# loading and installing required packages
library(Rcpp)
library(glmnet)
library(picasso)
# library(gcdnet)
# library(spams)

source("scripts.R")
sourceCpp("utils.cpp")

# Experiment parameters
# skip some comparison
skip = c("glmnet")
useRealData = FALSE
# for simulated data set
n = 2000
d = 10000


# Linear Regression
set.seed(111)

# Simulated data
sim_wc <- generate_sim(n=n, d=d, c=0.3, seed=111) 
sim_ic <- generate_sim(n=n, d=d, c=3.0, seed=112)
test_gausnet(sim_wc,skip=skip,trialN = 1,prec=1*1e-4)
test_gausnet(sim_ic,skip=skip,trialN = 1,prec=1*1e-4)

# Real Data
if(useRealData)
{
  load("dataset/linear reg/eyedata.RData")
  eyedata$X[which(is.na(eyedata$X))] <- 0 # missing values
  eyedata$X <- eyedata$X[ ,find_nonconstant_column(eyedata$X)]
  eyedata$X <- scale(eyedata$X)
  
  load("dataset/linear reg/DrivFace.RData")
  DrivFace$X[which(is.na(DrivFace$X))] <- 0 # missing values
  DrivFace$X <- DrivFace$X[ ,find_nonconstant_column(DrivFace$X)]
  DrivFace$X <- scale(DrivFace$X)
  test_gausnet(eyedata,skip=skip)
  test_gausnet(DrivFace,skip=skip)
}


#-----------------------------------------------------------------------------

# Logistic Regression
set.seed(111)

# Simulated data
sim_wc <- generate_sim_lognet(n=n, d=d, c=0.3, seed=111)
sim_ic <- generate_sim_lognet(n=n, d=d, c=3.0, seed=112)
test_lognet(sim_wc,skip=skip,trialN = 10,prec=1.0*1e-4)
test_lognet(sim_ic,skip=skip,trialN = 10,prec=1.0*1e-4)


# Real Data
if(useRealData)
{
  load("dataset/logistic reg/madelon.RData")
  madelon$X[which(is.na(madelon$X))] <- 0 # missing values
  madelon$X <- madelon$X[ ,find_nonconstant_column(madelon$X)]
  madelon$X <- scale(madelon$X)
  
  load("dataset/logistic reg/gisette.RData")
  gisette$X[which(is.na(gisette$X))] <- 0 # missing values
  gisette$X <- gisette$X[ ,find_nonconstant_column(gisette$X)]
  gisette$X <- scale(gisette$X)
  
  test_lognet(madelon,skip=skip)
  test_lognet(gisette,skip=skip)
}







#-----------------------------------------------------------------------------

# Poisson Regression
set.seed(111)

# Simulated data
sim_wc <- generate_sim_poi(n=n, d=d, c=0.03, seed=111)
sim_ic <- generate_sim_poi(n=n, d=d, c=0.1, seed=112)
test_poi(sim_wc,skip=skip,trialN = 10,prec=1*1e-4)
test_poi(sim_ic,skip=skip,trialN = 10,prec=1*1e-4)

# Real Data
if(useRealData)
{
  load("dataset/linear reg/DrivFace.RData")
  DrivFace$X[which(is.na(DrivFace$X))] <- 0 # missing values
  DrivFace$X <- DrivFace$X[ ,find_nonconstant_column(DrivFace$X)]
  DrivFace$X <- scale(DrivFace$X)
  DrivFace$Y <- DrivFace$Ang
  
  test_gausnet(DrivFace,skip=skip)
}