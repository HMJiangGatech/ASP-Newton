## Notice that the precision have to be fine tunned for getting comparable results
#setwd("~/Desktop/2018 MPB Prox-Newton/2018 mpb prox-newton/testcode")
# setwd("~/Desktop/ASP-Newton")
# loading and installing required packages
library(Rcpp)
library(glmnet)
library(picasso)
library(gcdnet)
# library(spams)

source("scripts.R")
sourceCpp("utils.cpp")

# Experiment parameters
# skip some comparison
skip = c("glmnet")
useRealData = FALSE
dolinreg = TRUE
dologreg = FALSE
trialN = 10
# for simulated data set
n = 2000
d = 10000


# Linear Regression
if(dolinreg)
{
set.seed(111)

# Simulated data
print('c=0.3')
sim_data <- generate_sim(n=n, d=d, c=0.3, seed=111) 
skip = c("glmnet","gcdnet")    # asp-newton
test_gausnet(sim_data,skip=skip,trialN = trialN,prec=1*1e-4,ratio=0.2, nlambda = 10)
skip = c("glmnet","gcdnet","oldpicasso")    # greedy
test_gausnet(sim_data,skip=skip,trialN = trialN,prec=4.5*1e-4,ratio=0.2, nlambda = 10)
skip = c("glmnet","gcdnet","oldpicasso","cyclic")    # cyclic
test_gausnet(sim_data,skip=skip,trialN = trialN,prec=4.3*1e-4,ratio=0.2, nlambda = 10)
skip = c("picasso","gcdnet")  # glmnet
test_gausnet(sim_data,skip=skip,trialN = trialN,prec=2*1e-6,ratio=0.2, nlambda = 10)
skip = c("picasso","glmnet")   # gcdnet
test_gausnet(sim_data,skip=skip,trialN = trialN,prec=2*1e-5,ratio=0.2, nlambda = 10)
rm(sim_data)
print('c=3.0')
sim_data <- generate_sim(n=n, d=d, c=3.0, seed=112)
skip = c("glmnet","gcdnet")    # asp-newton
test_gausnet(sim_data,skip=skip,trialN = trialN,prec=1*1e-6,ratio=0.2, nlambda = 10)
skip = c("glmnet","gcdnet","oldpicasso")    # greedy
test_gausnet(sim_data,skip=skip,trialN = trialN,prec=2*1e-5,ratio=0.2, nlambda = 10)
skip = c("glmnet","gcdnet","oldpicasso","cyclic")    # cyclic
test_gausnet(sim_data,skip=skip,trialN = trialN,prec=1*1e-5,ratio=0.2, nlambda = 10)
skip = c("picasso","gcdnet")  # glmnet
test_gausnet(sim_data,skip=skip,trialN = trialN,prec=1.5*1e-6,ratio=0.2, nlambda = 10)
skip = c("picasso","glmnet")   # gcdnet
test_gausnet(sim_data,skip=skip,trialN = trialN,prec=3.5*1e-5,ratio=0.2, nlambda = 10)
rm(sim_data)

# Real Data
if(useRealData)
{
  load("DrivFace.RData")
  x=as.matrix(x)
  y=as.matrix(y)
  x=scale(x)
  y=scale(y)
  skip = c("glmnet","gcdnet")    # asp-newton
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=1*1e-4,ratio=0.2, nlambda = 10)
  skip = c("glmnet","gcdnet","oldpicasso")    # greedy
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=5*1e-4,ratio=0.2, nlambda = 10)
  skip = c("glmnet","gcdnet","oldpicasso","cyclic")    # cyclic
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=5*1e-4,ratio=0.2, nlambda = 10)
  skip = c("picasso","gcdnet")  # glmnet
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=1*1e-4,ratio=0.2, nlambda = 10)
  skip = c("picasso","glmnet")   # gcdnet
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=2*1e-5,ratio=0.2, nlambda = 10)
  
  
  load("GHG.RData")
  x=scale(x)
  y=scale(y)
  skip = c("glmnet","gcdnet")    # asp-newton
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=1*1e-5,ratio=0.2, nlambda = 10)
  skip = c("glmnet","gcdnet","oldpicasso")    # greedy
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=8*1e-5,ratio=0.2, nlambda = 10)
  skip = c("glmnet","gcdnet","oldpicasso","cyclic")    # cyclic
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=8*1e-5,ratio=0.2, nlambda = 10)
  skip = c("picasso","gcdnet")  # glmnet
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=2*1e-5,ratio=0.2, nlambda = 10)
  skip = c("picasso","glmnet")   # gcdnet
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=5*1e-6,ratio=0.2, nlambda = 10)
  
  
  load("riboflavin.RData")
  x=scale(x)
  y=scale(y)
  skip = c("glmnet","gcdnet")    # asp-newton
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=1*1e-7,ratio=0.2, nlambda = 10)
  skip = c("glmnet","gcdnet","oldpicasso")    # greedy
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=1*1e-5,ratio=0.2, nlambda = 10)
  skip = c("glmnet","gcdnet","oldpicasso","cyclic")    # cyclic
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=1*1e-8,ratio=0.2, nlambda = 10)
  skip = c("picasso","gcdnet")  # glmnet
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=3*1e-8,ratio=0.2, nlambda = 10)
  skip = c("picasso","glmnet")   # gcdnet
  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=5*1e-9,ratio=0.2, nlambda = 10)
  
  test_gausnet(DrivFace,skip=skip)
}

}
#-----------------------------------------------------------------------------

# Logistic Regression
if(dologreg)
{
set.seed(111)

# Simulated data
print('c=0.3')
sim_data <- generate_sim_lognet(n=n, d=d, c=0.3, seed=111)
test_lognet(sim_data,skip=skip,trialN = trialN,prec=1.0*1e-4,ratio=0.2, nlambda = 10)
rm(sim_data)
print('c=3.0')
sim_data <- generate_sim_lognet(n=n, d=d, c=3.0, seed=112)
test_lognet(sim_data,skip=skip,trialN = trialN,prec=1.0*1e-4,ratio=0.2, nlambda = 10)
rm(sim_data)


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






}
#-----------------------------------------------------------------------------

# Poisson Regression
set.seed(111)

# Simulated data
# sim_data <- generate_sim_poi(n=n, d=d, c=0.03, seed=111)
# test_poi(sim_data,skip=skip,trialN = 1,prec=1*1e-4)
# sim_data <- generate_sim_poi(n=n, d=d, c=0.1, seed=112)
# test_poi(sim_data,skip=skip,trialN = 1,prec=1*1e-4)
# 
# # Real Data
# if(useRealData)
# {
#   load("dataset/linear reg/DrivFace.RData")
#   DrivFace$X[which(is.na(DrivFace$X))] <- 0 # missing values
#   DrivFace$X <- DrivFace$X[ ,find_nonconstant_column(DrivFace$X)]
#   DrivFace$X <- scale(DrivFace$X)
#   DrivFace$Y <- DrivFace$Ang
#   
#   test_gausnet(DrivFace,skip=skip)
# }