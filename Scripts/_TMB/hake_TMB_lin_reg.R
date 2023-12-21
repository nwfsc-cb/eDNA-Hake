library(here)
library(TMB)
TMB::compile(here("src","LinRegMat.cpp"))

dyn.load(dynlib(here("src","LinRegMat")))

data <- list(Y_i = Y_i, X_f=X_f)
parameters <- list(beta=rnorm(ncol_beta,0,0.1), lnSigma=0)
obj <- MakeADFun(data, parameters, DLL="LinRegMat")



obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj)


###########################


library(here)
library(TMB)
setwd("../src")
TMB::compile("RegSmooth3.cpp")
dyn.load(dynlib("RegSmooth3"))


n_bs     <- ncol(SM$Xs)
b_smooth_start <- SM$b_smooth_start
n_smooth <- length(b_smooth_start)
b_smooth <- if (SM$has_smooths) rep(0,sum(SM$sm_dims)) else array(0) 
has_smooths <- SM$has_smooths

tmb_data <- list(Y_i = Y_i, 
                 Xf =  Xf,
                 Xs = SM$Xs,
                 Zs = SM$Zs,
                 n_smooth = n_smooth,
                 has_smooths = as.integer(has_smooths),
                 b_smooth_start = b_smooth_start)

tmb_params <- list(# Regression terms
                   betaf=rnorm(ncol_beta,0,0.1), 
                   # smooth terms
                   bs=rep(0,n_bs),
                   b_smooth = b_smooth,
                   ln_smooth_sigma = rep(0,n_smooth),
                   lnSigma=0)
                   
                   
tmb_random <- c("bs","b_smooth")                   
                   
                 
obj <- MakeADFun(data=tmb_data, 
                 parameters=tmb_params, 
                 random=tmb_random, 
                 DLL="RegSmooth3")
dyn.unload(dynlib("RegSmooth3"))


obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj)

