#### Combine and Post-process stan results.


library(rstan)
library(dplyr)
library(abind)
library(loo)

results.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Stan Model Fits"
setwd(results.dir)

# Get the remainder of the information for the stanfit (priors, raw data, etc.)

load("qPCR 2019 hake lat.long.smooth 7_12_fix_nu_A Base_Var Fitted.RData")
stan_pars <- Output.qpcr$samp %>% names()
e1 <- rstan::extract(Output.qpcr$stanMod,pars=stan_pars, permuted = FALSE)
f1 <- rstan::extract(Output.qpcr$stanMod,pars=stan_pars, permuted = TRUE)
g1 <- get_sampler_params(Output.qpcr$stanMod)

load("qPCR 2019 hake lat.long.smooth 7_12_fix_nu_B Base_Var Fitted.RData")
e2 <- rstan::extract(Output.qpcr$stanMod,pars=stan_pars, permuted = FALSE)
f2 <- rstan::extract(Output.qpcr$stanMod,pars=stan_pars, permuted = TRUE)
g2 <- get_sampler_params(Output.qpcr$stanMod)

load("qPCR 2019 hake lat.long.smooth 7_12_fix_nu_C Base_Var Fitted.RData")
e3 <- rstan::extract(Output.qpcr$stanMod,pars=stan_pars, permuted = FALSE)
f3 <- rstan::extract(Output.qpcr$stanMod,pars=stan_pars, permuted = TRUE)
g3 <- get_sampler_params(Output.qpcr$stanMod)

load("qPCR 2019 hake lat.long.smooth 7_12_fix_nu_D Base_Var Fitted.RData")
e4 <- rstan::extract(Output.qpcr$stanMod,pars=stan_pars, permuted = FALSE)
f4 <- rstan::extract(Output.qpcr$stanMod,pars=stan_pars, permuted = TRUE)
g4 <- get_sampler_params(Output.qpcr$stanMod)

# load("qPCR 2019 hake lat.long.smooth 8_16_fix_nu_E Base_Var Fitted.RData")
# e5 <- rstan::extract(Output.qpcr$stanMod,pars=stan_pars, permuted = FALSE)
# f5 <- rstan::extract(Output.qpcr$stanMod,pars=stan_pars, permuted = TRUE)
# g5 <- get_sampler_params(Output.qpcr$stanMod)


### Identify the parameters that were supposed to be kept by the stanfit.

## Extract only the stan pars that are important 
# e1 <- rstan::extract(stanMod1,pars=stan_pars, permuted = FALSE)
# e2 <- rstan::extract(stanMod2,pars=stan_pars, permuted = FALSE)
# e3 <- rstan::extract(stanMod3,pars=stan_pars, permuted = FALSE)
# e4 <- rstan::extract(stanMod4,pars=stan_pars, permuted = FALSE)

out.array <- array(dim=c(dim(e1)[1],4,dim(e1)[3]))
out.array[,1,] <- e1
out.array[,2,] <- e2
out.array[,3,] <- e3
out.array[,4,] <- e4
# out.array[,5,] <- e5

## Convergence diagnostics
mon <- rstan::monitor(out.array, print = TRUE, warmup = 0)
rownames(mon) <- dimnames(e1)[3]$parameters
mon <- as.data.frame(mon)

pars <- list()

for(i in 1:length(stan_pars)){
  
  A <- f1[[stan_pars[i]]]
  B <- f2[[stan_pars[i]]]
  C <- f3[[stan_pars[i]]]
  D <- f4[[stan_pars[i]]]
  # E <- f5[[stan_pars[i]]]
  
  temp <- abind(A,B,along=1)
  temp <- abind(temp,C,along=1)
  temp <- abind(temp,D,along=1)
  # temp <- abind(temp,E,along=1)
  
  pars[[stan_pars[i]]] <- temp
}

# get_adaptation_info(stanMod)
samp_params <-  list(g1,g2,g3,g4)#,g5)
#get_sampler_params(stanMod4))

# Update LOO
n.chains <- 4
chain_id <- c(rep(1,nrow(pars$beta_0)/n.chains),
              rep(2,nrow(pars$beta_0)/n.chains),
              rep(3,nrow(pars$beta_0)/n.chains),
              rep(4,nrow(pars$beta_0)/n.chains))
              #rep(5,nrow(pars$beta_0)/n.chains))
log_lik_1 <- pars[["log_lik"]] %>% as.matrix()
r_eff <- relative_eff(exp(log_lik_1),chain_id = chain_id)
loo_val <- loo(log_lik_1, r_eff = r_eff)


Output.qpcr$loo_val_all <- loo_val
Output.qpcr$pars <- pars
Output.qpcr$converge = samp_params
Output.qpcr$stanMod_summary_all=mon

NAME <- "qPCR 2019 hake lat.long.smooth 7_12_fix_nu_ALL Base_Var Fitted"
setwd(results.dir)
save(Output.qpcr,file=paste(NAME,".RData",sep=""))


