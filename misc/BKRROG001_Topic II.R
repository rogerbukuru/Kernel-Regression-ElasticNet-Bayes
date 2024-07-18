# --------------------------------------- Gibbs Sampler -------------------------------------------------------

rm(list = ls())
library(mvtnorm)
library(truncnorm)
library(tidyverse)
rm(list = ls())
diabetes_data = read.csv("PimaIndiansDiabetes2.csv",sep=";", header = T)
N = nrow(diabetes_data)
D = 9
X = cbind(rep(1, N),diabetes_data[,1:8])
X = as.matrix(X)
colnames(X) = c("Intercept",colnames(diabetes_data[,1:8]))
y =ifelse(diabetes_data[,9]=="neg",0,1)
N1 = sum(y) 
N0 = N - N1


gibbs_sampler = function(N_sim){
  

  beta_0 <- rep(0, D)
  Q <- diag(4, D)
  
  
  beta_ <- rep(0, D)
  z <- rep(0, N)

  beta_chain <- matrix(0, nrow = N_sim, ncol = D)
  
  covariance <- solve(Q)
  V <- solve(covariance + crossprod(X, X))
  
  for (i in 2:N_sim) {
    mu_z <- X %*% beta_
    z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1, a = -Inf, b = 0)
    z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1, a = 0, b = Inf)
    M <- V %*% (covariance %*% beta_0 + crossprod(X, z))
    beta_ <- c(rmvnorm(1, M, V))
    beta_chain[i, ] <- beta_
  }
  return(beta_chain)
}


N_sim <- 10000 
beta_chain   = gibbs_sampler(10000)
beta_chain_2 = gibbs_sampler(10000)
beta_chain_3 = gibbs_sampler(10000)
burn_in     = 5000



#--------------------------------------------- Metropolis-Hastings Alogorithm----------------------------------


likelihood <- function(beta,X, m_i){

  Phi <- pnorm(X %*% beta)
  res <- sum(dbinom(x = m_i, size = 1, prob = Phi, log = TRUE))
  return(res)
}


# Compute prior distribution
prior <- function(beta, mu_beta_0 = 0, s_beta_0 = 4){
  return(sum(dnorm(beta, mean = mu_beta_0, sd = s_beta_0, log = T)))
}


posterior <- function(beta, X, m_i, mu_beta_0 = 0, s_beta_0 = 4){
  result = likelihood(beta,X,m_i) + prior(beta, mu_beta_0, s_beta_0)
  return(result)
}



proposal <- function(beta, Sigma_0){
  return(rmvnorm(1, mean = beta, sigma = Sigma_0))
}


metropolis <- function(beta = NULL, X, m_i, Sigma_0, mu_beta_0 = 0, 
                       s_beta_0 = 4, N_sim = 10000){
  
  # Number of parameters
  D <- NCOL(X)
  # Initialize regression coefficients
  if (is.null(beta)) { beta <- rep(0, D) }
  if (is.null(Sigma_0)) { Sigma_0 <- solve(diag(10, D)) }
  

  beta_chain <- matrix(0, nrow = N_sim, ncol = D)
  beta_chain[1, ] <- beta
  accepted = 0

  for (i in 2:N_sim) {
    beta_chain[i, ] <- proposal(beta_chain[i - 1, ], Sigma_0)
    
    log_acc_prob <- posterior(beta_chain[i, ], X, m_i, 
                              mu_beta_0, s_beta_0) - 
      posterior(beta_chain[i - 1, ], X, m_i, 
                mu_beta_0, s_beta_0)
    if (log(runif(1)) > log_acc_prob) {
      beta_chain[i, ] <- beta_chain[i - 1, ] # reject value
      
    }else{
      accepted = accepted +1
    }
    
  }
  data = list(accepted = accepted, chain=beta_chain )
  return(data)
}

beta <- rep(0, NCOL(X))
Sigma_0 <- solve(diag(0.01, NCOL(X)) + crossprod(X, X))
N_sim = 10000
result = metropolis(beta = beta, X = X, m_i = y, Sigma_0 = Sigma_0, N_sim = N_sim)
acc_1 = round((result$accepted/N_sim)*100,2)
mh_beta_chain = result$chain
result = metropolis(beta = beta, X = X, m_i = y, Sigma_0 = Sigma_0, N_sim = N_sim)
acc_2 = round((result$accepted/N_sim)*100,2)
mh_beta_chain_2 = result$chain
result = metropolis(beta = beta, X = X, m_i = y, Sigma_0 = Sigma_0, N_sim = N_sim)
acc_3 = round((result$accepted/N_sim)*100,2)
mh_beta_chain_3 = result$chain

#-------------------------- Traceplots-----------------------------------------------------------------------

library(latex2exp)
par(mfrow = c(3,3))
plot(beta_chain[, 1], type = "l", xlab="Iterations" , ylab="",
     main = TeX("Chain values of $\\beta_{0}$"))
lines(beta_chain_2[, 1], col="red")
lines(beta_chain_3[, 1], col="green")
lines(cumsum(beta_chain[, 1])/(1:N_sim), col="blue", lwd=2)

plot(beta_chain[, 2], type = "l", xlab="Iterations" ,ylab="",
     main = TeX("Chain values of $\\beta_{1}$"))
lines(beta_chain_2[, 2], col="red")
lines(beta_chain_3[, 2], col="green")
lines(cumsum(beta_chain[, 2])/(1:N_sim), col="blue", lwd=2)

plot(beta_chain[, 3], type = "l", xlab="Iterations" , ylab="",
     main = TeX("Chain values of $\\beta_{2}$"))
lines(beta_chain_2[, 3], col="red")
lines(beta_chain_3[, 3], col="green")
lines(cumsum(beta_chain[, 3])/(1:N_sim), col="blue", lwd=2)


plot(beta_chain[, 4], type = "l", xlab="" , ylab="",
     main = TeX("Chain values of $\\beta_{3}$"))
lines(beta_chain_2[, 4], col="red")
lines(beta_chain_3[, 4], col="green")
lines(cumsum(beta_chain[, 4])/(1:N_sim), col="blue", lwd=2)

plot(beta_chain[, 5], type = "l", xlab="Iterations" , ylab="",
     main = TeX("Chain values of $\\beta_{4}$"))
lines(beta_chain_2[, 5], col="red")
lines(beta_chain_3[, 5], col="green")
lines(cumsum(beta_chain[, 5])/(1:N_sim), col="blue", lwd=2)

plot(beta_chain[, 6], type = "l", xlab="" , ylab="",
     main = TeX("Chain values of $\\beta_{5}$"))
lines(beta_chain_2[, 6], col="red")
lines(beta_chain_3[, 6], col="green")
lines(cumsum(beta_chain[, 6])/(1:N_sim), col="blue", lwd=2)

plot(beta_chain[, 7], type = "l", xlab="Iterations" , ylab="",
     main = TeX("Chain values of $\\beta_{6}$"))
lines(beta_chain_2[, 7], col="red")
lines(beta_chain_3[, 7], col="green")
lines(cumsum(beta_chain[, 7])/(1:N_sim), col="blue", lwd=2)

plot(beta_chain[, 8], type = "l", xlab="" , ylab="",
     main = TeX("Chain values of $\\beta_{7}$"))
lines(beta_chain_2[, 8], col="red")
lines(beta_chain_3[, 8], col="green")
lines(cumsum(beta_chain[, 8])/(1:N_sim), col="blue", lwd=2)

plot(beta_chain[, 9], type = "l", xlab="Iterations" , ylab="",
     main = TeX("Chain values of $\\beta_{8}$"))
lines(beta_chain_2[, 9], col="red")
lines(beta_chain_3[, 9], col="green")
lines(cumsum(beta_chain[, 9])/(1:N_sim), col="blue", lwd=2)
#------------------------------------------ Convergence Diagnostics ------------------------------------------

# Geweke Test
require(coda)
library(ggmcmc)
library(knitr)
sim1.mcmc <- as.mcmc(beta_chain[-(1:burn_in), ])
sim2.mcmc <- as.mcmc(beta_chain_2[-(1:burn_in), ])
sim3.mcmc <- as.mcmc(beta_chain_3[-(1:burn_in), ])
sim.list <- mcmc.list(sim1.mcmc, sim2.mcmc, sim3.mcmc)
data = ggs(sim.list)

# Geweke Test
geweke.sim1=geweke.diag(sim1.mcmc, frac1=0.1, frac2=0.5)
geweke.sim2=geweke.diag(sim2.mcmc, frac1=0.1, frac2=0.5)
geweke.sim3=geweke.diag(sim3.mcmc, frac1=0.1, frac2=0.5)
p = pnorm(abs(geweke.sim1$z),lower.tail=FALSE)*2
p_2 = pnorm(abs(geweke.sim2$z),lower.tail=FALSE)*2
p_3 = pnorm(abs(geweke.sim3$z),lower.tail=FALSE)*2



# Gelman and Rubin's Statistics

# Rubin Statistics
gelman.sims.list <- gelman.diag(sim.list,confidence = 0.95)
gelman.sims.list


#-------------------------------------- Bayesian Variable Selection ------------------------------------------

library(cubature)
index.mods <- expand.grid(1,0:1, 0:1, 0:1, 0:1, 0:1,0:1,0:1,0:1)
posterior_probs = numeric(length = nrow(index.mods))

post_gamma = function(param,x){
  beta = rep(param[1],NCOL(x))
  result = dnorm(beta,mean=0,sd=4)*prod(pnorm(x %*% beta)^y-(1-pnorm(x %*% beta)^(1-y)))
  return(result)
}

for(i in 2:length(posterior_probs)){
  x = X[,which(index.mods[i,]==1)]
  prob = hcubature(f=post_gamma,lowerLimit=c(-Inf),upperLimit=c(Inf),x,maxEval=15)
  posterior_probs[i] = prob$integral
}