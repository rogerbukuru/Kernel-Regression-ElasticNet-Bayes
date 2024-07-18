rm(list = ls())
library(ggplot2)
dat=read.table("measurements.csv",header = T)
x1<- as.matrix(dat[1:100,])
x2<- as.matrix(dat[101:200,])
x<- as.matrix(dat[1:200,])

M = 5
beta_max = 2
alphas = seq(0.1, M,length=200)
betas  = seq(0.1, beta_max,length=200)


# distribution of x 
x_distribution = function(alpha,beta,x){ # for both alpha1 and alpha2
  part_1 = (beta/((beta^2)+(x-alpha)^2))
  part_2 = 1/(atan(alpha/beta)+atan(((M-alpha)/beta)))
  f = part_1*part_2
  return(f)
}


log_posterior_distribution = function(params,x){
  alpha = params[1]
  beta  = params[2]
  alpha2 = alpha+1
  p = params[3]
  f=1:200
  for (i in 1:200){
    f1 = x_distribution(alpha,beta,x[i])
    f2 = x_distribution(alpha2,beta,x[i])
    f[i] = (p*f1)+((1-p)*f2)
  }
  
  posterior = log(prod(f))
  return(posterior)
}


log_posterior_distribution_cont = function(alpha,beta,p,x){
  alpha2 = alpha+1
  #consts = k*(1/beta_max)*(1/beta)
  f=1:200
  for (i in 1:200){
    f1 = x_distribution(alpha,beta,x[i])
    f2 = x_distribution(alpha2,beta,x[i])
    f[i] = (p*f1)+((1-p)*f2)
    #f[i]=((0.5)*alpha_1(params,obs[i]))+((0.5)*alpha_2(params,obs[i])) 
  }
  
  posterior = log(prod(f))
  return(posterior)
}
optimal_results = optim(par= c(1,0.5,0.5), f=log_posterior_distribution,x=x,
                        lower=c(1,0.1,0.5),
                        upper=c(4,4,1),
                        control = list(fnscale=-1),
                        method="L-BFGS-B", hessian = TRUE)

p = optimal_results$par[3]
g <- function(alpha,beta) log_posterior_distribution_cont(alpha,beta,p,x)
z = outer(alphas,betas,Vectorize(g))
contour(alphas,betas,z, nlevels=20, xlab="Alpha", ylab="Beta")
points(optimal_results$par[1], optimal_results$par[2], pch = 4, col="red")
legend("topright",c("MLE Position"), col="red", pch=3, pt.cex=1.2, pt.lwd=2) 


#contour(alphas, betas, z, xlim=c(0,5), ylim=c(0,1), nlevels = 20,xlab="Alpha", ylab="Beta")
#points(optimal_results$par[1], optimal_results$par[2], pch=4, col="blue", lwd=2) # true value
#legend("topright",c("MLE Position"), col="blue", pch=3, pt.cex=1.2, pt.lwd=2) 

# Expectation and variance

optimal_results_2 = optim(par= c(1,0.5,0.5), f=log_posterior_distribution,x=x,
                          lower=c(1,0.1,0.5),
                          upper=c(4,2,1),
                          control = list(fnscale=-1),
                          method="L-BFGS-B", hessian = TRUE)

mu=optimal_results_2$par
cov_matrix=solve(-optimal_results_2$hessian)

#Sample from a mutlivariate normal distribtuion
library(mvtnorm)
target_dist= rmvnorm(1000, mean =mu, sigma = cov_matrix) # target distribution 

fvalues=1:1000
for (i in 1:1000){
  
  fvalues[i]=log_posterior_distribution(target_dist[i,],x)
}

#Accept/reject 
h=exp(log_posterior_distribution(mu,x))
set.seed(2)
unif_probs=runif(n=1000,0,h)
f_alpha_beta_p=exp(fvalues)
accpt=which(f_alpha_beta_p>unif_probs)
accepted_alpha_beta_p=target_dist[accpt,]
mean_res=apply(accepted_alpha_beta_p,2,mean)
var_res=apply(accepted_alpha_beta_p,2,var)

# Credibility Intervals
alpha_cred=c(mean_res[1]-(1.96*sqrt(var_res[1])),mean_res[1]+(1.96*sqrt(var_res[1])))
beta_cred=c(mean_res[2]-(1.96*sqrt(var_res[2])),mean_res[2]+(1.96*sqrt(var_res[2])))
p_cred=c(mean_res[3]-(1.96*sqrt(var_res[3])),mean_res[3]+(1.96*sqrt(var_res[3])))