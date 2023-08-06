#This is is a Bayesian implementation of the ZIBGLMM. 
####### Author: Lu Li 
####### Reference: https://pdhoff.github.io/book/

args <- commandArgs(trailingOnly = TRUE)
print(args)
# Ensure that start and end IDs are provided
if(length(args) < 2) {
  stop("You must supply a start and end ID.")
}

# Convert the arguments to integers
start_id <- as.integer(args[1])
end_id <- as.integer(args[2])
print(start_id)
print(end_id)
library(dplyr)
setwd(getwd())
sim.data = read.csv("/home/luli1/ZIBGLMM/cleaned.csv")
source("/home/luli1/ZIBGLMM/helper.R")
# sim.data=read.csv("/Users/luli1/Dropbox/_A_Research/ZIP/simulations/0509_1000_simulations.csv")
# 4781, 4790
expit<-function(x){
  return(exp(x)/(1+exp(x)))
}
probt<-function(x,df){
  return(1-pt(abs(x),df))
}
loglik_func <- function( fixed_effects_i_mean, Y, n, true_density) {
  logll = rep(0, length(Y[,1]))
  if(true_density){
    for (i in 1:length(Y[,1])){
          logll[i] <-  dbinom(Y[i,1], n[i,1], expit(fixed_effects_i_mean[i,1]), log=TRUE) + dbinom(Y[i,2], n[i,2], expit(fixed_effects_i_mean[i,2]), log=TRUE)
      }
  }
  else{
    # logll = rep(0, length(Y[,1]))
    for (i in 1:length(Y[,1])){
        logll[i] <- ( 
                    Y[i,1] * log(expit(fixed_effects_i_mean[i,1])) + (n[i,1] - Y[i,1]) * log(1 - expit(fixed_effects_i_mean[i,1])) + 
                    Y[i,2] * log(expit(fixed_effects_i_mean[i,2])) + (n[i,2] - Y[i,2]) * log(1 - expit(fixed_effects_i_mean[i,2])))[1]
    }
  }
  return (sum(logll))
}
 

BGLMM_Mh <- function(Y, n, num_iter = 100000, thinning = 100, burnin = 100){
  m<-dim(Y)[1]
  p <- 2 #dim(X)[2]
  ## start values
  fixed_effects = cbind(rep(-1, m), rep(-1, m))
  fixed_effects0<-apply(fixed_effects,2,mean)
  S0<- matrix(c(1,0,0,1), nrow=2, ncol=2) #cov(fixed_effects)
  eta0<-p
  iL0<-iSigma<-solve(S0)
  THETA.post<<-SIGMA.post<<-BETA.post<<-PI.post <<- NU.post<-NULL
  set.seed(1)
  accept = 0 
  
  estimates <- estimate_upper <- estimate_lower<- NULL
  for(s in 1: num_iter){
    ##updateTheta
    Lm<-solve(iL0+m*iSigma)
    fixed_effectsm<-Lm%*%(iL0%*%fixed_effects0+iSigma%*%apply(fixed_effects,2,sum))
    theta<-t(rmvnorm(1,fixed_effectsm,Lm)) 
    ##updateSigma
    mtheta<-matrix(theta,m,2,byrow=TRUE)
    iSigma<-rwish(1,eta0+m, solve(S0+t(fixed_effects-mtheta)%*%(fixed_effects-mtheta)/m))
    Sigma<-solve(iSigma)

    ##updatefixed_effects
    dSigma<-det(Sigma)
    for(j in 1:m){
      beta.p<-t(rmvnorm(1, (fixed_effects[j,]),.5*Sigma))
        lr <- 
          dbinom(Y[j,1], n[j,1], expit(beta.p[1]), log=TRUE) + dbinom(Y[j,2], n[j,2], expit(beta.p[2]), log=TRUE) - 
          dbinom(Y[j,1], n[j,1], expit(fixed_effects[j,1]), log=TRUE) - dbinom(Y[j,2], n[j,2], expit(fixed_effects[j,2]), log=TRUE) + 
          ldmvnorm(t(beta.p),theta,Sigma,iSigma=iSigma,dSigma=dSigma)-
          ldmvnorm(t(fixed_effects[j,]),theta,Sigma,iSigma=iSigma,dSigma=dSigma)

      if(log(runif(1))<lr){
        accept = accept + 1
        fixed_effects[j,]<-beta.p
      }
    }
    
    ##store every nth iteration
    if(s%%thinning==0){
      THETA.post<-rbind(THETA.post,t(theta))
      SIGMA.post<-rbind(SIGMA.post,c(Sigma))
      NU.post<-rbind(NU.post, c(fixed_effects))
    }   
  }
  
  c = 16 * sqrt(3) / (15 * pi)
  estimates = (1/(expit(THETA.post[burnin:1000,1]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,1])))/expit(THETA.post[burnin:1000,2]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,4])))))
  estimate_lower= quantile(((expit(THETA.post[burnin:1000,1]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,1])))/expit(THETA.post[burnin:1000,2]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,4]))))), c(0.025,0.975))[1]
  estimate_upper= quantile(((expit(THETA.post[burnin:1000,1]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,1])))/expit(THETA.post[burnin:1000,2]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,4]))))), c(0.025,0.975))[2]

  # calculate p value 
  relative_risk <- median(estimates)
  standard_deviation <- sd(estimates)
  sample_size <- dim(n)[1]
  

  # Calculate natural logarithm of RR
  ln_rr <- log(relative_risk)
  # Calculate standard error (if you don't have sample size, use standard_deviation as standard_error)
  standard_error <- standard_deviation / sqrt(sample_size)
  # Calculate Z-score
  z_score <- (ln_rr - 0) / standard_error
  # Calculate p-value (two-tailed test)
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  
  # pvalue 
  p_value <- (probt(abs(ln_rr/standard_error),1000))*2;
  # save the results in a list
  results <- list(THETA.post, SIGMA.post, NU.post)
  write.csv(results, file = paste0("/home/luli1/ZIBGLMM/0531_cdsr/BGLMM_results/specifics/",id, "_results.csv"))
  return (list(median(estimates),
                estimate_lower,
                estimate_upper,
                p_value, 
                standard_error, 
                NU.post, 
                accept))
}

estimate_lower <- NULL
estimate_upper <- NULL
estimates <- NULL
p_value <- NULL
standard_error<-NULL
DIC_array <- NULL
DIC_true_array <- NULL
# sim.data = read.csv("/Users/luli1/Dropbox/_A_Research/ZIP/restart/cdsr_data/cleaned.csv")
 
for(id in start_id:end_id){
  print(id)
  data <- sim.data%>%
    filter(ID==unique(sim.data$ID)[id])
  n = data[,c("n1", "n2")]
  Y = data[,c("r1", "r2")]
  result = BGLMM_Mh(Y,n)
  fixed_effects_i = result[[6]][100:1000]
  accept = result[[7]]
  # reshape to 10 * 2 * 1000
  fixed_effects_i_reshaped = array(fixed_effects_i, dim = c(dim(data)[1] ,2,1000))
  # calculate average of 1000
  fixed_effects_i_mean = apply(fixed_effects_i_reshaped, c(1,2), mean)
  # calculate log likelihood
  log_lik_bar = loglik_func( (fixed_effects_i_mean), Y,n, FALSE)
  # average over 1000 pairs of fixed effects
  Eloglik = mean(sapply(1:900, function(i) loglik_func( fixed_effects_i_reshaped[,,i], Y,n,FALSE)))
  pDIC = -2 * Eloglik + 2 * log_lik_bar
  DIC = -2 * log_lik_bar + 2 * pDIC
  log_lik_bar = loglik_func( (fixed_effects_i_mean), Y,n,TRUE)
  # average over 1000 pairs of fixed effects
  Eloglik = mean(sapply(1:900, function(i) loglik_func( fixed_effects_i_reshaped[,,i], Y,n,TRUE)))
  pDIC = -2 * Eloglik + 2 * log_lik_bar
  DIC_true = -2 * log_lik_bar + 2 * pDIC

  estimates = c(estimates, result[[1]])
  estimate_lower= c(estimate_lower, result[[2]])
  estimate_upper= c(estimate_upper, result[[3]])
  p_value = c(p_value, result[[4]])
  DIC_array = c(DIC_array, DIC)
  DIC_true_array = c(DIC_true_array, DIC_true)

}

a = cbind(estimates, estimate_lower, estimate_upper, p_value, DIC_array, DIC_true_array)
print("here")
print(paste0("/home/luli1/ZIBGLMM/0531_cdsr/BGLMM_results/sim_zibglmm_bayesian", start_id, "-" ,end_id, ".csv"))
write.csv(a, paste0("/home/luli1/ZIBGLMM/0531_cdsr/BGLMM_results/sim_bglmm_bayesian", start_id, "-" ,end_id, ".csv"))


# library(coda)
# autocorr.plot(as.mcmc(fixed_effects_i))
# thin_interval = 5
# thin_indx = seq(1,1000,thin_interval)
# traceplot(as.mcmc(fixed_effects_i))
# effectiveSize(as.mcmc(fixed_effects_i))


