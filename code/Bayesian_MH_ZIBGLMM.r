####### This is a Bayesian implementation of the ZIBGLMM.
####### Reference: https://pdhoff.github.io/book/


# clean environment
rm(list = ls())

setwd(getwd())
sim.data = read.csv("data/CDSR_data.csv")
source("code/helper.R")



# start and end id of the data to be processed,
# for example, if you want to process the first meta-analyses, you can set start_id = 1 and end_id = 1
start_id <- 1
end_id <- 1
num_iter = 100000
thinning = 100
burnin = 100
results_folder = "results_bayesian_zibglmm"


# log likelihood function for ZIBGLMM
loglik_func <- function(ZIR, fixed_effects_i_mean, Y, n, true_density) {
  logll = rep(0, length(Y[,1]))
  if(true_density){
    for (i in 1:length(Y[,1])){
        if (Y[i,1] == 0 && Y[i,2] == 0) {
          logll[i] <- log(
                dbinom(1, 1, ZIR ) + 
                exp( dbinom(0, 1, ZIR,log=TRUE) + 
                        dbinom(Y[i,1], n[i,1], expit(fixed_effects_i_mean[i,1]), log=TRUE) + dbinom(Y[i,2], n[i,2], expit(fixed_effects_i_mean[i,2]), log=TRUE)))   

        } else {
          logll[i] <- dbinom(0, 1, ZIR,log=TRUE) + dbinom(Y[i,1], n[i,1], expit(fixed_effects_i_mean[i,1]), log=TRUE) + dbinom(Y[i,2], n[i,2], expit(fixed_effects_i_mean[i,2]), log=TRUE)
        }
      }
  }
  else{ 
    for (i in 1:length(Y[,1])){
      if (Y[i,1] == 0 && Y[i,2] == 0) {
        logll[i] <- (log(ZIR + (1 - ZIR) * exp(n[i,1] * log(1 - expit(fixed_effects_i_mean[i,1])) + n[i,2] * log(1 - expit(fixed_effects_i_mean[i,2])))))[1]
      } else {
        logll[i] <- (log(1 - ZIR) + 
                    Y[i,1] * log(expit(fixed_effects_i_mean[i,1])) + (n[i,1] - Y[i,1]) * log(1 - expit(fixed_effects_i_mean[i,1])) + 
                    Y[i,2] * log(expit(fixed_effects_i_mean[i,2])) + (n[i,2] - Y[i,2]) * log(1 - expit(fixed_effects_i_mean[i,2])))[1]
      }
    }
  }
  return (sum(logll))
}

# Bayesian ZIBGLMM function
ZIBGLMM_Mh <- function(Y, n, num_iter = 100000, thinning = 100, burnin = 100){
  m<-dim(Y)[1]
  p <- 2  
  ##  initialize the fixed effects
  fixed_effects = cbind(logit(Y[,1]/n[,1]),logit(Y[,2]/n[,2]))
  fixed_effects[is.infinite(fixed_effects)] = NA
  for(i in 1:ncol(fixed_effects)){
    fixed_effects[is.na(fixed_effects[,i]), i] <- mean(fixed_effects[,i], na.rm = TRUE)
  }  
  fixed_effects0<-apply(fixed_effects,2,mean)
  S0 = cov(matrix(fixed_effects , ncol= 2))
  eta0<-p + 2
  iL0<-iSigma<-solve(S0)
  THETA.post<<-SIGMA.post<<-BETA.post<<-PI.post <<- NU.post<-NULL
  set.seed(1)
  accept = 0 
  a=1
  b=1
  ZIR<-rbeta(1, a, b)
  estimates <- estimate_upper <- estimate_lower<- NULL
  for(s in 1: num_iter){
    ##updateTheta
    Lm<-solve(iL0+m*iSigma)
    fixed_effectsm<-Lm%*%(iL0%*%fixed_effects0+iSigma%*%apply(fixed_effects,2,sum))
    theta<-t(rmvnorm(1,fixed_effectsm,Lm)) 
    ##updateSigma
    mtheta<-matrix(theta,dim(fixed_effects)[1],2,byrow=TRUE)
    iSigma<-rwish(1,eta0+m, solve(S0+t(fixed_effects-mtheta)%*%(fixed_effects-mtheta)/m))
    Sigma<-solve(iSigma)
    #updateGamma
    gamma <- NULL
    for(j in 1:dim(fixed_effects)[1]){
      nu = fixed_effects[j,]
      if(Y[j,1] == 0 & Y[j,2] == 0){
        prob_u = ZIR
        prob_l = exp((n[j,1] - Y[j,1])*log(1- expit(nu[1])) + (n[j,2] - Y[j,2])*log(1- expit(nu[2])) + log(1-ZIR)) + ZIR 
        prob = prob_u/prob_l
        gamma[j] = rbinom(1, 1, prob)
      }
      else{
        gamma[j]=0
      }
    }
    #updatePi
    a = sum(gamma)+a
    b = m-sum(gamma)+b
    ZIR <- rbeta(1, sum(gamma)+a, m-sum(gamma)+b)
    ##updatefixed_effects
    dSigma<-det(Sigma)
    for(j in 1:m){
      beta.p<-t(rmvnorm(1, (fixed_effects[j,]), Sigma))
      if(Y[j,1] == 0 & Y[j,2]==0){
        lr<- log(
          dbinom(1, 1, ZIR ) + 
            exp( dbinom(0, 1, ZIR,log=TRUE) + 
                   dbinom(Y[j,1], n[j,1], expit(beta.p[1]), log=TRUE) + dbinom(Y[j,2], n[j,2], expit(beta.p[2]), log=TRUE)))  - 
          log(
            dbinom(1, 1, ZIR ) + 
              exp( dbinom(0, 1, ZIR,log=TRUE) + 
                     dbinom(Y[j,1], n[j,1], expit(fixed_effects[j,1]), log=TRUE) + dbinom(Y[j,2], n[j,2], expit(fixed_effects[j,2]), log=TRUE))) + 
          ldmvnorm(t(beta.p),theta,Sigma,iSigma=iSigma,dSigma=dSigma)-
          ldmvnorm(t(fixed_effects[j,]),theta,Sigma,iSigma=iSigma,dSigma=dSigma)
      }
      else {
        lr <- 
          dbinom(Y[j,1], n[j,1], expit(beta.p[1]), log=TRUE) + dbinom(Y[j,2], n[j,2], expit(beta.p[2]), log=TRUE) - 
          dbinom(Y[j,1], n[j,1], expit(fixed_effects[j,1]), log=TRUE) - dbinom(Y[j,2], n[j,2], expit(fixed_effects[j,2]), log=TRUE) + 
          ldmvnorm(t(beta.p),theta,Sigma,iSigma=iSigma,dSigma=dSigma)-
          ldmvnorm(t(fixed_effects[j,]),theta,Sigma,iSigma=iSigma,dSigma=dSigma)
      }
      if(log(runif(1))<lr){
        accept = accept + 1
        fixed_effects[j,]<-beta.p
      }
    }
    
    ##store every nth iteration
    if(s%%thinning==0){
      THETA.post<-rbind(THETA.post,t(theta))
      SIGMA.post<-rbind(SIGMA.post,c(Sigma))
      PI.post<-rbind(PI.post,ZIR)
      NU.post<-rbind(NU.post, c(fixed_effects))
    }   
  }
  
  c = 16 * sqrt(3) / (15 * pi)
  estimates = ((expit(THETA.post[burnin:1000,1]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,1])))/expit(THETA.post[burnin:1000,2]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,4])))))
  estimate_lower= quantile(((expit(THETA.post[burnin:1000,1]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,1])))/expit(THETA.post[burnin:1000,2]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,4]))))), c(0.025,0.975))[1]
  estimate_upper= quantile(((expit(THETA.post[burnin:1000,1]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,1])))/expit(THETA.post[burnin:1000,2]/sqrt(1+c^2 * (SIGMA.post[burnin:1000,4]))))), c(0.025,0.975))[2]
  pi_estimate = median(PI.post[burnin:1000])
  pi_lower = quantile(PI.post[burnin:1000], c(0.025,0.975))[1]
  pi_upper = quantile(PI.post[burnin:1000], c(0.025,0.975))[2]
  
  # calculate p value 
  relative_risk <- median(estimates)
  standard_deviation <- sd(estimates)
  sample_size <- dim(n)[1]
  p_value_analog <- sum(estimates > 1) / length(estimates)
  
  # Two-sided test (relative risk is different than 1)
  p_value_analog <- 2 * min(sum(estimates > 1), sum(estimates < 1)) / length(estimates)

  # Calculate natural logarithm of RR
  ln_rr <- log(relative_risk)
  
  standard_deviation <- sd(log(estimates)) 
  # Calculate standard error (if you don't have sample size, use standard_deviation as standard_error)
  standard_error <- standard_deviation / sqrt(sample_size)
  # Calculate Z-score
  z_score <- (ln_rr - 0) / standard_error
  # Calculate p-value (two-tailed test)
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  
  # pvalue 
  p_value <- (probt(abs(ln_rr/standard_error),1000))*2;
  # save the results in a list
  results <- list(THETA.post, SIGMA.post, PI.post, NU.post)
  return (list(median(estimates),estimate_lower,estimate_upper, pi_estimate, pi_lower, pi_upper, p_value, standard_error, NU.post, PI.post, accept))
}

estimate_lower <- NULL
estimate_upper <- NULL
estimates <- NULL
pi_estimates <- NULL
pi_lower <- NULL
pi_upper <- NULL
p_value <- NULL
standard_error<-NULL
DIC_array <- NULL
DIC_true_array <- NULL


for(id in start_id:end_id){
  print(id)
  data <- sim.data%>%
    filter(ID==unique(sim.data$ID)[id])
  n = data[,c("n1", "n2")]
  Y = data[,c("r1", "r2")] # 
  result = ZIBGLMM_Mh(Y,n)
  fixed_effects_i = result[[9]][burnin: (num_iter/thinning)]
  accept = result[[11]]
  
  # reshape to 10 * 2 * 1000 
  fixed_effects_i_reshaped = array(fixed_effects_i, dim = c(dim(data)[1] ,2, (num_iter/thinning)))
  # calculate average of 1000
  fixed_effects_i_mean = apply(fixed_effects_i_reshaped, c(1,2), mean)
  ZI_rate = result[[10]][burnin:(num_iter/thinning)]
  # calculate log likelihood
  log_lik_bar = loglik_func(mean(ZI_rate), (fixed_effects_i_mean), Y,n, FALSE)
  # average over 1000 pairs of fixed effects
  Eloglik = mean(sapply(1:(num_iter/thinning-burnin), function(i) loglik_func(ZI_rate[i], fixed_effects_i_reshaped[,,i], Y,n,FALSE)))
  pDIC = -2 * Eloglik + 2 * log_lik_bar
  DIC = -2 * log_lik_bar + 2 * pDIC
  log_lik_bar = loglik_func(mean(ZI_rate), (fixed_effects_i_mean), Y,n,TRUE)
  # average over 1000 pairs of fixed effects
  Eloglik = mean(sapply(1:(num_iter/thinning-burnin), function(i) loglik_func(ZI_rate[i], fixed_effects_i_reshaped[,,i], Y,n,TRUE)))
  pDIC = -2 * Eloglik + 2 * log_lik_bar
  DIC_true = -2 * log_lik_bar + 2 * pDIC

  estimates = c(estimates, result[[1]])
  estimate_lower= c(estimate_lower, result[[2]])
  estimate_upper= c(estimate_upper, result[[3]])
  pi_estimates = c(pi_estimates, result[[4]])
  pi_lower = c(pi_lower, result[[5]])
  pi_upper = c(pi_upper, result[[6]])
  p_value = c(p_value, result[[7]])
  DIC_array = c(DIC_array, DIC)
  DIC_true_array = c(DIC_true_array, DIC_true)

}


 
results = cbind(estimates, estimate_lower, estimate_upper, pi_estimates, pi_lower, pi_upper, p_value, DIC_array, DIC_true_array)

if (!dir.exists(results_folder)){
  dir.create(results_folder)
}
write.csv(results, paste0(results_folder, "/sim_zibglmm_bayesian", start_id, "-" ,end_id, ".csv"))
 
