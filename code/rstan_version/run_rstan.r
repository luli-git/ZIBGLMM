# install.packages("rstan")
# install.packages("dplyr")
library("rstan")
library(dplyr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
logit <- function(x){
  return (log(x/(1-x)))
  }
cdsr.data = read.csv("CDSR_data.csv")
expit<-function(x){
  return(exp(x)/(1+exp(x)))
}
data <- cdsr.data%>%
  filter(ID==55574)

data_input = list(J=dim(data)[1], 
               zero = c(0,0), 
               y= matrix(data=cbind(data$r1, data$r2), ncol=2),
               sample=matrix(data=cbind(data$n1, data$n2), ncol=2)
               )


n=matrix(data=cbind(data$n1, data$n2), ncol=2)
Y= matrix(data=cbind(data$r1, data$r2), ncol=2)
fixed_effects = cbind(logit(Y[,1]/n[,1]),logit(Y[,2]/n[,2]))

fixed_effects[is.infinite(fixed_effects)] = NA
for(i in 1:ncol(fixed_effects)){
fixed_effects[is.na(fixed_effects[,i]), i] <- mean(fixed_effects[,i], na.rm = TRUE)
}  
fixed_effects0<-apply(fixed_effects,2,mean)

S0 = cov(matrix(fixed_effects , ncol= 2))


model2 = stan_model("stan_version.stan")

# mat_cov = solve(S0+t(fixed_effects-mtheta)%*%(fixed_effects-mtheta)/m)

fit = stan("stan_version.stan",data=data_input,iter=100000, warmup = 10000, chains=4, thin=100)
print(fit)
 
la <- extract(fit, permuted = TRUE) # return a list of arrays 
mu <- la$mu
a  <- extract(fit, permuted = FALSE) 
a2 <- as.array(fit)
m  <- as.matrix(fit)
d  <- as.data.frame(fit)
  
la  <- extract(fit, permuted = TRUE) 
dim(la$mu)
C = 16 * sqrt(3) / (15 * pi)

print(quantile(((expit((la$mu[,1])/sqrt(1+C^2 * (la$sigma_nu[,1])^2)))/(expit((la$mu[,2])/sqrt(1+C^2 * (la$sigma_nu[,2])^2)))), c(0.025,0.975)))

print(median(((expit((la$mu[,1])/sqrt(1+C^2 * (la$sigma_nu[,1])^2)))/(expit((la$mu[,2])/sqrt(1+C^2 * (la$sigma_nu[,2])^2))))))

