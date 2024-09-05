# Install and load required packages
if (!require(rstan)) install.packages("rstan")
if (!require(dplyr)) install.packages("dplyr")

library(rstan)
library(dplyr)

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

# Helper functions
logit <- function(x) log(x / (1 - x))
expit <- function(x) exp(x) / (1 + exp(x))

# Load and prepare data
cdsr.data = read.csv("data/CDSR_data.csv")
# select the meta-analysis you are interested in running the model for 
data <- cdsr.data%>%
  filter(ID==1) # 55574 for the case study 


data_input = list(J=dim(data)[1], 
               zero = c(0,0), 
               y= matrix(data=cbind(data$r1, data$r2), ncol=2),
               sample=matrix(data=cbind(data$n1, data$n2), ncol=2)
               )
# Prepare fixed effects
n=matrix(data=cbind(data$n1, data$n2), ncol=2)
Y = matrix(data = cbind(data$r1, data$r2), ncol = 2)
fixed_effects = cbind(logit(Y[,1]/n[,1]),logit(Y[,2]/n[,2]))
fixed_effects[is.infinite(fixed_effects)] = NA
# Handle infinite and NA values
for (i in 1:ncol(fixed_effects)) {
  fixed_effects[is.na(fixed_effects[, i]), i] <- mean(fixed_effects[, i], na.rm = TRUE)
}

fixed_effects0<-apply(fixed_effects,2,mean)

S0 = cov(matrix(fixed_effects , ncol= 2))

# Compile and run Stan model, replace here with bayesian_zibglmm.stan to run the zibglmm version
fit = stan("code/rstan_version/bayesian_bglmm.stan",data=data_input,iter=100000, warmup = 10000, chains=4, thin=100)

# Extract and analyze results
la <- extract(fit, permuted = TRUE) 
mu <- la$mu

C = 16 * sqrt(3) / (15 * pi)

# Print results
# 95% credible interval
print(quantile(((expit((la$mu[, 1]) / sqrt(1 + C^2 * (la$sigma_nu[, 1])^2))) / (expit((la$mu[, 2]) / sqrt(1 + C^2 * (la$sigma_nu[, 2])^2)))), c(0.025, 0.975)))
# median
print(median(((expit((la$mu[, 1]) / sqrt(1 + C^2 * (la$sigma_nu[, 1])^2))) / (expit((la$mu[, 2]) / sqrt(1 + C^2 * (la$sigma_nu[, 2])^2))))))

