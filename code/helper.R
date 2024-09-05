library(dplyr)

### log-density of the multivariate normal distribution
ldmvnorm<-function(X,mu,Sigma,iSigma=solve(Sigma),dSigma=det(Sigma)) 
{
  Y<-t( t(X)-mu)
  sum(diag(-.5*t(Y)%*%Y%*%iSigma)) - .5*(  prod(dim(X))*log(2*pi) + dim(X)[1]*log(dSigma) )
}
###


### sample from the multivariate normal distribution
rmvnorm<-function(n,mu,Sigma)
{
  p<-length(mu)
  res<-matrix(0,nrow=n,ncol=p)
  if( n>0 & p>0 )
  {
    E<-matrix(rnorm(n*p),n,p)
    res<-t(  t(E%*% chol(Sigma)) +c(mu))
  }
  res
}
###


### sample from the Wishart distribution
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}

# expit function
expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

# logit function
logit <- function(x) {
  return(log(x / (1 - x)))
}

# probit function
probt <- function(x, df) {
  return(1 - pt(abs(x), df))
}

# replace infinite values with mean
replace_inf_with_mean <- function(x) {
  finite_vals <- x[is.finite(x)]
  x[is.infinite(x)] <- mean(finite_vals, na.rm = TRUE)
  return(x)
}
