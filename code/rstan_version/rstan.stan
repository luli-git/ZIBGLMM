data {
  int<lower=0> J;                  // number of studies 
  vector[2] zero;                  // 0 0
  int y[J,2];                      // number of events, for each arm
  int sample[J,2];                 // number of sample, for each arm
  matrix[2,2] mat_cov; 
  real init;   
}
parameters {
  vector[2] mu;         // fixed effects for treatment and control
  vector[2] nu[J];         // nu[j] represents random effects associated with the j-th subject or group
  vector<lower=0>[2] sigma_nu; // 2-dimensional vector, standard deviations of random effects, each element must be greater than 0
  corr_matrix[2] omega_nu; // correlation matrix for random effects, positive definite and have diagonals of 1
  real<lower=0, upper=1> pi;  // ZI rate 
}
transformed parameters {
  cov_matrix[2] Sigma_nu;
}
model{
  
  Sigma_nu ~ inv_wishart(21.0, mat_cov); //cauchy(0, 2.5); //gamma(1.5, 0.1);
  nu ~ multi_normal(zero, Sigma_nu);
  omega_nu ~ lkj_corr(2.0);
  mu ~ normal(0, 5);
  pi ~ beta(1,1);
  for(n in 1:J){
    if(y[n,1] == 0 && y[n,2]==0){
      target += log_sum_exp(bernoulli_lpmf(1|pi),
                      bernoulli_logit_lpmf(0|pi)  + 
      binomial_logit_lpmf(y[n,1] | sample[n,1], mu[1] + nu[n,1]) + 
      binomial_logit_lpmf(y[n,2] | sample[n,2], mu[2] + nu[n,2] ));
    }
    else{
      target += (bernoulli_lpmf(0|pi) + 
      binomial_logit_lpmf(y[n,1] | sample[n,1], mu[1] + nu[n,1]) + 
      binomial_logit_lpmf(y[n,2] | sample[n,2], mu[2] + nu[n,2]));
    }
  }
}
