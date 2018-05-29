
data {
  int<lower=0> N;
  int<lower=0> N_spp;
  int<lower=0> spp[N];
  vector[N] x_mean;
  vector[N] x_se;
  vector[N] y_mean;
  vector[N] y_se;
}

parameters {
  vector[N_spp] z_alpha_spp;
  
  real mu_alpha;
  real mu_beta;
  
  vector[N] y;
  vector[N] x;
  
  real<lower=0> sigma_alpha_spp;
  real<lower=0> sigma_y;
}

transformed parameters {
  vector[N_spp] alpha_spp;
  
  vector[N] y_hat;
  
  alpha_spp = z_alpha_spp * sigma_alpha_spp;
  
  for (i in 1:N) {
    y_hat[i] = (mu_alpha + alpha_spp[spp[i]]) +
               (mu_beta * x[i]);
  }
}

model {
  z_alpha_spp ~ normal(0, 1);
  
  mu_alpha ~ normal(0, 10);
  mu_beta ~ normal(0, 10);
  
  sigma_alpha_spp ~ normal(0, 10);
  
  sigma_y ~ normal(0, 10);
  
  y_mean ~ normal(y, y_se);
  
  x ~ normal(0, 10);
  x_mean ~ normal(x, x_se);
  
  y ~ normal(y_hat, sigma_y);
}
