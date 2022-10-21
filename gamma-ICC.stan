// distribution: gamma
// truncation:   no
// doubly interval censored

data {
  int<lower = 0> N;                     // number of records
  vector<lower = 0>[N] a_minus;         // lower limit of event A
  vector<lower = 0>[N] a_plus;          // upper limit of event A
  vector<lower = 0>[N] b_minus;         // lower limit of event B
  vector<lower = 0>[N] b_plus;          // upper limit of event B
  int<lower = 0, upper = 1> incubation; // inclusion of time from O to A (incubation period)
}

parameters {
  real<lower = 0>mean_;                               // the distribution mean
  real<lower = 0>sd_;                                 // the distribution standard deviation
  vector<lower = 0, upper = 1>[N] a_window; // where time a lies in the event A window
  vector<lower = 0, upper = 1>[N] b_window; // where time b lies in the event B window
  vector<lower = 0>[N] t0;                  // time from O to A
}

transformed parameters {
  real<lower = 0> alpha = (mean_/sd_)^2;
  real<lower = 0> beta = mean_/(sd_^2);

  vector<lower = min(a_minus), upper = max(a_plus)>[N] a;
  vector<lower = min(b_minus), upper = max(b_plus)>[N] b;
  vector[N] ub;

  b = b_minus + (b_plus - b_minus) .* b_window;
  for (n in 1:N)
    ub[n] = min([a_plus[n], b[n]]');
  a = a_minus + (ub - a_minus) .* a_window;
}

model {
  mean_ ~ normal(5.0, 10.0);
  sd_ ~ cauchy(0, 5.0);

  if (incubation)
    t0 ~ lognormal(1.460337, 0.9135229);
  else
    t0 ~ normal(0, 1e-10);

  target += gamma_lpdf((b - a + t0) | alpha, beta);
}

generated quantities {
  vector[N] log_likelihood;
  for (n in 1:N)
    log_likelihood[n] = gamma_lpdf(b[n] - a[n] + t0[n] | alpha, beta);
}

