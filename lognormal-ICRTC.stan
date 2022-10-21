// distribution: lognormal
// truncation:   right
// doubly interval censored

data {
  int<lower = 0> N;                     // number of records
  vector<lower = 0>[N] a_minus;         // lower limit of event A
  vector<lower = 0>[N] a_plus;          // upper limit of event A
  vector<lower = 0>[N] b_minus;         // lower limit of event B
  vector<lower = 0>[N] b_plus;          // upper limit of event B
  int<lower = 0, upper = 1> incubation; // inclusion of time from O to A (incubation period)
  real<lower = 0> upper_bound;          // the latest time of observation
  real<lower = 0> r;                    // the exponential growth rate
}

transformed data {
  int X_i[0];                           // empty array
}

parameters {
  real logmean;                             // natural log of the mean
  real logsd;                               // natural log of the standard deviation
  vector<lower = 0, upper = 1>[N] a_window; // where time a lies in the event A window
  vector<lower = 0, upper = 1>[N] b_window; // where time b lies in the event B window
  vector<lower = 0, upper = 1>[N] a2_window; // where time a2 lies in the event A window
  vector<lower = 0>[N] t0;                  // time from O to A
}

transformed parameters {
  real<lower = 0> sigma = sqrt(log1p_exp(2 * (logsd - logmean)));
  real mu = logmean - (sigma^2)/2;

  vector<lower = min(a_minus), upper = max(a_plus)>[N] a;
  vector<lower = min(a_minus), upper = max(a_plus)>[N] a2;
  vector<lower = min(b_minus), upper = max(b_plus)>[N] b;
  vector[N] ub;
  vector<lower = 0>[N] tstar;

  b = b_minus + (b_plus - b_minus) .* b_window;
  for (n in 1:N)
    ub[n] = min([a_plus[n], b[n]]');
  a = a_minus + (ub - a_minus) .* a_window;
  a2 = a_minus + (ub - a_minus) .* a2_window;
  for (n in 1:N)
    ub[n] = min([upper_bound,a[n]+14]');
  tstar = upper_bound - a;
}

model {
  logmean ~ std_normal();
  logsd ~ std_normal();

  if (incubation)
    t0 ~ lognormal(1.63, 0.5);
  else
    t0 ~ normal(0, 1e-10);

  target += lognormal_lpdf((b - a + t0) | mu, sigma) - lognormal_lcdf((upper_bound - a2 + t0) | mu, sigma);
}

generated quantities {
  real<lower = 0> mean_ = exp(mu + (sigma^2)/2);
  real<lower = 0> sd_ = sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2));
  real<lower = 0> limit_val_ = exp(mu + sigma*(1.959963984540));
  vector[N] log_likelihood;
  for (n in 1:N)
    log_likelihood[n] = lognormal_lpdf(b[n] - a[n] + t0[n] | mu, sigma) - lognormal_lcdf((upper_bound - a2[n] + t0[n]) | mu, sigma);
}

