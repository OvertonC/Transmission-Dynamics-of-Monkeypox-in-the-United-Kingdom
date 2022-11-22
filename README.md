# Transmission Dynamics of Monkeypox in the United Kingdom
 Repository containing stan models for estimating delay distributions in the presence of right-truncation, epidemic phase, and interval-censoring biases


This repository contains the model files associated with the paper "Transmission Dynamics of Monkeypox in the United Kingdom" by Tom Ward et al (2022).

The stan models are easily generalisable to estimating other delays distributions.

The code requires input data of the form: EL, EU, SL, SU, where E is the primary event in the delay, S is the secondary event, L denotes the lower bound and U denotes the upper bound, for each observation.

To run the code, one needs to also specify the right-truncation date, T, which is the maximum date on which secondary events can be observed, i.e. T >= SU.  

The 4 folders contain the trace results for the results shown in the paper. 

For the Weibull and gamma distributions, we assume a normal distribution prior for the mean, with mean 5 and standard deviation 10. We assume an exponential prior with rate 0.0001 for the Weibull shape parameter, and a normal distribution with mean 5 and standard deviation 10 for the gamma standard deviation parameter. For the lognormal, we assume standard normal priors for both the log mean and log standard deviation.
