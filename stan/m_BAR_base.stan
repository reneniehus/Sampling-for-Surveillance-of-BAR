 data {
   // Define variables in data
   int<lower=0> N; // Number of observations (an integer)
   
   //
   real q_log[N]; // the measured gene quantities
   int<lower=0 , upper=1> gene_zero[N]; // when gene is zero
   //
   int n_gene_on;// how many have gene on
   int gene_on_loc[n_gene_on];// where gene is 
   //
 }
 parameters {
   // Define parameters to estimate
   real<lower=0, upper=1> f_gene_carr;
   real mu; // mean for gene distribution
   real sigma;// variance for gene distribution
   // effects
 }
 transformed parameters  {
   // 
   
 }
 model {
   // Prior part of Bayesian inference
   // (no need to specify if non-informative)
   f_gene_carr ~ uniform(0,1);
   mu ~ normal(4,50);
   sigma ~ cauchy(0,2);

   // Likelihood part of Bayesian inference
   // distribution of zeros by bernoulli
   for (i in 1:N) gene_zero[i] ~ bernoulli( 1-f_gene_carr );
   // distribution for where is some gene
   for (i in 1:n_gene_on) q_log[gene_on_loc[i]] ~ normal( mu , sigma ); 

 }
 generated quantities {
   // define predicted vector
     vector[N] q_pred_log;
     vector[N] log_lik;
     for (i in 1:N) q_pred_log[i] = bernoulli_rng(1-f_gene_carr)*normal_rng( mu, sigma ); 
     // make log likelihood
     for (i in 1:N) log_lik[i] = normal_lpdf( q_log[i] | mu, sigma );
 }
   