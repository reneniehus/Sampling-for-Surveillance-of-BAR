 data {
   // Define variables in data
   int<lower=0> N; // Number of observations (an integer)
   int<lower=0> N_c; // Number of obervations from gene carriers
   int<lower=1> pool_size; // number of samples in every pool
   //
   real q_log[N]; // the measured gene quantities
   real car_q_log[N_c]; // gene quantities from carriers
   int<lower=0 , upper=1> gene_zero[N]; // when gene is zero
   //
   int n_gene_on;// how many have gene on
   int gene_on_loc[n_gene_on];// where gene is 
   //
 }
 parameters {
   // Define parameters to estimate
   real<lower=0, upper=1> f_gene_carr;
   real mu_c; // the gene abundance mean in carriers
   real sigma_c; // gene variation in carriers
   real sigma; // the variation in pools
   int k; // latent parameter
   // of the real population
   //
   // effects
 }
 transformed parameters  {
   // 
 }
 // functions{
 //                 ZIpool = function( mu_c , sigma_c , f_gene_carr, pool_size, sigma ) {
 //                                 // when outcome q_log is -Inf, then all samples where non-carriers
 //                                 // otherwise outcome is a mix of non-carriers and carriers with mu_c, sigma_c
 //                                 
 //                 }
 // }
 model {
   // Priors
   // (no need to specify if non-informative)
   f_gene_carr ~ uniform(0,1);

   // Likelihoods
   // distribution of zeros by bernoulli
   //for (i in 1:N) gene_zero[i] ~ bernoulli( (1-f_gene_carr)^pool_size );
   
   // for pools, which are mixed carriers and non-carriers
   // for (i in 1:N) {
   //                 q_log[i] ~ ZIpool( mu_c , sigma_c , f_gene_carr, pool_size, sigma ); 
   // }
   for (i in 1:N) k ~ binomial(  pool_size , f_gene_carr );
   
   // about only carriers
   for (i in 1:N_c) car_q_log[i] ~ normal( mu_c , sigma_c );
 
 }
 // generated quantities {
 //   // define predicted vector
 //     vector[N] q_pred_log;
 //     vector[N] log_lik;
 //     for (i in 1:N) q_pred_log[i] = bernoulli_rng( (1-f_gene_carr)^pool_size ) * normal_rng( mu_pool , sigma_pool ); 
 //     // make log likelihood
 //     for (i in 1:N) log_lik[i] = normal_lpdf( q_log[i] | mu_pool, sigma_pool );
 // }
   