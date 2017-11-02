 data {
   // Define variables in data
   int<lower=0> N; // Number of observations (an integer)
   int<lower=1> pool_size; // number of samples in every pool
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
   real mu_pool;// mean from the pools
   real sigma_pool;// variance for gene distribution
   // of the real population
   //
   matrix[n_gene_on,pool_size] patients;
   // effects
 }
 transformed parameters  {
   // 
   real sigma_pop; // the gene sum in a patient pool
   real mu_pop; // the gene mean in the real population
   sigma_pop = sigma_pool*pool_size; // according to Caudill 2010
   mu_pop = mu_pool;
 }
 model {
   // Priors
   // (no need to specify if non-informative)
   f_gene_carr ~ uniform(0,1);
   sigma_pool ~ cauchy(0,2);

   // Likelihoods
   // distribution of zeros by bernoulli
   for (i in 1:N) gene_zero[i] ~ bernoulli( (1-f_gene_carr)^pool_size );
   // distribution for where is gene
   // the mean of the pools is the mean of the real population
   // this is not true when 
   for (i in 1:n_gene_on) {
                   q_log[gene_on_loc[i]] ~ normal( mu_pool , sigma_pool ); 
   }
   // imputation into the entire patient population
   // on the lowest level, individual patients follow mean from each pool, and a shared variance.
   // is this variance identifiable from the data about the means that we have?
   for (i in 1:n_gene_on) {
                   for (p in 1:pool_size) patients[i,p] ~ normal( q_log[gene_on_loc[i]] , sigma_pop );
   }
 }
 generated quantities {
   // define predicted vector
     vector[N] q_pred_log;
     vector[N] log_lik;
     for (i in 1:N) q_pred_log[i] = bernoulli_rng( (1-f_gene_carr)^pool_size ) * normal_rng( mu_pool , sigma_pool ); 
     // make log likelihood
     for (i in 1:N) log_lik[i] = normal_lpdf( q_log[i] | mu_pool, sigma_pool );
 }
   