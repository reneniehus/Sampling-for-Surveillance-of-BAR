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
   real mu; // mean for gene distribution
   real sigma;// variance for gene distribution
   // of the real population
   real mu_pop;
   real sigma_pop;
   //
   matrix[n_gene_on,pool_size] patients;
   // effects
 }
 transformed parameters  {
   // 
   real sum_patients[n_gene_on]; // the gene sum in a patient pool
   real mean_patients[n_gene_on];
   for (i in 1:n_gene_on) {
                   // for each pool with gene
                   sum_patients[i] = 0;
                   for (p in 1:pool_size) {
                                   sum_patients[i] = sum_patients[i]+patients[i,p] ;
                   }
                   mean_patients[i] = sum_patients[i]/pool_size;
   }
 }
 model {
   // Prior part of Bayesian inference
   // (no need to specify if non-informative)
   f_gene_carr ~ uniform(0,1);
   mu ~ normal(4,3);
   mu_pop ~ normal(4,3);
   sigma ~ cauchy(0,2);
   sigma_pop ~ cauchy(0,2);             
   // Likelihood part of Bayesian inference
   // distribution of zeros by bernoulli
   for (i in 1:N) gene_zero[i] ~ bernoulli( (1-f_gene_carr)^pool_size );
   
   // distribution for where is some gene
   for (i in 1:n_gene_on) {
                   q_log[gene_on_loc[i]] ~ normal( mu , sigma ); 
                   mean_patients[i] ~ normal( mu , sigma );
                    
   }
   // imputation
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
     for (i in 1:N) q_pred_log[i] = bernoulli_rng( (1-f_gene_carr)^pool_size ) * normal_rng( mu , sigma ); 
     // make log likelihood
     for (i in 1:N) log_lik[i] = normal_lpdf( q_log[i] | mu, sigma );
 }
   