 data {
   // Define variables in data
   int<lower=0> N; // Number of observations (an integer)
   int<lower=0> N_c; // Number of obervations from gene carriers
   int<lower=1> pool_size; // number of samples in every pool
   //
   real q[N]; // the measured gene quantities
   real car_q[N_c]; // gene quantities from carriers
   //
   int<lower=0 , upper=1> gene_on[N];// where gene is on
   real mu_c;
   real sigma_c;
   //
 }
 parameters {
   // Define parameters to estimate
   real<lower=0, upper=1> f_gene_carr;
   //real mu_c; // the gene abundance mean in carriers
   //real sigma_c; // gene variation in carriers
   // 
   //
   // 
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
   // Priors (no need to specify if non-informative)
   f_gene_carr ~ uniform(0,1);
   // Likelihoods
   for (i in 1:N) { // go through data points
                   real lp[pool_size]; // here you save the log likelihoods for all possible 
                   real sigma_p;
                   real mean_p;
                   // whether gene or not
                   if ( gene_on[i] == 0 ) { // when there is no gene
                                   target += ( binomial_lpmf( 0 | pool_size , f_gene_carr ) ); //could have calculated by hand
                   } else { // if there is gene
                          for (k in 1:pool_size){ // loop through how many are carriers
                          // the mean and sigma of the pooled samples
                          sigma_p = sigma_c*k;
                          mean_p = mu_c*k;
                          //
                                 lp[k] = binomial_lpmf( k | pool_size, f_gene_carr ) + 
                                 normal_lpdf( q[i] | mean_p , sigma_p );  
                                 //   wrong: (k/pool_size*(mu_c)) 
                                 // correct: log(k/pool_size*exp(mu_c))
                          }    
                          target += ( log_sum_exp(lp) ); //
                   }
   }
}
 // generated quantities {
 //   // define predicted vector
 //     vector[N] q_pred_log;
 //     vector[N] log_lik;
 //     for (i in 1:N) q_pred_log[i] = bernoulli_rng( (1-f_gene_carr)^pool_size ) * normal_rng( mu_pool , sigma_pool ); 
 //     // make log likelihood
 //     for (i in 1:N) log_lik[i] = normal_lpdf( q_log[i] | mu_pool, sigma_pool );
 // }