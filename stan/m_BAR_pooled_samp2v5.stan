 functions {
                 // a zero inflated likelihood
                 real zifunc(vector q ,real f_gene_carr ,int pool_size , real mu_c , real sigma_c ) {
                                vector[num_elements(q)] prob;
                                real lprob;
                                for (i in 1:num_elements(q)) {
                                                real ll; 
                                                if (q[i] == 0) { 
                                                                ll = log((1 - f_gene_carr)^pool_size);
                                                } else {
                                                                for (k in 1:pool_size) {
                                                                   ll = binomial_lpmf( k | pool_size, f_gene_carr ) + 
                                                                   normal_lpdf( q[i] | mu_c*k , sigma_c*k );
                                                                }
                                                }
                                                prob[i] = ll;
                                }
                                lprob = sum(prob);
                                return lprob;
                }
 }
 data {
   // Define variables in data
   int<lower=0> N; // Number of observations (an integer)
   int<lower=0> N_c; // Number of obervations from gene carriers
   int<lower=1> pool_size; // number of samples in every pool
   //
   vector[N] q; // the measured gene quantities
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
   //
 }
 transformed parameters  {
   // 
 }
 model {
   // Priors (no need to specify if non-informative)
   f_gene_carr ~ uniform(0,1);
   // Likelihoods
   q ~ zifunc( f_gene_carr , pool_size , mu_c , sigma_c );
}
 // generated quantities {
 //   // define predicted vector
 //     vector[N] q_pred_log;
 //     vector[N] log_lik;
 //     for (i in 1:N) q_pred_log[i] = bernoulli_rng( (1-f_gene_carr)^pool_size ) * normal_rng( mu_pool , sigma_pool ); 
 //     // make log likelihood
 //     for (i in 1:N) log_lik[i] = normal_lpdf( q_log[i] | mu_pool, sigma_pool );
 // }