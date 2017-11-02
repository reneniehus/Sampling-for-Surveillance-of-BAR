library(rethinking)
library(tidyverse)
# simutate a target population N
N <- 100000 # size of N
N_v <- 1:N # a vector holding N
f_resA <- 40/100 # fraction of people in N that carry ab-res gene A
# in the gene carriers, the gene follows a log-normal distribution
# the parameters are derived from the Saturn study that looks at carriage of ctx-m gene
# the parameters are
# mean-log = 4.5
# sdlog = 4.2
curve( dnorm(x,7,3), from=-6, to = 15) # the distribution when x is the gene abundance on log-scale
# simulate the carriers
n_resA = N*f_resA # absolute number of carriers 
carr_v <- c( rep(1,n_resA) , rep(0,N-n_resA) )
q_gene <- rnorm( n=n_resA , mean=7 , sd=3 )
simplehist(q_gene) # look at the distribution
dens(q_gene) # look at the distribution
gene_v <- c( q_gene , rep(0,N-n_resA) )
# put data-frame together
df.pop <- data_frame( N_v , carr_v , gene_v )
df.pop$gene_v[] %>% simplehist() # this will look zero inflated from non-carriers

# 1 sampling: individual
n <- 10000 # how many patients we sample from
samp1 <-  sample_n(df.pop , size=n , replace=F)
# 2 sampling: pooled
pool_size = 5
samp2 <- sample_n( df.pop , size=n*pool_size , replace=F )
pool_id <- rep(1:n, each=pool_size)
samp2$pool_id <- pool_id
# samples are pooled, this makes averages of the gene quantity
sampl2.pooled <- samp2 %>% 
                group_by(pool_id) %>% 
                summarise(pooled_q = sum(gene_v))
# 3 sampling: n_c carriers
n_c = 20
sampl3.c <- sample_n( df.pop[df.pop$carr_v==1,] , n_c , replace=F )

# get data into stan format
#sampl2.pooled$gene_v_log <- log(sampl2.pooled$pooled_q)
gene_on <- as.numeric(sampl2.pooled$pooled_q != 0) # where in the vector the gene is on
#
sampl3.c$gene_v.log <- log(sampl3.c$gene_v)
library(rstan)
stan.d <- list(
                N = nrow(sampl2.pooled),
                pool_size = pool_size,
                q = sampl2.pooled$pooled_q,
                n_gene_on = sum(gene_on),
                gene_on = gene_on,
                N_c = n_c,
                car_q = sampl3.c$gene_v,
                mu_c = 7,
                sigma_c = 3
)
# set up stan
rstan_options( auto_write = TRUE )
options(mc.cores = parallel::detectCores() )
# run stan
mod.1 <- stan( "stan/m_BAR_pooled_samp2v4.stan" , data=stan.d ,
              chains=4 , iter=4000 , warmup=1000 , thin=10 )
stan_trace(mod.1, inc_warmup = F)
library(rethinking)
precis(mod.1) %>% plot(pars=c('f_gene_carr'),xlim=c(0.2,0.6))
mod.1 %>% print(pars=c('f_gene_carr'))
mod.1 %>% plot(pars=c('f_gene_carr'))

## do the same using rethinking
model <- alist(
                q ~ dzin(f_gene_carr,pool_size,mean_p,sigma_p),
                logit(f_gene_carr) <- f_interc,
                f_interc ~ dnorm(0,1)
)
# define my own likelihood. I could turn q == -Inf where no gene and change that in if-statement
dzin <- function(q,f_gene_carr,pool_size,mean_p,sigma_p,log=TRUE){
                ll <- 0
                if (q == 0){
                                ll <- dbinom(0,pool_size,f_gene_carr)
                } else {
                                lp <- rep(NA,pool_size)
                                for(k in 1:pool_size){
                                                lp[k] <- dbinom(k,pool_size,f_gene_carr) * dnorm(q,mean_p*k,sigma_p*k)
                                }
                                ll <- sum(lp)
                }
                if (log == TRUE) ll <- log(ll)
                return(ll)
}
# 

#


