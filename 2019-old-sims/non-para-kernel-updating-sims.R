# non-parametric Bayesian updating by kernel densities
# simulation study; this script will run four chains on four cores
# and is intended to run in multiple parallel instances of R
# (I used digitalocean with 16 CPU cores to get 1000 simulations)

# Robert Grant, robert@bayescamp.com, @robertstats

# Presented at StanCon 2019, Cambridge, UK

library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
extract<-rstan::extract # avoid ambiguity


# change these for parallel instances of R
r_seed<-68963985
stan_master_seed<-92438196
setwd("~/instance_1")




###########################################################
# Section 1                                               #
# Define all Stan models and compile them                 #
###########################################################


stan_all_code <- '
data {
int n;
real x1[n];
real x2[n];
int y[n];
}
parameters {
real beta0;
real beta1;
real beta2;
}
model {
real mu[n];
real prob[n];
beta0 ~ normal(0,4);
beta1 ~ normal(0,4);
beta2 ~ normal(0,4);
for(i in 1:n) {
mu[i] = beta0 + beta1*x1[i] + beta2*x2[i];
prob[i] = exp(mu[i])/(1+exp(mu[i]));
target += bernoulli_lpmf(y[i] | prob[i]);
}
}
'
stan_all_fit_model <- stan_model(model_code=stan_all_code)

stan_next_code <- '
data {
int n;
real x1[n];
real x2[n];
int y[n];
vector[3] priormean;
matrix[3,3] priorcovar;
}
parameters {
vector[3] beta;
}
model {
real mu[n];
real prob[n];
beta ~ multi_normal(priormean,priorcovar);
for(i in 1:n) {
mu[i] = beta[1] + beta[2]*x1[i] + beta[3]*x2[i];
prob[i] = exp(mu[i])/(1+exp(mu[i]));
target += bernoulli_lpmf(y[i] | prob[i]);
}
}
'
stan_next_fit_model <- stan_model(model_code=stan_next_code)

stan_kernel_code <- '
data {
int n;
int k; // draws in the prior
int p; // columns in the prior; this should NOT include __lp, for example, which we will ignore
real x1[n];
real x2[n];
int y[n];
matrix[k,p] prior;
real bandwidth; // this is a common value for all parameters, but then so is stepsize!
}
parameters {
real beta[3]; // this is simpler as a vector
}
model {
real mu[n];
real prob[n];
real logpriorprob;
real priorprob;
// no prior specified here, so they are uniform
for(i in 1:n) {
mu[i] = beta[1] + beta[2]*x1[i] + beta[3]*x2[i];
prob[i] = exp(mu[i])/(1+exp(mu[i]));
// increment log-posterior by log-likelihood
target += bernoulli_lpmf(y[i] | prob[i]);
}
// increment log-posterior by log-prior
priorprob = 0.0;
for(i in 1:k) {
logpriorprob = 0.0;
for(j in 1:p) {
logpriorprob += normal_lpdf(beta[j] | prior[i,j], bandwidth); // multivariate distribution in parameter space
}
priorprob += exp(logpriorprob); // add all kernels together
}
priorprob = priorprob / (k+0.0); // average over k draws in the prior
target += log(priorprob); // add to log-likelihood!
}
'
stan_kernel_fit_model <- stan_model(model_code=stan_kernel_code)

stan_kernel_cube_code <- '
data {
int n;
int k; // number of non-zero-count hypercubes in the prior
int p; // columns in the prior; this should NOT include the count per hypercube
int<lower=1> count[k]; // number of draws in each hypercube
real x1[n];
real x2[n];
int y[n];
matrix[k,p] prior;
real bandwidth; // this is a common value for all parameters, but then so is stepsize!
}
transformed data {
int<lower=1> sum_count;
sum_count = sum(count);
}
parameters {
real beta[3]; // this is simpler as a vector
}
model {
real mu[n];
real prob[n];
real logcubeprob;
real priorprob;
// no prior specified here, so they are uniform
for(i in 1:n) {
mu[i] = beta[1] + beta[2]*x1[i] + beta[3]*x2[i];
prob[i] = exp(mu[i])/(1+exp(mu[i]));
// increment log-posterior by log-likelihood
target += bernoulli_lpmf(y[i] | prob[i]);
}
// increment log-posterior by log-prior
priorprob = 0.0;
for(i in 1:k) {
logcubeprob = 0.0;
for(j in 1:p) {
logcubeprob += normal_lpdf(beta[j] | prior[i,j], bandwidth); // multivariate distribution in parameter space
}
priorprob += exp(count[i]*logcubeprob); // add all kernels together, weighted by count
}
priorprob = priorprob / (sum_count+0.0); // weighted average over all draws in the prior
target += log(priorprob); // add to log-likelihood!
}
'
stan_kernel_cube_fit_model <- stan_model(model_code=stan_kernel_cube_code)






###########################################################
# Section 2                                               #
# Define function                                         #
###########################################################

# define all arguments here:
n<-200
batches<-10
n_of_chains<-4
true_beta0<-2
true_beta1<-0.5
true_beta2<-(-1.5)
r_seed<-370
stan_all_fit_seed<-97428
stan_first_fit_seed<-3383
stan_next_fit_seed<-473462355
stan_kernel_fit_seed<-599654765
stan_kernel_cube_fit_seed<-289110
n_of_draws_per_chain<-1000
bandwidth<-0.08
parameter_bins<-10
bw_scaling<-1
expit<-function(x) { exp(x)/(1+exp(x)) }





kernelsim<-function(n,
                    batches,
                    n_of_chains,
                    true_beta0,
                    true_beta1,
                    true_beta2,
                    r_seed,
                    stan_all_fit_seed,
                    stan_first_fit_seed,
                    stan_next_fit_seed,
                    stan_kernel_fit_seed,
                    stan_kernel_cube_fit_seed,
                    n_of_draws_per_chain,
                    bandwidth,
                    parameter_bins,
                    bw_scaling,
                    stan_all_fit_model,
                    stan_next_fit_model,
                    stan_kernel_fit_model,
                    stan_kernel_cube_fit_model) {
  
set.seed(r_seed)
# checks
if(n%%batches!=0) stop("n is not divisible by batches")

# generate batches of data
x1<-rnorm(n,0,1)
x2<-(x1/1.8)+rnorm(n,0,0.8)
y<-rbinom(n,size=1,prob=expit(true_beta0+(true_beta1*x1)+(true_beta2*x2)))
# the MLEs by iterative reweighted least squares:
mles<-summary(glm(y~x1+x2,family=binomial))$coefficients

# split the data into batches
data_batch<-NULL
for(i in 1:batches) {
  batch_start<-(i-1)*(n/batches)+1
  batch_end<-i*(n/batches)
  data_batch[[i]] <- list(x1=x1[batch_start:batch_end],
                          x2=x2[batch_start:batch_end],
                          y=y[batch_start:batch_end])
}

# what if we analysed all together (the "right" answer)?
stan_all_fit <- sampling(stan_all_fit_model,
                     data=list(n=n,
                               y=y,
                               x1=x1,
                               x2=x2),
                     chains=n_of_chains,
                     cores=4,
                     iter=5000,
                     warmup=2500,
                     seed=stan_all_fit_seed)
print(stan_all_fit)
stan_summary_all<-summary(stan_all_fit)$summary

# normal stats updating
stan_first_fit <- sampling(stan_all_fit_model,
                       data=list(n=n/batches,
                                 y=data_batch[[1]]$y,
                                 x1=data_batch[[1]]$x1,
                                 x2=data_batch[[1]]$x2),
                       chains=n_of_chains,
                       cores=4,
                       iter=5000,
                       warmup=2500,
                       seed=stan_first_fit_seed)
print(stan_first_fit)
stan_summary_stats<-array(NA,dim=c(3,2,batches))
stan_summary_stats[,,1]<-summary(stan_first_fit)$summary[1:3,c(1,3)]
next_prior_means<-stan_summary_stats[,1,1]
first_samples<-samples<-extract(stan_first_fit)
next_prior_covar<-cov(cbind(samples$beta0,samples$beta1,samples$beta2))
stan_normal_stats_time<-c(Sys.time(),"")
for(j in 2:batches) {
  stan_next_fit <- sampling(stan_next_fit_model,
                        data=list(n=n/batches,
                                  y=data_batch[[j]]$y,
                                  x1=data_batch[[j]]$x1,
                                  x2=data_batch[[j]]$x2,
                                  priormean=next_prior_means,
                                  priorcovar=next_prior_covar),
                        chains=n_of_chains,
                        cores=4,
                        iter=5000,
                        warmup=2500,
                        seed=stan_next_fit_seed+29*j)
  print(stan_next_fit)
  stan_summary_stats[,,j]<-summary(stan_next_fit)$summary[1:3,c(1,3)]
  next_prior_means<-stan_summary_stats[,1,j]
  samples<-extract(stan_next_fit)
  next_prior_covar<-cov(samples$beta)
}
stan_normal_stats_time[2]<-Sys.time()


# kernel densities for priors 
next_prior_draws<-cbind(first_samples$beta0,
                        first_samples$beta1,
                        first_samples$beta2)
next_prior_draws<-next_prior_draws[sample.int(n=NROW(next_prior_draws),
                                              size=n_of_chains*n_of_draws_per_chain),]
stan_summary_kernel<-array(NA,dim=c(3,2,batches))
stan_summary_kernel[,,1]<-summary(stan_first_fit)$summary[1:3,c(1,3)]
stan_kernel_time<-c(Sys.time(),"")
for(j in 2:batches) {
  stan_kernel_fit<-stan(model_code=stan_kernel_code,
                        data=list(n=n/batches,
                                  y=data_batch[[j]]$y,
                                  x1=data_batch[[j]]$x1,
                                  x2=data_batch[[j]]$x2,
                                  k=NROW(next_prior_draws),
                                  p=NCOL(next_prior_draws),
                                  prior=next_prior_draws,
                                  bandwidth=bandwidth),
                        chains=n_of_chains,
                        cores=4,
                        iter=2000+n_of_draws_per_chain,
                        warmup=2000,
                        seed=stan_kernel_fit_seed+991*j)
  print(stan_kernel_fit)
  stan_summary_kernel[,,j]<-summary(stan_kernel_fit)$summary[1:3,c(1,3)]
  next_prior_draws<-extract(stan_kernel_fit)$beta
}
stan_kernel_time[2]<-Sys.time()


# using hypercube counts of posterior draws 
bandwidth<-median(c((max(first_samples$beta0)-min(first_samples$beta0))/parameter_bins,
                    (max(first_samples$beta1)-min(first_samples$beta1))/parameter_bins,
                    (max(first_samples$beta2)-min(first_samples$beta2))/parameter_bins))*bw_scaling
parameter_bin_centres<-cbind(min(first_samples$beta0)+
                               ((max(first_samples$beta0)-min(first_samples$beta0))*
                                  ((1:parameter_bins)-0.5)/parameter_bins),
                             min(first_samples$beta1)+
                               ((max(first_samples$beta1)-min(first_samples$beta1))*
                                  ((1:parameter_bins)-0.5)/parameter_bins),
                             min(first_samples$beta2)+
                               ((max(first_samples$beta2)-min(first_samples$beta2))*
                                  ((1:parameter_bins)-0.5)/parameter_bins))
next_prior_draws<-cbind(cut(first_samples$beta0,breaks=parameter_bins,labels=FALSE),
                        cut(first_samples$beta1,breaks=parameter_bins,labels=FALSE),
                        cut(first_samples$beta2,breaks=parameter_bins,labels=FALSE))
next_prior_draws<-as.data.frame(next_prior_draws)
colnames(next_prior_draws)<-c("beta0_bins","beta1_bins","beta2_bins")
# adaptive bandwidth:
bandwidth<-median(c((max(next_prior_draws$beta0)-min(next_prior_draws$beta0))/parameter_bins,
                    (max(next_prior_draws$beta1)-min(next_prior_draws$beta1))/parameter_bins,
                    (max(next_prior_draws$beta2)-min(next_prior_draws$beta2))/parameter_bins))*bw_scaling
next_prior_draws %>% group_by(beta0_bins,beta1_bins,beta2_bins) %>% tally() -> next_prior_draws
cube_counts<-next_prior_draws$n
next_prior_draws<-next_prior_draws[,-4]
next_prior_draws[,1]<-parameter_bin_centres[as.numeric(next_prior_draws$beta0_bins),1]
next_prior_draws[,2]<-parameter_bin_centres[as.numeric(next_prior_draws$beta1_bins),2]
next_prior_draws[,3]<-parameter_bin_centres[as.numeric(next_prior_draws$beta2_bins),3]

stan_summary_kernel_cube<-array(NA,dim=c(3,2,batches))
stan_summary_kernel_cube[,,1]<-summary(stan_first_fit)$summary[1:3,c(1,3)]
stan_kernel_cube_time<-c(Sys.time(),"")
for(j in 2:batches) {
  stan_kernel_cube_fit<-sampling(stan_kernel_cube_fit_model,
                             data=list(n=n/batches,
                                       y=data_batch[[j]]$y,
                                       x1=data_batch[[j]]$x1,
                                       x2=data_batch[[j]]$x2,
                                       k=NROW(next_prior_draws),
                                       p=NCOL(next_prior_draws),
                                       prior=as.matrix(next_prior_draws),
                                       count=cube_counts,
                                       bandwidth=bandwidth),
                             chains=n_of_chains,
                             cores=4,
                             iter=2000+n_of_draws_per_chain,
                             warmup=2000,
                             seed=stan_kernel_cube_fit_seed+4577*j)
  print(stan_kernel_cube_fit)
  stan_summary_kernel_cube[,,j]<-summary(stan_kernel_cube_fit)$summary[1:3,c(1,3)]
  next_prior_draws<-extract(stan_kernel_cube_fit)$beta
  next_prior_draws<-as.data.frame(next_prior_draws)
  colnames(next_prior_draws)<-c("beta0_bins","beta1_bins","beta2_bins")
  # adaptive bandwidth:
  bandwidth<-median(c((max(next_prior_draws$beta0)-min(next_prior_draws$beta0))/parameter_bins,
                      (max(next_prior_draws$beta1)-min(next_prior_draws$beta1))/parameter_bins,
                      (max(next_prior_draws$beta2)-min(next_prior_draws$beta2))/parameter_bins))*bw_scaling
  next_prior_draws %>% group_by(beta0_bins,beta1_bins,beta2_bins) %>% tally() -> next_prior_draws
  cube_counts<-next_prior_draws$n
  next_prior_draws<-next_prior_draws[,-4]
}
stan_kernel_cube_time[2]<-Sys.time()

# compile all results
results<-list(mles=mles,
              stan_summary_all=stan_summary_all,
              stan_summary_stats=stan_summary_stats,
              stan_normal_stats_time=stan_normal_stats_time,
              stan_summary_kernel=stan_summary_kernel,
              stan_kernel_time=stan_kernel_time,
              stan_summary_kernel_cube=stan_summary_kernel_cube,
              stan_kernel_cube_time=stan_kernel_cube_time)
return(results)
}

iter<-200
for(k in 1:iter) {
  monitorfile<-file("monitorfile.txt","w")
  cat(paste("Now running iteration",k), file=monitorfile)
  close(monitorfile)
  r_seed<-round(r_seed*(1+(k/3579)))
  stan_master_seed<-round(stan_master_seed*(1+(k/4681)))
  all_results<-try( kernelsim(n=100,
                              batches=10,
                              n_of_chains=4,
                              true_beta0=2,
                              true_beta1=0.5,
                              true_beta2=(-1.5),
                              r_seed=r_seed,
                              stan_all_fit_seed=stan_master_seed,
                              stan_first_fit_seed=2*stan_master_seed,
                              stan_next_fit_seed=3*stan_master_seed,
                              stan_kernel_fit_seed=4*stan_master_seed,
                              stan_kernel_cube_fit_seed=5*stan_master_seed,
                              n_of_draws_per_chain=1000,
                              bandwidth=0.08,
                              parameter_bins=10,
                              bw_scaling=0.7,
                              stan_all_fit_model=stan_all_fit_model,
                              stan_next_fit_model=stan_next_fit_model,
                              stan_kernel_fit_model=stan_kernel_fit_model,
                              stan_kernel_cube_fit_model=stan_kernel_cube_fit_model))
  saveRDS(all_results,file=paste0("sim_results_",k,".rds"))
}
