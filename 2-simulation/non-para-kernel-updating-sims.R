# non-parametric Bayesian updating by kernel densities
# this script will run four chains on four cores
#

# Robert Grant, robert@bayescamp.com
# (c) BayesCamp Ltd, CC-ND-SA-BY-3.0

setwd('~/instance1')

library(dplyr)
library(diptest)
library(rstan)
rstan_options(auto_write = TRUE)
extract<-rstan::extract # avoid ambiguity
options(mc.cores = 4)


# change these for parallel instances of R
r_overall_seed<-   700000
stan_overall_seed<-900000


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
real bandwidth[p];
}
parameters {
real beta[p]; // this is simpler as a vector
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
logpriorprob += normal_lpdf(beta[j] | prior[i,j], bandwidth[j]);
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
int k; // number of non-zero-count hypercubes (bins) in the prior
int k_draws; // draws in the prior
int p; // columns in the prior; this should NOT include the count per hypercube
int<lower=1> count[k]; // number of draws in each hypercube
real x1[n];
real x2[n];
int y[n];
matrix[k,p] prior;
real bandwidth[p];
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
logcubeprob += normal_lpdf(beta[j] | prior[i,j], bandwidth[j]); // multivariate distribution in parameter space
}
priorprob += count[i]*exp(logcubeprob); // add all kernels together, weighted by count
}
priorprob = priorprob / (sum_count+0.0); // weighted average over all draws in the prior
target += log(priorprob); // add to log-likelihood!
}
'
stan_kernel_cube_fit_model <- stan_model(model_code=stan_kernel_cube_code)



stan_beta_code <- '
functions{
int count_neighbours(int k, int p, real[] proposal, matrix prior, real[] beta_widths) {
  // identify neighbour prior draws
  int n_neighbours=0;
  int neighbour_p=0;
  for(i in 1:k) {
    neighbour_p=0;
    for(j in 1:p) {
      if(fabs(prior[i,j]-proposal[j])<(0.99*beta_widths[j])) neighbour_p+=1;
      // 0.99 because we dont want rounding error to nudge a prior draw out of the domain of beta_ldpf()
    }
    if(neighbour_p==p) {
      n_neighbours+=1;
    }
  }
return(n_neighbours);
}

int[] find_neighbours(int k, int p, real[] proposal, matrix prior, real[] beta_widths, int n_neighbours) {
  // identify neighbour prior draws
  int loopcount=0;
  int neighbour_prior[n_neighbours];
  for(i in 1:k) {
    int neighbour_p=0;
    for(j in 1:p) {
      if(fabs(prior[i,j]-proposal[j])<(0.99*beta_widths[j])) neighbour_p+=1;
      // 0.99 because we dont want rounding error to nudge a prior draw out of the domain of beta_ldpf()
    }
    if(neighbour_p==p) {
      loopcount+=1;
      neighbour_prior[loopcount]=i;
    }
  }
  return(neighbour_prior);
}
}

data {
  int n;
  int k; // draws in the prior (in the hypercuboid version, including the phony one)
  int p; // columns in the prior; this should NOT include __lp
  real<lower=0> betapars; // shared alpha and beta parameters for the beta distribution
  real<lower=0> beta_widths[p];
  matrix[k,p] prior;
  real x1[n];
  real x2[n];
  int y[n];
}

parameters {
  real beta[3]; // this is simpler as a vector
}

model {
  real mu[n];
  real prob[n];
  real logpriorprob;
  real priorprob;
  int n_neighbours = count_neighbours(k, p, beta, prior, beta_widths);
  int neighbour_prior[n_neighbours] = find_neighbours(k, p, beta, prior, beta_widths, n_neighbours);

  for(i in 1:n) {
    mu[i] = beta[1] + beta[2]*x1[i] + beta[3]*x2[i];
    prob[i] = exp(mu[i])/(1+exp(mu[i]));
    // increment log-posterior by log-likelihood
    target += bernoulli_lpmf(y[i] | prob[i]);
  }
  // increment log-posterior by log-prior
  priorprob = 0.0;
  for(i in 1:n_neighbours) {
    logpriorprob = 0.0;
    for(j in 1:p) {
      logpriorprob += beta_lpdf(0.5+(beta[j]-prior[neighbour_prior[i],j])/(2*beta_widths[j]) | betapars, betapars);
    }
    priorprob += exp(logpriorprob); // add all kernels densities (not log-densities) together
  }
  priorprob = priorprob / (k+0.0); // average over k draws in the prior
  target += log(priorprob); // add log-prior to log-likelihood!
}
'
stan_beta_fit_model <- stan_model(model_code=stan_beta_code)







###########################################################
# Section 2                                               #
# Define function                                         #
###########################################################

# define all arguments here:
n<-2000 # n of data observations
batches<-10 # batches of data
n_of_chains<-4
true_beta0<-2
true_beta1<-0.5
true_beta2<-(-1.5)
r_seed<-r_overall_seed+10000
stan_first_fit_seed<-stan_overall_seed+10000
stan_next_fit_seed<-stan_overall_seed+20000
stan_kernel_fit_seed<-stan_overall_seed+30000
stan_kernel_cube_fit_seed<-stan_overall_seed+40000
stan_beta_fit_seed<-stan_overall_seed+50000
n_of_draws_per_chain<-1000 # sample this many for kernel method
n_of_draws_per_chain_beta<-1000 # sample for beta kernel method
bandwidth<-0.08 # kernel Gaussian bandwidth for hypercube centroids
parameter_bins<-10 # number of hypercubes across each dimension in parameter space
bw_scaling<-1 # controls hypercube kernel size relative to posterior
expit<-function(x) { exp(x)/(1+exp(x)) } # inverse logistic function
betapars<-8 # shared alpha and beta parameters for the beta distribution
currentiter<-1 # this is just used for live progress display

# # overwrite with seeds for a dodgy simulation
# badoldsim<-readRDS('bad_sim_results_3.rds')
# r_seed<-badoldsim$r_seed
# stan_all_fit_seed<-badoldsim$stan_all_fit_seed
# stan_first_fit_seed<-badoldsim$stan_first_fit_seed
# stan_next_fit_seed<-badoldsim$stan_next_fit_seed
# stan_kernel_fit_seed<-badoldsim$stan_kernel_fit_seed
# stan_kernel_cube_fit_seed<-badoldsim$stan_kernel_cube_fit_seed


# initial draw picking function
init_draw<-function(n_chains,prior_draws) {
  init<-list(list(beta=as.numeric(prior_draws[sample(1:NROW(prior_draws),1),])))
  if(n_chains>1) {
    for(i in 2:n_chains) {
      init[[i]]<-list(beta=as.numeric(prior_draws[sample(1:NROW(prior_draws),1),]))
    }
  }
  return(init)
}


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
                    stan_beta_fit_seed,
                    n_of_draws_per_chain,
                    n_of_draws_per_chain_beta,
                    parameter_bins,
                    bw_scaling,
                    stan_all_fit_model,
                    stan_next_fit_model,
                    stan_kernel_fit_model,
                    stan_kernel_cube_fit_model,
                    stan_beta_fit_model,
                    betapars,
                    currentiter) {

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


######################################################
# what if we analysed all together (the "right" answer)?
stan_all_time<-list(Sys.time())
# note that this counts all data, the others all but the first batch
cat(paste("now running simulation number",currentiter,"-- all data model\n"))
stan_all_fit <- sampling(stan_all_fit_model,
                     data=list(n=n,
                               y=y,
                               x1=x1,
                               x2=x2),
                     chains=n_of_chains,
                     cores=4,
                     iter=2500+n_of_draws_per_chain,
                     warmup=2500,
                     seed=stan_all_fit_seed)
print(stan_all_fit)
stan_summary_all<-list(summary=summary(stan_all_fit)$summary,
                       posterior=extract(stan_all_fit))
stan_all_time[[2]]<-Sys.time()

# we monitor these lists:
log_normal<-log_kernel<-log_hypercube<-log_beta<-NULL




######################################################
# normal stats updating
stan_first_fit <- sampling(stan_all_fit_model,
                       data=list(n=n/batches,
                                 y=data_batch[[1]]$y,
                                 x1=data_batch[[1]]$x1,
                                 x2=data_batch[[1]]$x2),
                       chains=n_of_chains,
                       cores=4,
                       iter=2500+n_of_draws_per_chain,
                       warmup=2500,
                       seed=stan_first_fit_seed)
print(stan_first_fit)
stan_summary_stats<-array(NA,dim=c(3,2,batches))
stan_summary_stats[,,1]<-summary(stan_first_fit)$summary[1:3,c(1,3)]
next_prior_means<-stan_summary_stats[,1,1]
# this is used as the starting sample for kernel, hypercube and wavelet:
first_samples<-samples<-extract(stan_first_fit)
next_prior_covar<-cov(cbind(samples$beta0,samples$beta1,samples$beta2))

# log for debug:
log_normal[[1]]<-list(next_prior_means=next_prior_means,
                      next_prior_covar=next_prior_covar,
                      posterior=samples)

stan_normal_stats_time<-list(Sys.time())
for(j in 2:batches) {
  cat(paste("now running simulation number",currentiter,
            "-- normal stats model -- batch",j,"of",batches,"\n"))
  stan_next_fit <- sampling(stan_next_fit_model,
                        data=list(n=n/batches,
                                  y=data_batch[[j]]$y,
                                  x1=data_batch[[j]]$x1,
                                  x2=data_batch[[j]]$x2,
                                  priormean=next_prior_means,
                                  priorcovar=next_prior_covar),
                        chains=n_of_chains,
                        cores=4,
                        iter=2500+n_of_draws_per_chain,
                        warmup=2500,
                        seed=stan_next_fit_seed+29*j)
  print(stan_next_fit)
  stan_summary_stats[,,j]<-summary(stan_next_fit)$summary[1:3,c(1,3)]
  next_prior_means<-stan_summary_stats[,1,j]
  samples<-extract(stan_next_fit)
  next_prior_covar<-cov(samples$beta)
  # log for debug:
  log_normal[[j]]<-list(next_prior_means=next_prior_means,
                        next_prior_covar=next_prior_covar,
                        posterior=samples)

}
stan_normal_stats_time[[2]]<-Sys.time()




######################################################
# kernel densities for priors
next_prior_draws<-cbind(first_samples$beta0,
                        first_samples$beta1,
                        first_samples$beta2)
next_prior_draws<-next_prior_draws[sample.int(n=NROW(next_prior_draws),
                                              size=n_of_chains*n_of_draws_per_chain),]

kernel_setup<-function(n_of_chains,
                       next_prior_draws,
                       hartigan_p_threshold=0.05,
                       multimodal_boost=1.5) {
  p<-ncol(next_prior_draws)
  # get Sheather-Jones bandwidth
  bandwidth<-apply(next_prior_draws,2,bw.SJ)
  #bandwidth<-median(parameter_bandwidth)
  # check for multimodality
  # loop over the parameters
  dipstats<-rep(NA,3)
  loopcount<-1
  dip_p<-rep(0.01,3)
  while(loopcount<100 & min(dip_p)<hartigan_p_threshold) {
    # get maximum Hartigans' dip statistic and increase bandwidth until non-sig
    for(i in 1:3) {
      tempdens<-density(next_prior_draws[,i],bw=bandwidth[i])
      tempscale<-n_of_chains*n_of_draws_per_chain/round(sum(tempdens$y))
      dip_p[i]<-unlist(diptest::dip.test(rep(tempdens$x,round(tempscale*tempdens$y)))$p.value)
      if(loopcount>1 & dip_p[i]<hartigan_p_threshold) { bandwidth[i]<-bandwidth[i]*multimodal_boost }
    }
    loopcount<-loopcount+1
  }
  if(loopcount>99) { stop("Unable to find unimodal univariate kernel bandwidth after 100 expansions") }

  # get initial values
  inits<-init_draw(n_of_chains,next_prior_draws)

  return(list(bandwidth=bandwidth,
              loopcount=loopcount,
              inits=inits))
}

# setup after first batch of data
kernel_vars<-kernel_setup(n_of_chains,next_prior_draws)

# log for debug:
log_kernel[[1]]<-list(next_prior_draws=next_prior_draws,
                      bandwidth=kernel_vars$bandwidth,
                      loopcount=kernel_vars$loopcount,
                      inits=kernel_vars$inits)

stan_summary_kernel<-array(NA,dim=c(3,2,batches))
stan_summary_kernel[,,1]<-summary(stan_first_fit)$summary[1:3,c(1,3)]
stan_kernel_time<-list(Sys.time())
for(j in 2:batches) {
  cat(paste("now running simulation number",currentiter,
            "-- kernel model -- batch",j,"of",batches,"\n"))

  stan_kernel_fit<-stan(model_code=stan_kernel_code,
                        data=list(n=n/batches,
                                  y=data_batch[[j]]$y,
                                  x1=data_batch[[j]]$x1,
                                  x2=data_batch[[j]]$x2,
                                  k=NROW(next_prior_draws),
                                  p=NCOL(next_prior_draws),
                                  prior=next_prior_draws,
                                  bandwidth=kernel_vars$bandwidth),
                        chains=n_of_chains,
                        cores=4,
                        iter=2500+n_of_draws_per_chain,
                        warmup=2500,
                        seed=stan_kernel_fit_seed+991*j,
                        init=kernel_vars$inits)
  print(stan_kernel_fit)
  stan_summary_kernel[,,j]<-summary(stan_kernel_fit)$summary[1:3,c(1,3)]
  next_prior_draws<-extract(stan_kernel_fit)$beta
  kernel_vars<-kernel_setup(n_of_chains,next_prior_draws)
  # log for debug:
  log_kernel[[j]]<-list(next_prior_draws=next_prior_draws,
                        bandwidth=kernel_vars$bandwidth,
                        loopcount=kernel_vars$loopcount,
                        inits=kernel_vars$inits)

}
stan_kernel_time[[2]]<-Sys.time()



######################################################

# using hypercube counts of posterior draws

## deprecated bandwidth calculation
# bandwidth<-median(c((max(first_samples$beta0)-min(first_samples$beta0))/parameter_bins,
#                     (max(first_samples$beta1)-min(first_samples$beta1))/parameter_bins,
#                     (max(first_samples$beta2)-min(first_samples$beta2))/parameter_bins))*bw_scaling



hypercube_setup<-function(n_of_chains,
                          parameter_bins,
                          next_prior_draws,
                          hartigan_p_threshold=0.05,
                          multimodal_boost=1.5) {

  hypercube_posterior<-as.data.frame(next_prior_draws)

  # divide into parameter bins
  parameter_bin_centres<-cbind(min(next_prior_draws[,1])+
                                 ((max(next_prior_draws[,1])-min(next_prior_draws[,1]))*
                                    ((1:parameter_bins)-0.5)/parameter_bins),
                               min(next_prior_draws[,2])+
                                 ((max(next_prior_draws[,2])-min(next_prior_draws[,2]))*
                                    ((1:parameter_bins)-0.5)/parameter_bins),
                               min(next_prior_draws[,3])+
                                 ((max(next_prior_draws[,3])-min(next_prior_draws[,3]))*
                                    ((1:parameter_bins)-0.5)/parameter_bins))
  next_prior_draws<-cbind(cut(next_prior_draws[,1],breaks=parameter_bins,labels=FALSE),
                          cut(next_prior_draws[,2],breaks=parameter_bins,labels=FALSE),
                          cut(next_prior_draws[,3],breaks=parameter_bins,labels=FALSE))
  p<-ncol(next_prior_draws)
  next_prior_draws<-cbind(parameter_bin_centres[as.numeric(next_prior_draws[,1]),1],
                         parameter_bin_centres[as.numeric(next_prior_draws[,2]),2],
                         parameter_bin_centres[as.numeric(next_prior_draws[,3]),3])
  next_prior_draws<-as.data.frame(next_prior_draws)
  colnames(next_prior_draws)<-c("beta0_bins","beta1_bins","beta2_bins")

  # the binned version
  next_prior_bins<-tally(group_by(next_prior_draws,beta0_bins,beta1_bins,beta2_bins))
  cube_counts<-next_prior_bins$n
  next_prior_bins<-select(next_prior_bins,-n)

  # get Sheather-Jones bandwidth
  bandwidth<-apply(next_prior_draws,2,bw.SJ)
  # check for multimodality
  # loop over the parameters
  dipstats<-rep(NA,3)
  loopcount<-1
  dip_p<-rep(0.01,3)
  while(loopcount<100 & min(dip_p)<hartigan_p_threshold) {
    # get maximum Hartigans' dip statistic and increase bandwidth until non-sig
    for(i in 1:3) {
      tempdens<-density(next_prior_draws[,i],bw=bandwidth[i])
      tempscale<-n_of_chains*n_of_draws_per_chain/round(sum(tempdens$y))
      dip_p[i]<-unlist(diptest::dip.test(rep(tempdens$x,round(tempscale*tempdens$y)))$p.value)
      if(loopcount>1 & dip_p[i]<hartigan_p_threshold) { bandwidth[i]<-bandwidth[i]*multimodal_boost }
    }
    loopcount<-loopcount+1
  }
  if(loopcount>99) { stop("Unable to find unimodal univariate kernel bandwidth after 100 expansions") }

  return(list(next_prior_draws=next_prior_draws,
              next_prior_bins=next_prior_bins,
              cube_counts=cube_counts,
              parameter_bin_centres=parameter_bin_centres,
              bandwidth=bandwidth,
              loopcount=loopcount,
              hypercube_posterior=hypercube_posterior))
}

# setup after first batch of data
next_prior_draws<-cbind(first_samples$beta0,
                        first_samples$beta1,
                        first_samples$beta2)
next_prior_draws<-next_prior_draws[sample.int(n=NROW(next_prior_draws),
                                              size=n_of_chains*n_of_draws_per_chain),]
hypercube_vars<-hypercube_setup(n_of_chains,
                                parameter_bins,
                                next_prior_draws)
next_prior_draws<-hypercube_vars$next_prior_bins # replace with binned version

# # visual check
# for(par in 1:3) {
#   plot(density(rep(unlist(next_prior_draws[,par]),hypercube_vars$cube_counts),
#                bw=hypercube_vars$bandwidth[par]))
# }

# get initial values
hypercube_inits<-init_draw(n_of_chains,next_prior_draws)

# log for debug:
log_hypercube[[1]]<-list(hypercube_vars=hypercube_vars,
                         inits=hypercube_inits)

stan_summary_kernel_cube<-array(NA,dim=c(3,2,batches))
stan_summary_kernel_cube[,,1]<-summary(stan_first_fit)$summary[1:3,c(1,3)]
stan_kernel_cube_time<-list(Sys.time())
for(j in 2:batches) {
  cat(paste("now running simulation number",currentiter,
            "-- hypercube model -- batch",j,"of",batches,"\n"))

  stan_kernel_cube_fit<-sampling(stan_kernel_cube_fit_model,
                             data=list(n=n/batches,
                                       y=data_batch[[j]]$y,
                                       x1=data_batch[[j]]$x1,
                                       x2=data_batch[[j]]$x2,
                                       k=NROW(next_prior_draws),
                                       k_draws=n_of_chains*n_of_draws_per_chain,
                                       p=NCOL(next_prior_draws),
                                       prior=as.matrix(next_prior_draws),
                                       count=hypercube_vars$cube_counts,
                                       bandwidth=hypercube_vars$bandwidth),
                             chains=n_of_chains,
                             cores=4,
                             iter=2500+n_of_draws_per_chain,
                             warmup=2500,
                             seed=stan_kernel_cube_fit_seed+4577*j,
                             init=hypercube_inits)
  print(stan_kernel_cube_fit)
  stan_summary_kernel_cube[,,j]<-summary(stan_kernel_cube_fit)$summary[1:3,c(1,3)]
  next_prior_draws<-extract(stan_kernel_cube_fit)$beta
  next_prior_draws<-next_prior_draws[sample.int(n=NROW(next_prior_draws),
                                                size=n_of_chains*n_of_draws_per_chain),]

  hypercube_vars<-hypercube_setup(n_of_chains,
                                  parameter_bins,
                                  next_prior_draws)
  next_prior_draws<-hypercube_vars$next_prior_bins # replace with binned version

  # get initial values
  hypercube_inits<-init_draw(n_of_chains,next_prior_draws)

  # log for debug:
  log_hypercube[[j]]<-list(hypercube_vars=hypercube_vars,
                           inits=hypercube_inits)

}
stan_kernel_cube_time[[2]]<-Sys.time()




######################################################
# multivariate beta kernel densities for priors

# get posterior from first batch (and possibly down-sample)
next_prior_draws<-cbind(first_samples$beta0,
                        first_samples$beta1,
                        first_samples$beta2)
next_prior_draws<-next_prior_draws[sample.int(n=NROW(next_prior_draws),
                                              size=n_of_chains*
                                                n_of_draws_per_chain_beta),]

beta_setup<-function(n_of_chains,
                     next_prior_draws,
                     hartigan_p_threshold=0.05,
                     multimodal_boost=1.5) {
  p<-ncol(next_prior_draws)
  # get Sheather-Jones bandwidths
  bandwidth<-apply(next_prior_draws,2,bw.SJ)

  # check for multimodality
  # loop over the parameters
  dipstats<-rep(NA,p)
  loopcount<-0
  dip_p<-rep(0.01,p)
  while(loopcount<100 & min(dip_p)<hartigan_p_threshold) {
    # get maximum Hartigans' dip statistic and increase bandwidth by 50% until non-sig
    for(i in 1:p) {
      tempdens<-density(next_prior_draws[,i],bw=bandwidth[i])
      tempscale<-n_of_chains*n_of_draws_per_chain/round(sum(tempdens$y))
      dip_p[i]<-unlist(diptest::dip.test(rep(tempdens$x,round(tempscale*tempdens$y)))$p.value)
    }
    if(loopcount>0) { bandwidth<-bandwidth*multimodal_boost }
    loopcount<-loopcount+1
  }
  if(loopcount>99) { stop("Unable to find unimodal univariate kernel bandwidth after 100 expansions") }

  # convert to beta width
  beta_widths<-bandwidth*6 # ±3 SDs

  # get initial values
  inits<-init_draw(n_of_chains,next_prior_draws)



  return(list(next_prior_draws=next_prior_draws,
              bandwidth=bandwidth,
              loopcount=loopcount,
              beta_widths=beta_widths,
              inits=inits))
}

# setup after first batch of data
beta_vars<-beta_setup(n_of_chains,next_prior_draws)

# log for debug:
log_beta[[1]]<-beta_vars

stan_summary_beta<-array(NA,dim=c(3,2,batches))
stan_summary_beta[,,1]<-summary(stan_first_fit)$summary[1:3,c(1,3)]
stan_beta_time<-list(Sys.time())

# loop over updates
for(j in 2:batches) {
  cat(paste("now running simulation number",currentiter,
            "-- beta model -- batch",j,"of",batches,"\n"))

  stan_beta_fit<-sampling(stan_beta_fit_model,
                         data=list(n=n/batches,
                                   y=data_batch[[j]]$y,
                                   x1=data_batch[[j]]$x1,
                                   x2=data_batch[[j]]$x2,
                                   k=NROW(beta_vars$next_prior_draws),
                                   p=NCOL(beta_vars$next_prior_draws),
                                   beta_widths=beta_vars$beta_widths,
                                   prior=beta_vars$next_prior_draws,
                                   betapars=betapars),
                          chains=n_of_chains,
                          cores=4,
                          iter=2500+n_of_draws_per_chain,
                          warmup=2500,
                          seed=stan_beta_fit_seed+968547*j,
                          init=beta_vars$inits)

  stan_summary_beta[,,j]<-summary(stan_beta_fit)$summary[1:3,c(1,3)]
  next_prior_draws<-extract(stan_beta_fit)$beta
  next_prior_draws<-as.data.frame(next_prior_draws)

  beta_vars<-beta_setup(n_of_chains,next_prior_draws)

  # log for debug:
  log_beta[[j]]<-beta_vars

}
stan_beta_time[[2]]<-Sys.time()






######################################################
# compile all results
results<-list(mles=mles,
              stan_summary_all=stan_summary_all,
              stan_all_time=stan_all_time,
              stan_summary_stats=stan_summary_stats,
              stan_normal_stats_time=stan_normal_stats_time,
              stan_summary_kernel=stan_summary_kernel,
              stan_kernel_time=stan_kernel_time,
              stan_summary_kernel_cube=stan_summary_kernel_cube,
              stan_kernel_cube_time=stan_kernel_cube_time,
              stan_summary_beta=stan_summary_beta,
              stan_beta_time=stan_beta_time,
              r_seed=r_seed,
              stan_all_fit_seed=stan_all_fit_seed,
              stan_first_fit_seed=stan_first_fit_seed,
              stan_next_fit_seed=stan_next_fit_seed,
              stan_kernel_fit_seed=stan_kernel_fit_seed,
              stan_kernel_cube_fit_seed=stan_kernel_cube_fit_seed,
              data_batches=data_batch,
              log_normal=log_normal,
              log_kernel=log_kernel,
              log_hypercube=log_hypercube,
              log_beta=log_beta)
return(results)
}



###########################################################
# Section 3                                               #
# Run the simulation                                      #
###########################################################


iter<-334
for(k in 1:iter) {
  monitorfile<-file("monitorfile.txt","w")
  cat(paste("Now running iteration",k), file=monitorfile)
  close(monitorfile)
  all_results<-try( kernelsim(n=2000,
                              batches=10,
                              n_of_chains=4,
                              true_beta0=2,
                              true_beta1=0.5,
                              true_beta2=(-1.5),
                              r_seed=r_overall_seed+10000+(3*k),
                              stan_all_fit_seed=stan_overall_seed+10000+(3*k),
                              stan_first_fit_seed=stan_overall_seed+20000+(3*k),
                              stan_next_fit_seed=stan_overall_seed+30000+(3*k),
                              stan_kernel_fit_seed=stan_overall_seed+40000+(3*k),
                              stan_kernel_cube_fit_seed=stan_overall_seed+50000+(3*k),
                              stan_beta_fit_seed=stan_overall_seed+60000+(3*k),
                              n_of_draws_per_chain=1000,
                              n_of_draws_per_chain_beta=1000,
                              parameter_bins=10,
                              bw_scaling=0.7,
                              stan_all_fit_model=stan_all_fit_model,
                              stan_next_fit_model=stan_next_fit_model,
                              stan_kernel_fit_model=stan_kernel_fit_model,
                              stan_kernel_cube_fit_model=stan_kernel_cube_fit_model,
                              stan_beta_fit_model=stan_beta_fit_model,
                              betapars=betapars,
                              currentiter=k))
  saveRDS(all_results,file=paste0("sim_results_",k,".rds"))
}
