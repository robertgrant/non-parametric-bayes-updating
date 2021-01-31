# Non-parametric Bayesian updating
# demonstrations with prior-only models

setwd('~/Documents/Shop-window-NOC/non-parametric-updating/2020-sims/prior-only')
library(dplyr)
library(ggplot2)
library(rstan)
options(mc.cores=4)
rstan_options(auto_write = TRUE)
overall_seed<-80955000
sj_coef_option<-0 # change this for different bandwidths (see findbw() below)


# running order:
#   three distributions, all in 2 parameters:
#     well-behaved MVN
#     banana
#     multimodal (LASSO-style)
#   for each distribution:
#     set up first "prior" sample by drawing from distributions
#     for each batch:
#       supply prior draws in the usual way
#       add either a constant, or no, log-likelihood
#       store posterior draws for the next batch


# define distribution functions
logd_mvn <- function(theta) {
  if(length(theta)!=2) { stop("parameter vector is not of length 2")}
  mnormt::dmnorm(x=theta,
                 mean=rep(0,2),
                 varcov=matrix(c(1.0,0.6,
                                 0.6,1.0),nrow=2)^2,
                 log=TRUE)
}
r_mvn <- function(n) {
  mnormt::rmnorm(n=n,
                 mean=rep(0,2),
                 varcov=matrix(c(1.0,0.6,
                                 0.6,1.0),nrow=2)^2)
}

logd_banana <- function(theta) {
  if(length(theta)!=2) { stop("parameter vector is not of length 2")}
  mtheta<-t(matrix(theta)) # convert to row matrix
  rotangle<-(-0.25)
  rot<-matrix(c(cos(rotangle),
                sin(rotangle),
                -sin(rotangle),
                cos(rotangle)),nrow=2)
  rtheta<-t(rot%*%theta)
  res<-rtheta[2]+rtheta[1]^2
  return(dnorm(rtheta[1],0,1,log=TRUE)+
           dnorm(res,0,0.7,log=TRUE))
}
r_banana <- function(n) {
  dx<-rnorm(n,0,1)
  dy<-(-1)*dx^2+rnorm(n,0,0.7)
  rotangle<-0.25
  rot<-matrix(c(cos(rotangle),
                sin(rotangle),
                -sin(rotangle),
                cos(rotangle)),nrow=2)
  return(t(rot%*%t(cbind(dx,dy))))
}

logd_lasso <- function(theta) {
  if(length(theta)!=2) { stop("parameter vector is not of length 2")}
  return(log((0.5*mnormt::dmnorm(x=theta,
                                 mean=c(1,0),
                                 varcov=matrix(c(0.3,0.0,
                                                 0.0,0.06),nrow=2)^2)) +
               (0.5*mnormt::dmnorm(x=theta,
                                  mean=c(0,1),
                                  varcov=matrix(c(0.06,0.0,
                                                  0.0,0.3),nrow=2)^2))))
}
r_lasso <- function(n) {
  temp<-matrix(NA,nrow=n,ncol=2)
  for(i in 1:n) {
    yaxisdist<-rbinom(1,1,0.5)
    if(yaxisdist==1) {
      temp[i,]<-mnormt::rmnorm(n=1,
                            mean=c(1,0),
                            varcov=matrix(c(0.3,0.0,
                                            0.0,0.06),nrow=2)^2)
    }
    else {
      temp[i,]<-mnormt::rmnorm(n=1,
                            mean=c(0,1),
                            varcov=matrix(c(0.06,0.0,
                                            0.0,0.3),nrow=2)^2)
    }
  }
  return(temp)
}

# first prior samples
first_prior_mvn<-r_mvn(4000)
first_prior_banana<-r_banana(4000)
first_prior_lasso<-r_lasso(4000)
saveRDS(list(first_prior_mvn=first_prior_mvn,
             first_prior_banana=first_prior_banana,
             first_prior_lasso=first_prior_lasso),
        "first_prior_samples.rds")

# initial draw picking function
init_draw<-function(n_chains,prior_draws) {
  init<-list(list(theta=as.numeric(prior_draws[sample(1:NROW(prior_draws),1),])))
  if(n_chains>1) {
    for(i in 2:n_chains) {
      init[[i]]<-list(theta=as.numeric(prior_draws[sample(1:NROW(prior_draws),1),]))
    }
  }
  return(init)
}


# illustrative densities
illus<-function(xfrom,xto,yfrom,yto,densfunc,trim,filename) {
  gridx<-rep(seq(from=xfrom,to=xto,length.out=100),100)
  gridy<-rep(seq(from=yfrom,to=yto,length.out=100),each=100)
  lp<-rep(NA,10000)
  for(i in 1:10000) {
    lp[i]<-densfunc(c(gridx[i],gridy[i]))
    if(lp[i]<trim) { lp[i]<-trim } # to focus contours on maximum region
  }
  grid<-data.frame(gridx,gridy,lp)
  ggplot(grid,aes(gridx,gridy))+
    geom_contour(aes(z=lp),colour='black')+
    theme_classic()+
    labs(x="x",y="y")
  ggsave(filename,width=5,height=5)
}
illus(-6,6,-6,6,logd_mvn,-10,"illustrative_mvn.png")
illus(-4,7,-12,4,logd_banana,-10,"illustrative_banana.png")
illus(-1,2.5,-1,2.5,logd_lasso,-8,"illustrative_lasso.png")

# find bandwidths (2-parameter)
findbw<-function(x,sj.coef=0) {
  # get Sheather-Jones bandwidth
  bandwidth<-apply(x,2,bw.SJ)
  # sj.coef is a fixed coeffcicient for SJ, unless it is zero:
  if(sj.coef>0) {
    return(list(bandwidth=bandwidth*sj.coef,loopcount=NA))
  }
  else if(sj.coef==0) {
   # check for multimodality
    # loop over the parameters
    dipstats<-rep(NA,2)
    loopcount<-1
    dip_p<-rep(0.01,2)
    while(loopcount<100 & min(dip_p)<0.05) {
      # get maximum Hartigans' dip statistic and increase bandwidth until non-sig
      for(i in 1:2) {
        tempdens<-density(x[,i],bw=bandwidth[i])
        tempscale<-4000/round(sum(tempdens$y))
        dip_p[i]<-unlist(diptest::dip.test(rep(tempdens$x,round(tempscale*tempdens$y)))$p.value)
        if(loopcount>1 & dip_p[i]<0.05) { bandwidth[i]<-bandwidth[i]*1.5 }
      }
      loopcount<-loopcount+1
    }
    if(loopcount>99) { stop("Unable to find unimodal univariate kernel bandwidth after 100 expansions") }
    return(list(bandwidth=bandwidth,loopcount=loopcount-1))
  }
  else {
    stop("argument sj.coef in function findbw was negative")
  }
}




######################################################
##################   MVN kernels   ###################
######################################################

# lists for posteriors
posteriors_mvn<-list(first_prior_mvn)
posteriors_banana<-list(first_prior_banana)
posteriors_lasso<-list(first_prior_lasso)

# empty arrays for summaries
summary_mvn<-summary_banana<-summary_lasso<-array(NA,dim=c(2,2,11))

# summaries of first prior samples
summary_mvn[,1,1]<-apply(posteriors_mvn[[1]],2,mean)
summary_mvn[,2,1]<-apply(posteriors_mvn[[1]],2,sd)
summary_banana[,1,1]<-apply(posteriors_banana[[1]],2,mean)
summary_banana[,2,1]<-apply(posteriors_banana[[1]],2,sd)
summary_lasso[,1,1]<-apply(posteriors_lasso[[1]],2,mean)
summary_lasso[,2,1]<-apply(posteriors_lasso[[1]],2,sd)

# lists for bandwidths
bandwidth_mvn<-list(findbw(posteriors_mvn[[1]],sj_coef_option))
bandwidth_banana<-list(findbw(posteriors_banana[[1]],sj_coef_option))
bandwidth_lasso<-list(findbw(posteriors_lasso[[1]],sj_coef_option))

# no-likelihood Stan model for MVN kernels
prior_only_mvn_stan<-
  'data {
          int n; // how many prior draws
          real px[n]; // vector of x locations of prior draws
          real py[n]; // vector of y locations of prior draws
          real bwx; // bandwidth (SD) for x dimension
          real bwy; // bandwidth (SD) for y dimension
  }
  parameters {
          real theta[2];
  }
  model {
          real accuprior=0;
          for(i in 1:n) {
            accuprior += exp(normal_lpdf(theta[1] | px[i], bwx) +
                             normal_lpdf(theta[2] | py[i], bwy));
          }
          target += log(accuprior/(n+0.0));
  }
'
prior_only_mvn_model <- stan_model(model_code=prior_only_mvn_stan)


# MVN distribution, MVN method
for(i in 1:10) {
  # fit Stan model
  prior_only_mvn_fit <- sampling(prior_only_mvn_model,
                                 data=list(n=NROW(posteriors_mvn[[i]]),
                                           px=posteriors_mvn[[i]][,1],
                                           py=posteriors_mvn[[i]][,2],
                                           bwx=bandwidth_mvn[[i]]$bandwidth[1],
                                           bwy=bandwidth_mvn[[i]]$bandwidth[2]),
                                 chains=4,
                                 cores=4,
                                 iter=2500+1000,
                                 warmup=2500,
                                 seed=overall_seed+i,
                                 init=init_draw(4,posteriors_mvn[[i]]))
  summary_mvn[,,i+1]<-summary(prior_only_mvn_fit)$summary[1:2,c(1,3)]
  temp<-extract(prior_only_mvn_fit,permuted=TRUE)
  posteriors_mvn[[i+1]]<-temp$theta
  bandwidth_mvn[[i+1]]<-findbw(posteriors_mvn[[i+1]],sj_coef_option)
}  

# banana distribution, MVN method
for(i in 1:10) {
  # fit Stan model
  prior_only_mvn_fit <- sampling(prior_only_mvn_model,
                                 data=list(n=NROW(posteriors_banana[[i]]),
                                           px=posteriors_banana[[i]][,1],
                                           py=posteriors_banana[[i]][,2],
                                           bwx=bandwidth_banana[[i]]$bandwidth[1],
                                           bwy=bandwidth_banana[[i]]$bandwidth[2]),
                                 chains=4,
                                 cores=4,
                                 iter=2500+1000,
                                 warmup=2500,
                                 seed=overall_seed+20+i,
                                 init=init_draw(4,posteriors_banana[[i]]))
  summary_banana[,,i+1]<-summary(prior_only_mvn_fit)$summary[1:2,c(1,3)]
  temp<-extract(prior_only_mvn_fit,permuted=TRUE)
  posteriors_banana[[i+1]]<-temp$theta
  bandwidth_banana[[i+1]]<-findbw(posteriors_banana[[i+1]],sj_coef_option)
}  

# lasso distribution, MVN method
for(i in 1:10) {
  # fit Stan model
  prior_only_mvn_fit <- sampling(prior_only_mvn_model,
                                 data=list(n=NROW(posteriors_lasso[[i]]),
                                           px=posteriors_lasso[[i]][,1],
                                           py=posteriors_lasso[[i]][,2],
                                           bwx=bandwidth_lasso[[i]]$bandwidth[1],
                                           bwy=bandwidth_lasso[[i]]$bandwidth[2]),
                                 chains=4,
                                 cores=4,
                                 iter=2500+1000,
                                 warmup=2500,
                                 seed=overall_seed+40+i,
                                 init=init_draw(4,posteriors_lasso[[i]]))
  summary_lasso[,,i+1]<-summary(prior_only_mvn_fit)$summary[1:2,c(1,3)]
  temp<-extract(prior_only_mvn_fit,permuted=TRUE)
  posteriors_lasso[[i+1]]<-temp$theta
  bandwidth_lasso[[i+1]]<-findbw(posteriors_lasso[[i+1]],sj_coef_option)
}  

mvn_kernel_outputs<-list(posteriors_mvn=posteriors_mvn,
                         posteriors_banana=posteriors_banana,
                         posteriors_lasso=posteriors_lasso,
                         bandwidth_mvn=bandwidth_mvn,
                         bandwidth_banana=bandwidth_banana,
                         bandwidth_lasso=bandwidth_lasso,
                         summary_mvn=summary_mvn,
                         summary_banana=summary_banana,
                         summary_lasso=summary_lasso,
                         sj_coef_option=sj_coef_option)







######################################################
##############   hypercuboid MVN kernels   ###########
######################################################

# lists for posteriors
posteriors_mvn<-list(first_prior_mvn)
posteriors_banana<-list(first_prior_banana)
posteriors_lasso<-list(first_prior_lasso)

# empty arrays for summaries
summary_mvn<-summary_banana<-summary_lasso<-array(NA,dim=c(2,2,11))

# summaries of first prior samples
summary_mvn[,1,1]<-apply(posteriors_mvn[[1]],2,mean)
summary_mvn[,2,1]<-apply(posteriors_mvn[[1]],2,sd)
summary_banana[,1,1]<-apply(posteriors_banana[[1]],2,mean)
summary_banana[,2,1]<-apply(posteriors_banana[[1]],2,sd)
summary_lasso[,1,1]<-apply(posteriors_lasso[[1]],2,mean)
summary_lasso[,2,1]<-apply(posteriors_lasso[[1]],2,sd)

# lists for bandwidths
bandwidth_mvn<-list(findbw(posteriors_mvn[[1]],sj_coef_option))
bandwidth_banana<-list(findbw(posteriors_banana[[1]],sj_coef_option))
bandwidth_lasso<-list(findbw(posteriors_lasso[[1]],sj_coef_option))

# function to find hypercuboid centres and their counts
gethypercuboids<-function(x,parameter_bins) {
  # divide into parameter bins
  parameter_bin_centres<-cbind(min(x[,1])+
                                 ((max(x[,1])-min(x[,1]))*
                                    ((1:parameter_bins)-0.5)/parameter_bins),
                               min(x[,2])+
                                 ((max(x[,2])-min(x[,2]))*
                                    ((1:parameter_bins)-0.5)/parameter_bins))
  x<-cbind(cut(x[,1],breaks=parameter_bins,labels=FALSE),
           cut(x[,2],breaks=parameter_bins,labels=FALSE))
  p<-ncol(x)
  x<-cbind(parameter_bin_centres[as.numeric(x[,1]),1],
           parameter_bin_centres[as.numeric(x[,2]),2])
  x<-as.data.frame(x)
  colnames(x)<-c("x_bins","y_bins")
  # the binned version
  xbins<-tally(group_by(x,x_bins,y_bins))
  cube_counts<-xbins$n
  xbins<-dplyr::select(xbins,-n)
  return(list(cube_counts=cube_counts,
              prior_bins=as.matrix(xbins),
              prior_draws=as.matrix(x)))
}

# lists of hypercuboids
hypercuboids_mvn<-list(gethypercuboids(posteriors_mvn[[1]],10))
hypercuboids_banana<-list(gethypercuboids(posteriors_banana[[1]],10))
hypercuboids_lasso<-list(gethypercuboids(posteriors_lasso[[1]],10))

prior_only_mvncube_stan <- '
data {
  int n; // number of non-zero-count hypercubes (bins) in the prior
  int<lower=1> count[n]; // number of draws in each hypercube
  matrix[n,2] prior;
  real bandwidth[2];
}
transformed data {
  int<lower=1> sum_count;
  sum_count = sum(count);
}
parameters {
  real theta[2];
}
model {
  real logcubeprob;
  real priorprob;
  priorprob = 0.0;
  for(i in 1:n) {
    logcubeprob = 0.0;
    for(j in 1:2) {
      logcubeprob += normal_lpdf(theta[j] | 
                                prior[i,j], bandwidth[j]); 
    }
    priorprob += count[i]*exp(logcubeprob); 
  }
  priorprob = priorprob / (sum_count+0.0); 
  target += log(priorprob);
}
'
prior_only_mvncube_model <- stan_model(model_code=prior_only_mvncube_stan)

# MVN distribution, hypercuboid MVN method
for(i in 1:10) {
  # fit Stan model
  prior_only_mvncube_fit <- sampling(prior_only_mvncube_model,
                                     data=list(n=NROW(hypercuboids_mvn[[i]]$prior_bins),
                                               prior=hypercuboids_mvn[[i]]$prior_bins,
                                               count=hypercuboids_mvn[[i]]$cube_counts,
                                               bandwidth=bandwidth_mvn[[i]]$bandwidth),
                                     chains=4,
                                     cores=4,
                                     iter=2500+1000,
                                     warmup=2500,
                                     seed=overall_seed+100+i,
                                     init=init_draw(4,hypercuboids_mvn[[i]]$prior_bins))
  summary_mvn[,,i+1]<-summary(prior_only_mvncube_fit)$summary[1:2,c(1,3)]
  temp<-extract(prior_only_mvncube_fit,permuted=TRUE)
  posteriors_mvn[[i+1]]<-temp$theta
  bandwidth_mvn[[i+1]]<-findbw(posteriors_mvn[[i+1]],sj_coef_option)
  hypercuboids_mvn[[i+1]]<-gethypercuboids(posteriors_mvn[[i+1]],10)
}  

# banana distribution, hypercuboid MVN method
for(i in 1:10) {
  # fit Stan model
  prior_only_mvncube_fit <- sampling(prior_only_mvncube_model,
                                     data=list(n=NROW(hypercuboids_banana[[i]]$prior_bins),
                                               prior=hypercuboids_banana[[i]]$prior_bins,
                                               count=hypercuboids_banana[[i]]$cube_counts,
                                               bandwidth=bandwidth_banana[[i]]$bandwidth),
                                     chains=4,
                                     cores=4,
                                     iter=2500+1000,
                                     warmup=2500,
                                     seed=overall_seed+120+i,
                                     init=init_draw(4,hypercuboids_banana[[i]]$prior_bins))
  summary_banana[,,i+1]<-summary(prior_only_mvncube_fit)$summary[1:2,c(1,3)]
  temp<-extract(prior_only_mvncube_fit,permuted=TRUE)
  posteriors_banana[[i+1]]<-temp$theta
  bandwidth_banana[[i+1]]<-findbw(posteriors_banana[[i+1]],sj_coef_option)
  hypercuboids_banana[[i+1]]<-gethypercuboids(posteriors_banana[[i+1]],10)
}  

# lasso distribution, hypercuboid MVN method
for(i in 1:10) {
  # fit Stan model
  prior_only_mvncube_fit <- sampling(prior_only_mvncube_model,
                                     data=list(n=NROW(hypercuboids_lasso[[i]]$prior_bins),
                                               prior=hypercuboids_lasso[[i]]$prior_bins,
                                               count=hypercuboids_lasso[[i]]$cube_counts,
                                               bandwidth=bandwidth_lasso[[i]]$bandwidth),
                                     chains=4,
                                     cores=4,
                                     iter=2500+1000,
                                     warmup=2500,
                                     seed=overall_seed+140+i,
                                     init=init_draw(4,hypercuboids_lasso[[i]]$prior_bins))
  summary_lasso[,,i+1]<-summary(prior_only_mvncube_fit)$summary[1:2,c(1,3)]
  temp<-extract(prior_only_mvncube_fit,permuted=TRUE)
  posteriors_lasso[[i+1]]<-temp$theta
  bandwidth_lasso[[i+1]]<-findbw(posteriors_lasso[[i+1]],sj_coef_option)
  hypercuboids_lasso[[i+1]]<-gethypercuboids(posteriors_lasso[[i+1]],10)
}  

mvncube_kernel_outputs<-list(posteriors_mvn=posteriors_mvn,
                             posteriors_banana=posteriors_banana,
                             posteriors_lasso=posteriors_lasso,
                             bandwidth_mvn=bandwidth_mvn,
                             bandwidth_banana=bandwidth_banana,
                             bandwidth_lasso=bandwidth_lasso,
                             summary_mvn=summary_mvn,
                             summary_banana=summary_banana,
                             summary_lasso=summary_lasso,
                             sj_coef_option=sj_coef_option)



######################################################
###################   beta kernels   #################
######################################################


# lists for posteriors
posteriors_mvn<-list(first_prior_mvn)
posteriors_banana<-list(first_prior_banana)
posteriors_lasso<-list(first_prior_lasso)

# empty arrays for summaries
summary_mvn<-summary_banana<-summary_lasso<-array(NA,dim=c(2,2,11))

# summaries of first prior samples
summary_mvn[,1,1]<-apply(posteriors_mvn[[1]],2,mean)
summary_mvn[,2,1]<-apply(posteriors_mvn[[1]],2,sd)
summary_banana[,1,1]<-apply(posteriors_banana[[1]],2,mean)
summary_banana[,2,1]<-apply(posteriors_banana[[1]],2,sd)
summary_lasso[,1,1]<-apply(posteriors_lasso[[1]],2,mean)
summary_lasso[,2,1]<-apply(posteriors_lasso[[1]],2,sd)

# lists for bandwidths
bandwidth_mvn<-list(findbw(posteriors_mvn[[1]],sj_coef_option))
bandwidth_banana<-list(findbw(posteriors_banana[[1]],sj_coef_option))
bandwidth_lasso<-list(findbw(posteriors_lasso[[1]],sj_coef_option))

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
  int k; // draws in the prior 
  int p; // columns in the prior; this should NOT include __lp
  real<lower=0> betapars; // shared alpha and beta parameters for the beta distribution
  real<lower=0> beta_widths[p];
  matrix[k,p] prior;
}

parameters {
  real theta[2]; 
}

model {
  real logpriorprob;
  real priorprob;
  int n_neighbours = count_neighbours(k, p, theta, prior, beta_widths);
  int neighbour_prior[n_neighbours] = find_neighbours(k, p, theta, prior, beta_widths, n_neighbours);
  
  priorprob = 0.0;
  for(i in 1:n_neighbours) {
    logpriorprob = 0.0;
    for(j in 1:p) {
      logpriorprob += beta_lpdf(0.5+(theta[j]-prior[neighbour_prior[i],j])/(2*beta_widths[j]) | betapars, betapars);
    }
    priorprob += exp(logpriorprob); // add all kernels densities (not log-densities) together
  }
  priorprob = priorprob / (k+0.0); // average over k draws in the prior
  target += log(priorprob); // add log-prior to log-likelihood!
}
'
prior_only_beta_model <- stan_model(model_code=stan_beta_code)


# MVN distribution, beta method
for(i in 1:10) {
  # fit Stan model
  prior_only_beta_fit <- sampling(prior_only_beta_model,
                                  data=list(k=NROW(posteriors_mvn[[i]]),
                                            p=NCOL(posteriors_mvn[[i]]),
                                           prior=posteriors_mvn[[i]],
                                           betapars=8.0,
                                           beta_widths=6*bandwidth_mvn[[i]]$bandwidth),
                                 chains=4,
                                 cores=4,
                                 iter=2500+1000,
                                 warmup=2500,
                                 seed=overall_seed+200+i,
                                 init=init_draw(4,posteriors_mvn[[i]]))
  summary_mvn[,,i+1]<-summary(prior_only_beta_fit)$summary[1:2,c(1,3)]
  temp<-extract(prior_only_beta_fit,permuted=TRUE)
  posteriors_mvn[[i+1]]<-temp$theta
  bandwidth_mvn[[i+1]]<-findbw(posteriors_mvn[[i+1]],sj_coef_option)
}  


# banana distribution, beta method
for(i in 1:10) {
  # fit Stan model
  prior_only_beta_fit <- sampling(prior_only_beta_model,
                                  data=list(k=NROW(posteriors_banana[[i]]),
                                            p=NCOL(posteriors_banana[[i]]),
                                            prior=posteriors_banana[[i]],
                                            betapars=8.0,
                                            beta_widths=6*bandwidth_banana[[i]]$bandwidth),
                                  chains=4,
                                  cores=4,
                                  iter=2500+1000,
                                  warmup=2500,
                                  seed=overall_seed+220+i,
                                  init=init_draw(4,posteriors_banana[[i]]))
  summary_banana[,,i+1]<-summary(prior_only_beta_fit)$summary[1:2,c(1,3)]
  temp<-extract(prior_only_beta_fit,permuted=TRUE)
  posteriors_banana[[i+1]]<-temp$theta
  bandwidth_banana[[i+1]]<-findbw(posteriors_banana[[i+1]],sj_coef_option)
}  


# lasso distribution, beta method
for(i in 1:10) {
  # fit Stan model
  prior_only_beta_fit <- sampling(prior_only_beta_model,
                                  data=list(k=NROW(posteriors_lasso[[i]]),
                                            p=NCOL(posteriors_lasso[[i]]),
                                            prior=posteriors_lasso[[i]],
                                            betapars=8.0,
                                            beta_widths=6*bandwidth_lasso[[i]]$bandwidth),
                                  chains=4,
                                  cores=4,
                                  iter=2500+1000,
                                  warmup=2500,
                                  seed=overall_seed+240+i,
                                  init=init_draw(4,posteriors_lasso[[i]]))
  summary_lasso[,,i+1]<-summary(prior_only_beta_fit)$summary[1:2,c(1,3)]
  temp<-extract(prior_only_beta_fit,permuted=TRUE)
  posteriors_lasso[[i+1]]<-temp$theta
  bandwidth_lasso[[i+1]]<-findbw(posteriors_lasso[[i+1]],sj_coef_option)
}  

beta_kernel_outputs<-list(posteriors_mvn=posteriors_mvn,
                          posteriors_banana=posteriors_banana,
                         posteriors_lasso=posteriors_lasso,
                         bandwidth_mvn=bandwidth_mvn,
                         bandwidth_banana=bandwidth_banana,
                         bandwidth_lasso=bandwidth_lasso,
                         summary_mvn=summary_mvn,
                         summary_banana=summary_banana,
                         summary_lasso=summary_lasso,
                         sj_coef_option=sj_coef_option)





######################################################
#########   sampling importance resampling   #########
######################################################

# lists for posteriors
posteriors_mvn<-list(first_prior_mvn)
posteriors_banana<-list(first_prior_banana)
posteriors_lasso<-list(first_prior_lasso)

# empty arrays for summaries
summary_mvn<-summary_banana<-summary_lasso<-array(NA,dim=c(2,2,11))

# summaries of first prior samples
summary_mvn[,1,1]<-apply(posteriors_mvn[[1]],2,mean)
summary_mvn[,2,1]<-apply(posteriors_mvn[[1]],2,sd)
summary_banana[,1,1]<-apply(posteriors_banana[[1]],2,mean)
summary_banana[,2,1]<-apply(posteriors_banana[[1]],2,sd)
summary_lasso[,1,1]<-apply(posteriors_lasso[[1]],2,mean)
summary_lasso[,2,1]<-apply(posteriors_lasso[[1]],2,sd)

# lists for bandwidths
bandwidth_mvn<-list(findbw(posteriors_mvn[[1]],sj_coef_option))
bandwidth_banana<-list(findbw(posteriors_banana[[1]],sj_coef_option))
bandwidth_lasso<-list(findbw(posteriors_lasso[[1]],sj_coef_option))




# MVN distribution, SIR method
for(i in 1:10) {
  # set seed
  set.seed(overall_seed+300+i)
  # use all prior draws, no sampling
  prior_draws<-posteriors_mvn[[i]]

  # find bandwidths
  bandwidth_mvn[[i]]<-findbw(prior_draws,sj_coef_option)
  
  # constant log-likelihood means one posterior draw per prior draw
  posteriors_mvn[[i+1]]<-mnormt::rmnorm(1,
                                        as.matrix(prior_draws[1,]),
                                        diag((bandwidth_mvn[[i]]$bandwidth)^2))
  for(k in 2:NROW(prior_draws)) {
    posteriors_mvn[[i+1]]<-rbind(posteriors_mvn[[i+1]],
                                 mnormt::rmnorm(1,
                                                as.matrix(prior_draws[k,]),
                                                diag((bandwidth_mvn[[i]]$bandwidth)^2)))
  }
  # summary stats
  summary_mvn[,1,i+1]<-apply(posteriors_mvn[[i+1]],2,mean)
  summary_mvn[,2,i+1]<-apply(posteriors_mvn[[i+1]],2,sd)
}




# banana distribution, SIR method
for(i in 1:10) {
  # set seed
  set.seed(overall_seed+320+i)
  # use all prior draws, no sampling
  prior_draws<-posteriors_banana[[i]]
  
  # find bandwidths
  bandwidth_banana[[i]]<-findbw(prior_draws,sj_coef_option)
  
  # constant log-likelihood means one posterior draw per prior draw
  posteriors_banana[[i+1]]<-mnormt::rmnorm(1,
                                           as.matrix(prior_draws[1,]),
                                           diag((bandwidth_banana[[i]]$bandwidth)^2))
  for(k in 2:NROW(prior_draws)) {
    posteriors_banana[[i+1]]<-rbind(posteriors_banana[[i+1]],
                                    mnormt::rmnorm(1,
                                                   as.matrix(prior_draws[k,]),
                                                   diag((bandwidth_banana[[i]]$bandwidth)^2)))
  }
  # summary stats
  summary_banana[,1,i+1]<-apply(posteriors_banana[[i+1]],2,mean)
  summary_banana[,2,i+1]<-apply(posteriors_banana[[i+1]],2,sd)
}




# lasso distribution, SIR method
for(i in 1:10) {
  # set seed
  set.seed(overall_seed+340+i)
  # use all prior draws, no sampling
  prior_draws<-posteriors_lasso[[i]]
  
  # find bandwidths
  bandwidth_lasso[[i]]<-findbw(prior_draws,sj_coef_option)
  
  # constant log-likelihood means one posterior draw per prior draw
  posteriors_lasso[[i+1]]<-mnormt::rmnorm(1,
                                          as.matrix(prior_draws[1,]),
                                          diag((bandwidth_lasso[[i]]$bandwidth)^2))
  for(k in 2:NROW(prior_draws)) {
    posteriors_lasso[[i+1]]<-rbind(posteriors_lasso[[i+1]],
                                   mnormt::rmnorm(1,
                                                  as.matrix(prior_draws[k,]),
                                                  diag((bandwidth_lasso[[i]]$bandwidth)^2)))
  }
  # summary stats
  summary_lasso[,1,i+1]<-apply(posteriors_lasso[[i+1]],2,mean)
  summary_lasso[,2,i+1]<-apply(posteriors_lasso[[i+1]],2,sd)
}


sir_kernel_outputs<-list(posteriors_mvn=posteriors_mvn,
                          posteriors_banana=posteriors_banana,
                          posteriors_lasso=posteriors_lasso,
                          bandwidth_mvn=bandwidth_mvn,
                          bandwidth_banana=bandwidth_banana,
                          bandwidth_lasso=bandwidth_lasso,
                          summary_mvn=summary_mvn,
                          summary_banana=summary_banana,
                          summary_lasso=summary_lasso,
                          sj_coef_option=sj_coef_option)


######################################################
##############   graphics and stats   ################
######################################################

# small multiple scatterplots of posteriors

# convex hull function (not used)
gethull<-function(x,n=1) {
  hull<-chull(x)
  if(n>1) {
    for(i in 1:(n-1)) {
      x<-x[-hull,]
      hull<-chull(x)
    }
  }
  return(x[c(hull,hull[1]),])
}

# hard-coding ahoy!
x_lims<-list(c(-5,5),
             c(-10,12),
             c(-2,3))
y_lims<-list(c(-5,5),
             c(-16,6),
             c(-2,3))

# loop over methods
for(meth in c("mvn","mvncube","beta","sir")) {
# loop over distributions
  dist_names<-c("mvn","banana","lasso")
  for(dist in 1:3) {
    temp<-get(paste0(meth,"_kernel_outputs"))
    png(paste0("posteriors_method-",meth,"_dist-",dist_names[dist],".png"),
        width=1000,height=3000)
    par(mfrow=c(6,2))
    #hull<-gethull(temp[[dist]][[1]]) # hull of initial prior draws
    for(update in 1:11) {
      plot(temp[[dist]][[update]],
           cex=1,col='#00000040',
           xaxt='n',yaxt='n',
           xlab='',ylab='',
           xlim=x_lims[[dist]],ylim=y_lims[[dist]],
           asp=1,
           main=paste0("Update ",update-1),
           cex.main=4)  
      #lines(hull)
    }
    dev.off()
  }
}

######################################################
##################   save files   ####################
######################################################

saveRDS(list(mvn_kernel_outputs=mvn_kernel_outputs,
             mvncube_kernel_outputs=mvncube_kernel_outputs,
             beta_kernel_outputs=beta_kernel_outputs,
             sir_kernel_outputs=sir_kernel_outputs),
        "prior-only-outputs.rds")