# Calculation of sliced 2-Wasserstein distances for
# samples from prior-only examples

# change this directory for different bandwidth calculations
setwd('~/Documents/Shop-window-NOC/non-parametric-updating/2020-sims/prior-only/double-SJ-bandwidth')
n_slices<-20 # number of random rotations and slices
set.seed(8155)


# load the results of the prior-only updates
results<-readRDS("prior-only-outputs.rds")

# empty array to store 2-Ws:
wassersteins<-array(NA,dim=c(4,3,10))

pb<-txtProgressBar(min=0,max=120,style=3)
# loop over 4 methods (in order mvn, mvncube, beta, sir)
for(meth in 1:4) {
  # loop over 3 distributions (in order mvn, banana, lasso)
  for(dist in 1:3) {
    # get original prior sample
    prior_sample<-results[[meth]][[dist]][[1]]
    # loop over 10 updates
    for(update in 1:10) {
      setTxtProgressBar(pb,update + (dist-1)*10 + (meth-1)*30)
      # get posterior sample
      posterior_sample<-results[[meth]][[dist]][[update+1]]
      wasserstein_slices<-rep(NA,n_slices)
      # loop over n_slices slices
      for(slice in 1:n_slices) {
        # generate random rotation matrix
        rotangle<-runif(1,0,2*pi)
        rotmat<-matrix(c(cos(rotangle),
                         sin(rotangle),
                         -sin(rotangle),
                         cos(rotangle)),
                       nrow=2)
        # get 1-D marginal distributions for this update's posterior and
        #   the original prior sample
        prior_marginal<-sort(t(rotmat%*%t(prior_sample))[,1])
        posterior_marginal<-sort(t(rotmat%*%t(posterior_sample))[,1])
        # estimate 2-Wasserstein
        wasserstein_slices[slice]<-sqrt(sum((posterior_marginal-prior_marginal)^2))
      }
      # average over the slices
      wassersteins[meth,dist,update]<-mean(wasserstein_slices)
    }
  }
}
close(pb)
saveRDS(wassersteins,"wassersteins.rds")
# In another script, we will compile 2-Wassersteins across bandwidths, and
# plot them against updates
