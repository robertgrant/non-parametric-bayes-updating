# obtain distribution of log10 20-sliced 2-Wasserstein distances
# arising from sampling error in our 3 distributions

setwd('~/Documents/Shop-window-NOC/non-parametric-updating/2020-sims/prior-only')
set.seed(40151)
n_slices<-20 # number of "slices" 
n_sims<-1000 # number of simulations


# same function code as in prior-only.R:
r_mvn <- function(n) {
  mnormt::rmnorm(n=n,
                 mean=rep(0,2),
                 varcov=matrix(c(1.0,0.6,
                                 0.6,1.0),nrow=2)^2)
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

getwass<-function(x0,x,n_slices) {
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
    # get 1-D marginal distributions 
    x0_marginal<-sort(t(rotmat%*%t(x0))[,1])
    x_marginal<-sort(t(rotmat%*%t(x))[,1])
    # estimate 2-Wasserstein
    wasserstein_slices[slice]<-sqrt(sum((x0_marginal-x_marginal)^2))
  }
  # average over the slices
  return(mean(wasserstein_slices))
}

# empty vectors
wasserstein_mvn<-wasserstein_banana<-wasserstein_lasso<-rep(NA,n_sims)

# "original" samples to compare others to:
x0_mvn<-r_mvn(4000)
x0_banana<-r_banana(4000)
x0_lasso<-r_lasso(4000)

# get Wassers
pb<-txtProgressBar(min=0,max=n_sims,style=3)
for(j in 1:n_sims) {
  setTxtProgressBar(pb,j)
  x_mvn<-r_mvn(4000)
  wasserstein_mvn[j]<-getwass(x0_mvn,x_mvn,n_slices)
  x_banana<-r_banana(4000)
  wasserstein_banana[j]<-getwass(x0_banana,x_banana,n_slices)
  x_lasso<-r_lasso(4000)
  wasserstein_lasso[j]<-getwass(x0_lasso,x_lasso,n_slices)
}
close(pb)

png("wasserstein_mvn_sampling.png")
hist(log10(wasserstein_mvn),breaks=50)
dev.off()
png("wasserstein_banana_sampling.png")
hist(log10(wasserstein_banana),breaks=50)
dev.off()
png("wasserstein_lasso_sampling.png")
hist(log10(wasserstein_lasso),breaks=50)
dev.off()

quantile(log10(wasserstein_mvn),probs=c(0.025,0.25,0.5,0.75,0.975))
quantile(log10(wasserstein_banana),probs=c(0.025,0.25,0.5,0.75,0.975))
quantile(log10(wasserstein_lasso),probs=c(0.025,0.25,0.5,0.75,0.975))


