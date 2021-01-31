# recombine simulations for non-parametric kernel Bayesian updating

library(lubridate)
setwd('~/Documents/Shop-window-NOC/non-parametric-updating/2020-sims')

# -------------- this block to be changed for each import ----------------
subfolder<-"cloud-sims-4000-draws"
k_draws<-4000
instances<-(1:4)
max_sims_per_instance<-100
row_craft_chart<-10
col_craft_chart<-3
sir_folder<-"importance-sampling-1k-from-4k"
# ------------------------------------------------------------------------

n_instances<-max(instances)
total_sims<-n_instances*max_sims_per_instance
n_craft_chart<-row_craft_chart*col_craft_chart

# function returning average count of beta neighbours (L-inf norm)
countneighbours<-function(draws,widths,n_sample) {
  drawsample<-draws[sample(1:NROW(draws),n_sample,replace=FALSE),]
  samplecounts<-rep(NA,n_sample)
  for(i in 1:n_sample) {
    univar_neighbours<-abs(draws[,1]-drawsample[i,1])<widths[1]
    if(NCOL(draws)>2) {
      for(j in 2:NCOL(draws)) {
        univar_neighbours<-cbind(univar_neighbours,
                                 abs(draws[,j]-drawsample[i,j])<widths[j])
      }
    }
    samplecounts[i]<-sum(apply(univar_neighbours,1,all))
  }
  return(mean(samplecounts))
}

# empty objects
mydims<-list(parameters=c('beta0','beta1','beta2'),
             statistics=c('estimate','SE/SD'),
             batches=as.character(1:10),
             simulations=as.character(1:total_sims))
mle<-summary_all<-summary_para<-summary_mvn<-summary_mvncube<-summary_beta<-
  array(NA,dim=c(3,2,10,total_sims),dimnames=mydims)
time_all<-time_para<-time_mvn<-time_mvncube<-time_beta<-rep(NA,total_sims)
posterior_all<-posterior_para<-posterior_mvn<-posterior_mvncube<-posterior_beta<-NULL
bandwidth_mvn<-bandwidth_mvncube<-bandwidth_beta<-NULL
occupied_hypercubes<-beta_neighbours<-matrix(NA,k_draws,10)

pb<-txtProgressBar(min=0,max=total_sims,style=3)


# loop over instances
for(i in instances) {
  # loop over simulations
  for(j in 1:max_sims_per_instance) {
    simcount<-(i-1)*max_sims_per_instance + j
    setTxtProgressBar(pb,simcount)
    currentsim<-readRDS(paste0(subfolder,"/instance",i,"/sim_results_",j,".rds"))
    if(simcount<=n_craft_chart) {
      posterior_all[[simcount]]<-cbind(currentsim$stan_summary_all$posterior$beta0,
                                       currentsim$stan_summary_all$posterior$beta1,
                                       currentsim$stan_summary_all$posterior$beta2)
      posterior_para[[simcount]]<-
        posterior_mvn[[simcount]]<-
        posterior_mvncube[[simcount]]<-
        posterior_beta[[simcount]]<-
        bandwidth_mvn[[simcount]]<-
        bandwidth_mvncube[[simcount]]<-
        bandwidth_beta[[simcount]]<-list(42) # temporary placeholder
      for(k in 1:10) {
        if(k==1) {
          posterior_para[[simcount]][[k]]<-cbind(currentsim$log_normal[[k]]$posterior$beta0,
                                                 currentsim$log_normal[[k]]$posterior$beta1,
                                                 currentsim$log_normal[[k]]$posterior$beta2)
        }
        else {
          posterior_para[[simcount]][[k]]<-currentsim$log_normal[[k]]$posterior$beta
          
        }
        posterior_mvn[[simcount]][[k]]<-currentsim$log_kernel[[k]]$next_prior_draws
        posterior_mvncube[[simcount]][[k]]<-as.matrix(
                              currentsim$log_hypercube[[k]]$hypercube_vars$next_prior_draws)
        posterior_beta[[simcount]][[k]]<-as.matrix(currentsim$log_beta[[k]]$next_prior_draws)
        dimnames(posterior_para[[simcount]][[k]])<-NULL
        dimnames(posterior_mvn[[simcount]][[k]])<-NULL
        dimnames(posterior_mvncube[[simcount]][[k]])<-NULL
        dimnames(posterior_beta[[simcount]][[k]])<-NULL
        
        bandwidth_mvn[[simcount]][[k]]<-currentsim$log_kernel[[k]]$bandwidth
        bandwidth_mvncube[[simcount]][[k]]<-currentsim$log_hypercube[[k]]$hypercube_vars$bandwidth
        bandwidth_beta[[simcount]][[k]]<-currentsim$log_beta[[k]]$bandwidth
        names(bandwidth_mvn[[simcount]][[k]])<-NULL
        names(bandwidth_mvncube[[simcount]][[k]])<-NULL
        names(bandwidth_beta[[simcount]][[k]])<-NULL
      }
    }
    mle[,,10,simcount]<-currentsim$mles[,1:2]
    summary_all[,,10,simcount]<-currentsim$stan_summary_all$summary[1:3,c(1,3)]
    summary_mvn[,,,simcount]<-currentsim$stan_summary_kernel
    summary_mvncube[,,,simcount]<-as.matrix(currentsim$stan_summary_kernel_cube)
    summary_beta[,,,simcount]<-currentsim$stan_summary_beta
    
    # this part takes about 6 hours
    if(subfolder=="cloud-sims-4000-draws") {
      for(k in 1:10) {
        occupied_hypercubes[simcount,k]<-length(currentsim$log_hypercube[[k]]$hypercube_vars$cube_counts)
        beta_neighbours[simcount,k]<-countneighbours(draws=currentsim$log_beta[[k]]$next_prior_draws,
                                                     widths=currentsim$log_beta[[k]]$beta_widths,
                                                     n_sample=1000)
      }
    }
    
    time_para[simcount]<-dminutes(as_datetime(currentsim$stan_normal_stats_time[[1]])%--%
                                  as_datetime(currentsim$stan_normal_stats_time[[2]]))/60  
    time_mvn[simcount]<-dminutes(as_datetime(currentsim$stan_kernel_time[[1]])%--%
                                 as_datetime(currentsim$stan_kernel_time[[2]]))/60
    time_mvncube[simcount]<-dminutes(as_datetime(currentsim$stan_kernel_cube_time[[1]])%--%
                                     as_datetime(currentsim$stan_kernel_cube_time[[2]]))/60
    time_beta[simcount]<-dminutes(as_datetime(currentsim$stan_beta_time[[1]])%--%
                                  as_datetime(currentsim$stan_beta_time[[2]]))/60
  }
}
close(pb)
#rm(currentsim)

save.image(paste0('combined-results-',subfolder,'.RData'))



#########################################################
#################   get data batches   ##################
#########################################################
pb<-txtProgressBar(min=0,max=total_sims,style=3)
# loop over instances
for(i in instances) {
  # loop over simulations
  for(j in 1:max_sims_per_instance) {
    simcount<-(i-1)*max_sims_per_instance + j
    setTxtProgressBar(pb,simcount)
    currentsim<-readRDS(paste0(subfolder,"/instance",i,"/sim_results_",j,".rds"))
    db<-currentsim$data_batches
    saveRDS(db,paste0('data_batches/data_batches_',simcount,'.rds'))
  }
}
close(pb)
rm(db)
rm(currentsim)

###############################################################
#################   get batch 1 posteriors   ##################
###############################################################
pb<-txtProgressBar(min=0,max=total_sims,style=3)
# loop over instances
for(i in instances) {
  # loop over simulations
  for(j in 1:max_sims_per_instance) {
    simcount<-(i-1)*max_sims_per_instance + j
    setTxtProgressBar(pb,simcount)
    currentsim<-readRDS(paste0(subfolder,"/instance",i,"/sim_results_",j,".rds"))
    fp<-currentsim$log_kernel[[1]]$next_prior_draws
    saveRDS(fp,paste0('first_posteriors/first_posterior_',simcount,'.rds'))
  }
}
close(pb)
rm(fp)
rm(currentsim)

########################################################
#################   get SIR results   ##################
########################################################

# empty objects
posterior_sir<-bandwidth_sir<-NULL
summary_sir<-array(NA,dim=c(3,2,10,total_sims),dimnames=mydims)
time_sir<-rep(NA,total_sims)
# loop over sims
pb<-txtProgressBar(min=0,max=total_sims,style=3)
for(i in 1:total_sims) {
  setTxtProgressBar(pb,i)
  current_results<-readRDS(paste0(sir_folder,"/results_impsamp_",i,".rds"))
  bandwidth_sir[[i]]<-list(NA)
  if(i<=n_craft_chart) {
    posterior_sir[[i]]<-list(NA)
    for(k in 1:10) {
      posterior_sir[[i]][[k]]<-current_results[[k]]$posterior
    }
  }
  for(k in 1:10) {
    summary_sir[,1,k,i]<-current_results[[k]]$posterior_mean
    summary_sir[,2,k,i]<-current_results[[k]]$posterior_sd
  }
  for(k in 2:10) {
    bandwidth_sir[[i]][[k]]<-current_results[[k]]$bandwidth
  }
  current_time<-readRDS(paste0(sir_folder,"/time_impsamp_",i,".rds"))
  time_sir[i]<-dminutes(as_datetime(current_time[[1]])%--%
                        as_datetime(current_time[[2]]))/60  
  deleteme<-i
}
close(pb)
rm(current_results)
combined_results_sir<-list(posterior_sir=posterior_sir,
                           summary_sir=summary_sir,
                           bandwidth_sir=bandwidth_sir,
                           time_sir=time_sir)
saveRDS(combined_results_sir,'combined-results-sir.rds')










plotest<-function(name,rownum,colnum,xmin,xmax) {
  png(paste0("sim_",name,".png"))
  par(mfrow=c(3,2))
  plot(1,1,bty='n',col='white',xaxt='n',yaxt='n',xlab='',ylab='')
  text(1,1,labels=name,cex=1.5)
  hist(mle[rownum,colnum,],breaks=seq(from=xmin,to=xmax,length.out=100),
       xlim=c(xmin,xmax),main="",ylab="n of simulations",xlab="MLE (IRLS)")
  hist(stan_all[rownum,colnum,],breaks=seq(from=xmin,to=xmax,length.out=100),
       xlim=c(xmin,xmax),main="",ylab="n of simulations",xlab="Stan: all data")
  hist(stan_stats[rownum,colnum,],breaks=seq(from=xmin,to=xmax,length.out=100),
       xlim=c(xmin,xmax),main="",ylab="n of simulations",xlab="Stan: normal stats")
  hist(stan_kernel[rownum,colnum,],breaks=seq(from=xmin,to=xmax,length.out=100),
       xlim=c(xmin,xmax),main="",ylab="n of simulations",xlab="Stan: kernels")
  hist(stan_cube[rownum,colnum,],breaks=seq(from=xmin,to=xmax,length.out=100),
       xlim=c(xmin,xmax),main="",ylab="n of simulations",xlab="Stan: hypercubes")
  dev.off()
}
plotest(name="beta0_estimate",rownum=1,colnum=1,xmin=(-0.5),xmax=10)
plotest(name="beta0_SE",rownum=1,colnum=2,xmin=0,xmax=3)

plotest(name="beta1_estimate",rownum=2,colnum=1,xmin=(-2),xmax=6)
plotest(name="beta1_SE",rownum=2,colnum=2,xmin=0,xmax=1.6)

plotest(name="beta2_estimate",rownum=3,colnum=1,xmin=(-9),xmax=3)
plotest(name="beta2_SE",rownum=3,colnum=2,xmin=0,xmax=3)

# convex hull function
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

# scatterplots of beta0 vs beta1
hull<-gethull(t(stan_all[1:2,1,]),3)
png("sim_beta0_v_beta1.png",width=600,height=900)
par(mfrow=c(3,2))
plot(1,1,bty='n',col='white',xaxt='n',yaxt='n',xlab='',ylab='')
text(1,1,labels="Posterior means",cex=1.5)
plot(mle[1,1,],mle[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0,8),ylim=c(-2,5),xlab="beta0",ylab="beta1",
     main="Method: IRLS")
  lines(hull,col='#579a04a0')
plot(stan_all[1,1,],stan_all[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0,8),ylim=c(-2,5),xlab="beta0",ylab="beta1",
     main="Method: all data")
  lines(hull,col='#579a04a0')
plot(stan_stats[1,1,],stan_stats[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0,8),ylim=c(-2,5),xlab="beta0",ylab="beta1",
     main="Method: parametric")
  lines(hull,col='#579a04a0')
plot(stan_kernel[1,1,],stan_kernel[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0,8),ylim=c(-2,5),xlab="beta0",ylab="beta1",
     main="Method: kernels")
lines(hull,col='#579a04a0')
plot(stan_cube[1,1,],stan_cube[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0,8),ylim=c(-2,5),xlab="beta0",ylab="beta1",
     main="Method: hypercubes")
lines(hull,col='#579a04a0')
dev.off()
par(mfrow=c(1,1))

# beta1 vs beta2
hull<-gethull(t(stan_all[c(3,2),1,]),3)
png("sim_beta2_v_beta1.png",width=600,height=900)
par(mfrow=c(3,2))
plot(1,1,bty='n',col='white',xaxt='n',yaxt='n',xlab='',ylab='')
text(1,1,labels="Posterior means",cex=1.5)
plot(mle[3,1,],mle[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(-8,2),ylim=c(-2,5),xlab="beta2",ylab="beta1",
     main="Method: IRLS")
  lines(hull,col='#579a04a0')
plot(stan_all[3,1,],stan_all[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(-8,2),ylim=c(-2,5),xlab="beta2",ylab="beta1",
     main="Method: all data")
lines(hull,col='#579a04a0')
plot(stan_stats[3,1,],stan_stats[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(-8,2),ylim=c(-2,5),xlab="beta2",ylab="beta1",
     main="Method: parametric")
lines(hull,col='#579a04a0')
plot(stan_kernel[3,1,],stan_kernel[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(-8,2),ylim=c(-2,5),xlab="beta2",ylab="beta1",
     main="Method: kernels")
lines(hull,col='#579a04a0')
plot(stan_cube[3,1,],stan_cube[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(-8,2),ylim=c(-2,5),xlab="beta2",ylab="beta1",
     main="Method: hypercubes")
lines(hull,col='#579a04a0')
dev.off()
par(mfrow=c(1,1))

# SEs: beta0 vs beta1
hull<-gethull(t(stan_all[c(1,2),2,]),3)
png("sim_beta0_v_beta1_SEs.png",width=600,height=900)
par(mfrow=c(3,2))
plot(1,1,bty='n',col='white',xaxt='n',yaxt='n',xlab='',ylab='')
text(1,1,labels="Posterior SDs",cex=1.5)
plot(mle[1,2,],mle[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2),xlab="beta0",ylab="beta1",
     main="Method: IRLS")
lines(hull,col='#579a04a0')
plot(stan_all[1,2,],stan_all[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2),xlab="beta0",ylab="beta1",
     main="Method: all data")
lines(hull,col='#579a04a0')
plot(stan_stats[1,2,],stan_stats[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2),xlab="beta0",ylab="beta1",
     main="Method: parametric")
lines(hull,col='#579a04a0')
plot(stan_kernel[1,2,],stan_kernel[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2),xlab="beta0",ylab="beta1",
     main="Method: kernels")
lines(hull,col='#579a04a0')
plot(stan_cube[1,2,],stan_cube[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2),xlab="beta0",ylab="beta1",
     main="Method: hypercubes")
lines(hull,col='#579a04a0')
dev.off()
par(mfrow=c(1,1))

# SEs: beta2 vs beta1
hull<-gethull(t(stan_all[c(3,2),2,]),3)
png("sim_beta2_v_beta1_SEs.png",width=600,height=900)
par(mfrow=c(3,2))
plot(1,1,bty='n',col='white',xaxt='n',yaxt='n',xlab='',ylab='')
text(1,1,labels="Posterior SDs",cex=1.5)
plot(mle[3,2,],mle[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2),xlab="beta2",ylab="beta1",
     main="Method: IRLS")
lines(hull,col='#579a04a0')
plot(stan_all[3,2,],stan_all[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2),xlab="beta2",ylab="beta1",
     main="Method: all data")
lines(hull,col='#579a04a0')
plot(stan_stats[3,2,],stan_stats[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2),xlab="beta2",ylab="beta1",
     main="Method: parametric")
lines(hull,col='#579a04a0')
plot(stan_kernel[3,2,],stan_kernel[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2),xlab="beta2",ylab="beta1",
     main="Method: kernels")
lines(hull,col='#579a04a0')
plot(stan_cube[3,2,],stan_cube[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2),xlab="beta2",ylab="beta1",
     main="Method: hypercubes")
lines(hull,col='#579a04a0')
dev.off()
par(mfrow=c(1,1))


# investigate spike in slope SEs
spike<-which(stan_kernel[2,2,]>0.18 & stan_kernel[2,2,]<0.23 &
               stan_kernel[3,2,]>0.18 & stan_kernel[3,2,]<0.23)
length(spike) # how many are in there?
plot(stan_kernel[3,2,spike],stan_kernel[2,2,spike]) # zoom in
# find one
spike[1]
stan_kernel[,,8]