# recombine simulations for non-parametric kernel Bayesian updating

library(lubridate)
setwd('~/Dropbox/Shop window NOC/non-parametric-updating')

n_instances<-4
sims_per_instance<-260
total_sims<-n_instances*sims_per_instance

mydims<-list(parameters=c('beta0','beta1','beta2'),
             statistics=c('estimate','SE/SD'),
             simulations=as.character(1:total_sims))
mle<-stan_all<-stan_stats<-stan_kernel<-stan_cube<-array(NA,dim=c(3,2,total_sims),dimnames=mydims)
time_stats<-time_kernel<-time_cube<-rep(NA,total_sims)
pb<-txtProgressBar(min=0,max=total_sims,style=3)
for(inst in 1:n_instances) {
  for(sim in 1:sims_per_instance) {
    i<-((inst-1)*sims_per_instance)+sim
    setTxtProgressBar(pb,i)
    temp<-readRDS(paste0('instance_',inst,
                         '/sim_results_',sim,
                         '.rds'))
    if(class(temp)!="try-error") {
      mle[,,i]<-temp$mles[,1:2]
      stan_all[,,i]<-temp$stan_summary_all[1:3,c(1,3)]
      stan_stats[,,i]<-temp$stan_summary_stats[,,10]
      stan_kernel[,,i]<-temp$stan_summary_kernel[,,10]
      stan_cube[,,i]<-temp$stan_summary_kernel_cube[,,10]
      time_stats[i]<-dminutes(as_datetime(temp$stan_normal_stats_time[1])%--%
                                as_datetime(temp$stan_normal_stats_time[2]))  
      time_kernel[i]<-dminutes(as_datetime(temp$stan_kernel_time[1])%--%
                                as_datetime(temp$stan_kernel_time[2]))  
      time_cube[i]<-dminutes(as_datetime(temp$stan_kernel_cube_time[1])%--%
                                as_datetime(temp$stan_kernel_cube_time[2]))  
    }
  }
}
close(pb)

# remove faulty sim(s) and reduce to 1000
selecta<-c(1:487,489:1001)
mle<-mle[,,selecta]
stan_all<-stan_all[,,selecta]
stan_stats<-stan_stats[,,selecta]
stan_kernel<-stan_kernel[,,selecta]
stan_cube<-stan_cube[,,selecta]
time_stats<-time_stats[selecta]
time_kernel<-time_stats[selecta]
time_cube<-time_stats[selecta]

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

plot(mle[1,1,],mle[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0,8),ylim=c(-2,5))
plot(stan_all[1,1,],stan_all[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0,8),ylim=c(-2,5))
plot(stan_stats[1,1,],stan_stats[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0,8),ylim=c(-2,5))
plot(stan_kernel[1,1,],stan_kernel[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0,8),ylim=c(-2,5))
plot(stan_cube[1,1,],stan_cube[2,1,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0,8),ylim=c(-2,5))


plot(mle[1,2,],mle[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2))
plot(stan_all[1,2,],stan_all[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2))
plot(stan_stats[1,2,],stan_stats[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2))
plot(stan_kernel[1,2,],stan_kernel[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2))
plot(stan_cube[1,2,],stan_cube[2,2,],cex=0.3,col='#00000010',pch=19,
     xlim=c(0.1,1.2),ylim=c(0.1,1.2))
