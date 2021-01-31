# simulation study stats

setwd("~/Documents/Shop-window-NOC/non-parametric-updating/2020-sims")
load("combined-results-cloud-sims-4000-draws.RData")
results_sir<-readRDS("combined-results-sir.rds")

# bias and mean posterior SD (latter not used)
apply(summary_mvn[,,10,],c(1,2),mean)
apply(summary_mvncube[,,10,],c(1,2),mean)
apply(summary_beta[,,10,],c(1,2),mean)
apply(results_sir$summary_sir[,,10,],c(1,2),mean)

# coverage
deleteme<-(summary_mvn[,1,10,]-matrix(rep(c(2,0.5,-1.5),400),nrow=3)) / 
            summary_mvn[,2,10,]
apply(deleteme,1,function(x){sum(abs(x)<1.96)})/400

deleteme<-(summary_mvncube[,1,10,]-matrix(rep(c(2,0.5,-1.5),400),nrow=3)) / 
  summary_mvncube[,2,10,]
apply(deleteme,1,function(x){sum(abs(x)<1.96)})/400

deleteme<-(summary_beta[,1,10,]-matrix(rep(c(2,0.5,-1.5),400),nrow=3)) / 
  summary_beta[,2,10,]
apply(deleteme,1,function(x){sum(abs(x)<1.96)})/400

deleteme<-(results_sir$summary_sir[,1,10,]-
             matrix(rep(c(2,0.5,-1.5),400),nrow=3)) / 
  results_sir$summary_sir[,2,10,]
apply(deleteme,1,function(x){sum(abs(x)<1.96)})/400


# time taken
mean(time_mvn)
mean(time_mvncube)
mean(time_beta)
mean(results_sir$time_sir)

# proportion of occupied hypercubes
apply(occupied_hypercubes[1:400,],2,mean)/400


# mean means and mean CIs
makemeans<-function(x) {
  mm<-array(NA,dim=c(3,3,10))
  mm[,1,]<-apply(x[,1,,],c(1,2),mean)
  mm[,2,]<-mm[,1,] - 1.96*apply(x[,2,,],c(1,2),mean)
  mm[,3,]<-mm[,1,] + 1.96*apply(x[,2,,],c(1,2),mean)
  return(mm)
}
mm_mvn<-makemeans(summary_mvn)
mm_mvncube<-makemeans(summary_mvncube)
mm_beta<-makemeans(summary_beta)
mm_sir<-makemeans(results_sir$summary_sir)
plotmeans<-function(x,title) {
  # plot means first
  plot(1:10,x[1,1,],
       type='l',
       lwd=1.8,
       ylim=c(min(x),max(x)),
       xlim=c(1,10.8),
       col='#b62525',
       ylab="Parameter values",
       xlab="Update number",
       main=title)
  lines(1:10,x[2,1,],
        lwd=1.8,
        col='#25b66e')
  lines(1:10,x[3,1,],
        lwd=1.8,
        col='#266eb6')
  # lower CIs
  lines(1:10,x[1,2,],
        lwd=1,
        col='#b62525')
  lines(1:10,x[2,2,],
        lwd=1,
        col='#25b66e')
  lines(1:10,x[3,2,],
        lwd=1,
        col='#266eb6')
  # upper CIs
  lines(1:10,x[1,3,],
        lwd=1,
        col='#b62525')
  lines(1:10,x[2,3,],
        lwd=1,
        col='#25b66e')
  lines(1:10,x[3,3,],
        lwd=1,
        col='#266eb6')
  # true values
  lines(c(1,10.7),rep(2,2),
        lty=3,col='#b62525')
  lines(c(1,10.7),rep(0.5,2),
        lty=3,col='#25b66e')
  lines(c(1,10.7),rep(-1.5,2),
        lty=3,col='#266eb6')
}
plotmeans(mm_mvn,"Algorithm 2 (MVN)")
plotmeans(mm_mvncube,"Algorithm 3 (Hypercuboid)")
plotmeans(mm_beta,"Algorithm 4 (beta)")
plotmeans(mm_sir,"Algorithm 5 (SIR)")



# plot of all posterior means
plotall<-function(x,title) {
  plot(1:10,x[1,1,,1],type='l',
       ylim=c(min(x[1,1,,]),
              max(x[1,1,,])),
       col='#00000020',
       ylab="Beta_0 values",
       xlab="Update number",
       main=title)
  for(i in 2:400) {
    lines(1:10,x[1,1,,i],
          col='#00000020')
  }
}
plotall(summary_mvn,"Algorithm 2 (MVN)")
plotall(summary_mvncube,"Algorithm 3 (Hypercuboid)")
plotall(summary_beta,"Algorithm 4 (beta)")
plotall(results_sir$summary_sir,"Algorithm 5 (SIR)")
