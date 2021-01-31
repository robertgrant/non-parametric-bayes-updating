# Wasserstein plot across bandwidths
setwd('~/Documents/Shop-window-NOC/non-parametric-updating/2020-sims/prior-only')

folders<-c("quarter-SJ-bandwidth",
           "half-SJ-bandwidth",
           "SJ-bandwidth",
           "double-SJ-bandwidth",
           "quadruple-SJ-bandwidth")

# collect stats
wassers<-list(readRDS("wassersteins.rds"))
fldrcount<-2
for(fldr in folders) {
  wassers[[fldrcount]]<-readRDS(paste0(fldr,"/wassersteins.rds"))
  fldrcount<-fldrcount+1
}
names(wassers)<-c("multimodality-expansion",
                  "quarter-Sheather-Jones",
                  "half-Sheather-Jones",
                  "Sheather-Jones",
                  "double-Sheather-Jones",
                  "quadruple-Sheather-Jones")

# draw plots, one for each combination of method & distribution
linelabs<-c("Multimod. Exp.",
            "Quarter-SJ",
            "Half-SJ",
            "Sheather-Jones",
            "Double-SJ",
            "Quadruple-SJ")
methlabs<-c("MVN method, ",
            "Hypercuboid method, ",
            "Beta method, ",
            "SIR method, ")
distlabs<-c("MVN distribution",
            "banana distribution",
            "lasso distribution")
samplestats<-matrix(c(0.319,0.637,-0.280,
                      0.459,0.806,0.039,
                      0.619,0.979,0.431),
                    nrow=3)
rownames(samplestats)<-c("mvn","banana","lasso")
colnames(samplestats)<-c("2.5%","50%","97.5%")
svglite::svglite("wassersteins.svg")
par(mfrow=c(4,3))
for(meth in 1:4) {
  for(dist in 1:3) {
    wassplot<-matrix(NA,nrow=6,ncol=10)
    for(i in 1:6) {
      wassplot[i,]<-wassers[[i]][meth,dist,]
    }
    plot(1:10,log10(wassplot[1,]),type='n',
         ylim=c(min(c(min(log10(wassplot)),samplestats[dist,1]))-0.2,
                max(c(max(log10(wassplot)),samplestats[dist,3]))+0.2),
         xlim=c(1,13),
         yaxs='i',xaxt='n',cex.axis=0.8,
         ylab='log_10 sliced 2-Wasserstein',
         xlab='Number of updates',
         main=paste0(methlabs[meth],distlabs[dist]),cex.main=1)
    polygon(x=c(1,10,10,1),
            y=rep(samplestats[dist,c(1,3)],each=2),
            col='#a0a0a0')
    lines(x=c(1,10),
          y=rep(samplestats[dist,2],2),
          lty=3,lwd=1.5)
    for(i in 1:6) {
      lines(1:10,log10(wassplot[i,]))
      text(10.2,log10(wassplot[i,10]),
           labels=linelabs[i],
           adj=0,cex=0.6)
    }
    # axis(2,
    #      at=seq(from=0.4,to=2.8,by=0.4),
    #      labels=as.character(seq(from=0.4,to=2.8,by=0.4)),
    #      cex.axis=0.8)
    axis(1,
         at=1:10,
         labels=as.character(1:10),
         cex.axis=0.8)
  }
}
dev.off()