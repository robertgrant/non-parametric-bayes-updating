# Generate summary stats and graphics for simulation
# study of non-parametric Bayesian updating
setwd('~/Documents/Shop-window-NOC/non-parametric-updating/2020-sims/outputs')

# load compiled results
load('../combined-results-cloud-sims-4000-draws.RData')
comb_sir<-readRDS('../combined-results-sir.rds')
truepars<-c(2.0,0.5,-1.5)





############################################################
##########  craft charts of first simulations  #############
############################################################

# 'posteriors' expects a list of length nrows*ncols, each
  # element of which is a list of length 10, each
  # element of which is a data frame or matrix of the posterior draws
  #
  # 'all_posterior' is a list of length nrows*ncols, each
  # element of which is a data frame or matrix containing the posterior
  # draws from analysing all data together
  #
  # 'bw' is either:
  #     a list of length nrows*ncols, each
  #     element of which is a list of length 10, each
  #     element of which is a vector of length 3
  # or:
  #     "SJ", in which case it is calculated ad hoc (without 
  #     multimodality expansion) -- this is really just for the
  #     parametric updates.
craft_chart<-function(posteriors,bw,nrows,ncols,all_posterior,
                      savename,truepars,vspacing) {
    n_charts<-nrows*ncols
    if(length(posteriors)!=n_charts) { warning("nrows*ncols != length(posteriors)")}
    for(p in 1:3) {
      png(paste0("updates_",savename,"_beta",p,".png"),width=1000,height=2000)
      par(mfrow=c(nrows,ncols))
      for(sim in 1:n_charts) {
        # draw posterior from first batch
        if(is.character(bw)) {
          this_bw<-bw.SJ(posteriors[[sim]][[1]][,p])
        }
        else {
          this_bw<-bw[[sim]][[2]][p]
        }
        plot(density(posteriors[[sim]][[1]][,p],bw=this_bw)$x,
             density(posteriors[[sim]][[1]][,p],bw=this_bw)$y+(10-1)*vspacing,
             ylim=c(0,10*vspacing),yaxt='n',type='l',
             main=paste0("Simulation ",sim),
             xlab=paste0("beta",p),ylab="Update batches")
        axis(side=2,at=seq(from=0,by=vspacing,length.out=10),
             labels=10:1,cex.axis=0.8)
        # draw posteriors for subsequent batches
        for(i in 2:9) {
          if(is.character(bw)) {
            this_bw<-bw.SJ(posteriors[[sim]][[i]][,p])
          }
          else {
            this_bw<-bw[[sim]][[i+1]][p]
          }
          lines(density(posteriors[[sim]][[i]][,p],bw=this_bw)$x,
                density(posteriors[[sim]][[i]][,p],bw=this_bw)$y+(10-i)*vspacing)
          polygon(c(density(posteriors[[sim]][[i]][,p],bw=this_bw)$x,
                    density(posteriors[[sim]][[i]][,p],bw=this_bw)$x[1]),
                  c(density(posteriors[[sim]][[i]][,p],bw=this_bw)$y+(10-i)*vspacing,
                    density(posteriors[[sim]][[i]][,p],bw=this_bw)$y[1]+(10-i)*vspacing),
                  col='#ffffff90',border=NA)
        }
        # draw posterior from final batch
        if(is.character(bw)) {
          this_bw<-bw.SJ(posteriors[[sim]][[10]][,p])
        }
        else {
          # re-use previous batch's bandwidths (for hypercuboid & SIR)
          this_bw<-bw[[sim]][[10]][p] 
        }
        lines(density(posteriors[[sim]][[10]][,p],bw=this_bw)$x,
              density(posteriors[[sim]][[10]][,p],bw=this_bw)$y)
        polygon(c(density(posteriors[[sim]][[10]][,p],bw=this_bw)$x,
                  density(posteriors[[sim]][[10]][,p],bw=this_bw)$x[1]),
                c(density(posteriors[[sim]][[10]][,p],bw=this_bw)$y,
                  density(posteriors[[sim]][[10]][,p],bw=this_bw)$y[1]),
                col='#ffffff90',border=NA)
        # draw posterior from analysing all data together
        bw_all<-bw.SJ(all_posterior[[sim]][,p])
        lines(density(all_posterior[[sim]][,p],bw=bw_all)$x,
              density(all_posterior[[sim]][,p],bw=bw_all)$y,
              lwd=2)
        # draw true parameter value
        #abline(v=truepars[p],lty=3)
      }
      dev.off()
    }
}

# parametric updating
craft_chart(posteriors=posterior_para,
            bw="SJ",
            nrows=10,
            ncols=3,
            all_posterior=posterior_all,
            savename="para",
            truepars=truepars,
            vspacing=1.4)

# unbounded kernel (MVN) updating
craft_chart(posteriors=posterior_mvn,
            bw=bandwidth_mvn,
            nrows=10,
            ncols=3,
            all_posterior=posterior_all,
            savename="mvn",
            truepars=truepars,
            vspacing=1.4)


# hypercuboid unbounded kernel (MVN) updating
craft_chart(posteriors=posterior_mvncube,
            bw=bandwidth_mvncube,
            nrows=10,
            ncols=3,
            all_posterior=posterior_all,
            savename="mvncube",
            truepars=truepars,
            vspacing=1.7)


# bounded kernel (beta) updating
craft_chart(posteriors=posterior_beta,
            bw=bandwidth_beta,
            nrows=10,
            ncols=3,
            all_posterior=posterior_all,
            savename="beta",
            truepars=truepars,
            vspacing=1.4)


# sampling importance resampling MVN updating
craft_chart(posteriors=comb_sir$posterior_sir,
            bw=comb_sir$bandwidth_sir,
            nrows=10,
            ncols=3,
            all_posterior=posterior_all,
            savename="SIR",
            truepars=truepars,
            vspacing=1.4)



############################################################
############  mean and 95% reference range of  #############
##########  simulation estimates at each batch  ############
############################################################

summary_chart<-function(means,sds) {
  
}


######################################################################
###############  table of of time, bias, coverage  ###################
##########  K-L divergence and efficiency (posterior SD)  ############
######################################################################

sim_summ<-function(means,sds,all_mean,all_sd) {
  
}


##############################################################################
###############  theoretical chart of beta neighbours and  ###################
##############  occupied hypercubes as dimension increases  ##################
##############################################################################







######################################################################
###############  convergence of all means and SDs  ###################
######################################################################



