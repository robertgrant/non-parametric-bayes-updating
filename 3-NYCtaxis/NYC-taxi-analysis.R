# Feasibility case study on non-parametric Bayesian updating
# Daily SIR updates on New York taxi data

setwd("~/Documents/Shop-window-NOC/non-parametric-updating/NYCtaxis")
library(readr)
library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
extract<-rstan::extract # avoid ambiguity
options(mc.cores = 4)
rmvn<-mnormt::rmnorm
rtmvn<-mnormt::rmtruncnorm
mclapply<-parallel::mclapply # getll() uses 4 cores
# how many journeys to take from each day for an initial
# full model in Stan
n_pump_prime<-30000 

# function to turn taxidata file into matrices for analysis
data2mats<-function(filename,weekday) {
  # read in file
  taxidata<-read_csv(filename,col_types=cols())
  
  # make rotated coordinates
  rotangle<-pi/2-atan(1.366604) # anticlockwise to make 5th avenue vertical
  taxidata<-mutate(taxidata,
                   updown=sin(rotangle)*pickup_longitude + 
                     cos(rotangle)*pickup_latitude,
                   westeast=cos(rotangle)*pickup_longitude - 
                     sin(rotangle)*pickup_latitude)
  
  # classify neighbourhoods
  bottom<- sin(rotangle)*-74.016382 + 
    cos(rotangle)*40.700062 # updown for off Battery Point
  top<-sin(rotangle)*-73.975864 + 
    cos(rotangle)*40.809038 # updown for 110th street
  fifth<-cos(rotangle)*-73.949209 - 
    sin(rotangle)*40.796788 # eastwest for 5th Ave
  taxidata<-mutate(taxidata,
                   eastof5th=as.numeric(westeast>fifth),
                   updown10=floor((updown-bottom)/(top-bottom)*10))
  # restrict to Manhattan rectangle         
  taxidata<-filter(taxidata,
                   between(updown10,0,9))
  taxidata<-filter(taxidata,
                   between(westeast,-83.8,-83.74))
  # neighbourhood number
  taxidata<-mutate(taxidata,
                   area=10*eastof5th+updown10)
  
  # predictor matrices
  taxidata$area<-as.factor(taxidata$area)
  taxidata$wkday<-as.factor(taxidata$wkday)
  taxidata$hour<-as.factor(taxidata$hour)
  areamat<-model.matrix(~0+area,
                        data=taxidata)
  weekdaymat<-matrix(0,nrow=NROW(taxidata),ncol=7)
  weekdaymat[,weekday]<-1
  colnames(weekdaymat)<-c("Monday",
                          "Tuesday",
                          "Wednesday",
                          "Thursday",
                          "Friday",
                          "Saturday",
                          "Sunday")
  hourmat<-model.matrix(~0+hour,
                        data=taxidata)
  # combined design matrix uses one weekday column like an intercept
  predmat<-cbind(areamat[,-1],weekdaymat,hourmat[,-1])
  # qr_predmat<-qr(predmat)
  # q_predmat<-qr.Q(qr_predmat)
  # r_predmat<-qr.R(qr_predmat)
  logtotal<-taxidata$logtotal
  
  return(list(logtotal=logtotal,
              predmat=predmat,
              n=NROW(predmat),
              p=NCOL(predmat)))
}

mats<-data2mats("taxi_data_2013_1_1.csv",1)
logtotal7k<-mats$logtotal[1:n_pump_prime]
predmat7k<-mats$predmat[1:n_pump_prime,]
for(i in 2:7) {
  mats<-data2mats(paste0("taxi_data_2013_1_",i,".csv"),i)
  logtotal7k<-c(logtotal7k,mats$logtotal[1:n_pump_prime])
  predmat7k<-rbind(predmat7k,mats$predmat[1:n_pump_prime,])
}

# 
# stanQRcode<-'
# data {
#   int n;
#   real logtotal[n];
#   matrix[n,51] q_pred;
# }
# parameters {
#   vector[51] theta;
#   real<lower=0> sigma;
# } 
# model {
#   logtotal ~ normal((q_pred*sqrt(n-1))*theta,sigma);
# }
# '
# 
# stanQRmod<-stan_model(model_code=stanQRcode)
# 
# # experimental fit with 10,000 journeys
# qr_predmat<-qr(predmat[1:10000,])
# q_predmat<-qr.Q(qr_predmat)
# r_predmat<-qr.R(qr_predmat)
# stan_QR1E4_fit <- sampling(stanQRmod,
#                          data=list(logtotal=logtotal[1:10000],
#                                    q_pred=q_predmat,
#                                    n=10000),
#                          chains=2,
#                          cores=4,
#                          iter=2000,
#                          warmup=1000,
#                          seed=477)
# stan_summary_QR1E4<-list(summary=summary(stan_QR1E4_fit)$summary,
#                        posterior=extract(stan_QR1E4_fit))
# beta_posterior<-solve(r_predmat)%*%stan_summary_QR1E4$posterior$theta

stancode<-'
data {
  int n;
  int p;
  real logtotal[n];
  matrix[n,p] predictors;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma;
} 
model {
  beta ~ normal(0,2); // prior
  logtotal ~ normal(predictors*beta, sigma); // likelihood
}
'
stanmod<-stan_model(model_code=stancode)
# stan_1E4_fit <- sampling(stanmod,
#                            data=list(logtotal=logtotal[1:10000],
#                                      predictors=predmat[1:10000,],
#                                      n=10000),
#                            chains=2,
#                            cores=4,
#                            iter=2000,
#                            warmup=1000,
#                            seed=943)
# 
# stan_summary_1E4<-list(summary=summary(stan_1E4_fit)$summary,
#                          posterior=extract(stan_1E4_fit))

stan_7k_fit <- sampling(stanmod,
                         data=list(logtotal=logtotal7k,
                                   predictors=predmat7k,
                                   n=length(logtotal7k),
                                   p=NCOL(predmat7k)),
                         chains=2,
                         cores=4,
                         iter=3000,
                         warmup=1000,
                         seed=9623)

stan_summary_7k<-list(summary=summary(stan_7k_fit)$summary,
                       posterior=extract(stan_7k_fit))



############################################################
##################   SIR updating 1-7 Jan   ################
############################################################
min_posteriors<-1
max_posteriors<-200000
bw_scale<-1 # multiply SJ bandwidth


# bandwidth function
getbw<-function(prior_draws,scale) { scale*apply(prior_draws,2,bw.SJ) }

# # function to find nearest neighbours in each dimension from prior draws
# getnn<-function(prior_draws,index) {
#   edges<-matrix(NA,nrow=NCOL(prior_draws),ncol=2)
#   for(i in 1:NCOL(prior_draws)) {
#     sorted_prior<-sort(prior_draws[,i])
#     target_draw<-prior_draws[index,i]
#     target_index<-which(sorted_prior==target_draw)
#     edges[i,]<-c(sorted_prior[target_index-1],
#                  sorted_prior[target_index+1])
#   }
#   return(edges)
# }

# function to sample from truncated kernel density within hypercuboid
# 3 SEs either side of the MLE
# mle and se are both p-vectors
# n_ses is the number of SEs to extend out from the MLE
# posi is a vector of column numbers that must be positive
oversample<-function(prior_draws,mle,se,n_sample,bw,n_ses,posi) {
  p<-NCOL(prior_draws)
  k<-NROW(prior_draws)
  osamp<-matrix(NA,nrow=n_sample,ncol=p)
  trunc_cdf<-rep(NA,k)
  # this loop is an obvious candidate for parallelising
  for(i in 1:k) {
    # in each dimension, find the probability between +/- n_ses SEs
    trunc_lcdfs<-rep(NA,p)
    for(j in 1:p) {
      #cat(paste("evaluating dimension",j,"\n"))
      # we are using MVN (uncorrelated) kernels
      if(any(posi==j)) {
        if(mle[j]+n_ses*se[j]>0) {
          trunc_lcdfs[j]<-log(pnorm(mle[j]+n_ses*se[j],
                                    mean=prior_draws[i,j],
                                    sd=bw[j],
                                    log.p=FALSE) -
                              pnorm(0,
                                    mean=prior_draws[i,j],
                                    sd=bw[j],
                                    log.p=FALSE))
        }
        else if(mle[j]+n_ses*se[j]<=0) {
          # it would be more efficient not to evaluate any dimension
          trunc_lcdfs[j]<-(-Inf)
        }
      }
      else {
        trunc_lcdfs[j]<-log(pnorm(mle[j]+n_ses*se[j],
                                 mean=prior_draws[i,j],
                                 sd=bw[j],
                                 log.p=FALSE) -
                            pnorm(mle[j]-n_ses*se[j],
                                  mean=prior_draws[i,j],
                                  sd=bw[j],
                                  log.p=FALSE))
      }
    }
    # joint prob of lying in hypercuboid:
    trunc_cdf[i]<-exp(sum(trunc_lcdfs)) # it's OK if this rounds to zero
  }
  # if likelihood and prior overlap poorly (as with a first
  # batch), then use a uniform distribution instead
  if(sum(trunc_cdf)==0) {
    for(j in 1:p) {
      osamp[,j]<-runif(n_sample,
                       mle[j]-n_ses*se[j],
                       mle[j]+n_ses*se[j])
    }
  }
  else {
    # find number of prior local oversamples per prior draw:
    n_samples_per_draw<-n_sample*trunc_cdf/sum(trunc_cdf)
    # sample from truncated MVNs:
    osamp_count<-1
    for(i in 1:k) {
      if(n_samples_per_draw[i]>0) {
        if(any(posi==j)) {
          osamp[osamp_count:(osamp_count+n_samples_per_draw[i]-1)]<-
            rtmvn(n=n_samples_per_draw[i],
                  mean=prior_draws[i,],
                  varcov=diag(bw^2),
                  lower=rep(0,p),
                  upper=mle+n_ses*se)
          osamp_count<-osamp_count+n_samples_per_draw[i]
        }
        else {
          osamp[osamp_count:(osamp_count+n_samples_per_draw[i]-1)]<-
            rtmvn(n=n_samples_per_draw[i],
                  mean=prior_draws[i,],
                  varcov=diag(bw^2),
                  lower=mle-n_ses*se,
                  upper=mle+n_ses*se)
          osamp_count<-osamp_count+n_samples_per_draw[i]
        }
      }
    }
  }
  return(osamp)
}

# add 7k posterior sample to list 
priors<-list(cbind(stan_summary_7k$posterior$beta,
                   stan_summary_7k$posterior$sigma))

# # oversample the 7k posterior sample
# # find MLE
# lmfit<-lm(logtotal~0+predmat)
# maxloglik<-logLik(lmfit)
# lmsumm<-summary(lmfit)
# # note that we have to plug in something for the SE of the residual SD
# # and that any parameters uninformed by the new batch must be carried 
# # forward from the prior sample
# delete_mle<-c(lmsumm$coefficients[1:19,1],
#               apply(priors[[1]][,20:26],2,mean),
#               lmsumm$coefficients[21:43,1],
#               lmsumm$sigma)
# delete_mle[20]<-lmsumm$coefficients[20,1]
# delete_se<-c(lmsumm$coefficients[1:19,2],
#              apply(priors[[1]][,20:26],2,sd),
#              lmsumm$coefficients[21:43,2],
#              sd(priors[[1]][,50]))
# delete_se[20]<-lmsumm$coefficients[20,2]
# 
# deleteme<-oversample(prior_draws=priors[[1]],
#                      mle=delete_mle,
#                      se=delete_se,
#                      n_sample=10000,
#                      bw=getbw(priors[[1]],bw_scale),
#                      n_ses=2,
#                      posi=50)

# empty list for summaries
summaries<-NULL
# empty list for n_posteriors (per prior draw)
n_posteriors<-NULL
# empty list for ML inference
ml_summaries<-NULL
# empty list for degeneracy flags
degen<-NULL

# function to get likelihoods
getll<-function(prior_draws,predmat,logtotal) {
  k<-NROW(prior_draws)
  p<-NCOL(prior_draws)
  llik<-rep(NA,k)
  cat("Getting likelihoods\n")
  pb<-txtProgressBar(min=0,max=NROW(prior_draws),style=3)
  for(i in 1:k) {
    setTxtProgressBar(pb,i)
    pred<-predmat%*%prior_draws[i,1:(p-1)] # not sigma
    resid<-logtotal-pred
    llik[i]<-sum(dnorm(resid,0,prior_draws[i,p],
                       log=TRUE)) 
  }
  close(pb)
  return(llik)
}

# function to get n of posterior draws per prior draw
getnpost<-function(prior_draws,max_posteriors,predmat,logtotal){
  k<-NROW(prior_draws)
  kk<-floor(k/4)
  list_llik<-mclapply(list(prior_draws[1:kk,],
                           prior_draws[(kk+1):(2*kk),],
                           prior_draws[(2*kk+1):(3*kk),],
                           prior_draws[(3*kk+1):k,]),
                      FUN=getll,
                      mc.cores=4,
                      predmat=predmat,
                      logtotal=logtotal)
  llik<-do.call(rbind,list_llik)
  return(floor(exp(llik+(log(max_posteriors)-max(llik)))))
}


# loop over days
for(update in 1:7) {
  print(paste("Now updating for",update,"January"))
  print(Sys.time())
  # weekday variable
  if(update%%7==0) {
    wkday<-7
  }
  else {
    wkday<-update%%7
  }
  # get prior
  prior_draws<-priors[[update]]
  # get data (not 1st 1000 journeys)
  mats<-data2mats(paste0("taxi_data_2013_1_",update,".csv"),wkday)
  logtotal<-mats$logtotal[(n_pump_prime+1):mats$n]
  predmat<-mats$predmat[(n_pump_prime+1):mats$n,]
  print(Sys.time())
  
  # bandwidths
  bandwidths<-getbw(prior_draws,bw_scale)
  # find MLE
  lmfit<-lm(logtotal~0+predmat)
  maxloglik<-logLik(lmfit)
  lmsumm<-summary(lmfit)
  ml_summaries[[update]]<-lmsumm
  # get likelihoods and number of posterior draws at each prior draw
  n_posteriors[[update]]<-getnpost(prior_draws,
                                   max_posteriors,
                                   predmat,
                                   logtotal)
  
  # if degenerate, oversample and try again
  if(sum(n_posteriors[[update]]>0)==1) {
    degen[[update]]<-TRUE
    # assemble MLE & SE, dealing with weekdays not present in batch
    delete_mle<-c(lmsumm$coefficients[1:19,1],
                  apply(priors[[1]][,20:26],2,mean),
                  lmsumm$coefficients[21:43,1],
                  lmsumm$sigma)
    delete_mle[19+wkday]<-lmsumm$coefficients[19+wkday,1]
    delete_se<-c(lmsumm$coefficients[1:19,2],
                 apply(priors[[1]][,20:26],2,sd),
                 lmsumm$coefficients[21:43,2],
                 sd(priors[[1]][,50]))
    delete_se[19+wkday]<-lmsumm$coefficients[19+wkday,2]
    
    prior_draws<-oversample(prior_draws,
                            mle=delete_mle,
                            se=delete_se,
                            n_sample=NROW(prior_draws),
                            bw=bandwidths,
                            n_ses=2,
                            posi=50)
    n_posteriors[[update]]<-getnpost(prior_draws,
                                     max_posteriors,
                                     predmat,
                                     logtotal)
  }
  else {
    degen[[update]]<-FALSE
  }
  
  # in practice, we had to reduce the 7th January (and
  # really, the max_posteriors variable should be lowered
  # automatically when n_posteriors gets big)
  # n_posteriors[[7]]<-floor(n_posteriors[[7]]/100)
  
  # sample from kernels
  pb<-txtProgressBar(min=0,max=NROW(prior_draws),style=3)
  priors[[update+1]]<-matrix(NA,
                             nrow=sum(n_posteriors[[update]]),
                             ncol=NCOL(prior_draws))
  lastrow<-c(0,cumsum(n_posteriors[[update]]))
  for(i in 1:NROW(prior_draws)) {
    setTxtProgressBar(pb,i)
    if(n_posteriors[[update]][i]>0) {
      priors[[update+1]][((lastrow[i]+1):lastrow[i+1]),]<-
                                  rmvn(n_posteriors[[update]][i],
                                       as.matrix(prior_draws[i,]),
                                       diag(bandwidths^2))
    }
  }
  close(pb)
  # down-sample and store in posteriors list
  downindex<-sample(1:sum(n_posteriors[[update]]),NROW(priors[[1]]),
                      replace=(sum(n_posteriors[[update]])<NROW(prior_draws)))
  priors[[update+1]]<-priors[[update+1]][downindex,]
  print(Sys.time())
  summaries[[update]]<-cbind(apply(priors[[update+1]],
                                   2,
                                   mean),
                             apply(priors[[update+1]],
                                   2,
                                   sd))
}
# end loop

# how many prior draws are being sampled in the posterior?
unlist(lapply(n_posteriors,function(x){sum(x>0)}))

# save
saveRDS(list(priors=priors,
             summaries=summaries,
             n_posteriors=n_posteriors,
             degen=degen,
             ml_summaries=ml_summaries),
        "SIR-output.rds")

# some summary stats
summs<-array(NA,dim=c(50,8,4))
for(paranumber in 1:50) {
  summs[paranumber,1,1]<-stan_summary_7k$summary[paranumber,1]
  summs[paranumber,1,2]<-stan_summary_7k$summary[paranumber,2]
  summs[paranumber,2:8,1]<-unlist(lapply(summaries,function(x){x[paranumber,1]}))
  summs[paranumber,2:8,2]<-unlist(lapply(summaries,function(x){x[paranumber,2]}))
  summs[paranumber,,3]<-summs[paranumber,,1]-1.96*summs[paranumber,,2]
  summs[paranumber,,4]<-summs[paranumber,,1]+1.96*summs[paranumber,,2]
}

paranames<-c("Area 2 parameter",
             "Area 3 parameter",
             "Area 4 parameter",
             "Area 5 parameter",
             "Area 6 parameter",
             "Area 7 parameter",
             "Area 8 parameter",
             "Area 9 parameter",
             "Area 10 parameter",
             "Area 11 parameter",
             "Area 12 parameter",
             "Area 13 parameter",
             "Area 14 parameter",
             "Area 15 parameter",
             "Area 16 parameter",
             "Area 17 parameter",
             "Area 18 parameter",
             "Area 19 parameter",
             "Area 20 parameter",
             "Monday parameter",
             "Tuesday parameter",
             "Wednesday parameter",
             "Thursday parameter",
             "Friday parameter",
             "Saturday parameter",
             "Sunday parameter",
             "1am parameter",
             "2am parameter",
             "3am parameter",
             "4am parameter",
             "5am parameter",
             "6am parameter",
             "7am parameter",
             "8am parameter",
             "9am parameter",
             "10am parameter",
             "11am parameter",
             "12am parameter",
             "1pm parameter",
             "2pm parameter",
             "3pm parameter",
             "4pm parameter",
             "5pm parameter",
             "6pm parameter",
             "7pm parameter",
             "8pm parameter",
             "9pm parameter",
             "10pm parameter",
             "11pm parameter",
             "Residual SD")

png("means_update_0_1.png",height=360,width=480)
plot(1:50,summs[,1,1],
     xaxt='n',xlab='',ylab="Posterior means")
points(1:50,summs[,2,1],pch=19,cex=0.6)
#points(1:50,summs[,8,1],pch=19,cex=0.2,col='#b62525')
axis(1,
     at=c(1,10,20,27,32),
     labels=c("Area 2",
              "Area 11",
              "Monday",
              "1am",
              "6am"),
     cex.axis=0.7)
dev.off()

# plots


summplots<-function(summs,paranumber,paraname) {
  plot(0:7,
       summs[paranumber,,1],
       type='b',
       cex=0.5,
       pch=19,
       main=paraname,
       ylab="posterior mean +/- 1.96SD",
       xlab="Day / update",
       ylim=c(min(summs[paranumber,,3]),max(summs[paranumber,,4])))
  for(i in 1:8) {
    lines(rep(i-1,2),
          summs[paranumber,i,3:4],
          lwd=1.6)
  }
}

for(p in 1:50) {
  png(paste0("summary_plot_",p,".png"),
      width=800,
      height=600)
  summplots(summs,p,paranames[p])
  dev.off()
}

# SD evolution
png("SD_evolution_area.png",width=800,height=600)
plot(0:7,summs[1,,2],type='l',col='#00000060',
     ylim=c(min(summs[1:19,,2]),max(summs[1:19,,2])),
     xlab="Updates",ylab="Posterior SD",main="Area parameters")
for(j in 2:19) {
  lines(0:7,summs[j,,2],col='#00000060')
}
dev.off()

png("SD_evolution_weekdays.png",width=800,height=600)
plot(0:7,summs[20,,2],type='l',col='#00000060',
     ylim=c(min(summs[20:26,,2]),max(summs[20:26,,2])),
     xlab="Updates",ylab="Posterior SD",main="Weekday parameters")
for(j in 21:26) {
  lines(0:7,summs[j,,2],col='#00000060')
}
dev.off()

png("SD_evolution_hours.png",width=800,height=600)
plot(0:7,summs[27,,2],type='l',col='#00000060',
     ylim=c(min(summs[27:49,,2]),max(summs[27:49,,2])),
     xlab="Updates",ylab="Posterior SD",main="Hour parameters")
for(j in 28:49) {
  lines(0:7,summs[j,,2],col='#00000060')
}
dev.off()

# mean evolution
png("mean_evolution_areas.png",width=800,height=600)
plot(0:7,summs[1,,1],type='l',col='#00000060',
     ylim=c(min(summs[1:19,,1]),max(summs[1:19,,1])),
     xlab="Updates",ylab="Posterior mean",main="Area parameters")
for(j in 2:19) {
  lines(0:7,summs[j,,1],col='#00000060')
}
dev.off()

png("mean_evolution_weekdays.png",width=800,height=600)
plot(0:7,summs[20,,1],type='l',col='#00000060',
     ylim=c(min(summs[20:26,,1]),max(summs[20:26,,1])),
     xlab="Updates",ylab="Posterior mean",main="Weekday parameters")
for(j in 21:26) {
  lines(0:7,summs[j,,1],col='#00000060')
}
dev.off()

png("mean_evolution_hours.png",width=800,height=600)
plot(0:7,summs[27,,1],type='l',col='#00000060',
     ylim=c(min(summs[27:49,,1]),max(summs[27:49,,1])),
     xlab="Updates",ylab="Posterior mean",main="Hour parameters")
for(j in 28:49) {
  lines(0:7,summs[j,,1],col='#00000060')
}
dev.off()

# bivariate priors
hullpoints<-function(x) {
  temp<-chull(x)
  return(rbind(x[temp,],x[temp[1],]))
}
priorplot<-function(para1,para2,update) {
  if(update==1) {maintext<-"Stan 210k model"}
  else { maintext<-paste("After",update-1,"Jan SIR")}
  plot(priors[[update]][,c(para1,para2)],cex=0.2,col='#00000060',
       xlab=paranames[para1],ylab=paranames[para2],
       main=maintext)
  if(update<8) {
    lines(hullpoints(priors[[update+1]][,c(para1,para2)]),
          col='#b62525')
  }
}
png("bivariate_priors_areas.png",width=960,height=1440)
par(mfrow=c(4,2))
for(j in 1:8) {
  priorplot(1,2,j)
}
dev.off()

png("bivariate_priors_weekdays.png",width=960,height=1440)
par(mfrow=c(4,2))
for(j in 1:8) {
  priorplot(20,24,j)
}
dev.off()

png("bivariate_priors_hours.png",width=960,height=1440)
par(mfrow=c(4,2))
for(j in 1:8) {
  priorplot(31,39,j)
}
dev.off()

