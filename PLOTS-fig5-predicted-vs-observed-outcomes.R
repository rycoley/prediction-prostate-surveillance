
### Rebecca Yates Coley rycoley@gmail.com
### Code for "Bayesian Joint Hierarchical Model for Prediction of Latent Health States with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in biopsy, reclassification, and surgery data for individuals, MCMC predictions of these outcomes, and make calibration plots in Figure 5


### LOAD PACKAGES
library("splines")

### DEFINE FUNCTIONS
expit<-function(x){return(exp(x)/(1+exp(x)))}


### LOAD DATA
pt.data<-read.csv("simulation-data/pt-data-sim.csv")
data.use<-read.csv("simulation-data/bx-data-sim.csv")

##remove observations (since this model was run with some observations removed in order to demonstrate individualized predictions)
pred.ids<-c(250:261) 
data.use<-data.use[!(data.use$subj%in%pred.ids & data.use$time>5),]


## biopsy data
bx.data<-data.use[!is.na(data.use$bx.here),]
(n_bx<-dim(bx.data)[1])
BX<-as.numeric(bx.data$bx.here)

## reclassification data
rc.data<-data.use[data.use$bx.here==1 & !is.na(data.use$bx.here),]
(n_rc<-dim(rc.data)[1])
RC<-as.numeric(rc.data$rc)

## surgery data
SURG<-as.numeric(data.use$rrp)
(n_surg<-dim(data.use)[1])


### LOAD MCMC RESULTS

## biopsy predictions
p.bx<-read.csv(paste("results/jags-prediction-iop-p_bx-1.csv",sep=""))
p.bx<-as.matrix(p.bx[,2:(n_bx+1)])

## reclassification predictions
p.rc<-read.csv(paste("results/jags-prediction-iop-p_rc-1.csv",sep=""))
p.rc<-as.matrix(p.rc[,2:(n_rc+1)])

##sugery predictions
p.surg<-read.csv(paste("results/jags-prediction-iop-p_surg-1.csv",sep=""))
p.surg<-as.matrix(p.surg[,2:(n_surg+1)])



### FIND NS FIT BETWEEN PREDICTIONS AND OBSERVATIONS

## biopsy
p.bx.mean<-apply(p.bx,2,mean)
mod.bx.mean<-glm(BX~ns(p.bx.mean,4),family="binomial")

pred.p.bx<-predict.glm(mod.bx.mean, ns(p.bx.mean,4), se.fit=TRUE)

order.p.bx<-order(p.bx.mean)
post.p.bx<-p.bx.mean[order.p.bx]

fit.p.bx<-expit(pred.p.bx$fit[order.p.bx])
se.p.bx<-pred.p.bx$se.fit[order.p.bx]
fit.up.p.bx<-expit(pred.p.bx$fit[order.p.bx] + qnorm(0.975)*se.p.bx)
fit.low.p.bx<-expit(pred.p.bx$fit[order.p.bx] + qnorm(0.025)*se.p.bx)

##reclassifcation
p.rc.mean<-apply(p.rc,2,mean)
mod.rc.mean<-glm(RC~ns(p.rc.mean,4),family="binomial")
pred.p.rc<-predict.glm(mod.rc.mean, ns(p.rc.mean,4), se.fit=TRUE)



order.p.rc<-order(p.rc.mean)
post.p.rc<-p.rc.mean[order.p.rc]

fit.p.rc<-expit(pred.p.rc$fit[order.p.rc])
se.p.rc<-pred.p.rc$se.fit[order.p.rc]
fit.up.p.rc<-expit(pred.p.rc$fit[order.p.rc] + qnorm(0.975)*se.p.rc)
fit.low.p.rc<-expit(pred.p.rc$fit[order.p.rc] + qnorm(0.025)*se.p.rc)


##surgery
p.surg.mean<-apply(p.surg,2,mean)
mod.rrp.mean<-glm(RRP~ns(p.surg.mean,4),family="binomial")

pred.p.surg<-predict.glm(mod.rrp.mean, ns(p.surg.mean,4), se.fit=TRUE)

order.p.surg<-order(p.surg.mean)
post.p.surg<-p.surg.mean[order.p.surg]

fit.p.surg<-pred.p.surg$fit[order.p.surg]
se.p.surg<-pred.p.surg$se.fit[order.p.surg]
fit.up.p.surg<-pred.p.surg$fit[order.p.surg] + qnorm(0.975)*se.p.surg
fit.low.p.surg<-pred.p.surg$fit[order.p.surg] + qnorm(0.025)*se.p.surg




### MAKE CALIBRATION PLOTS

### biopsy

pdf("plots/figure5-calibration-biopsy.pdf", height=7, width=9)
par(mar=c(5,5.5,2,5))

plot(BX~p.bx.mean, pch=3,  ylab="", xlab="", xaxt="n", yaxt="n") 
mtext("Posterior P(Biopsy)",1,line=3, cex=1.5)
mtext("Observed P(Biopsy)",2,line=2.75, cex=1.5)

Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)


polygon(x=c(post.p.bx, rev(post.p.bx)), c(fit.low.p.bx, rev(fit.up.p.bx)) , col="lightblue", border="lightblue")
abline(0,1, lwd=2, lty="dotted")

lines(fitted(mod.bx.mean)[order(p.bx.mean)]~p.bx.mean[order(p.bx.mean)], col="blue", lwd=4)
dev.off()



### reclassifcation
pdf("plots/figure5-calibration-reclassification.pdf", height=7, width=9)
par(mar=c(5,5.5,2,5))

plot(RC~p.rc.mean, pch=3,  ylab="", xlab="", xaxt="n", yaxt="n") 
mtext("Posterior P(Reclassification)",1,line=3, cex=1.5)
mtext("Observed P(Reclassification)",2,line=2.75, cex=1.5)

Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.p.rc, rev(post.p.rc)), c(fit.low.p.rc, rev(fit.up.p.rc)) , col="lightblue", border="lightblue")
points(RC~p.rc.mean, pch=3)
abline(0,1, lwd=2, lty="dotted")

lines(fitted(mod.rc.mean)[order(p.rc.mean)]~p.rc.mean[order(p.rc.mean)], col="blue", lwd=4)
dev.off()


### surgery
pdf("plots/figure5-calibration-surgery.pdf", height=7, width=9)
par(mar=c(5,5.5,2,5))

plot(SURG~p.surg.mean, pch=3,  ylab="", xlab="", xaxt="n", yaxt="n") 
mtext("Posterior P(Surgery)",1,line=3, cex=1.5)
mtext("Observed P(Surgery)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)


polygon(x=c(post.p.surg, rev(post.p.surg)), c(expit(fit.low.p.surg), rev(expit(fit.up.p.surg))) , col="lightblue", border="lightblue")

abline(0,1, lwd=2, lty="dotted")

lines(fitted(mod.rrp.mean)[order(p.surg.mean)]~p.surg.mean[order(p.surg.mean)], col="blue", lwd=4)

dev.off()


