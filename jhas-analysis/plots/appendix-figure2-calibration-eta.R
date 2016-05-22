### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
#This code reproduces Figure 2 in the appendix
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.

library("pROC")
library("ROCR")

pt.data<-read.csv("jhas-analysis/simulation-data/pt-data-sim.csv")

eta.data<-pt.data$obs_eta
obs.eta<-eta.data[!is.na(eta.data)]
(N<-length(obs.eta))


##get unadjusted results
id.test<-c(1:N)
fitted.eta<-vector(length=N)
for(i in 1:N){
	fitted.eta[i]<-read.csv(paste("jhas-analysis/UNADJ/results/eta-fitted-unadj-",id.test[i],".csv",sep=""))[2]}
fitted.unadj<-unlist(fitted.eta)



##get IOP-BX results
id.test<-c(1:N)
fitted.eta<-vector(length=N)
for(i in 1:N){
	fitted.eta[i]<-read.csv(paste("jhas-analysis/IOP-BX/results/eta-fitted-IOP-BX-",id.test[i],".csv",sep=""))[2]}
fitted.bx<-unlist(fitted.eta)


##get IOP-SURG results
id.test<-c(1:N)
fitted.eta<-vector(length=N)
for(i in 1:N){
	fitted.eta[i]<-read.csv(paste("jhas-analysis/IOP-SURG/results/eta-fitted-IOP-SURG-",id.test[i],".csv",sep=""))[2]}
fitted.surg<-unlist(fitted.eta)


## get LOGISTIC results
pred.ek<-read.csv("jhas-analysis/LOGISTIC/results/pred-eta-logistic-ek.csv")
fitted.logistic<-pred.ek$pred.eta



expit<-function(x){return(exp(x)/(1+exp(x)))}





##CALIBRATION PLOTS FOR ALL MODELS


ns.unadj<-glm(obs.eta~ns(fitted.unadj,2), family="binomial")
pred.unadj<-predict.glm(ns.unadj, ns(fitted.unadj,2), se.fit=TRUE)

order.unadj<-order(fitted.unadj)
post.unadj<-fitted.unadj[order.unadj] #these are my posterior predictions
fit.unadj<-pred.unadj$fit[order.unadj] #this is ns model fit
se.unadj<-pred.unadj$se.fit[order.unadj]
fit.up.unadj<-pred.unadj$fit[order.unadj] + qnorm(0.975)*se.unadj
fit.low.unadj<-pred.unadj$fit[order.unadj] + qnorm(0.025)*se.unadj


pdf("jhas-analysis/plots/appendix-figure2a-calibration-ek-unadj.pdf", width=9, height=7)
par(mar=c(5,4.5,2,2))

plot(obs.eta~fitted.unadj, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.unadj[post.unadj>0.1], rev(post.unadj[post.unadj>0.1])), c(expit(fit.low.unadj[post.unadj>0.1]), rev(expit(fit.up.unadj[post.unadj>0.1]))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.unadj[post.unadj>0.1])~post.unadj[post.unadj>0.1], col="dark blue", lwd=3)

dev.off()




###IOP-Biopsy
ns.bx<-glm(obs.eta~ns(fitted.bx,2), family="binomial")
pred.bx<-predict.glm(ns.bx, ns(fitted.bx,2), se.fit=TRUE)

order.bx<-order(fitted.bx)
post.bx<-fitted.bx[order.bx] #these are my posterior predictions
fit.bx<-pred.bx$fit[order.bx] #this is ns model fit
se.bx<-pred.bx$se.fit[order.bx]
fit.up.bx<-pred.bx$fit[order.bx] + qnorm(0.975)*se.bx
fit.low.bx<-pred.bx$fit[order.bx] + qnorm(0.025)*se.bx


pdf("jhas-analysis/plots/appendix-figure2b-calibration-ek-iop-bx.pdf", width=9, height=7)
par(mar=c(5,4.5,2,2))

plot(obs.eta~fitted.bx, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.bx[post.bx>0.1], rev(post.bx[post.bx>0.1])), c(expit(fit.low.bx[post.bx>0.1]), rev(expit(fit.up.bx[post.bx>0.1]))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.bx[post.bx>0.1])~post.bx[post.bx>0.1], col="dark blue", lwd=3)

dev.off()





##IOP-Surgery
ns.surg<-glm(obs.eta~ns(fitted.surg,2), family="binomial")
pred.surg<-predict.glm(ns.surg, ns(fitted.surg,2), se.fit=TRUE)

order.surg<-order(fitted.surg)
post.surg<-fitted.surg[order.surg] #these are my posterior predictions
fit.surg<-pred.surg$fit[order.surg] #this is ns model fit
se.surg<-pred.surg$se.fit[order.surg]
fit.up.surg<-pred.surg$fit[order.surg] + qnorm(0.975)*se.surg
fit.low.surg<-pred.surg$fit[order.surg] + qnorm(0.025)*se.surg


pdf("jhas-analysis/plots/appendix-figure2c-calibration-ek-iop-surg.pdf", width=9, height=7)
par(mar=c(5,4.5,2,2))

plot(obs.eta~fitted.surg, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.surg[post.surg>0.1], rev(post.surg[post.surg>0.1])), c(expit(fit.low.surg[post.surg>0.1]), rev(expit(fit.up.surg[post.surg>0.1]))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.surg[post.surg>0.1])~post.surg[post.surg>0.1], col="dark blue", lwd=3)

dev.off()




##LOGISTIC REGRESSION
ns.logistic<-glm(obs.eta~ns(fitted.logistic,2), family="binomial")
pred.logistic<-predict.glm(ns.logistic, ns(fitted.logistic,2), se.fit=TRUE)

order.logistic<-order(fitted.logistic)
post.logistic<-fitted.logistic[order.logistic] #these are my posterior predictions
fit.logistic<-pred.logistic$fit[order.logistic] #this is ns model fit
se.logistic<-pred.logistic$se.fit[order.logistic]
fit.up.logistic<-pred.logistic$fit[order.logistic] + qnorm(0.975)*se.logistic
fit.low.logistic<-pred.logistic$fit[order.logistic] + qnorm(0.025)*se.logistic



pdf("jhas-analysis/plots/appendix-figure2d-calibration-ek-logistic.pdf", width=9, height=7)
par(mar=c(5,4.5,2,2))

plot(obs.eta~fitted.logistic, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.logistic[post.logistic>0.1], rev(post.logistic[post.logistic>0.1])), c(expit(fit.low.logistic[post.logistic>0.1]), rev(expit(fit.up.logistic[post.logistic>0.1]))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.logistic[post.logistic>0.1])~post.logistic[post.logistic>0.1], col="dark blue", lwd=3)

dev.off()

