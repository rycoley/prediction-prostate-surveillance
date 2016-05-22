### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
#This code reproduces Figure 21 in the appendix


library("pROC")

SEED<-1
set.seed(SEED)

source("simulations/simulation-data/simulate-pca-data.R")

true.eta<-pt_data$eta_true
n_eta_known<-sum(pt_data$surg)


pred.ek<-read.csv("simulations/UNADJ/results/eta-known-true-vs-pred-unadj-1.csv")
pred.unadj.ek<-pred.ek$predicted

pred.ek<-read.csv("simulations/IOP-BX/results/eta-known-true-vs-pred-IOP-BX-1.csv")
pred.bx.ek<-pred.ek$predicted


pred.ek<-read.csv("simulations/IOP-SURG/results/eta-known-true-vs-pred-IOP-SURG-1.csv")
pred.surg.ek<-pred.ek$predicted


pred.ek<-read.csv("simulations/IOP/results/eta-known-true-vs-pred-iop-1.csv")
pred.iop.ek<-pred.ek$predicted

res<-read.csv("simulations/LOGISTIC/results/eta-predictions-LOGISTIC.csv")
names(res)
pred.log.ek<-res$predicted[1:n_eta_known]



expit<-function(x){return(exp(x)/(1+exp(x)))}



### CALIBRATION - ETA KNOWN
obs.eta<-true.eta[1:n_eta_known]


#UNADJ
fitted.unadj<-pred.unadj.ek

ns.unadj<-glm(obs.eta~ns(fitted.unadj,2), family="binomial")
pred.unadj<-predict.glm(ns.unadj, ns(fitted.unadj,2), se.fit=TRUE)

order.unadj<-order(fitted.unadj)
post.unadj<-fitted.unadj[order.unadj] #these are my posterior predictions
fit.unadj<-pred.unadj$fit[order.unadj] #this is ns model fit
se.unadj<-pred.unadj$se.fit[order.unadj]
fit.up.unadj<-pred.unadj$fit[order.unadj] + qnorm(0.975)*se.unadj
fit.low.unadj<-pred.unadj$fit[order.unadj] + qnorm(0.025)*se.unadj

pdf("simulations/plots/appendix-figure21a-SIM-calibration-eta-known-UNADJ.pdf", width=9, height=7)
par(mar=c(5,4.5,2,2))

plot(obs.eta~fitted.unadj, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.unadj, rev(post.unadj)), c(expit(fit.low.unadj), rev(expit(fit.up.unadj))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.unadj)~post.unadj, col="dark blue", lwd=3)
dev.off()



##IOP-BX
fitted.bx<-pred.bx.ek

ns.bx<-glm(obs.eta~ns(fitted.bx,2), family="binomial")
pred.bx<-predict.glm(ns.bx, ns(fitted.bx,2), se.fit=TRUE)

order.bx<-order(fitted.bx)
post.bx<-fitted.bx[order.bx] #these are my posterior predictions
fit.bx<-pred.bx$fit[order.bx] #this is ns model fit
se.bx<-pred.bx$se.fit[order.bx]
fit.up.bx<-pred.bx$fit[order.bx] + qnorm(0.975)*se.bx
fit.low.bx<-pred.bx$fit[order.bx] + qnorm(0.025)*se.bx

pdf("simulations/plots/appendix-figure21b-SIM-calibration-eta-known-IOP-BX.pdf", width=9, height=7)
par(mar=c(5,4.5,2,2))

plot(obs.eta~fitted.bx, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.bx, rev(post.bx)), c(expit(fit.low.bx), rev(expit(fit.up.bx))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.bx)~post.bx, col="dark blue", lwd=3)
dev.off()





##IOP-SURG
fitted.surg<-pred.surg.ek

ns.surg<-glm(obs.eta~ns(fitted.surg,2), family="binomial")
pred.surg<-predict.glm(ns.surg, ns(fitted.surg,2), se.fit=TRUE)

order.surg<-order(fitted.surg)
post.surg<-fitted.surg[order.surg] #these are my posterior predictions
fit.surg<-pred.surg$fit[order.surg] #this is ns model fit
se.surg<-pred.surg$se.fit[order.surg]
fit.up.surg<-pred.surg$fit[order.surg] + qnorm(0.975)*se.surg
fit.low.surg<-pred.surg$fit[order.surg] + qnorm(0.025)*se.surg

pdf("simulations/plots/appendix-figure21c-SIM-calibration-eta-known-IOP-SURG.pdf", width=9, height=7)
par(mar=c(5,4.5,2,2))

plot(obs.eta~fitted.surg, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.surg, rev(post.surg)), c(expit(fit.low.surg), rev(expit(fit.up.surg))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.surg)~post.surg, col="dark blue", lwd=3)
dev.off()




##IOP 
fitted.iop<-pred.iop.ek

ns.iop<-glm(obs.eta~ns(fitted.iop,2), family="binomial")
pred.iop<-predict.glm(ns.iop, ns(fitted.iop,2), se.fit=TRUE)

order.iop<-order(fitted.iop)
post.iop<-fitted.iop[order.iop] #these are my posterior predictions
fit.iop<-pred.iop$fit[order.iop] #this is ns model fit
se.iop<-pred.iop$se.fit[order.iop]
fit.up.iop<-pred.iop$fit[order.iop] + qnorm(0.975)*se.iop
fit.low.iop<-pred.iop$fit[order.iop] + qnorm(0.025)*se.iop

pdf("simulations/plots/appendix-figure21d-SIM-calibration-eta-known-IOP.pdf", width=9, height=7)
par(mar=c(5,4.5,2,2))

plot(obs.eta~fitted.iop, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.iop, rev(post.iop)), c(expit(fit.low.iop), rev(expit(fit.up.iop))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.iop)~post.iop, col="dark blue", lwd=3)
dev.off()



##logistic
fitted.log<-pred.log.ek

ns.log<-glm(obs.eta~ns(fitted.log,2), family="binomial")
pred.log<-predict.glm(ns.log, ns(fitted.log,2), se.fit=TRUE)

order.log<-order(fitted.log)
post.log<-fitted.log[order.log] #these are my posterior predictions
fit.log<-pred.log$fit[order.log] #this is ns model fit
se.log<-pred.log$se.fit[order.log]
fit.up.log<-pred.log$fit[order.log] + qnorm(0.975)*se.log
fit.low.log<-pred.log$fit[order.log] + qnorm(0.025)*se.log

pdf("simulations/plots/appendix-figure21e-SIM-calibration-eta-known-LOGISTIC.pdf", width=9, height=7)
par(mar=c(5,4.5,2,2))

plot(obs.eta~fitted.log, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.log, rev(post.log)), c(expit(fit.low.log), rev(expit(fit.up.log))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.log)~post.log, col="dark blue", lwd=3)
dev.off()


