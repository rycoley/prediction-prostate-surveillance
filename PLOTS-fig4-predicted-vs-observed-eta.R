
### Rebecca Yates Coley rycoley@gmail.com
### Code for "Bayesian Joint Hierarchical Model for Prediction of Latent Health States with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in results and make plots for CROSS-VALIDATION on the latent state (Figure 4 in Coley et al.)



### LOAD PACKAGES
library("ROCR")

### DEFINE FUNCTIONS
expit<-function(x){return(exp(x)/(1+exp(x)))}


### LOAD DATA
pt.data<-read.csv("simulation-data/pt-data-sim.csv")

eta.data<-pt.data$obs.eta
obs.eta<-eta.data[!is.na(eta.data)]
(N<-length(obs.eta))


### GET RESULTS- UNADJUSTED ANALYSIS

id.test<-c(1:N)

fitted.eta<-vector(length=N)

for(i in 1:N){
	fitted.eta[i]<-read.csv(paste("results/eta-fitted-unadj-",id.test[i],".csv",sep=""))[2]}

fitted.unadj<-unlist(fitted.eta)

(auc.unadj<-performance(prediction(fitted.unadj,obs.eta),"auc")@y.values[[1]]) #0.85
roc.unadj<-performance(prediction(fitted.unadj,obs.eta),"tpr","fpr")

### GET RESULTS- IOP ANALYSIS

fitted.eta<-vector(length=N)

for(i in 1:N){
	fitted.eta[i]<-read.csv(paste("results/eta-fitted-iop-",id.test[i],".csv",sep=""))[2]}

fitted.iop<-unlist(fitted.eta)


(auc.iop<-performance(prediction(fitted.iop,obs.eta),"auc")@y.values[[1]]) #0.87
roc.iop<-performance(prediction(fitted.iop,obs.eta),"tpr","fpr")




### ROC PLOT
pdf("plots/figure4-left-ROC.pdf", height=7, width=9)

par(mar=c(5,4.5,2,2))
plot(c(0,1),c(0,1), type="n", xlab="", ylab="", xaxt="n", yaxt="n")
mtext("False Positive Rate",1,line=3, cex=1.5)
mtext("True Positive Rate",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)
lines(roc.unadj@y.values[[1]]~roc.unadj@x.values[[1]], col="red", lwd=3, lty="dotted")
lines(roc.iop@y.values[[1]]~roc.iop@x.values[[1]], col="blue", lwd=3)
legend("bottomright", legend=c("Unadjusted", "I.O.P."), col=c("red","blue"), lwd=rep(2), cex=1.5, bty="n", lty=c("dotted","solid"))

text(0.1,auc.iop, paste("AUC=", round(auc.iop,2), sep=""), col="blue", cex=1.5)
text(0.3,auc.unadj-0.2, paste("AUC=", round(auc.unadj,2), sep=""), col="red", cex=1.5)

dev.off()

### CALIBRATION PLOT

ns.iop<-glm(obs.eta~ns(fitted.iop,2), family="binomial")
pred.iop<-predict.glm(ns.iop, ns(fitted.iop,2), se.fit=TRUE)

order.iop<-order(fitted.iop)
post.iop<-fitted.iop[order.iop] #these are my posterior predictions
fit.iop<-pred.iop$fit[order.iop] #this is ns model fit
se.iop<-pred.iop$se.fit[order.iop]
fit.up.iop<-pred.iop$fit[order.iop] + qnorm(0.975)*se.iop
fit.low.iop<-pred.iop$fit[order.iop] + qnorm(0.025)*se.iop


pdf("plots/figure4-right-calibration-eta.pdf", width=9, height=7)
par(mar=c(5,5.5,2,5))

plot(obs.eta~fitted.iop, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.iop, rev(post.iop)), c(expit(fit.low.iop), rev(expit(fit.up.iop))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.iop)~post.iop, col="blue", lwd=3)

dev.off()

