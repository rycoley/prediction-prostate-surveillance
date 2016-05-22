### NAME AND EMAIL REDACTED
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in results and make plots for out-of-sample predictive accuracy of the latent state (Figure 4 in Coley et al.)
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.



### LOAD PACKAGES
list.of.packages <- c("pROC","ROCR","splines")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=T)

library("pROC")
library("ROCR")
library("splines")



### DEFINE FUNCTIONS
expit<-function(x){return(exp(x)/(1+exp(x)))}


### LOAD DATA
pt_data<-read.csv("jhas-analysis/simulation-data/pt-data-sim.csv")

eta_data<-pt_data$obs_eta
obs_eta<-eta_data[!is.na(eta_data)]
(N<-length(obs_eta))


### GET RESULTS- UNADJUSTED ANALYSIS

id_test<-c(1:N)

fitted_eta<-vector(length=N)

for(i in 1:N){
	fitted_eta[i]<-read.csv(paste("jhas-analysis/UNADJ/results/eta-fitted-unadj-",id_test[i],".csv",sep=""))[2]}

fitted_unadj<-unlist(fitted_eta)

(auc.unadj<-performance(prediction(fitted_unadj,obs_eta),"auc")@y.values[[1]]) #0.81
roc.unadj<-performance(prediction(fitted_unadj,obs_eta),"tpr","fpr")

str(roc.unadj) #x is FPR, y is TPR
cbind(roc.unadj@x.values[[1]][roc.unadj@y.values[[1]]>0.60 & roc.unadj@y.values[[1]]<0.63], roc.unadj@y.values[[1]][roc.unadj@y.values[[1]]>0.60 & roc.unadj@y.values[[1]]<0.63])
#FPR=0.089 
1-0.089 #specificity=0.911 #This is much better than in the actual application because we generated data under the proposed model

ci.unadj<-ci.auc(response=obs_eta, predictor=fitted_unadj, method="bootstrap") 
(ci.unadj<-c(ci.unadj[1], ci.unadj[3])) #(0.75, 0.86)


### GET RESULTS- IOP ANALYSIS
fitted_eta<-vector(length=N)

for(i in 1:N){
	fitted_eta[i]<-read.csv(paste("jhas-analysis/IOP/results/eta-fitted-iop-",id_test[i],".csv",sep=""))[2]}


fitted_iop<-unlist(fitted_eta)

(auc.iop<-performance(prediction(fitted_iop,obs_eta),"auc")@y.values[[1]]) #0.79
roc.iop<-performance(prediction(fitted_iop,obs_eta),"tpr","fpr")


cbind(roc.iop@x.values[[1]][roc.iop@y.values[[1]]>0.60 & roc.iop@y.values[[1]]<0.63], roc.iop@y.values[[1]][roc.iop@y.values[[1]]>0.60 & roc.iop@y.values[[1]]<0.63])
#FPR=0.0889 at TPR=0.615
#specificity=0.911 #Again, much better than in the actual application because we generated data under the proposed model

ci.io<-ci.auc(response=obs_eta, predictor=fitted_iop, method="bootstrap")
(ci.io<-c(ci.io[1], ci.io[3])) #0.73, 0.85

#Also, not surprosing that IOP model does not dominate non-IOP. Since non-IOP assumes that those with surgery are representative of all patients, we would expect to see improvement in accuracy for those who didn't get surgery



### ROC PLOT
pdf("jhas-analysis/plots/figure6a-ROC.pdf", height=7, width=9)

par(mar=c(5,4.5,2,2))
plot(c(0,1),c(0,1), type="n", xlab="", ylab="", xaxt="n", yaxt="n")
mtext("False Positive Rate",1,line=3, cex=1.5)
mtext("True Positive Rate",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

lines(y=c(0.615,0.615), x=c(-1, 0.0889), lwd=1, lty="dashed")
lines(y=c(-1, 0.615), x=rep(0.0889,2), lwd=1, lty="dashed")
lines(y=c(-1, 0.615), x=rep(0.0889,2), lwd=1, lty="dashed")


lines(roc.unadj@y.values[[1]]~roc.unadj@x.values[[1]], col="coral1", lwd=3, lty="dotted")
lines(roc.iop@y.values[[1]]~roc.iop@x.values[[1]], col="coral3", lwd=3)
#legend("bottomright", legend=c("Unadjusted", "I.O.P."), col=c("coral1","coral"), lwd=rep(2), cex=1.5, bty="n", lty=c("dotted","solid"))

points(y=0.615, x=(1-0.795), pch="*", cex=3)

#text(0.1,auc.iop, paste("AUC=", round(auc.iop,2), sep=""), col="blue", cex=1.5)
#text(0.3,auc.unadj-0.2, paste("AUC=", round(auc.unadj,2), sep=""), col="red", cex=1.5)

dev.off()

## Labels, AUC can be put on plot in Powerpoint/Keynote or Adobe Illustrator



### CALIBRATION PLOT

ns.iop<-glm(obs_eta~ns(fitted_iop,2), family="binomial")
pred.iop<-predict.glm(ns.iop, ns(fitted_iop,2), se.fit=TRUE)

order.iop<-order(fitted_iop)
post.iop<-fitted_iop[order.iop] #these are my posterior predictions
fit.iop<-pred.iop$fit[order.iop] #this is ns model fit
se.iop<-pred.iop$se.fit[order.iop]
fit.up.iop<-pred.iop$fit[order.iop] + qnorm(0.975)*se.iop
fit.low.iop<-pred.iop$fit[order.iop] + qnorm(0.025)*se.iop


pdf("plots/figure6b-calibration-eta.pdf", width=9, height=7)
par(mar=c(5,5.5,2,5))

plot(obs_eta~fitted_iop, pch=3, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1))
mtext("Posterior P(Aggressive PCa)",1,line=3, cex=1.5)
mtext("Observed P(Aggressive PCa)",2,line=2.75, cex=1.5)
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)

polygon(x=c(post.iop, rev(post.iop)), c(expit(fit.low.iop), rev(expit(fit.up.iop))) , col="lightblue", border="lightblue")

abline(0,1, lwd=3, lty="dotted")

lines(expit(fit.iop)~post.iop, col="blue", lwd=3)

dev.off()

##can write longer labels for x- and y-axis and add annotations on plot in Powerpoint/Keynote or Adobe Illustrator