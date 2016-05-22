### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
#This code reproduced Figure 19 in the appendix



SEED<-1
set.seed(SEED)

source("simulations/simulation-data/simulate-pca-data.R")

true.eta<-pt_data$eta_true

sum(pt_data$surg==0 & pt_data$rc==0) #672
sum(pt_data$surg==0 & pt_data$rc==1) #125
sum(pt_data$surg==1 & pt_data$rc==0) #123
sum(pt_data$surg==1 & pt_data$rc==1) #80



pred.ek<-read.csv("simulations/UNADJ/results/eta-known-true-vs-pred-unadj-1.csv")
pred.nek<-read.csv("simulations/UNADJ/results/eta-true-vs-pred-unadj-1.csv")
pred.unadj<-c(pred.ek$predicted, pred.nek$predicted)


pred.ek<-read.csv("simulations/IOP-BX/results/eta-known-true-vs-pred-IOP-BX-1.csv")
pred.nek<-read.csv("simulations/IOP-BX/results/eta-true-vs-pred-IOP-BX-1.csv")
pred.bx<-c(pred.ek$predicted, pred.nek$predicted)


pred.ek<-read.csv("simulations/IOP-SURG/results/eta-known-true-vs-pred-IOP-SURG-1.csv")
pred.nek<-read.csv("simulations/IOP-SURG/results/eta-true-vs-pred-IOP-SURG-1.csv")
pred.rrp<-c(pred.ek$predicted, pred.nek$predicted)


pred.ek<-read.csv("simulations/IOP/results/eta-known-true-vs-pred-iop-1.csv")
pred.nek<-read.csv("simulations/IOP/results/eta-true-vs-pred-iop-1.csv")
pred.iop<-c(pred.ek$predicted, pred.nek$predicted)

res<-read.csv("simulations/LOGISTIC/results/eta-predictions-LOGISTIC.csv")
pred.log<-res$predicted






d.nek.un.0<-density(pred.unadj[pt_data$surg==0 & pt_data$rc==0])
d.nek.bx.0<-density(pred.bx[pt_data$surg==0 & pt_data$rc==0])
d.nek.rrp.0<-density(pred.rrp[pt_data$surg==0 & pt_data$rc==0])
d.nek.iop.0<-density(pred.iop[pt_data$surg==0 & pt_data$rc==0])
d.nek.log.0<-density(pred.log[pt_data$surg==0 & pt_data$rc==0])

d.nek.un.1<-density(pred.unadj[pt_data$surg==0 & pt_data$rc==1])
d.nek.bx.1<-density(pred.bx[pt_data$surg==0 & pt_data$rc==1])
d.nek.rrp.1<-density(pred.rrp[pt_data$surg==0 & pt_data$rc==1])
d.nek.iop.1<-density(pred.iop[pt_data$surg==0 & pt_data$rc==1])
d.nek.log.1<-density(pred.log[pt_data$surg==0 & pt_data$rc==1])

d.ek.un.0<-density(pred.unadj[pt_data$surg==1 & pt_data$rc==0])
d.ek.bx.0<-density(pred.bx[pt_data$surg==1 & pt_data$rc==0])
d.ek.rrp.0<-density(pred.rrp[pt_data$surg==1 & pt_data$rc==0])
d.ek.iop.0<-density(pred.iop[pt_data$surg==1 & pt_data$rc==0])
d.ek.log.0<-density(pred.log[pt_data$surg==1 & pt_data$rc==0])

d.ek.un.1<-density(pred.unadj[pt_data$surg==1 & pt_data$rc==1])
d.ek.bx.1<-density(pred.bx[pt_data$surg==1 & pt_data$rc==1])
d.ek.rrp.1<-density(pred.rrp[pt_data$surg==1 & pt_data$rc==1])
d.ek.iop.1<-density(pred.iop[pt_data$surg==1 & pt_data$rc==1])
d.ek.log.1<-density(pred.log[pt_data$surg==1 & pt_data$rc==1])







pdf("simulations/plots/appendix-figure19-SIM-posterior-pred-eta.pdf", width=7, height=5)

layout(mat=matrix(c(1,2,3,4,5,5), byrow=T, nrow=3, ncol=2), width=c(3.5, 3.5), height=c(3,3,1))

par(mar=c(4,3,2,2))
plot(0, type="n", xlim=c(0,1), ylim=c(0,8.7), ylab="", xlab="")
lines(d.nek.un.0$y~d.nek.un.0$x, col="gray50", lty="dotdash", lwd=2)
lines(d.nek.bx.0$y~d.nek.bx.0$x, col="red", lty="dotted", lwd=2)
lines(d.nek.rrp.0$y~d.nek.rrp.0$x, col="green3", lty="dashed", lwd=2)
lines(d.nek.iop.0$y~d.nek.iop.0$x, col="blue", lty="solid", lwd=2)
lines(d.nek.log.0$y~d.nek.log.0$x, col="black", lty="longdash", lwd=2)
abline(v=0.33, lwd=1, lty="dotted")

par(mar=c(4,3,2,2))
plot(0, type="n", xlim=c(0,1), ylim=c(0,8.7), ylab="", xlab="")
lines(d.nek.un.1$y~d.nek.un.1$x, col="gray50", lty="dotdash", lwd=2)
lines(d.nek.bx.1$y~d.nek.bx.1$x, col="red", lty="dotted", lwd=2)
lines(d.nek.rrp.1$y~d.nek.rrp.1$x, col="green3", lty="dashed", lwd=2)
lines(d.nek.iop.1$y~d.nek.iop.1$x, col="blue", lty="solid", lwd=2)
lines(d.nek.log.1$y~d.nek.log.1$x, col="black", lty="longdash", lwd=2)
abline(v=0.9, lwd=1, lty="dotted")


par(mar=c(4,3,2,2))
plot(0, type="n", xlim=c(0,1), ylim=c(0,8.7), ylab="", xlab="")
lines(d.ek.un.0$y~d.ek.un.0$x, col="gray50", lty="dotdash", lwd=2)
lines(d.ek.bx.0$y~d.ek.bx.0$x, col="red", lty="dotted", lwd=2)
lines(d.ek.rrp.0$y~d.ek.rrp.0$x, col="green3", lty="dashed", lwd=2)
lines(d.ek.iop.0$y~d.ek.iop.0$x, col="blue", lty="solid", lwd=2)
lines(d.ek.log.0$y~d.ek.log.0$x, col="black", lty="longdash", lwd=2)
abline(v=0.33, lwd=1, lty="dotted")


par(mar=c(4,3,2,2))
plot(0, type="n", xlim=c(0,1), ylim=c(0,8.7), ylab="", xlab="")
lines(d.ek.un.1$y~d.ek.un.1$x, col="gray50", lty="dotdash", lwd=2)
lines(d.ek.bx.1$y~d.ek.bx.1$x, col="red", lty="dotted", lwd=2)
lines(d.ek.rrp.1$y~d.ek.rrp.1$x, col="green3", lty="dashed", lwd=2)
lines(d.ek.iop.1$y~d.ek.iop.1$x, col="blue", lty="solid", lwd=2)
lines(d.ek.log.1$y~d.ek.log.1$x, col="black", lty="longdash", lwd=2)
abline(v=0.9, lwd=1, lty="dotted")


par(mar=c(0,0,0,0))
plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n" )
legend("center", legend=c("No IOP", "Biopsy IOP", "Surgery IOP", "Biopsy, Surgery IOP", "Logistic"), lty=c("dotdash", "dotted", "dashed", "solid", "longdash"), col=c("gray50", "red", "green3", "blue", "black"), lwd=rep(2,5), bty="n", ncol=3)

dev.off()

#Annotation can be added in keynote/powerpoint or adobe illustrator