### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
#This code reproduces Figure 1 in the appendix
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.




pt.data<-read.csv("jhas-analysis/simulation-data/pt-data-sim.csv")
bx.data<-read.csv("jhas-analysis/simulation-data/bx-data-sim.csv")
psa.data<-read.csv("jhas-analysis/simulation-data/psa-data-sim.csv")


eta.data<-pt.data$obs_eta
(N<-length(eta.data))

obs.eta<-eta.data[!is.na(eta.data)]
(n_eta_known<-length(obs.eta))


pt.ek<-pt.data[!is.na(eta.data),]
pt.nek<-pt.data[is.na(eta.data),]


sum(pt.data$surg==0 & pt.data$rc==0) #672
sum(pt.data$surg==0 & pt.data$rc==1) #125
sum(pt.data$surg==1 & pt.data$rc==0) #123
sum(pt.data$surg==1 & pt.data$rc==1) #80



##unadjusted model
#eta known
id.test<-c(1:n_eta_known)
fitted.eta<-vector(length=n_eta_known)
for(i in 1:n_eta_known){
	fitted.eta[i]<-read.csv(paste("jhas-analysis/UNADJ/results/eta-fitted-unadj-",id.test[i],".csv",sep=""))[2]}
pred.ek.unadj<-unlist(fitted.eta)

#eta unknown
etahat<-read.csv(paste("jhas-analysis/UNADJ/results/jags-prediction-unadj-eta_hat-1.csv",sep=""))
dim(etahat)
etahat<-as.matrix(etahat[,2:798])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/UNADJ/results/jags-prediction-unadj-eta_hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:798])}
pred.nek.unadj<-apply(etahat,2,mean)



##IOP-BX model
#eta known
id.test<-c(1:n_eta_known)
fitted.eta<-vector(length=n_eta_known)
for(i in 1:n_eta_known){
	fitted.eta[i]<-read.csv(paste("jhas-analysis/IOP-BX/results/eta-fitted-IOP-BX-", id.test[i], ".csv", sep=""))[2]}
pred.ek.bx<-unlist(fitted.eta)

#eta unknown
etahat<-read.csv(paste("jhas-analysis/IOP-BX/results/jags-prediction-IOP-BX-eta_hat-1.csv",sep=""))
dim(etahat)
etahat<-as.matrix(etahat[,2:798])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP-BX/results/jags-prediction-IOP-BX-eta_hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:798])}
pred.nek.bx<-apply(etahat,2,mean)


##IOP-SURG model
#eta known
id.test<-c(1:n_eta_known)
fitted.eta<-vector(length=n_eta_known)
for(i in 1:n_eta_known){
	fitted.eta[i]<-read.csv(paste("jhas-analysis/IOP-SURG/results/eta-fitted-IOP-SURG-", id.test[i], ".csv", sep=""))[2]}
pred.ek.surg<-unlist(fitted.eta)
summary(pred.ek.surg)

#eta unknown
etahat<-read.csv(paste("jhas-analysis/IOP-SURG/results/jags-prediction-IOP-SURG-eta_hat-1.csv",sep=""))
dim(etahat)
etahat<-as.matrix(etahat[,2:798])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP-SURG/results/jags-prediction-IOP-SURG-eta_hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:798])}
pred.nek.surg<-apply(etahat,2,mean)
summary(pred.nek.surg)


##IOP model
#eta known
id.test<-c(1:n_eta_known)
fitted.eta<-vector(length=n_eta_known)
for(i in 1:n_eta_known){
	fitted.eta[i]<-read.csv(paste("jhas-analysis/IOP/results/eta-fitted-iop-", id.test[i], ".csv", sep=""))[2]}
pred.ek.iop<-unlist(fitted.eta)
summary(pred.ek.iop)

#eta unknown
etahat<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-eta_hat-1.csv",sep=""))
dim(etahat)
etahat<-as.matrix(etahat[,2:798])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-eta_hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:798])}
pred.nek.iop<-apply(etahat,2,mean)
summary(pred.nek.iop)




#logistic regression
res<-read.csv("jhas-analysis/LOGISTIC/results/pred-eta-logistic-ek.csv")
pred.ek.log<-res$pred.eta
summary(pred.ek.log)

res<-read.csv("jhas-analysis/LOGISTIC/results/pred-eta-logistic-nek.csv")
pred.nek.log<-res$pred.eta
summary(pred.nek.log)

rm(res)





##vertical line at accuracy of biopsies

sum(pt.data$surg==1 & pt.data$obs_eta==1 & pt.data$rc==0)/sum(pt.data$surg==1& pt.data$rc==0) #0.33

sum(pt.data$surg==1 & pt.data$obs_eta==1 & pt.data$rc==1)/sum(pt.data$surg==1& pt.data$rc==1) #0.9






d.nek.un.0<-density(pred.nek.unadj[pt.nek$rc==0])
d.nek.bx.0<-density(pred.nek.bx[pt.nek$rc==0])
d.nek.surg.0<-density(pred.nek.surg[pt.nek$rc==0])
d.nek.iop.0<-density(pred.nek.iop[pt.nek$rc==0])
d.nek.log.0<-density(pred.nek.log[pt.nek$rc==0])

d.nek.un.1<-density(pred.nek.unadj[pt.nek$rc==1])
d.nek.bx.1<-density(pred.nek.bx[pt.nek$rc==1])
d.nek.surg.1<-density(pred.nek.surg[pt.nek$rc==1])
d.nek.iop.1<-density(pred.nek.iop[pt.nek$rc==1])
d.nek.log.1<-density(pred.nek.log[pt.nek$rc==1])

d.ek.un.0<-density(pred.ek.unadj[pt.ek$rc==0])
d.ek.bx.0<-density(pred.ek.bx[pt.ek$rc==0])
d.ek.surg.0<-density(pred.ek.surg[pt.ek$rc==0])
d.ek.iop.0<-density(pred.ek.iop[pt.ek$rc==0])
d.ek.log.0<-density(pred.ek.log[pt.ek$rc==0])

d.ek.un.1<-density(pred.ek.unadj[pt.ek$rc==1])
d.ek.bx.1<-density(pred.ek.bx[pt.ek$rc==1])
d.ek.surg.1<-density(pred.ek.surg[pt.ek$rc==1])
d.ek.iop.1<-density(pred.ek.iop[pt.ek$rc==1])
d.ek.log.1<-density(pred.ek.log[pt.ek$rc==1])





pdf("jhas-analysis/plots/appendix-figure1-jhasposterior-pred-eta.pdf", width=7, height=5)

#par(mfrow=c(2,2), mar=c(2,2,2,2))
layout(mat=matrix(c(1,2,3,4,5,5), byrow=T, nrow=3, ncol=2), width=c(3.5, 3.5), height=c(3,3,1))


par(mar=c(4,3,2,2))
plot(0, type="n", xlim=c(0,1), ylim=c(0,10), ylab="", xlab="")
lines(d.nek.un.0$y~d.nek.un.0$x, col="gray50", lty="dotdash", lwd=2)
lines(d.nek.bx.0$y~d.nek.bx.0$x, col="red", lty="dotted", lwd=2)
lines(d.nek.surg.0$y~d.nek.surg.0$x, col="green3", lty="dashed", lwd=2)
lines(d.nek.iop.0$y~d.nek.iop.0$x, col="blue", lty="solid", lwd=2)
lines(d.nek.log.0$y~d.nek.log.0$x, col="black", lty="longdash", lwd=2)
abline(v=0.33, lwd=1, lty="dotted")


par(mar=c(4,3,2,2))
plot(0, type="n", xlim=c(0,1), ylim=c(0,10), ylab="", xlab="")
lines(d.nek.un.1$y~d.nek.un.1$x, col="gray50", lty="dotdash", lwd=2)
lines(d.nek.bx.1$y~d.nek.bx.1$x, col="red", lty="dotted", lwd=2)
lines(d.nek.surg.1$y~d.nek.surg.1$x, col="green3", lty="dashed", lwd=2)
lines(d.nek.iop.1$y~d.nek.iop.1$x, col="blue", lty="solid", lwd=2)
lines(d.nek.log.1$y~d.nek.log.1$x, col="black", lty="longdash", lwd=2)
abline(v=0.9, lwd=1, lty="dotted")


par(mar=c(4,3,2,2))
plot(0, type="n", xlim=c(0,1), ylim=c(0,10), ylab="", xlab="")
lines(d.ek.un.0$y~d.ek.un.0$x, col="gray50", lty="dotdash", lwd=2)
lines(d.ek.bx.0$y~d.ek.bx.0$x, col="red", lty="dotted", lwd=2)
lines(d.ek.surg.0$y~d.ek.surg.0$x, col="green3", lty="dashed", lwd=2)
lines(d.ek.iop.0$y~d.ek.iop.0$x, col="blue", lty="solid", lwd=2)
lines(d.ek.log.0$y~d.ek.log.0$x, col="black", lty="longdash", lwd=2)
abline(v=0.33, lwd=1, lty="dotted")



par(mar=c(4,3,2,2))
plot(0, type="n", xlim=c(0,1), ylim=c(0,10), ylab="", xlab="")
lines(d.ek.un.1$y~d.ek.un.1$x, col="gray50", lty="dotdash", lwd=2)
lines(d.ek.bx.1$y~d.ek.bx.1$x, col="red", lty="dotted", lwd=2)
lines(d.ek.surg.1$y~d.ek.surg.1$x, col="green3", lty="dashed", lwd=2)
lines(d.ek.iop.1$y~d.ek.iop.1$x, col="blue", lty="solid", lwd=2)
lines(d.ek.log.1$y~d.ek.log.1$x, col="black", lty="longdash", lwd=2)
abline(v=0.9, lwd=1, lty="dotted")

par(mar=c(0,0,0,0))
plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n" )
legend("center", legend=c("No IOP", "Biopsy IOP", "Surgery IOP", "Biopsy, Surgery IOP", "Logistic"), lty=c("dotdash", "dotted", "dashed", "solid", "longdash"), col=c("gray50", "red", "green3", "blue", "black"), lwd=rep(2,5), bty="n", ncol=3)


dev.off()



#Annotation can be added in keynote/powerpoint or adobe illustrator.

