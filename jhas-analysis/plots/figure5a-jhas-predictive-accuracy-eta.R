### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
#This code calculates bootstrap intervals for predictive accuracy for Figure 5a and appendix table 9
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.


library("pROC")



#true values
pt.data<-read.csv("jhas-analysis/simulation-data/pt-data-sim.csv")
eta.data<-pt.data$obs_eta
obs.eta<-eta.data[!is.na(eta.data)]
(N<-length(obs.eta))
#obs.eta  #all 0s than all 1s

#binary classifier
rc<-pt.data$rc[!is.na(pt.data$obs_eta)]

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



##get inf-obs results
fitted.eta<-vector(length=N)
for(i in 1:N){
	fitted.eta[i]<-read.csv(paste("jhas-analysis/IOP/results/eta-fitted-iop-",id.test[i],".csv",sep=""))[2]}
fitted.iop<-unlist(fitted.eta)


##get logistic results
res<-read.csv("jhas-analysis/LOGISTIC/results/pred-eta-logistic-ek.csv")
fitted.log<-res$pred.eta






B<- 5000 #10000 #10k used in analysis of JHAS data

auc.unadj<-auc.bx<-auc.surg<-auc.iop<-auc.log<-vector(length=B)
mse.unadj<-mse.bx<-mse.surg<-mse.iop<-mse.log<-mse.rc<-vector(length=B)
fpr.unadj<-fpr.bx<-fpr.surg<-fpr.iop<-fpr.log<-fpr.rc<-vector(length=B)
#b<-1

for(b in 1:B){

set.seed(b)
boot.samp<-sample(c(1:N), N, replace=TRUE)

roc.rc<-performance(prediction(rc[boot.samp], obs.eta[boot.samp]), "tpr", "fpr")
sens<-roc.rc@y.values[[1]][2]
fpr.rc[b]<-roc.rc@x.values[[1]][2]
mse.rc[b]<-mean((rc[boot.samp]-obs.eta[boot.samp])^2)

auc.unadj[b]<-performance(prediction(fitted.unadj[boot.samp], obs.eta[boot.samp]), "auc")@y.values[[1]]
roc.unadj<-performance(prediction(fitted.unadj[boot.samp], obs.eta[boot.samp]), "tpr", "fpr")
if(sum(roc.unadj@y.values[[1]]==sens)>0){
	fpr.unadj[b]<-min(roc.unadj@x.values[[1]][roc.unadj@y.values[[1]]==sens])}
if(sum(roc.unadj@y.values[[1]]==sens)==0){
	use.sens <- roc.unadj@y.values[[1]] [ abs(roc.unadj@y.values[[1]]-sens) == min(abs(roc.unadj@y.values[[1]]-sens))]
	fpr.unadj[b] <- min(roc.unadj@x.values[[1]][roc.unadj@y.values[[1]] %in% use.sens])}
mse.unadj[b]<-mean((fitted.unadj[boot.samp]-obs.eta[boot.samp])^2)


auc.bx[b]<-performance(prediction(fitted.bx[boot.samp], obs.eta[boot.samp]), "auc")@y.values[[1]]
roc.bx<-performance(prediction(fitted.bx[boot.samp], obs.eta[boot.samp]), "tpr", "fpr")
if(sum(roc.bx@y.values[[1]]==sens)>0){
	fpr.bx[b]<-min(roc.bx@x.values[[1]][roc.bx@y.values[[1]]==sens])}
if(sum(roc.bx@y.values[[1]]==sens)==0){
	use.sens <- roc.bx@y.values[[1]] [ abs(roc.bx@y.values[[1]]-sens) == min(abs(roc.bx@y.values[[1]]-sens))]
	fpr.bx[b] <- min(roc.bx@x.values[[1]][roc.bx@y.values[[1]] %in% use.sens])}
mse.bx[b]<-mean((fitted.bx[boot.samp]-obs.eta[boot.samp])^2)



auc.surg[b]<-performance(prediction(fitted.surg[boot.samp], obs.eta[boot.samp]), "auc")@y.values[[1]]
roc.surg<-performance(prediction(fitted.surg[boot.samp], obs.eta[boot.samp]), "tpr", "fpr")
if(sum(roc.surg@y.values[[1]]==sens)>0){
	fpr.surg[b]<-min(roc.surg@x.values[[1]][roc.surg@y.values[[1]]==sens])}
if(sum(roc.surg@y.values[[1]]==sens)==0){
	use.sens <- roc.surg@y.values[[1]] [ abs(roc.surg@y.values[[1]]-sens) == min(abs(roc.surg@y.values[[1]]-sens))]
	fpr.surg[b] <- min(roc.surg@x.values[[1]][roc.surg@y.values[[1]] %in% use.sens])}
mse.surg[b]<-mean((fitted.surg[boot.samp]-obs.eta[boot.samp])^2)



auc.iop[b]<-performance(prediction(fitted.iop[boot.samp], obs.eta[boot.samp]), "auc")@y.values[[1]]
roc.iop<-performance(prediction(fitted.iop[boot.samp], obs.eta[boot.samp]), "tpr", "fpr")
if(sum(roc.iop@y.values[[1]]==sens)>0){
	fpr.iop[b]<-min(roc.iop@x.values[[1]][roc.iop@y.values[[1]]==sens])}
if(sum(roc.iop@y.values[[1]]==sens)==0){
	use.sens <- roc.iop@y.values[[1]] [ abs(roc.iop@y.values[[1]]-sens) == min(abs(roc.iop@y.values[[1]]-sens))]
	fpr.iop[b] <- min(roc.iop@x.values[[1]][roc.iop@y.values[[1]] %in% use.sens])}
mse.iop[b]<-mean((fitted.iop[boot.samp]-obs.eta[boot.samp])^2)



auc.log[b]<-performance(prediction(fitted.log[boot.samp], obs.eta[boot.samp]), "auc")@y.values[[1]]
roc.log<-performance(prediction(fitted.log[boot.samp], obs.eta[boot.samp]), "tpr", "fpr")
if(sum(roc.log@y.values[[1]]==sens)>0){
	fpr.log[b]<-min(roc.log@x.values[[1]][roc.log@y.values[[1]]==sens])}
if(sum(roc.log@y.values[[1]]==sens)==0){
	use.sens <- roc.log@y.values[[1]] [ abs(roc.log@y.values[[1]]-sens) == min(abs(roc.log@y.values[[1]]-sens))]
	fpr.log[b] <- min(roc.log@x.values[[1]][roc.log@y.values[[1]] %in% use.sens])}
mse.log[b]<-mean((fitted.log[boot.samp]-obs.eta[boot.samp])^2)

if(b %in% seq(100,B,100)){print(b)}

}


##These results are reported in Table 9
quantile(auc.unadj, p=c(0.025, 0.975))
quantile(auc.bx, p=c(0.025, 0.975))
quantile(auc.surg, p=c(0.025, 0.975))
quantile(auc.iop, p=c(0.025, 0.975))
quantile(auc.log, p=c(0.025, 0.975))



quantile(mse.unadj, p=c(0.025, 0.975))
quantile(mse.bx, p=c(0.025, 0.975))
quantile(mse.surg, p=c(0.025, 0.975))
quantile(mse.iop, p=c(0.025, 0.975))
quantile(mse.log, p=c(0.025, 0.975))
quantile(mse.rc, p=c(0.025, 0.975))


quantile(fpr.unadj, p=c(0.025, 0.975))
quantile(fpr.bx, p=c(0.025, 0.975))
quantile(fpr.surg, p=c(0.025, 0.975))
quantile(fpr.iop, p=c(0.025, 0.975))
quantile(fpr.log, p=c(0.025, 0.975))
quantile(fpr.rc, p=c(0.025, 0.975))

sum(fpr.iop<fpr.rc)/B

quantile( (fpr.rc - fpr.iop)/ fpr.rc, p=c(0.025, 0.975)) 

quantile(fpr.rc - fpr.iop, p=c(0.025, 0.975))




##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.

###

pdf("jhas-analysis/plots/figure5a-JHAS-pred-accuracy-eta.pdf", width=7, height=2.5)
#layout(matrix(c(1,2), nrow=1, ncol=2), width=c(4.15,2.85), height=3)

#par(mfrow=c(1,2))

par(mar=c(3,8,2,1), las=1)
plot(0, type="n", xlim=c(0.5,1), ylim=c(0.5,5.5), xlab="", ylab="", yaxt="n", main="JHAS: True State Observed Post-Surgery")
mtext("AUC (95% Interval)", side=1, line=2 )
text(x=rep(0.48, 5), y=c(5:1), labels=c("Biopsy, Surgery IOP", "Biopsy IOP only", "Surgery IOP only", "Unadjusted", "Logistic"), srt=0, pos=2, xpd=TRUE)

lines(x=c(0.5, mean(auc.iop)), y=c(5,5), lty="dotted")
points(x=mean(auc.iop), y=5, pch=19, cex=2)
lines(x=quantile(auc.iop, p=c(0.025, 0.975)), y=c(5,5), lty="solid", lwd=3)

lines(x=c(0.5, mean(auc.bx)), y=c(4,4), lty="dotted")
points(x=mean(auc.bx), y=4, pch=19, cex=2)
lines(x=quantile(auc.bx, p=c(0.025, 0.975)), y=c(4,4), lty="solid", lwd=3)

lines(x=c(0.5, mean(auc.surg)), y=c(3,3), lty="dotted")
points(x=mean(auc.surg), y=3, pch=19, cex=2)
lines(x=quantile(auc.surg, p=c(0.025, 0.975)), y=c(3,3), lty="solid", lwd=3)

lines(x=c(0.5, mean(auc.unadj)), y=c(2,2), lty="dotted")
points(x=mean(auc.unadj), y=2, pch=19, cex=2)
lines(x=quantile(auc.unadj, p=c(0.025, 0.975)), y=c(2,2), lty="solid", lwd=3)

lines(x=c(0.5, mean(auc.log)), y=c(1,1), lty="dotted")
points(x=mean(auc.log), y=1, pch=19, cex=2)
lines(x=quantile(auc.log, p=c(0.025, 0.975)), y=c(1,1), lty="solid", lwd=3)


dev.off()












