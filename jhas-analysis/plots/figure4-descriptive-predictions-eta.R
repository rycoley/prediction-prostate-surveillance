### NAME AND EMAIL REDACTED
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in results and make descriptive plots for out-of-sample predictions of the latent state (Figure 3 in Coley et al.)
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.



### LOAD PACKAGES
list.of.packages <- c("pROC","ROCR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=T)

library("pROC")
library("ROCR")





### LOAD DATA
pt_data<-read.csv("jhas-analysis/simulation-data/pt-data-sim.csv")
bx_data<-read.csv("jhas-analysis/simulation-data/bx-data-sim.csv")
psa_data<-read.csv("jhas-analysis/simulation-data/psa-data-sim.csv")


eta_data<-pt_data$obs_eta
(N<-length(eta_data))

obs_eta<-eta_data[!is.na(eta_data)]
(n_eta_known<-length(obs_eta))






### GET PREDICTIONS- UNADJUSTED MODEL

id_test<-c(1:n_eta_known)
fitted_eta<-vector(length=n_eta_known)
for(i in 1:n_eta_known){
	fitted_eta[i]<-read.csv(paste("jhas-analysis/UNADJ/results/eta-fitted-unadj-",id_test[i],".csv",sep=""))[2]}
fitted_unadj<-unlist(fitted_eta)

etahat<-read.csv(paste("jhas-analysis/UNADJ/results/jags-prediction-unadj-eta_hat-1.csv",sep=""))
dim(etahat)
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/UNADJ/results/jags-prediction-unadj-eta_hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:dim(res)[2]])}
	
post_eta_unadj<-c(fitted_unadj, apply(etahat,2,mean))





### GET PREDICTIONS- IOP MODEL
id_test<-c(1:n_eta_known)
fitted_eta<-vector(length=n_eta_known)
for(i in 1:n_eta_known){
	fitted_eta[i]<-read.csv(paste("jhas-analysis/IOP/results/eta-fitted-iop-",id_test[i],".csv",sep=""))[2]}
fitted_iop<-unlist(fitted_eta)

etahat<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-eta_hat-1.csv",sep=""))
dim(etahat)
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-eta_hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:dim(res)[2]])}
	
post_eta_iop<-c(fitted_iop, apply(etahat,2,mean))




## HISTOGRAM OF IOP PREDICTIONS

#Draw stacked histogram
#Groups with different colors and shading: eta unknown (without and with RC), eta=0 (w and wo RC), eta=1 (w and wo RC)
counts<-matrix(nrow=6, ncol=10)
breaks<-seq(0,1,0.1)
shading<-rep(c(10,NA),3)
coloring<-c("black","black", rep("green3",2), rep("red",2))

for(j in 1:10){
	counts[1,j]<-sum(is.na(pt_data$obs_eta) & pt_data$rc==0 & post_eta_iop>=breaks[j] & post_eta_iop<breaks[(j+1)]) 
	counts[2,j]<-sum(is.na(pt_data$obs_eta) & pt_data$rc==1 & post_eta_iop>=breaks[j] & post_eta_iop<breaks[(j+1)])
	counts[3,j]<-sum(!is.na(pt_data$obs_eta) & pt_data$obs_eta==0 & pt_data$rc==0 & post_eta_iop>=breaks[j] & post_eta_iop<breaks[(j+1)])
	counts[4,j]<-sum(!is.na(pt_data$obs_eta) & pt_data$obs_eta==0 & pt_data$rc==1 & post_eta_iop>=breaks[j] & post_eta_iop<breaks[(j+1)])
	counts[5,j]<-sum(!is.na(pt_data$obs_eta) & pt_data$obs_eta==1 & pt_data$rc==0 & post_eta_iop>=breaks[j] & post_eta_iop<breaks[(j+1)])
	counts[6,j]<-sum(!is.na(pt_data$obs_eta) & pt_data$obs_eta==1 & pt_data$rc==1 & post_eta_iop>=breaks[j] & post_eta_iop<breaks[(j+1)])
	}


pdf("jhas-analysis/plots/figure4a-histogram.pdf", width=9, height=7)

plot(1, type="n", xlim=c(0,1), ylim=c(0,325), frame.plot=F, xlab="", ylab="", xaxt="n", yaxt="n")
mtext("Posterior P(Aggressive Prostate Cancer)",1,line=3, cex=1.5)
mtext("Frequency",2,line=2.75, cex=1.5)
Axis(side=1, at=seq(0,1,0.2), labels=seq(0,1,0.2), cex=1.25)
Axis(side=2, at=seq(0,300,50), labels=seq(0,300,50), cex=1.25)

for(j in 1:10){
polygon(x=c(rep(breaks[j],2), rep(breaks[(j+1)],2)), y=c(0,counts[1,j], counts[1,j], 0), col=coloring[1], density=shading[1] )
	for(k in 2:6){
polygon(x=c(rep(breaks[j],2), rep(breaks[(j+1)],2)), y=c(sum(counts[1:(k-1),j]), sum(counts[1:k,j]), sum(counts[1:k,j]), sum(counts[1:(k-1),j])), col=coloring[k], density=shading[k] ) } }

legend(x=0.65, y=340, title="Observed Cancer State", legend=c("Unobserved", "Indolent", "Aggressive"), pt.bg=c("black", "green3", "red"), pch=rep(22,3), bty="n", cex=1.5 )
legend(x=0.375, y=340, title="Final Biopsy", legend=c("Gleason 6", "Gleason 7+"), density=c(10, NA), col=c("black", "black"), bty="n", cex=1.5)

dev.off()




## SCATTERPLOT OF IOP VS UNADJ PREDICTIONS

color_eta<-c(rep("green3", sum(obs_eta==0)), rep("red", sum(obs_eta==1)), rep("black", N-n_eta_known))

num_bx<-yr_bx<-vector(length=N)
for(i in 1:N){num_bx[i]<-sum(bx_data$subj==i & bx_data$bx_here==1 & !is.na(bx_data$bx_here))
	yr_bx[i]<-sum(!is.na(bx_data$bx_here) & bx_data$subj==i)}
summary(num_bx)
summary(num_bx/yr_bx)


pdf("jhas-analysis/plots/figure4b-scatter-eta-iop-vs-unadj.pdf", width=9, height=7)
par(mar=c(5,5.5,2,5))
plot(post_eta_iop~post_eta_unadj, col=color_eta, cex=(num_bx/yr_bx)*2, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n")
mtext("Non-IOP Model P(Aggressive Prostate Cancer)",1,line=3, cex=1.5)
mtext("IOP Model P(Aggressive Prostate Cancer)",2,line=2.75, cex=1.5)
#add title in post-processing
Axis(side=1, at=seq(0,1,0.2), labels=seq(0,1,0.2), cex=1.25)
Axis(side=2, at=seq(0,1,0.2), labels=seq(0,1,0.2), cex=1.25)
points(post_eta_iop[1:n_eta_known]~post_eta_unadj[1:n_eta_known], col=color_eta[1:n_eta_known], cex=(num_bx[1:n_eta_known]/yr_bx[1:n_eta_known])*2)
legend("topleft", title="Biopsy:Year", legend=c("1:1", "1:2"),  pch=c(21,21), pt.cex=c(1,0.5)*1.5, cex=1.5)
dev.off()

