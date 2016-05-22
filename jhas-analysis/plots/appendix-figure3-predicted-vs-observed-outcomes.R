### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in biopsy, reclassification, and surgery data for individuals, MCMC predictions of these outcomes, and make calibration plots in the appendix, figure 3
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.


### LOAD PACKAGES
list.of.packages <- c("splines")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=T)

library("splines")

### DEFINE FUNCTIONS
expit<-function(x){return(exp(x)/(1+exp(x)))}


### LOAD DATA
pt_data<-read.csv("jhas-analysis/simulation-data/pt-data-sim.csv")
data_use<-read.csv("jhas-analysis/simulation-data/bx-data-sim.csv")

(n0<-sum(pt_data$obs_eta==0,na.rm=T))
(n1<-sum(pt_data$obs_eta==1,na.rm=T))
eta_data<-pt_data$obs_eta
(n_eta_known<-sum(!is.na(eta_data))) 
n<-dim(pt_data)[1]

##remove observations (since this model was run with some observations removed in order to demonstrate individualized predictions)
set.seed(4)
pred_ids<-sample(c((n_eta_known+1):n), 12)
(pred_ids<-sort(pred_ids))
data_use<-data_use[!(data_use$subj%in%pred_ids & data_use$time>5),]


## biopsy data
bx_data<-data_use[!is.na(data_use$bx_here),]
(n_bx<-dim(bx_data)[1])
BX<-as.numeric(bx_data$bx_here)
subj_bx<-bx_data$subj

## reclassification data
rc_data<-data_use[data_use$bx_here==1 & !is.na(data_use$bx_here),]
(n_rc<-dim(rc_data)[1])
RC<-as.numeric(rc_data$rc)
subj_rc<-rc_data$subj

## surgery data
SURG<-as.numeric(data_use$surg)
(n_surg<-dim(data_use)[1])
subj_surg<-data_use$subj


### LOAD MCMC RESULTS

#We use "reduced" files here so that they are a manageable size to share, post on github, etc. Use RESULTS-reduce-size in IOP folder to obtain.

## biopsy predictions
p.bx<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-p_bx-1-reduced.csv",sep=""))
dim(p.bx)
p.bx<-as.matrix(p.bx[,2:(n_bx+1)])
p.bx.mean<-apply(p.bx,2,mean)

## reclassification predictions
p.rc<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-p_rc-1-reduced.csv",sep=""))
dim(p.rc)
p.rc<-as.matrix(p.rc[,2:(n_rc+1)])
p.rc.mean<-apply(p.rc,2,mean)

##sugery predictions
p.surg<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-p_surg-1-reduced.csv",sep=""))
dim(p.surg)
p.surg<-as.matrix(p.surg[,2:(n_surg+1)])
p.surg.mean<-apply(p.surg,2,mean)

B<-dim(p.surg)[1]



b.use<-seq(1,B,1)

pdf("jhas-analysis/plots/appendix-figure3-predicted-vs-observed-outcomes.pdf", width=12, height=6)
layout(matrix(c(1,2,3,4,4,4),2,3,byrow=TRUE), widths=rep(4,3), heights=c(5,1))


par(mar=c(2.5,2,3,2))

##BX
plot(BX~p.bx.mean, type="n",  ylab="", xlab="", xaxt="n", yaxt="n") 
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)
mtext("Biopsy Performed",3,line=1, cex=1.5)

for(b in 1:length(b.use)){
	p.bx.mean<-p.bx[b.use[b],]
	mod.bx.mean<-glm(BX~ns(p.bx.mean,4),family="binomial")
	lines(fitted(mod.bx.mean)[order(p.bx.mean)]~p.bx.mean[order(p.bx.mean)], col="lightblue", lwd=1)}

abline(0,1, lwd=3, lty="dotted")
points(BX[subj_bx>(n0+n1)]~p.bx.mean[subj_bx>(n0+n1)], pch=3)
points(BX[subj_bx%in%c(1:n0)]-0.02~p.bx.mean[subj_bx%in%c(1:n0)], pch=0, col="green3")
points(BX[subj_bx%in%c((n0+1):(n0+n1))]+0.02~p.bx.mean[subj_bx%in%c((n0+1):(n0+n1))], pch=5, col="red")



##RC
plot(RC~p.rc.mean, type="n",  ylab="", xlab="", xaxt="n", yaxt="n") 
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)
mtext("Grade Reclassification",3,line=1, cex=1.5)

for(b in 1:length(b.use)){
	p.rc.mean<-p.rc[b.use[b],]
	mod.rc.mean<-glm(RC~ns(p.rc.mean,4),family="binomial")
	lines(fitted(mod.rc.mean)[order(p.rc.mean)]~p.rc.mean[order(p.rc.mean)], col="lightblue", lwd=1)}

abline(0,1, lwd=3, lty="dotted")
points(RC[subj_rc>(n0+n1)]~p.rc.mean[subj_rc>(n0+n1)], pch=3)
points(RC[subj_rc%in%c(1:n0)]-0.02~p.rc.mean[subj_rc%in%c(1:n0)], pch=0, col="green3")
points(RC[subj_rc%in%c((n0+1):(n0+n1))]+0.02~p.rc.mean[subj_rc%in%c((n0+1):(n0+n1))], pch=5, col="red")

##SURG

plot(SURG~p.surg.mean, type="n",  ylab="", xlab="", xaxt="n", yaxt="n") 
Axis(side=1,labels=TRUE, cex.axis=1.25)
Axis(side=2,labels=TRUE, cex.axis=1.25)
mtext("Surgical Removal of Prostate",3,line=1, cex=1.5)
for(b in 1:length(b.use)){
p.surg.mean<-p.surg[b.use[b],]
mod.surg.mean<-glm(SURG~ns(p.surg.mean,4),family="binomial")
lines(fitted(mod.surg.mean)[order(p.surg.mean)]~p.surg.mean[order(p.surg.mean)], col="lightblue", lwd=1)}

abline(0,1, lwd=3, lty="dotted")
points(SURG[subj_surg>(n0+n1)]~p.surg.mean[subj_surg>(n0+n1)], pch=3)
points(SURG[subj_surg%in%c(1:n0)]-0.02~p.surg.mean[subj_surg%in%c(1:n0)], pch=0, col="green3")
points(SURG[subj_surg%in%c((n0+1):(n0+n1))]+0.02~p.surg.mean[subj_surg%in%c((n0+1):(n0+n1))], pch=5, col="red")

par(mar=c(0,0,0,0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(x="center", inset=0, title="Observed True Cancer State", legend=c("Aggressive (Gleason>6)", "Indolent (Gleason=6)", "Unobserved/No Surgery"), horiz=TRUE, pch=c(5,0,3), col=c("red","green3","black"), cex=1.5) 

dev.off()

#Annotation can be added in keynote/powerpoint or adobe illustrator