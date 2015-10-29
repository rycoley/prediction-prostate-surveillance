
### Rebecca Yates Coley rycoley@gmail.com
### Code for "Bayesian Joint Hierarchical Model for Prediction of Latent Health States with Application to Active Surveillance of Prostate Cancer"
### This code is used to check posterior results from IOP results (convergence, etc.)


### WORKFLOW: Load libraries, define functions for data viewing, for each model parameter- pull in posterior, examine results, and save plots in plots/results-check


### LOAD LIBRARIES
list.of.packages <- c("splines")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=T)

library("splines")


### DEFINE FUNCTIONS
expit<-function(x){return(exp(x)/(1+exp(x)))}

color.post<-c("blue","red","green","orange","navy","yellow","purple","dark green","sky blue","magenta")

plot.post.trace<-function(par,nc,n.row, n.col){
	num.ea<-dim(par)[1]/nc
	par(mfrow=c(n.row,n.col), mar=c(2,2,1,1))
	for(i in 1:dim(par)[2]){
		plot(c(min(par[,i]),max(par[,i]))~c(1,num.ea), type="n", xlab="", ylab=paste(i))
		for(j in 1:nc){
			lines(par[((j-1)*num.ea+1):(j*num.ea),i], col=color.post[j])}}}
		
plot.post.dens<-function(par, nc, n.row, n.col){
	num.ea<-dim(par)[1]/nc
	par(mfrow=c(n.row,n.col), mar=c(2,2,1,1))
	for(i in 1:dim(par)[2]){
		plot(c(0, max(density(par[,i])$y)*1.1) ~ c(min(par[,i]),max(par[,i])), type="n", xlab=paste(i), ylab="Density")
		for(j in 1:nc){
			lines(density(par[((j-1)*num.ea+1):(j*num.ea),i]), col=color.post[j])}}}

### EXAMINE POSTERIORS

## LATENT CLASS
etahat<-read.csv(paste("results/jags-prediction-iop-eta_hat-1.csv",sep=""))
dim(etahat)
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-eta_hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:dim(res)[2]])}
	
eta.mean<-apply(etahat,2,mean)
summary(eta.mean)

(a<-dim(etahat)[1]/5*c(0:4)+1)
(b<-dim(etahat)[1]/5*c(1:5))

eta.mean.mat<-as.matrix(cbind(apply(etahat[a[1]:b[1],],2,mean), apply(etahat[a[2]:b[2],],2,mean), apply(etahat[a[3]:b[3],],2,mean), apply(etahat[a[4]:b[4],],2,mean), apply(etahat[a[5]:b[5],],2,mean)))

pdf("plots/results-check/post-density-eta.pdf")
plot(density(eta.mean.mat[,1]),type="l",col="blue",xlim=c(0,1), ylim=c(0,4.1), lwd=2, main="Probablity Lethal PCa")
lines(density(eta.mean.mat[,2]),col="red",lwd=2)
lines(density(eta.mean.mat[,3]),col="green",lwd=2)
lines(density(eta.mean.mat[,4]),col="orange",lwd=2)
lines(density(eta.mean.mat[,5]),col="navy",lwd=2)
dev.off()

## here, we see that the posterior predictions for latent state are similar across 5 separate starting chains





## RHO (shared P(eta=1) across all subjects)
rho<-read.csv(paste("results/jags-prediction-iop-p_eta-1.csv",sep=""))
rho<-as.matrix(rho[,2])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-p_eta-",i,".csv",sep=""))
	rho<-c(rho,res[,2])}
summary(rho)

pdf("plots/results-check/trace-rho.pdf")
plot.post.trace(as.matrix(rho),5,1,1)
dev.off()

pdf("plots/results-check/dens-rho.pdf")
plot.post.dens(as.matrix(rho),5,1,1)
dev.off()
#here, we see that the trace plots for each chain show convergence and that results are similar across sampling chains


##biopsy model
nu.bx<-read.csv(paste("results/jags-prediction-iop-nu_BX-1.csv",sep=""))
dim(nu.bx)
nu.bx<-as.matrix(nu.bx[,2:dim(nu.bx)[2]])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-nu_BX-",i,".csv",sep=""))
	nu.bx<-rbind(nu.bx,res[,2:dim(res)[2]])}

apply(nu.bx,2,summary)

pdf("plots/results-check/trace-nu_bx.pdf")
plot.post.trace(as.matrix(nu.bx),5,5,3)
dev.off()

pdf("plots/results-check/dens-nu_bx.pdf")
plot.post.dens(as.matrix(nu.bx),5,5,3)
dev.off()
#here, we see that the trace plots for each chain show convergence and that results are similar across sampling chains



### reclassification model
gamma.rc<-read.csv(paste("results/jags-prediction-iop-gamma_RC-1.csv",sep=""))
dim(gamma.rc)
gamma.rc<-as.matrix(gamma.rc[,2:dim(gamma.rc)[2]])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-gamma_RC-",i,".csv",sep=""))
	gamma.rc<-rbind(gamma.rc,res[,2:dim(res)[2]])}

apply(gamma.rc,2,summary)


pdf("plots/results-check/trace-gamma_rc.pdf")
plot.post.trace(as.matrix(gamma.rc),5,4,2)
dev.off()

pdf("plots/results-check/dens-gamma_rc.pdf")
plot.post.dens(as.matrix(gamma.rc),5,4,2)
dev.off()
#here, we see that the trace plots for each chain show convergence and that results are similar across sampling chains


##surgery model
omega.surg<-read.csv(paste("results/jags-prediction-iop-omega_SURG-1.csv",sep=""))
dim(omega.surg)
omega.surg<-as.matrix(omega.surg[,2:dim(omega.surg)[2]])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-omega_SURG-",i,".csv",sep=""))
	omega.surg<-rbind(omega.surg,res[,2:dim(res)[2]])}

apply(omega.surg,2,summary)


pdf("plots/results-check/trace-omega_surg.pdf")
plot.post.trace(as.matrix(omega.surg),5,5,3)
dev.off()

pdf("plots/results-check/dens-omega_surg.pdf")
plot.post.dens(as.matrix(omega.surg),5,5,3)
dev.off()
#here, we see that the trace plots for each chain show convergence and that results are similar across sampling chains



###PSA model
#fixed effect for prostate volume
beta<-read.csv(paste("results/jags-prediction-iop-beta-1.csv",sep=""))
dim(beta)
beta<-as.matrix(beta[,2])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-beta-",i,".csv",sep=""))
	beta<-c(beta,res[,2])}
summary(beta)

pdf("plots/results-check/trace-beta.pdf")
plot.post.trace(as.matrix(beta),5,1,1)
dev.off()

pdf("plots/results-check/dens-beta.pdf")
plot.post.dens(as.matrix(beta),5,1,1)
dev.off()
#here, we see that the trace plots for each chain show convergence and that results are similar across sampling chains


#mean random intercept for each class
mu_int<-read.csv(paste("results/jags-prediction-iop-mu_int-1.csv",sep=""))
mu_int<-as.matrix(mu_int[,2:3])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-mu_int-",i,".csv",sep=""))
	mu_int<-rbind(mu_int,res[,2:3])}
apply(mu_int,2,summary)

#mean random slope for each class
mu_slope<-read.csv(paste("results/jags-prediction-iop-mu_slope-1.csv",sep=""))
mu_slope<-as.matrix(mu_slope[,2:3])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-mu_slope-",i,".csv",sep=""))
	mu_slope<-rbind(mu_slope,res[,2:3])}
apply(mu_slope,2,summary)

pdf("plots/results-check/trace-mu.pdf")
plot.post.trace(cbind(mu_int,mu_slope),5,2,2)
dev.off()

pdf("plots/results-check/dens-mu.pdf")
plot.post.dens(cbind(mu_int,mu_slope),5,2,2)
dev.off()
#here, we see that the trace plots for each chain show convergence and that results are similar across sampling chains



#standard deviation for random intercepts
sigma_int<-read.csv(paste("results/jags-prediction-iop-sigma_int-1.csv",sep=""))
sigma_int<-as.matrix(sigma_int[,2])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-sigma_int-",i,".csv",sep=""))
	sigma_int<-c(sigma_int,res[,2])}
summary(sigma_int)

#sandard deviation for random slope
sigma_slope<-read.csv(paste("results/jags-prediction-iop-sigma_slope-1.csv",sep=""))
sigma_slope<-as.matrix(sigma_slope[,2])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-sigma_slope-",i,".csv",sep=""))
	sigma_slope<-c(sigma_slope,res[,2])}
summary(sigma_slope)

#covariance of random effects
cov_int_slope<-read.csv(paste("results/jags-prediction-iop-cov_int_slope-1.csv",sep=""))
cov_int_slope<-as.matrix(cov_int_slope[,2])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-cov_int_slope-",i,".csv",sep=""))
	cov_int_slope<-c(cov_int_slope,res[,2])}
summary(cov_int_slope)

#residual variance for random effects model
sigma_res<-read.csv(paste("results/jags-prediction-iop-sigma_res-1.csv",sep=""))
sigma_res<-as.matrix(sigma_res[,2])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-sigma_res-",i,".csv",sep=""))
	sigma_res<-c(sigma_res,res[,2])}
summary(sigma_res)



pdf("plots/results-check/trace-cov-and-sigma.pdf")
plot.post.trace(cbind(sigma_int, sigma_slope, cov_int_slope, sigma_res),5,2,2)
dev.off()

pdf("plots/results-check/dens-cov-and-simga.pdf")
plot.post.dens(cbind(sigma_int, sigma_slope, cov_int_slope, sigma_res),5,2,2)
dev.off()
#here, we see that the trace plots for each chain show convergence and that results are similar across sampling chains

