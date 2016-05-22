### NAME AND EMAIL REDACTED
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
# This code reproduces the plot in Appendix Figure 4.
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.


library("mixtools")

pt.data<-read.csv("jhas-analysis/simulation-data/pt-data-sim.csv")
psa.data<-read.csv("jhas-analysis/simulation-data/psa-data-sim.csv")

(n<-dim(pt.data)[1]) #1000

(n_obs_psa<-dim(psa.data)[1])
Y<-psa.data$log_psa
subj_psa<-psa.data$subj


names(psa.data)

#covariates with random effects
Z.data<-as.matrix(cbind(rep(1,n_obs_psa), psa.data$age_std)) 

#covariates with only fixed effects
X.data<-as.matrix(cbind(psa.data$vol_std)) 



## PSA model estimates
mu_int<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-mu_int-1.csv",sep=""))
mu_int<-as.matrix(mu_int[,2:3])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-mu_int-",i,".csv",sep=""))
	mu_int<-rbind(mu_int,res[,2:3])}
apply(mu_int,2,summary)
mu_int<-apply(mu_int,2,mean)

mu_slope<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-mu_slope-1.csv",sep=""))
mu_slope<-as.matrix(mu_slope[,2:3])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-mu_slope-",i,".csv",sep=""))
	mu_slope<-rbind(mu_slope,res[,2:3])}
apply(mu_slope,2,summary)
mu_slope<-apply(mu_slope,2,mean)

mu0<-c(mu_int[1], mu_slope[1])
mu1<-c(mu_int[2], mu_slope[2])

sigma_int<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-sigma_int-1.csv",sep=""))
sigma_int<-as.matrix(sigma_int[,2])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-sigma_int-",i,".csv",sep=""))
	sigma_int<-c(sigma_int,res[,2])}
summary(sigma_int)
var_int<-mean(sigma_int)^2

sigma_slope<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-sigma_slope-1.csv",sep=""))
sigma_slope<-as.matrix(sigma_slope[,2])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-sigma_slope-",i,".csv",sep=""))
	sigma_slope<-c(sigma_slope,res[,2])}
summary(sigma_slope)
var_slope<-mean(sigma_slope)^2

cov_int_slope<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-cov_int_slope-1.csv",sep=""))
cov_int_slope<-as.matrix(cov_int_slope[,2])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-cov_int_slope-",i,".csv",sep=""))
	cov_int_slope<-c(cov_int_slope,res[,2])}
summary(cov_int_slope)
cov_int_slope<-mean(cov_int_slope)

(Sigma<-matrix(c(var_int, cov_int_slope, cov_int_slope, var_slope), nrow=2, ncol=2))

bet<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-beta-1.csv",sep=""))
bet<-as.matrix(bet[,2])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-beta-",i,".csv",sep=""))
	bet<-c(bet,res[,2])}
bet<-mean(bet)



#individual-level estimates
etahat<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-eta_hat-1.csv",sep=""))
dim(etahat)
etahat<-as.matrix(etahat[,2:798])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-eta_hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:798])}
eta.mean<-apply(etahat,2,mean)



b.vec<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-b_vec-1-reduced.csv",sep=""))
#dim(b.vec)
b.int<-as.matrix(b.vec[,2:(n+1)])
b.slope<-as.matrix(b.vec[,(n+2):(2*n+1)])
#for(i in 2:5){
#	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-b.vec-",i,".csv",sep=""))
#	b.int<-rbind(b.int,res[,2:(n+1)])
#	b.slope<-rbind(b.slope,res[,(n+2):(2*n+1)])}

b.int<-apply(b.int,2,mean)
b.slope<-apply(b.slope,2,mean)





pdf("jhas-analysis/plots/appendix-figure4-individual-effects.pdf", width=9, height=7)


plot(c(min(b.slope)-0.1, max(b.slope)+0.1)~c(min(b.int), max(b.int)), type="n", xlab="Intercept", ylab="Slope" )

ellipse(mu0, Sigma, alpha=0.05, lwd=3, col="green3", lty="dotted")
ellipse(mu0, Sigma, alpha=0.5, lwd=3, col="green3", lty="dashed")
ellipse(mu1, Sigma, alpha=0.05, lwd=3, col="red", lty="dotted")
ellipse(mu1, Sigma, alpha=0.5, lwd=3, col="red", lty="dashed")

#sum(pt.data$obs_eta==0 & !is.na(pt.data$obs_eta))

for(i in 1:90){
	points(b.slope[i]~b.int[i], pch=19, col="green3", cex=0.4)}
for(i in 91:203){
	points(b.slope[i]~b.int[i], pch=19, col="red", cex=0.4)}


for(i in 204:n){
	if(eta.mean[(i-161)]<=quantile(eta.mean,0.25)){points(b.slope[i]~b.int[i], pch=21, col="green3", cex=0.4)}
	else{ if(eta.mean[(i-161)]>quantile(eta.mean,0.25) & eta.mean[(i-161)]<=median(eta.mean)){points(b.slope[i]~b.int[i], pch=21, col="greenyellow", cex=0.4)}
		else{
	if(eta.mean[(i-161)]<=quantile(eta.mean,0.75) & eta.mean[(i-161)]>median(eta.mean)){points(b.slope[i]~b.int[i], pch=21, col="orange", cex=0.4)}
		else{
	if(eta.mean[(i-161)]>quantile(eta.mean,0.75)){points(b.slope[i]~b.int[i], pch=21, col="red", cex=0.4)}}
	}}}


legend("topleft", title="P(Aggressive)", legend=c("0-25%", "26-50%", "51-75%", "75-100%"), col=c("green3", "greenyellow", "orange", "red"), pch=rep(21,4), bty="n")
legend(0.2,1.75, title="State Observed", legend=c("Indolent", "Aggressive"), col=c("green3","red"), pch=rep(19,2), bty="n")

legend("bottomleft", title="Credible Intervals", legend=c("50%", "95%"), lty=c("dashed","dotted"), lwd=rep(3,2), bty="n")


dev.off()

