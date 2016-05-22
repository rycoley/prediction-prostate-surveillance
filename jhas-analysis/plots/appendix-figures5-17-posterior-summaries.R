### NAME AND EMAIL REDACTED
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
#This code reproduces plots in Figures 5-17 in the appendix and posterior summaries reported in tables 3-7 in the applendix/
#This code obtains results for the IOP model only. Change file names to get results for the unadjusted, biopsy IOP, an surgery IOP models.
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.


library(coda)

color.post<-c("blue","red","green","orange","navy")


#this function produces a trace plot, cumulative quantile plot, and density plot for all posteriors
plot.post<-function(par, nc=5, var.names, prior.samp, pdf.name){
	P<-dim(par)[2]
	pdf.height<-ifelse(P<5,2*P,9)
	pdf(file=pdf.name, width=7, height=2*P)
	num.ea<-dim(par)[1]/nc
	par(mfrow=c(dim(par)[2], 3), mar=c(3.5,3.5,1,1))
	for(i in 1:dim(par)[2]){
#i<-1
		#trace plot
		plot(c(min(par[,i]),max(par[,i]))~c(1,num.ea), type="n", xlab="", ylab="")
		mtext("Iteration", 1, line=2.25)
		mtext(var.names[i], 2, line=2)
		
		for(j in 1:nc){lines(par[((j-1)*num.ea+1):(j*num.ea),i], col=color.post[j])}
			
		#cumulative quantile plot
		cumuplot(par[1:num.ea,i], ylab="", xlab="", auto.layout=FALSE)
		mtext("Iteration", 1, line=2.25)
		mtext(var.names[i], 2, line=2)
				
		#prior vs posterior
		plot(0, type="n", ylim=c(0,1.1*max(density(par[,i])$y)), xlim=range(par[,i]), ylab="", xlab="")
		mtext("Density", 2, line=2)
		mtext(var.names[i], 1, line=2.25)

		lines(density(prior.samp), lty="dotted", lwd=2)
		for(j in 1:nc){lines(density(par[((j-1)*num.ea+1):(j*num.ea),i]), col=color.post[j])}	
		} 
		dev.off()}



p_eta<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-p_eta-1.csv",sep=""))
p_eta<-as.matrix(p_eta[,2])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-p_eta-",i,".csv",sep=""))
	p_eta<-c(p_eta,res[,2])}

quantile(p_eta,c(0.025, 0.5, 0.975))

plot.post(par=as.matrix(p_eta), var.names="rho", prior.samp=rbeta(10000,1,1), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-rho.pdf")


nu.bx<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-nu_BX-1.csv",sep=""))
nu.bx<-as.matrix(nu.bx[,2:dim(nu.bx)[2]])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-nu_BX-",i,".csv",sep=""))
	nu.bx<-rbind(nu.bx,res[,2:dim(res)[2]])}

round(apply(nu.bx,2,quantile, p=c(0.5, 0.025, 0.975)),3)

dim(nu.bx)

plot.post(par=nu.bx[,1:5], var.names=c("Intercept", rep("Time (NS)",4)), prior.samp=rnorm(100000,0,10), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-nu-1.pdf")
plot.post(par=nu.bx[,6:9], var.names=c(rep("Date (NS)",4)), prior.samp=rnorm(100000,0,10), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-nu-2.pdf")
plot.post(par=nu.bx[,10:14], var.names=c(rep("Age (NS)",2), rep("# Prev Biopsies",2), "True State"), prior.samp=rnorm(100000,0,10), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-nu-3.pdf")


gam<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-gamma_RC-1.csv",sep=""))
gam<-as.matrix(gam[,2:dim(gam)[2]])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-gamma_RC-",i,".csv",sep=""))
	gam<-rbind(gam,res[,2:dim(res)[2]])}
round(apply(gam,2,quantile, p=c(0.5, 0.025, 0.975)),3)
dim(gam) #7

plot.post(par=gam[,1:3], var.names=c("Intercept", rep("Time (NS)",2)), prior.samp=rnorm(100000,0,10), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-gamma-1.pdf")
plot.post(par=gam[,4:7], var.names=c(rep("Date (NS)",2), "Age", "True State"), prior.samp=rnorm(100000,0,10), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-gamma-2.pdf")

omega.surg<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-omega_SURG-1.csv",sep=""))
dim(omega.surg)
omega.surg<-as.matrix(omega.surg[,2:15])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-omega_SURG-",i,".csv",sep=""))
	omega.surg<-rbind(omega.surg,res[,2:15])}
round(apply(omega.surg,2,quantile, p=c(0.5, 0.025, 0.975)),3)

plot.post(par=omega.surg[,1:5], var.names=c("Intercept", rep("Time (NS)",4)), prior.samp=rnorm(100000,0,10), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-omega-1.pdf")
plot.post(par=omega.surg[,6:10], var.names=c(rep("Date (NS)",3), rep("Age (NS)", 2)), prior.samp=rnorm(100000,0,10), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-omega-2.pdf")
plot.post(par=omega.surg[,11:14], var.names=c("# Prev Biopsies","Previous R=1", "True State", "Previous R=1 x True State"), prior.samp=rnorm(100000,0,10), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-omega-3.pdf")



##PSA Model

mu_int<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-mu_int-1.csv",sep=""))
mu_int<-as.matrix(mu_int[,2:3])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-mu_int-",i,".csv",sep=""))
	mu_int<-rbind(mu_int,res[,2:3])}
apply(mu_int,2,quantile, p=c(0.5,0.025, 0.975))

mu_slope<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-mu_slope-1.csv",sep=""))
mu_slope<-as.matrix(mu_slope[,2:3])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-mu_slope-",i,".csv",sep=""))
	mu_slope<-rbind(mu_slope,res[,2:3])}
apply(mu_slope,2,quantile, p=c(0.5,0.025, 0.975))

plot.post(par=cbind(mu_int, mu_slope), var.names=c("Mean Intercept, eta=0", "Mean Intercept, eta=1", "Mean Slope, eta=0", "Mean Slope, eta=1"), prior.samp=rnorm(100000,0,10), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-mu.pdf")



bet<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-beta-1.csv",sep=""))
dim(bet)
bet<-as.matrix(bet[,2])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-beta-",i,".csv",sep=""))
	bet<-c(bet,res[,2])}
quantile(bet,p=c(0.025,0.5,0.975))

plot.post(par=as.matrix(bet), var.names="Fixed Effect, Volume", prior.samp=rnorm(10000,0,10), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-beta.pdf")



sigma_int<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-sigma_int-1.csv",sep=""))
sigma_int<-as.matrix(sigma_int[,2])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-sigma_int-",i,".csv",sep=""))
	sigma_int<-c(sigma_int,res[,2])}
quantile(sigma_int, p=c(0.5, 0.025, 0.975))

sigma_slope<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-sigma_slope-1.csv",sep=""))
sigma_slope<-as.matrix(sigma_slope[,2])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-sigma_slope-",i,".csv",sep=""))
	sigma_slope<-c(sigma_slope,res[,2])}
quantile(sigma_slope, p=c(0.5, 0.025, 0.975))

cov_int_slope<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-cov_int_slope-1.csv",sep=""))
cov_int_slope<-as.matrix(cov_int_slope[,2])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-cov_int_slope-",i,".csv",sep=""))
	cov_int_slope<-c(cov_int_slope,res[,2])}
quantile(cov_int_slope, p=c(0.5, 0.025, 0.975))


sigma_res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-sigma_res-1.csv",sep=""))
sigma_res<-as.matrix(sigma_res[,2])
for(i in 2:5){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-sigma_res-",i,".csv",sep=""))
	sigma_res<-c(sigma_res,res[,2])}
quantile(sigma_res, p=c(0.5, 0.025, 0.975))

plot.post(par=cbind(sigma_int, sigma_slope, cov_int_slope, sigma_res), var.names=c("Variance, Intercepts", "Variance, Slopes", "Covariance, Int & Slope", "Variance, Residuals"), prior.samp=runif(1000000,0,1), pdf.name="jhas-analysis/plots/appendix-figures5-17/post-sigma.pdf")