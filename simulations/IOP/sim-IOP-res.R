### NAME AND EMAIL REDACTED
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to summarize posterior estimates and coverage of the joint posterior from the biopsy and surgery IOP model

load("simulations/simulation-data/pars-for-data-sim.RData")


params <- c("p_eta",  "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "cov_int_slope", "beta", "nu_BX", "gamma_RC", "omega_SURG")  


nsims<-200


cov.fcn<-function(post, true){
	cov.vec<-vector(length=length(true))
	for(d in 1:length(true)){
		cov.vec[d] <- as.numeric( quantile( post[,d], p=0.025 ) < true[d] & quantile( post[,d], p=0.975) > true[d] ) }
	return(cov.vec)} 

p_eta.cov<-p_eta.est<-vector(length=nsims)
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-p_eta-",j,".csv") )
	p_eta.est[j]<-median(res[,2])
	p_eta.cov[j]<-as.numeric( quantile( res[,2], p=0.025 ) < p_eta & quantile( res[,2], p=0.975) > p_eta )}


mu_int.cov<-mu_int.est<-matrix(nrow=nsims, ncol=2)
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-mu_int-",j,".csv") )
	mu_int.est[j,]<-apply(res[,2:dim(res)[2]],2,median)
	mu_int.cov[j,]<-cov.fcn(post=res[,2:dim(res)[2]], true=mu_mat[1,]) }	

mu_slope.cov<-mu_slope.est<-matrix(nrow=nsims, ncol=2)
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-mu_slope-",j,".csv") )
	mu_slope.est[j,]<-apply(res[,2:dim(res)[2]],2,median)
	mu_slope.cov[j,]<-cov.fcn(post=res[,2:dim(res)[2]], true=mu_mat[2,]) }	


sigma_int.cov<-sigma_int.est<-vector(length=nsims)
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-sigma_int-",j,".csv") )
	sigma_int.est[j]<-median(res[,2])
	sigma_int.cov[j]<-as.numeric( quantile( res[,2], p=0.025 ) < sqrt(Sigma[1,1]) & quantile( res[,2], p=0.975) > sqrt(Sigma[1,1]) )}

sigma_slope.cov<-sigma_slope.est<-vector(length=nsims)
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-sigma_slope-",j,".csv") )
	sigma_slope.est[j]<-median(res[,2])
	sigma_slope.cov[j]<-as.numeric( quantile( res[,2], p=0.025 ) < sqrt(Sigma[2,2]) & quantile( res[,2], p=0.975) > sqrt(Sigma[2,2]) )}

cov_int_slope.cov<-cov_int_slope.est<-vector(length=nsims)
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-cov_int_slope-",j,".csv") )
	cov_int_slope.est[j]<-median(res[,2])
	cov_int_slope.cov[j]<-as.numeric( quantile( res[,2], p=0.025 ) < Sigma[1,2] & quantile( res[,2], p=0.975) > Sigma[1,2] )}

beta.cov<-beta.est<-vector(length=nsims)
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-beta-",j,".csv") )
	beta.est[j]<-median(res[,2])
	beta.cov[j]<-as.numeric( quantile( res[,2], p=0.025 ) < beta & quantile( res[,2], p=0.975) > beta )}

sigma_res.cov<-sigma_res.est<-vector(length=nsims)
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-sigma_res-",j,".csv") )
	sigma_res.est[j]<-median(res[,2])
	sigma_res.cov[j]<-as.numeric( quantile( res[,2], p=0.025 ) < sigma_res & quantile( res[,2], p=0.975) > sigma_res )}



nu_BX.cov<-nu_BX.est<-matrix(nrow=nsims, ncol=length(nu_bx))
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-nu_BX-",j,".csv") )
	nu_BX.est[j,]<-apply(res[,2:dim(res)[2]],2,median)
	nu_BX.cov[j,]<-cov.fcn(post=res[,2:dim(res)[2]], true=nu_bx) }

gamma_RC.cov<-gamma_RC.est<-matrix(nrow=nsims, ncol=length(gamma_rc))
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-gamma_RC-",j,".csv") )
	gamma_RC.est[j,]<-apply(res[,2:dim(res)[2]],2,median)
	gamma_RC.cov[j,]<-cov.fcn(post=res[,2:dim(res)[2]], true=gamma_rc) }

omega_SURG.cov<-omega_SURG.est<-matrix(nrow=nsims, ncol=length(omega_surg))
for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/jags-prediction-iop-omega_SURG-",j,".csv") )
	omega_SURG.est[j,]<-apply(res[,2:dim(res)[2]],2,median)
	omega_SURG.cov[j,]<-cov.fcn(post=res[,2:dim(res)[2]], true=omega_surg) }



save(p_eta.est=p_eta.est, p_eta.cov=p_eta.cov, mu_int.est=mu_int.est, mu_int.cov=mu_int.cov, mu_slope.est=mu_slope.est, mu_slope.cov=mu_slope.cov, sigma_int.est=sigma_int.est, sigma_int.cov=sigma_int.cov, sigma_slope.est=sigma_slope.est, sigma_slope.cov=sigma_slope.cov, cov_int_slope.est=cov_int_slope.est, cov_int_slope.cov=cov_int_slope.cov, sigma_res.est=sigma_res.est, sigma_res.cov=sigma_res.cov, beta.cov=beta.cov, beta.est=beta.est, nu_BX.est=nu_BX.est, nu_BX.cov=nu_BX.cov, gamma_RC.est=gamma_RC.est, gamma_RC.cov=gamma_RC.cov, omega_SURG.est=omega_SURG.est, omega_SURG.cov=omega_SURG.cov, file="simulations/IOP/results/sim-IOP-results.RData")
