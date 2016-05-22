### NAME AND EMAIL REDACTED
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to call the JAGS model that allows for a possible surgery informative observation process (IOP)
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.


### WORKFLOW: load packages, define data, initialize model parameters, define MCMC settings, run JAGS and save output in results folder
#***create folder for results before running***

rm(list=ls())

### LOAD NECESSARY PACKAGES
list.of.packages <- c("lme4", "bayesm", "R2jags")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=T)

library("lme4")
library("bayesm")
library("R2jags")



### DEFINE DATA

#get and format data for JAGS model
source("jhas-analysis/IOP-SURG/IOP-SURG-prep-data-for-jags.R")

K<-2 #number of latent classes

jags_data<-list(K=K, n=n, eta_data=eta_data, n_eta_known=n_eta_known, n_obs_psa=n_obs_psa, Y=Y, subj_psa=subj_psa, Z=Z_data, X=X_data, d_Z=d_Z, d_X=d_X, I_d_Z=diag(d_Z), RC=RC, n_rc=n_rc, subj_rc=subj_rc, V_RC=V_RC_data, d_V_RC=d_V_RC, SURG=SURG, n_surg=n_surg, subj_surg=subj_surg, W_SURG=W_SURG_data, d_W_SURG=d_W_SURG)


### INITIALIZE MODEL PARAMETERS
#note that not all "parameters" need to be initialized. specifically, don't initialize random effects, but do need to set initial values for mean and covariance of random effects
# also need to set initial values for latent class when unobserved (here, eta)

inits <- function(){
		
	p_eta<-rbeta(1,1,1)

	eta_hat<-pt_data$rc[is.na(eta_data)]

	xi<-c(min(rlnorm(1),100), min(rlnorm(1),100))
	mu_raw<-as.matrix(cbind(rnorm(d_Z),rnorm(d_Z)))
	Tau_B_raw<-rwishart((d_Z+1),diag(d_Z)*var_vec)$W
	sigma_res<-min(rlnorm(1),1)

	beta<-rnorm(d_X)


	gamma_RC<-rnorm((d_V_RC+1), mean=0, sd=0.1) #ditto
	omega_SURG<-c(rnorm((d_W_SURG+2), mean=0, sd=0.01))  #here, include interaction with last prediction and eta=1

	list(p_eta=p_eta, eta_hat=eta_hat, xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta, gamma_RC=gamma_RC, omega_SURG=omega_SURG)
}



### MCMC SETTINGS

# parameters to track
params <- c("p_eta", "eta_hat", "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "rho_int_slope", "cov_int_slope", "b_vec", "beta",  "gamma_RC", "omega_SURG", "p_rc", "p_surg")  


# change length; burn-in; number thinned; number of chains
#ni <- 100; nb <- 20; nt <- 5; nc <- 1 
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1 



### RUN JAGS MODEL, SAVE RESULTS
#this function fits the model in JAGS with a given starting seed, and saves the output in .csv files named based on seed and parameter name
do_one<-function(seed){
	set.seed(seed)	
	outj<-jags(jags_data, inits=inits, parameters.to.save=params, model.file="jhas-analysis/IOP-SURG/IOP-SURG-jags-model.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

out<-outj$BUGSoutput

for(j in 1:length(out$sims.list)){
	write.csv(out$sims.list[[j]], paste("jhas-analysis/IOP-SURG/results/jags-prediction-iop-surg-", names(out$sims.list)[j],"-",seed,".csv",sep=""))}
}

### this code was written to run in parallel (note- only 1 chain specified above)
### specifically written for a SGE cluster where task ids can be sent with -t option
### users can edit this part of the code for their own system
(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))
do_one(seed=SEED)







