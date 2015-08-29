
### Rebecca Yates Coley rycoley@gmail.com
### Code for "Bayesian Joint Hierarchical Model for Prediction of Latent Health States with Application to Active Surveillance of Prostate Cancer"
### This code is used to call the JAGS model that allows for an informative observation process (IOP)


### WORKFLOW: load packages, define data, initialize model parameters, define MCMC settings, run JAGS and save output


### LOAD NECESSARY PACKAGES
library("lme4")
library("bayesm")
library("R2jags")



### DEFINE DATA

#get and format data for JAGS model
source("IOP-prep-data-for-jags.R")

K<-2 #number of latent classes

jags_data<-list(K=K, n=n, eta.data=eta.data, n_eta_known=n_eta_known, n_obs_psa=n_obs_psa, Y=Y, subj_psa=subj_psa, Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z), BX=BX, n_bx=n_bx, subj_bx=subj_bx, U.BX=U.BX.data, d.U.BX=d.U.BX, RC=RC, n_rc=n_rc, subj_rc=subj_rc, V.RC=V.RC.data, d.V.RC=d.V.RC, SURG=SURG, n_surg=n_surg, subj_surg=subj_surg, W.SURG=W.SURG.data, d.W.SURG=d.W.SURG)


### INITIALIZE MODEL PARAMETERS
#note that not all "parameters" need to be initialized. specifically, don't initialize random effects, but do need to set initial values for mean and covariance of random effects
# also need to set initial values for latent class when unobserved (here, eta)

inits <- function(){
		
	p_eta<-rbeta(1,1,1)

	eta.hat<-pt.data$rc[is.na(eta.data)]

	xi<-c(min(rlnorm(1),100), min(rlnorm(1),100))
	mu_raw<-as.matrix(cbind(rnorm(d.Z),rnorm(d.Z)))
	Tau_B_raw<-rwishart((d.Z+1),diag(d.Z)*var_vec)$W
	sigma_res<-min(rlnorm(1),1)

	beta<-rnorm(d.X)

	nu.BX<-rnorm((d.U.BX+1), mean=0, sd=0.1) #last coefficient is effect of eta=1
	gamma.RC<-rnorm((d.V.RC+1), mean=0, sd=0.1) #ditto
	omega.SURG<-c(rnorm((d.W.SURG+2), mean=0, sd=0.01))  #here, include interaction with last prediction and eta=1

	list(p_eta=p_eta, eta.hat=eta.hat, xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta, nu.BX=nu.BX, gamma.RC=gamma.RC, omega.SURG=omega.SURG)
}



### MCMC SETTINGS

# parameters to track
params <- c("p_eta", "eta.hat", "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "rho_int_slope", "cov_int_slope", "b.vec", "beta", "nu.BX", "gamma.RC", "omega.SURG", "p_bx", "p_rc", "p_surg")  


# change length; burn-in; number thinned; number of chains
#ni <- 100; nb <- 20; nt <- 5; nc <- 1 
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1 



### RUN JAGS MODEL, SAVE RESULTS
#this function fits the model in JAGS with a given starting seed, and saves the output in .csv files named based on seed and parameter name
do.one<-function(seed){
	set.seed(seed)	
	outj<-jags(jags_data, inits=inits, parameters.to.save=params, model.file="IOP-jags-model.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

out<-outj$BUGSoutput

for(j in 1:length(out$sims.list)){
	write.csv(out$sims.list[[j]], paste("results/jags-prediction-iop-", names(out$sims.list)[j],"-",seed,".csv",sep=""))}
}

### this code was written to run in parallel (note- only 1 chain specified above)
### specifically written for a SGE cluster where task ids can be sent with -t option
### users can edit this part of the code for their own system
(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))
do.one(seed=SEED)







