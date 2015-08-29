
### Rebecca Yates Coley rycoley@gmail.com
### Code for "Bayesian Joint Hierarchical Model for Prediction of Latent Health States with Application to Active Surveillance of Prostate Cancer"
### This code is used to call the UNADJUSTED JAGS model, that is, it does not allow for an informative observation process (IOP) when doing CROSS-VALIDATION for the latent state


### WORKFLOW: identify observation to be masked, load packages, define data, initialize model parameters, define MCMC settings, run JAGS and save output

### MASK LATENT CLASS OBSERVATION FOR CV
#this code was designed for a SGE cluster system, where task ids can be sent with -t option
#this task id is the subj that will have latent state data excluded for this analysis out of 214 with observed observed state
(to.mask<-as.numeric(Sys.getenv("SGE_TASK_ID")))



### LOAD NECESSARY PACKAGES
library("bayesm")
library("R2jags")


### DEFINE DATA

#get and format data for JAGS model
source("UNADJ-CV-prep-data-for-jags.R")

#bundle data for call to JAGS
#this is observed data and constant variables that have already been assigned values (e.g. number of class K=2, number of subjects n, etc.)

K<-2 #number of latent classes

jags_data<-list(K=K, n=n, eta.data=eta.data, n_eta_known=n_eta_known, n_obs_psa=n_obs_psa, Y=Y, subj_psa=subj_psa, Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z), RC=RC, n_rc=n_rc, subj_rc=subj_rc, V.RC=V.RC.data, d.V.RC=d.V.RC)


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

	gamma.RC<-rnorm((d.V.RC+1), mean=0, sd=0.1) #ditto

	list(p_eta=p_eta, eta.hat=eta.hat, xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta, gamma.RC=gamma.RC)
}


### MCMC SETTINGS

# parameters to track
params <- c("eta.hat")  


# change length; burn-in; number thinned; number of chains
#ni <- 100; nb <- 20; nt <- 5; nc <- 1 
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1 



### RUN JAGS MODEL
set.seed(to.mask)
outj<-jags(jags_data, inits=inits, parameters.to.save=params, model.file="UNADJ-jags-model.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)
out<-outj$BUGSoutput


### IDENTIFY POSTERIOR FOR LATENT STATE, SAVE
(eta.fitted<-out$mean$eta.hat[1])

write.csv(eta.fitted,paste("results/eta-fitted-unadj-",to.mask,".csv",sep=""))
