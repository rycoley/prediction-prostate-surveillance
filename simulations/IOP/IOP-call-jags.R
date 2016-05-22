### NAME AND EMAIL REDACTED
### Simulation code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to call the JAGS model that allows for possible biopsy and surgery informative observation processes (IOP)


### WORKFLOW: load packages, define data, initialize model parameters, define MCMC settings, run JAGS and save output



### LOAD NECESSARY PACKAGES
library("bayesm")
library("R2jags")



### DEFINE DATA


K<-2 #number of latent classes

jags_data<-list(K=K, n=n, eta_data=eta_data, n_eta_known=n_eta_known, n_obs_psa=n_obs_psa, Y=Y, subj_psa=subj_psa, Z=Z_data, X=X_data, d_Z=d_Z, d_X=d_X, I_d_Z=diag(d_Z), BX=BX, n_bx=n_bx, subj_bx=subj_bx, U_BX=U_BX_data, d_U_BX=d_U_BX, RC=RC, n_rc=n_rc, subj_rc=subj_rc, V_RC=V_RC_data, d_V_RC=d_V_RC, SURG=SURG, n_surg=n_surg, subj_surg=subj_surg, W_SURG=W_SURG_data, d_W_SURG=d_W_SURG)


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

	nu_BX<-rnorm((d_U_BX+1), mean=0, sd=0.1) #last coefficient is effect of eta=1
	gamma_RC<-rnorm((d_V_RC+1), mean=0, sd=0.1) #ditto
	omega_SURG<-c(rnorm((d_W_SURG+2), mean=0, sd=0.01))  #here, include interaction with last prediction and eta=1

	list(p_eta=p_eta, eta_hat=eta_hat, xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta, nu_BX=nu_BX, gamma_RC=gamma_RC, omega_SURG=omega_SURG)
}



### MCMC SETTINGS

# parameters to track
params <- c("eta_hat", "p_eta",  "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "cov_int_slope", "beta", "nu_BX", "gamma_RC", "omega_SURG")  


# change length; burn-in; number thinned; number of chains
#ni <- 100; nb <- 20; nt <- 5; nc <- 2 
ni <- 50000; nb <- 25000; nt <- 20; nc <- 2 



### RUN JAGS MODEL
set.seed(SEED)	
outj<-jags(jags_data, inits=inits, parameters.to.save=params, model.file="simulations/IOP/IOP-jags-model.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

out<-outj$BUGSoutput


