### NAME AND EMAIL REDACTED
### Simulation code for "Bayesian Joint Hierarchical Model for Prediciton of Latent Health States with Application to Active Surveillance of Prostate Cancer"
### This code is used to run a single simulation for the jags model that accounts for possible biopsy and surgery informative observation processes (IOP).


##SIMULATIONS WORKFLOW
# (1) Define, set random seed
# (2) Simulate data
# (3) Shape data in prep for JAGS run
# (4) Run JAGS
# (5) Save posterior samples
# (6) Save true and predicted etas



# (1) Define, set random seed
(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))


set.seed(SEED)

# (2) Simulate data
source("simulations/simulation-data/simulate-pca-data.R")

#Remove true coefficient values
rm(beta, gamma_rc, mu_mat, nu_bx, omega_surg, p_eta, Sigma, sigma_res)


# (3) Shape data in prep for JAGS run
source("simulations/IOP/IOP-prep-data-for-jags.R")

# (4) Run JAGS
source("simulations/IOP/IOP-call-jags.R")

# (5) Save posterior samples
for(j in 1:length(out$sims.list)){
	if(!names(out$sims.list)[j]=="eta_hat"){
	write.csv(out$sims.list[[j]], paste("simulations/IOP/results/jags-prediction-iop-", names(out$sims.list)[j],"-",SEED,".csv",sep=""))}}

##Note- These results are not saved in the github folder. (It is too cumbersome). Use sim-IOP-res.R, sim-IOP-res-eta.R, and IS-for-eta-known-IOP.R to obtain summaries of posterior results and predictions

# (6) Save true and predicted etas
index<-c(1:length(names(out$sims.list)))[names(out$sims.list)=="eta_hat"]
eta_hat.mean<-apply(out$sims.list[[index]],2,mean)

eta_to_save<-data.frame(cbind(pt_data$eta_true[pt_data$surg==0], eta_hat.mean))
names(eta_to_save)<-c("true","predicted")

write.csv(eta_to_save, paste0("simulations/IOP/results/eta-true-vs-pred-iop-",SEED,".csv"))

