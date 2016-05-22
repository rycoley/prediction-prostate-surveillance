### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
#This code performs the following:
#1. Get posterior predictions of eta for patients who had surgery through importance sampling algorithm and
#2. Compare these to known eta. 
#3. Repeat for all 200 simulated datasets

#More information on the IS algorithm can be found in Fisher et al. (2015) "Fast out-of-sample predictions for Bayesian hierarchical models of latent health states"

library("dplyr")
library("ROCR")

nsims<-200

auc.res<-mse.res<-mae.res<-vector(length=nsims)

for(SEED in 1:nsims){
##loop through the following for all simulations 

#SEED<-1
set.seed(SEED)


#Simulate and shape data
source("simulations/simulation-data/simulate-pca-data.R")
rm(beta, gamma_rc, mu_mat, nu_bx, omega_surg, p_eta, Sigma, sigma_res)

source("simulations/UNADJ/UNADJ-prep-data-for-jags.R")
K<-2


#get posterior samples
params <- c("p_eta",  "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "cov_int_slope", "beta", "gamma_RC")  


nrep<-10

for(p in 1:length(params)){
	res<-NA
	res<-read.csv(paste0("simulations/UNADJ/results/jags-prediction-unadj-", params[p], "-", SEED, ".csv"))
	D<-dim(res)[2]
	if(D==2){assign(params[p], rep(res[,2], nrep) )}
	else{ assign(params[p], matrix(rep(t(res[,2:D]), nrep), ncol=D-1, byrow=T) ) }
		}

(B<-length(p_eta))


#generate candidate values for eta, b_vec, and epsilon
cands_eta<-rbinom(B,1,p_eta)

cands_bvec<-matrix(nrow=B, ncol=2) #is there a non-for loop way to do this?
for(b in 1:B){
	cov.mat<-diag(c(sigma_int[b]^2, sigma_slope[b]^2))
	cov.mat[1,2]<-cov.mat[2,1]<-cov_int_slope[b]
	cands_bvec[b,]<-mvrnorm(1, mu=c(mu_int[b,(cands_eta[b]+1)], mu_slope[b,(cands_eta[b]+1)]), Sigma=cov.mat)}


#I think I only generated this for the prediction intervals
#cands_ep<-matrix(nrow=B, ncol=30) #I assume there are at most 30 PSA obs per person
#for(b in 1:B){cands_ep[b,]<-rnorm(30,mean=0, sd=sigma_res[b])}


#clean up workspace a bit
rm(mu_int, mu_slope, p_eta, params, cov_int_slope, sigma_int, sigma_slope, cov.mat, res)


##loop through the following for each patient with eta known
#should I also loop through generating candidate values separately for each patient? It is probably not necessary with such large posterior samples

eta_known <- eta_data[!is.na(eta_data)]
eta_pred <- vector(length=n_eta_known)


for(i in 1:n_eta_known){
#i<-1

## define data for that patient
#psa data
Yi <- Y[subj_psa==i]
n_psai <- length(Yi)
Zi <- Z_data[subj_psa==i,]
Xi <- X_data[subj_psa==i,]

#reclassification outcomes
if(sum(subj_rc==i)>0){
RCi <- RC[subj_rc==i]
n_rci <- length(RCi)
Vi <- V_RC_data[subj_rc==i,]
}



## get weights for patient

#psa
Zb <- tcrossprod(cands_bvec, Zi)
Xbeta <- tcrossprod(beta,Xi)
mu_psa <- c( t( Zb+Xbeta ) )
log_lik_psa_j <- log(dnorm(rep(Yi,B), mean=mu_psa, sd=rep(sigma_res, each=n_psai)))
log_lik_psa <- (data.frame('log_lik_psa_all'=log_lik_psa_j, 'p_ind'=as.factor(rep(1:B, each=n_psai))) %>% group_by(p_ind) %>% dplyr::summarize(sum=sum(log_lik_psa_all)))$sum


#reclassification
if(sum(subj_rc==i)>0){
etagamma <- rep(gamma_RC[,(d_V_RC+1)]*cands_eta, each=n_rci)
Vgamma <- c( tcrossprod( Vi, gamma_RC[,1:d_V_RC] ) )
logit_p <- etagamma + Vgamma
p_rc <- c( t( expit( logit_p ) ) )
log_lik_rc_j <- log( dbinom(x=rep(RCi,times=B), size=1, prob=p_rc) )
log_lik_rc <- (data.frame('log_lik_rc_all'=c(t(log_lik_rc_j)), 'p_ind'=as.factor(rep(1:B, each=n_rci)) ) %>% group_by(p_ind) %>% dplyr::summarize(sum=sum(log_lik_rc_all)) )$sum 
}


if(sum(subj_rc==i)>0){lik <- exp(log_lik_psa + log_lik_rc)}
else{lik <- exp(log_lik_psa)}
#un <- lik #unscaled weights
wt.std <- lik/sum(lik) #scaled weights


# weight candidate etas to get posterior prediction of eta
eta_pred[i] <- crossprod(cands_eta, wt.std)
#print(i)
}


#save true and predicted eta for each simulated dataset
eta_to_save <- data.frame(cbind(eta_known, eta_pred))
names(eta_to_save) <- c("true", "predicted")

write.csv(eta_to_save, paste0("simulations/UNADJ/results/eta-known-true-vs-pred-unadj-",SEED,".csv"))


#calculate preditive accuracy

auc.res[SEED] <- performance(prediction(eta_pred, eta_known), "auc")@y.values[[1]]
mse.res[SEED] <- mean((eta_known-eta_pred)^2)
mae.res[SEED] <- mean(eta_known-eta_pred)

print(SEED)

}

save(auc.res=auc.res, mse.res=mse.res, mae.res=mae.res, file="simulations/UNADJ/results/sim-UNADJ-res-eta-known.RData")

