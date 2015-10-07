### Rebecca Yates Coley rycoley@gmail.com
### Code for "Bayesian Joint Hierarchical Model for Prediction of Latent Health States with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in and format all the data needed to run the unadjusted jags model.
### Simulated data is ammended to remove last 5 years data for 12 patients, in order to demonstrate making individualized predictions (Figure 3). 


#get data
pt.data<-read.csv("simulation-data/pt-data-sim.csv")
#this contains one record per patient
#patients are ordered so that those with surgery, i.e. eta observed, occur first. 

psa.data<-read.csv("simulation-data/psa-data-sim.csv")
#this contains one record per psa observations per patient

data.use<-read.csv("simulation-data/bx-data-sim.csv")
#this contains one record per annual interval for each patient until surgery or censoring


#Before call to JAGS, get the data into simple matrices and vectors to send to JAGS
(n<-dim(pt.data)[1]) #there are 1000 patients

#get observed latent class for those with surgery
eta.data<-pt.data$obs.eta
(n_eta_known<-sum(!is.na(eta.data))) #214


#remove last 5 years of data for selected patients, will then use these to demonstrate predictions
set.seed(4)
pred_ids<-sample(c((n_eta_known+1):n), 12)
(pred_ids<-sort(pred_ids))

psa_data<-psa_data[!(psa_data$subj%in%pred_ids & psa_data$psa_time>5),] 
data_use<-data_use[!(data_use$subj%in%pred_ids & data_use$time>5),]


#PSA model
(n_obs_psa<-dim(psa_data)[1])
Y<-psa_data$log_psa
subj_psa<-psa_data$subj

#covariates with random effects 
#here, intercept and age (standardized)
Z_data<-as.matrix(cbind(rep(1,n_obs_psa), psa_data$age_std))
(d_Z<-dim(Z_data)[2])

#covariates with only fixed effects
#here, prostate volume (standardized)
X_data<-as.matrix(cbind(psa_data$vol_std))
(d_X<-dim(X_data)[2])


#lmer fit to get starting value for covariance parameter in JAGS
mod_lmer<-lmer(log_psa~ vol_std + (1+ age_std |id), data=psa_data)
(var_vec <- apply(coef(mod_lmer)$id, 2, var)[1:d_Z])
(var_vec <- c(var_vec[2], var_vec[1]))


###bx data

#outcome model (logistic regression for reclassification)
#only use records where a biopsy occurred
rc_data<-data_use[data_use$bx_here==1 & !is.na(data_use$bx_here),] 
(n_rc<-dim(rc_data)[1])
RC<-as.numeric(rc_data$rc) #indicator of reclassificaiton in this interval
subj_rc<-rc_data$subj

#covariate matrix V for pooled logistic regression predicting reclassification
#here, age (standardized), time since diagnosis (ns with 2 df), calendar time
V_RC_data<-as.matrix(cbind(rep(1,n_rc),  rc_data[,grep("rc_time_ns", names(rc_data))], rc_data[,grep("rc_date_ns", names(rc_data))], rc_data$rc_age_std ))
apply(V_RC_data,2,summary)
(d_V_RC<-dim(V_RC_data)[2]) 

