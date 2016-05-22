### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in and format all the data needed to run the jags model that accounts for possible biopsy and surgery informative observation processes (IOP).
### Simulated data is ammended to remove last 5 years data for 12 patients, in order to demonstrate making individualized predictions (Figure 3). 
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.



##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.



#get data
pt_data<-read.csv("jhas-analysis/simulation-data/pt-data-sim.csv")
#this contains one record per patient
#patients are ordered so that those with surgery, i.e. eta observed, occur first. 

psa_data<-read.csv("jhas-analysis/simulation-data/psa-data-sim.csv")
#this contains one record per psa observations per patient

data_use<-read.csv("jhas-analysis/simulation-data/bx-data-sim.csv")
#this contains one record per annual interval for each patient until surgery or censoring

#Before call to JAGS, get the data into simple matrices and vectors to send to JAGS
(n<-dim(pt_data)[1]) #there are 1000 patients


#get observed latent class for those with surgery
eta_data<-pt_data$obs_eta
(n_eta_known<-sum(!is.na(eta_data))) 


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

#IOP- biopsy performed model
#remove patients who have already had RC observed but haven't had surgery or been censored
bx_data<-data_use[!is.na(data_use$bx_here),] 
(n_bx<-dim(bx_data)[1])
BX<-as.numeric(bx_data$bx_here) #indicator of biopsy in this interval
subj_bx<-bx_data$subj

#covariate matrix U for logistic regression predicting biopsy
U_BX_data<-as.matrix(cbind(rep(1,n_bx), bx_data[,grep("bx_time_ns",names(bx_data))], bx_data[,grep("bx_date_ns",names(bx_data))], bx_data[,grep("bx_age_ns",names(bx_data))], bx_data[,grep("bx_num_prev_bx_ns",names(bx_data))]  ))
apply(U_BX_data,2,summary)
(d_U_BX<-dim(U_BX_data)[2]) 


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


#logistic regression for surgery (radical retropubic prostatectomy)
#this uses all records in data_use, because patients always at risk of choosing surgery
SURG<-as.numeric(data_use$surg)
(n_surg<-dim(data_use)[1])
subj_surg<-data_use$subj

#covariate matrix W for pooled logistic regression predicting surgery
#here, age (standardized and ns with df=2), time since diagnosis (ns with df=4), calendar time (standardized around Jan 1, 2005 and ns with 3df), number of previous biopsies, previous grade reclassification; interaction with eta and previous RC
W_SURG_data<-as.matrix(cbind(rep(1,n_surg),  data_use[,grep("surg_time_ns",names(data_use))], data_use[,grep("surg_date_ns",names(data_use))], data_use[,grep("surg_age_ns",names(data_use))], data_use$surg_num_prev_bx_std, data_use$prev_G7)) #
apply(W_SURG_data,2,summary)
(d_W_SURG<-dim(W_SURG_data)[2]) 