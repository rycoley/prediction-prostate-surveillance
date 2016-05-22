### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in and format all the data needed to run the jags model that accounts for a possible surgery informative observation process (IOP) when doing CROSS-VALIDATION, meaning that some data on observed latent state will be excluded
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


#get observed latent class for those with surgery, excluding subj for CV

(eta_true<-pt_data$obs_eta[pt_data$subj==to_mask]) #save

pt_data$obs_eta[pt_data$subj==to_mask]<-NA

pt_data<-pt_data[order(pt_data$obs_eta),]
eta_data<-pt_data$obs_eta
(n_eta_known<-sum(!is.na(eta_data)))

#re-number subjects
pt_data$subj2<-c(1:n)
#(subj.test<-pt.data$subj2[pt.data$subj==to.mask]-n_eta_known) #will be 1


psa_data$subj2<-vector(length=dim(psa_data)[1])
for(i in 1:n){psa_data$subj2[psa_data$subj==i]<-pt_data$subj2[pt_data$subj==i]}
data_use$subj2<-vector(length=dim(data_use)[1])
for(i in 1:n){data_use$subj2[data_use$subj==i]<-pt_data$subj2[pt_data$subj==i]}


#PSA model
(n_obs_psa<-dim(psa_data)[1])
Y<-psa_data$log_psa
subj_psa<-psa_data$subj2

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
subj_rc<-rc_data$subj2

#covariate matrix V for pooled logistic regression predicting reclassification
#here, age (standardized), time since diagnosis (ns with 2 df), calendar time
V_RC_data<-as.matrix(cbind(rep(1,n_rc),  rc_data[,grep("rc_time_ns", names(rc_data))], rc_data[,grep("rc_date_ns", names(rc_data))], rc_data$rc_age_std ))
apply(V_RC_data,2,summary)
(d_V_RC<-dim(V_RC_data)[2]) 




#logistic regression for surgery (radical retropubic prostatectomy)
#this uses all records in data.use, because patients always at risk of choosing surgery
SURG<-as.numeric(data_use$surg)
(n_surg<-dim(data_use)[1])
subj_surg<-data_use$subj2


#covariate matrix W for pooled logistic regression predicting surgery
#here, age (standardized and ns with df=2), time since diagnosis (ns with df=4), calendar time (standardized around Jan 1, 2005 and ns with 3df), number of previous biopsies, previous grade reclassification; interaction with eta and previous RC
W_SURG_data<-as.matrix(cbind(rep(1,n_surg),  data_use[,grep("surg_time_ns",names(data_use))], data_use[,grep("surg_date_ns",names(data_use))], data_use[,grep("surg_age_ns",names(data_use))], data_use$surg_num_prev_bx_std, data_use$prev_G7)) #
apply(W_SURG_data,2,summary)
(d_W_SURG<-dim(W_SURG_data)[2]) #

