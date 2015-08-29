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

#remove last 5 years of data for selected patients
pred.ids<-c(250:261) 
psa.data<-psa.data[!(psa.data$subj%in%pred.ids & psa.data$psa.time>5),] 
data.use<-data.use[!(data.use$subj%in%pred.ids & data.use$time>5),]



#Before call to JAGS, get the data into simple matrices and vectors to send to JAGS
(n<-dim(pt.data)[1]) #there are 1000 patients

#get observed latent class for those with surgery
eta.data<-pt.data$obs.eta
(n_eta_known<-sum(!is.na(eta.data))) #214


#PSA model
(n_obs_psa<-dim(psa.data)[1])
Y<-psa.data$log.psa
subj_psa<-psa.data$subj

#covariates with random effects 
#here, intercept and age (standardized)
Z.data<-as.matrix(cbind(rep(1,n_obs_psa), psa.data$age.std))
(d.Z<-dim(Z.data)[2])

#covariates with only fixed effects
#here, prostate volume (standardized)
X.data<-as.matrix(cbind(psa.data$std.vol))
(d.X<-dim(X.data)[2])


#lmer fit to get starting value for covariance parameter in JAGS
mod.lmer<-lmer(log.psa~ std.vol + (1+ age.std |id), data=psa.data)
(var_vec <- apply(coef(mod.lmer)$id, 2, var)[1:d.Z])
(var_vec <- c(var_vec[2], var_vec[1]))



###bx data

#outcome model (logistic regression for reclassification)
#only use records where a biopsy occurred
rc.data<-data.use[data.use$bx.here==1 & !is.na(data.use$bx.here),] 
(n_rc<-dim(rc.data)[1])
RC<-as.numeric(rc.data$rc) #indicator of reclassificaiton in this interval
subj_rc<-rc.data$subj

#covariate matrix V for pooled logistic regression predicting reclassification
#here, age (standardized), time since diagnosis (ns with 2 df), calendar time
V.RC.data<-as.matrix(cbind(rep(1,n_rc),  rc.data$age.std, rc.data$time, rc.data$time.ns, rc.data$sec.time.std ))
(d.V.RC<-dim(V.RC.data)[2]) #5


