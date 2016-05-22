### NAME AND EMAIL REDACTED
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to run a logistic regression model on the subset of patients with true state observed. 
#OR estimates reported in Appendix Table 8
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.


#get data
pt_data<-read.csv("jhas-analysis/simulation-data/pt-data-sim.csv")
#this contains one record per patient
#patients are ordered so that those with surgery, i.e. eta observed, occur first. 

psa_data<-read.csv("jhas-analysis/simulation-data/psa-data-sim.csv")
#this contains one record per psa observations per patient

data_use<-read.csv("jhas-analysis/simulation-data/bx-data-sim.csv")
#this contains one record per annual interval for each patient until surgery or censoring



(n<-dim(pt_data)[1]) #874
eta.data<-pt_data$obs_eta
(n_eta_known<-sum(!is.na(eta.data)))

# Shape data in prep for logistic regression, predictions
# Want a dataset with one observation per patient with following variables: most recent bx, # bx w/o RC, age, time in AS, PSAD at diagnosis, slope of log(PSA) 
names(pt_data)

# final bx
table(pt_data$rc)

# # bx
pt_data$num.bx<-vector(length=n)
for(i in 1:n){
	pt_data$num.bx[i]<-sum(data_use$subj==pt_data$subj[i] & !is.na(data_use$bx_here) & data_use$bx_here==1 & data_use$time>0 & data_use$rc==0)}
table(pt_data$num.bx)


# age
pt_data$age<-vector(length=n)
for(i in 1:n){
	if(pt_data$num.bx[i]>0){
	pt_data$age[i]<-max( c( max(psa_data$age[psa_data$subj==pt_data$subj[i]]), max(data_use$age[data_use$subj==pt_data$subj[i] & data_use$bx_here==1 & !is.na(data_use$bx_here)]) )) }
	if(pt_data$num.bx[i]==0){
	pt_data$age[i]<-max(psa_data$age[psa_data$subj==pt_data$subj[i]])}}
summary(pt_data$age)

# time in AS
pt_data$time.fup<-vector(length=n)
for(i in 1:n){
	pt_data$time.fup[i]<-max(c(max(data_use$time[data_use$subj==i]), max(psa_data$psa_time[psa_data$subj==i]))) }
summary(pt_data$time.fup)

# PSAD at dx
#get volume at time of diagnosis
dx_vol_mean<-57.4983
dx_vol_sd<-24.87141
pt_data$vol_avg <- (pt_data$vol_std*dx_vol_sd)+dx_vol_mean
summary(pt_data$vol_avg)
#these volumes are unreasonable (i.e. some are negative) bc data is simulated. we will truncate volume
pt_data$vol_avg[pt_data$vol_avg > quantile(pt_data$vol_avg, 0.975)] <- quantile(pt_data$vol_avg, 0.975)
pt_data$vol_avg[pt_data$vol_avg < max(c(quantile(pt_data$vol_avg, 0.025), 5))] <- max(c(quantile(pt_data$vol_avg, 0.025), 5))
summary(pt_data$vol_avg)


pt_data$psad.dx<-vector(length=n)
for(i in 1:n){
	min.time.i<-min(psa_data$psa_time[psa_data$subj==i])
	pt_data$psad.dx[i]<-exp(psa_data$log_psa[psa_data$subj==i & psa_data$psa_time==min.time.i])/pt_data$vol_avg[i]}


for(i in 1:n){
	psa_data$time.since.dx[psa_data$subj==pt_data$subj[i]] <- (psa_data$psa.dt.num[psa_data$subj==pt_data$subj[i]] - pt_data$dx.dt.num[i])/365 }
	

# slope of log-PSA
pt_data$psa.slope<-vector(length=n)
for(i in 1:n){
	if(sum(psa_data$subj==i)<2){pt_data$psa.slope[i]<-0}
	else{pt_data$psa.slope[i]<-as.numeric(lm(psa_data$log_psa[psa_data$subj==i]~psa_data$psa_time[psa_data$subj==i])$coef[2]) } }
summary(pt_data$psa.slope)
hist(pt_data$psa.slope)


# Run logistic regression with patients with surgery/eta known

pt.ek<-pt_data[!is.na(pt_data$obs_eta),]
pt.nek<-pt_data[is.na(pt_data$obs_eta),]

library(splines)



pt.ek$age10<-pt.ek$age*10
pt.ek$psad.dx10<-pt.ek$psad.dx*10
pt.ek$psa.slope10<-pt.ek$psa.slope*10
summary(pt.ek$psa.slope)

mod <- glm(obs_eta ~ rc + age + num.bx + time.fup + psad.dx10 + psa.slope10, family="binomial", data=pt.ek) #  psad.dx
summary(mod)


##odds ratio and 95% interval below, appendix Table 8
exp(mod$coef)

round(exp(as.vector(mod$coef) + qnorm(0.025)*as.vector(sqrt(diag(summary(mod)$cov.scaled)))),4)

round(exp(as.vector(mod$coef) + qnorm(0.975)*as.vector(sqrt(diag(summary(mod)$cov.scaled)))),4)

#These are different than the values reported in the appendix Table 8 because this analysis is based on simulated data, not actual data.




# Make predictions for patients without surgery
pt.nek$age10<-pt.nek$age*10
pt.nek$psad.dx10<-pt.nek$psad.dx*10
pt.nek$psa.slope10<-pt.nek$psa.slope*10

pred.eta<-predict.glm(mod, newdata=pt.nek, type="response")
summary(pred.eta)
hist(pred.eta)

pred.nek<-as.matrix(cbind(pt.nek$subj,pred.eta))
colnames(pred.nek)<-c("subj", "pred.eta")
pred.nek
write.csv(pred.nek, "jhas-analysis/LOGISTIC/results/pred-eta-logistic-nek.csv")


# Make out-of-sample predictions for patients with surgery/eta known, calculate AUC and MSE, save

pred.eta.in<-predict.glm(mod, newdata=pt.ek, type="response")
summary(pred.eta.in)

library(ROCR)

pred.eta.out<-vector(length=n_eta_known)
for(i in 1:n_eta_known){
	modi<-glm(obs_eta ~ rc + num.bx + age + time.fup + psad.dx + psa.slope, family="binomial", data=pt.ek[-i,])
	pred.eta.out[i]<-predict.glm(modi, newdata=pt.ek[i,], type="response")}
hist(pred.eta.out)

(auc.ek.out <- performance(prediction(pred.eta.out, pt.ek$obs_eta), "auc")@y.values[[1]]) #0.78
(mse.ek.out <- mean((pt.ek$obs_eta-pred.eta.out)^2)) #0.18

library(pROC)
(auc.ek.ci<-ci.auc(response=pt.ek$obs_eta, predictor=pred.eta.out, method="bootstrap", boot.n=5000)) #0.72, 0.85


pred.ek<-as.matrix(cbind(pt.ek$subj,pred.eta.out))
colnames(pred.ek)<-c("subj", "pred.eta")
pred.ek
write.csv(pred.ek, "jhas-analysis/LOGISTIC/results/pred-eta-logistic-ek.csv")





###COME BACK HERE TO FINISH FOR PLOTS< TABLES



######

pred.ek<-read.csv("/Users/ryc/Documents/inhealth/prediction-model/logistic-regression/pred-eta-logistic-ek.csv")

pred.eta.out<-pred.ek$pred.eta



roc.logistic<-performance(prediction(pred.eta.out,pt.ek$eta),"tpr","fpr")

cbind(roc.logistic@x.values[[1]][roc.logistic@y.values[[1]]>0.60 & roc.logistic@y.values[[1]]<0.63], roc.logistic@y.values[[1]][roc.logistic@y.values[[1]]>0.60 & roc.logistic@y.values[[1]]<0.63])
#FPR=0.1928 at TRP=0.615
#1-0.1928 = specificity=0.8072, rounds to 81%

