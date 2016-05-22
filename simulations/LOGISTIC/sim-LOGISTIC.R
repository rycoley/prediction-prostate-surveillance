### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"

##SIMULATIONS WORKFLOW for comparison to logistic regression
# (1) Define, set random seed
# (2) Simulate data
# (3) Shape data in prep for logistic regression, predictions
# (4) Run logistic regression with patients with surgery/eta known
# (5) Make predictions for patients without surgery, calculate AUC and MSE
# (6) Make in- and out-of-sample predictions for patients with surgery/eta known, calculate AUC and MSE


library(ROCR)

nsims<-200
auc.res<-mse.res<-mae.res<-vector(length=nsims)
auc.ek.in<-mse.ek.in<-mae.ek.in<-vector(length=nsims)
auc.ek.out<-mse.ek.out<-mae.ek.out<-vector(length=nsims)


for(SEED in 1:nsims){

# (1) Set seed
set.seed(SEED)

# (2) Simulate data
source("simulations/simulation-data/simulate-pca-data.R")

bx_sim<-bx_sim[!is.na(bx_sim$bx_here),]


# (3) Shape data in prep for logistic regression, predictions
# Want a dataset with one observation per patient with following variables: most recent bx, # bx w/o RC, age, time in AS, PSAD at diagnosis, slope of log(PSA) 

names(pt_data)

#table(pt_data$rc) #this is final BX
pt_data$num.bx<-vector(length=n)
for(i in 1:n){
	pt_data$num.bx[i]<-sum(bx_sim$subj==i & bx_sim$bx_here==1 & !is.na(bx_sim$bx_here) & bx_sim$rc==0) }

pt_data$age<-vector(length=n)
for(i in 1:n){	
	pt_data$age[i]<-max(c(max(bx_sim$age[bx_sim$subj==i]), max(psa_data$age[psa_data$subj==i]))) }


pt_data$time.fup<-vector(length=n)
for(i in 1:n){
	pt_data$time.fup[i]<-max(c(max(bx_sim$time[bx_sim$subj==i]), max(psa_data$psa_time[psa_data$subj==i]))) }


pt_data$vol_avg[pt_data$vol_avg > quantile(pt_data$vol_avg, 0.975)] <- quantile(pt_data$vol_avg, 0.975)
pt_data$vol_avg[pt_data$vol_avg < max(c(quantile(pt_data$vol_avg, 0.025), 5))] <- max(c(quantile(pt_data$vol_avg, 0.025), 5))

pt_data$psad.dx<-vector(length=n)
for(i in 1:n){
	min.time.i<-min(psa_data$psa_time[psa_data$subj==i])
	pt_data$psad.dx[i]<-exp(psa_data$log_psa[psa_data$subj==i & psa_data$psa_time==min.time.i])/pt_data$vol_avg[i]}

pt_data$psa.slope<-vector(length=n)
for(i in 1:n){
	if(sum(psa_data$subj==i)<2){pt_data$psa.slope[i]<-0}
	else{pt_data$psa.slope[i]<-as.numeric(lm(psa_data$log_psa[psa_data$subj==i]~psa_data$psa_time[psa_data$subj==i])$coef[2]) } }
summary(pt_data$psa.slope)


# (4) Run logistic regression with patients with surgery/eta known

pt.ek<-pt_data[!is.na(pt_data$obs_eta),]
pt.nek<-pt_data[is.na(pt_data$obs_eta),]

mod <- glm(obs_eta ~ rc + num.bx + age + time.fup + psad.dx + psa.slope, family="binomial", data=pt.ek)
#summary(mod)


# (5) Make predictions for patients without surgery, calculate AUC and MSE save

pred.eta<-predict.glm(mod, newdata=pt.nek, type="response")

auc.res[SEED] <- performance(prediction(pred.eta, pt.nek$eta_true), "auc")@y.values[[1]]
mse.res[SEED] <- mean((pt.nek$eta_true-pred.eta)^2)
mae.res[SEED] <- mean(pt.nek$eta_true-pred.eta)


# (6) Make in- and out-of-sample predictions for patients with surgery/eta known, calculate AUC and MSE, save

pred.eta.in<-predict.glm(mod, newdata=pt.ek, type="response")

auc.ek.in[SEED] <- performance(prediction(pred.eta.in, pt.ek$eta_true), "auc")@y.values[[1]]
mse.ek.in[SEED] <- mean((pt.ek$eta_true-pred.eta.in)^2)
mae.ek.in[SEED] <- mean(pt.ek$eta_true-pred.eta.in)

n_ek<-dim(pt.ek)[1]
pred.eta.out<-vector(length=n_ek)
for(i in 1:n_ek){
	modi<-glm(obs_eta ~ rc + num.bx + age + time.fup + psad.dx + psa.slope, family="binomial", data=pt.ek[-i,])
	pred.eta.out[i]<-predict.glm(modi, newdata=pt.ek[i,], type="response")}

auc.ek.out[SEED] <- performance(prediction(pred.eta.out, pt.ek$eta_true), "auc")@y.values[[1]]
mse.ek.out[SEED] <- mean((pt.ek$eta_true-pred.eta.out)^2)
mae.ek.out[SEED] <- mean(pt.ek$eta_true-pred.eta.out)

print(SEED)

}

save(auc.res=auc.res, mse.res=mse.res, mae.ek.in=mae.ek.in, auc.ek.in=auc.ek.in, mse.ek.in=mse.ek.in, mae.ek.in=mae.ek.in, auc.ek.out=auc.ek.out, mse.ek.out=mse.ek.out, mae.ek.out=mae.ek.out, file="simulations/LOGISTIC/results/sim-LOGISTIC-res-eta.RData")


res<-cbind(c(pt.ek$subj, pt.nek$subj), c(pred.eta.out, pred.eta))
colnames(res)<-c("subj", "predicted")
write.csv(res, "simulations/LOGISTIC/results/eta-predictions-LOGISTIC.csv")