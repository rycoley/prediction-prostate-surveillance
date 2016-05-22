### NAME AND EMAIL REDACTED
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to summarize predictive accuracy of true state predictions from the biopsy and surgery IOP model

library("ROCR")

nsims<-200

auc.res<-mse.res<-mae.res<-vector(length=nsims)

for(j in 1:nsims){
	res<-read.csv(paste0("simulations/IOP/results/eta-true-vs-pred-iop-",j,".csv"))
	auc.res[j] <- performance(prediction(res$predicted, res$true), "auc")@y.values[[1]]
	mse.res[j] <- mean((res$true-res$predicted)^2)
	mae.res[j] <- mean(res$true-res$predicted)}
	

save(auc.res=auc.res, mse.res=mse.res, mae.res=mae.res, file="simulations/IOP/sim-IOP-res-eta.RData")



