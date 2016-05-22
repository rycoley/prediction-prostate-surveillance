### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in very large results files and reduce number of digits saved so that (1) the file is small enough to post on github and (2) producing figures is faster
## Note: only reduced result files are posted on github. JAGS scripts on github can be used to produce larger files
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.


### WORKFLOW: load result files, reduce number of significant digits, save smaller files

#use all seeds for random effects
for(seed in 1:5){
b_vec<-read.csv(paste0("results/jags-prediction-iop-b_vec-",seed,".csv"))
write.csv(round(b_vec[,2:dim(b_vec)[2]],3), paste0("results/jags-prediction-iop-b_vec-",seed,"-reduced.csv") )}




#only used seed=1 posterior fitted probabilities for biopsy, reclassication, and surgery
p_bx<-read.csv("results/jags-prediction-iop-p_bx-1.csv")
write.csv(round(p_bx[,2:dim(p_bx)[2]],3), "results/jags-prediction-iop-p_bx-1-reduced.csv" )

p_rc<-read.csv("results/jags-prediction-iop-p_rc-1.csv")
write.csv(round(p_rc[,2:dim(p_rc)[2]],3), "results/jags-prediction-iop-p_rc-1-reduced.csv" )


p_surg<-read.csv("results/jags-prediction-iop-p_surg-1.csv")
write.csv(round(p_surg[,2:dim(p_surg)[2]],3), "results/jags-prediction-iop-p_surg-1-reduced.csv" )

