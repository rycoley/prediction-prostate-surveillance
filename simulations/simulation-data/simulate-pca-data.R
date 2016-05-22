### Rebecca Yates Coley rycoley@gmail.com
### Simulation code for "Bayesian Joint Hierarchical Model for Prediciton of Latent Health States with Application to Active Surveillance of Prostate Cancer"
### This code is used to simulate a single dataset. 

### SIMULATE DATA FOR PROSTATE CANCER PREDICTION MODEL
### WORKFLOW: 
# (1) Load necessary libraries, define functions
# (2) Load characteristcs of data (mean, sd, knots for natural splines) and model parameters
# (3) Define sample size and unique ids
# (4) Simulate data for diagnosis (dx) and true state (eta)
# (5) Simulate data for biopsies (bx) received, biopsy results (reclassification- rc), and surgery (surg)
# (6) Simulate PSA data
# (7) Order all datasets based on observed eta
# (8) Save all datasets

### (1) Load libraries, define functions
library(MASS)
library(splines)

#inverse logit function
expit<-function(x){return(exp(x)/(1+exp(x)))}


### (2) Load data characteristics, parameters
load("simulations/simulation-data/pars-for-data-sim.RData")

### (3) Define sample size, unqiue ids
n <- 1000
id <- c(1:1000)

### (4) Simulate dx and eta data
#set.seed(SEED)

age_dx <- rnorm(n, mean=dx_age_mean, sd=dx_age_sd) 
date_dx <- rnorm(n, mean=dx_dt_mean, sd=dx_dt_sd)  

pt_data<-as.data.frame(cbind(id, age_dx, date_dx))
names(pt_data) <- c("id","age_dx","date_dx")

#latent class (eta)
pt_data$eta_true <- eta_true <- rbinom(n,1,p_eta)
#table(pt_data$eta_true)



### (5) Simulate data for biopsies (bx) received, biopsy results (reclassification- rc), and surgery (surg)

#define annual intervals
times<-seq(1,10,1)

ids<-c(rep(1,10))
for(i in 2:n){ids<-c(ids,rep(i,10))}

bx_sim<-as.data.frame(cbind(ids, rep(times,n)))
names(bx_sim)<-c("id","time")

(N<-dim(bx_sim)[1])

#bring in true state data
bx_sim$eta<-rep(0,N)
for(i in 1:n){
	bx_sim$eta[bx_sim$id==i]<-pt_data$eta_true[pt_data$id==i]}

#define age and date at each time
bx_sim$age<-bx_sim$date<-rep(0,N)
for(i in 1:n){
	bx_sim$date[bx_sim$id==i]<-pt_data$date_dx[i] + 365*(bx_sim$time[bx_sim$id==1] + 0.5)
	bx_sim$age[bx_sim$id==i]<-pt_data$age_dx[i] + bx_sim$time[bx_sim$id==1] + 0.5}


#natural spline (ns) representations for time, date, and age for each outcome (bx, rc, surg)
bx_sim$bx_time_ns <- ns(bx_sim$time, knots=bx_time_knots, Boundary.knots=bx_time_bknots)
bx_sim$bx_date_ns <- ns(bx_sim$date, knots=bx_date_knots, Boundary.knots=bx_date_bknots)
bx_sim$bx_age_ns <- ns(bx_sim$age, knots=bx_age_knots, Boundary.knots=bx_age_bknots)


bx_sim$rc_time_ns <- ns(bx_sim$time, knots=rc_time_knots, Boundary.knots=rc_time_bknots)
bx_sim$rc_date_ns <- ns(bx_sim$date, knots=rc_date_knots, Boundary.knots=rc_date_bknots)

bx_sim$surg_time_ns <- ns(bx_sim$time, knots=surg_time_knots, Boundary.knots=surg_time_bknots)
bx_sim$surg_date_ns <- ns(bx_sim$date, knots=surg_date_knots, Boundary.knots=surg_date_bknots)
bx_sim$surg_age_ns <- ns(bx_sim$age, knots=surg_age_knots, Boundary.knots=surg_age_bknots)

#standardize age for RC outcome
bx_sim$rc_age_std <- scale(bx_sim$age, center=rc_age_mean, scale=rc_age_sd)

#will use this indicator variable to remove records after censoring
bx_sim$rm<-rep(0,N)


##simulate outcome data for biopsies received (bx)
bx_sim$bx_here<-rep(0,N)
bx_sim$num_prev_bx<-rep(1,N)

bx_sub<-bx_sim[bx_sim$time==1,]
(n_bx<-dim(bx_sub)[1])

U_BX <- as.matrix(cbind( rep(1,n_bx), bx_sub$bx_time_ns, bx_sub$bx_date_ns, bx_sub$bx_age_ns, ns(bx_sub$num_prev_bx, knots=bx_num_prev_bx_knots, Boundary.knots=bx_num_prev_bx_bknots), bx_sub$eta ))
#summary(as.vector(expit(U_BX%*%nu_bx)))

bx_sim$bx_here[bx_sim$time==1]<-rbinom(n,1,as.vector(expit(U_BX%*%nu_bx)))
#table(bx_sim$bx_here[bx_sim$time==1])

for(j in 2:10){
	for(i in 1:n){
		bx_sim$num_prev_bx[bx_sim$id==i & bx_sim$time==j]<-sum(bx_sim$bx_here[bx_sim$id==i & bx_sim$time<j]) + 1}

	bx_sub<-bx_sim[bx_sim$time==j,]
	(n_bx<-dim(bx_sub)[1])
	U_BX<-as.matrix(cbind( rep(1,n_bx), bx_sub$bx_time_ns, bx_sub$bx_date_ns, bx_sub$bx_age_ns, ns(bx_sub$num_prev_bx, knots=bx_num_prev_bx_knots, Boundary.knots=bx_num_prev_bx_bknots), bx_sub$eta ))
	
	bx_sim$bx_here[bx_sim$time==j]<-rbinom(n,1,as.vector(expit(U_BX%*%nu_bx)))}
#table(bx_sim$bx_here)	
#table(bx_sim$num_prev_bx)

bx_sim$bx_num_prev_bx_ns <- ns(bx_sim$num_prev_bx, knots=bx_num_prev_bx_knots, Boundary.knots=bx_num_prev_bx_bknots)

### simulate outcomes for biopsies- reclassification (rc)
bx_sim$rc<-bx_sim$prev_G7<-rep(0,N)
rc_sub<-bx_sim[bx_sim$bx_here==1,]
(n_rc<-dim(rc_sub)[1])

V_RC<-as.matrix(cbind(rep(1,n_rc), rc_sub$rc_time_ns, rc_sub$rc_date_ns, rc_sub$rc_age_std, rc_sub$eta))

bx_sim$rc[bx_sim$bx_here==1]<-rbinom(n_rc,1,as.vector(expit(V_RC%*%gamma_rc)))

for(i in 1:n){
	if(sum(bx_sim$rc[bx_sim$id==i]==1)>0){
		rc_time<-min(bx_sim$time[bx_sim$rc==1 & bx_sim$id==i])
		bx_sim$rc[bx_sim$id==i & bx_sim$time>rc_time]<-0
		bx_sim$bx_here[bx_sim$id==i & bx_sim$time>rc_time]<-0
		bx_sim$num_prev_bx[bx_sim$id==i & bx_sim$time>rc_time]<-(bx_sim$num_prev_bx[bx_sim$id==i & bx_sim$time==rc_time] + 1)
		bx_sim$prev_G7[bx_sim$id==i & bx_sim$time>=rc_time]<-1
		bx_sim$rm[bx_sim$id==i & bx_sim$time>(rc_time+2)]<-1}}

#table(bx_sim$rc)

# surgery
bx_sim$surg<-rep(0,N) 
bx_sim$surg_num_prev_bx_std <- scale(x=(bx_sim$num_prev_bx + bx_sim$bx_here), center=surg_num_prev_bx_mean, scale=surg_num_prev_bx_sd)
#omega_surg

omega_surg[1] <- -5 #adjusting a bit to get higher number of surgeries (closer to true data)


W_SURG<-as.matrix(cbind(rep(1,N), bx_sim$surg_time_ns, bx_sim$surg_date_ns, bx_sim$surg_age_ns, bx_sim$surg_num_prev_bx_std, bx_sim$prev_G7, bx_sim$eta, (bx_sim$prev_G7*bx_sim$eta)))
bx_sim$surg<-rbinom(N,1,as.vector(expit(W_SURG%*%omega_surg)))

#messes up design matrices to delete columns earlier
bx_sim<-bx_sim[bx_sim$rm==0,]
(N<-dim(bx_sim)[1])

pt_data$rc<-pt_data$surg<-rep(0,n)

for(i in 1:n){
	if(sum(bx_sim$surg[bx_sim$id==i])>0){
		surg_time<-min(bx_sim$time[bx_sim$id==i & bx_sim$surg==1])
		bx_sim$rm[bx_sim$id==i & bx_sim$time>surg_time]<-1	
		pt_data$surg[pt_data$id==i]<-1}	}
#table(pt_data$surg)


bx_sim<-bx_sim[bx_sim$rm==0,]
(N<-dim(bx_sim)[1])

for(i in 1:n){
	pt_data$rc[i]<-sum(bx_sim$rc[bx_sim$id==pt_data$id[i]])}
#table(pt_data$rc) 

pt_data$obs_eta<-rep(NA,n)
pt_data$obs_eta[pt_data$surg==1]<-pt_data$eta_true[pt_data$surg==1]
#table(pt_data$obs_eta)

for(i in 1:n){
	if(max(bx_sim$rc[bx_sim$id==i])==1){
		rc_time<-bx_sim$time[bx_sim$rc==1 & bx_sim$id==i]
		bx_sim$bx_here[bx_sim$id==i & bx_sim$time>rc_time]<-NA	} }

#table(bx_sim$bx_here)
#summary(bx_sim$bx_here)


#length(unique(bx_sim$id[bx_sim$rc==1])) #205
#length(bx_sim$id[bx_sim$rc==1]) #205


##psa data

psa_time<-seq(-1, max(bx_sim$time[bx_sim$id==1 & !is.na(bx_sim$bx_here) ]),0.5)
psa_id<-rep(1, length(psa_time))
for(i in 2:n){
	psa_add<-seq(-1, max(bx_sim$time[bx_sim$id==i & !is.na(bx_sim$bx_here)]), 0.5)
	psa_time<-c(psa_time,psa_add)
	psa_id<-c(psa_id, rep(i, length(psa_add)))}	
	
psa_data<-as.data.frame(cbind(psa_id, psa_time))
names(psa_data)<-c("id","psa_time")
(n_obs_psa<-dim(psa_data)[1])


psa_data$psa_time<-psa_data$psa_time + runif(n_obs_psa, min=-0.25, max=0.25)
psa_data$age<-vector(length=n_obs_psa)
for(j in 1:n_obs_psa){
	psa_data$age[j] <- psa_data$psa_time[j] + pt_data$age[pt_data$id==psa_data$id[j]]}
	

psa_data$age_std<-(psa_data$age-psa_age_mean)/psa_age_sd

pt_data$vol_avg<-rnorm(n,dx_vol_mean, dx_vol_sd)
pt_data$vol_std<-(pt_data$vol_avg-dx_vol_mean)/dx_vol_sd

psa_data$vol_std<-vector(length=n_obs_psa)
for(i in 1:n){
	psa_data$vol_std[psa_data$id==i] <- pt_data$vol_std[i]}

b.vec <- matrix(nrow=n, ncol=2)
for(i in 1:n){
	b.vec[i,] <- mvrnorm(n=1, mu=mu_mat[,(pt_data$eta_true[pt_data$id==i]+1)], Sigma=Sigma)}


psa_data$log_psa <- vector(length=n_obs_psa)
for(j in 1:n_obs_psa){
	lin_pred <- NULL
	lin_pred <- sum(b.vec[pt_data$id==psa_data$id[j],] * c(1, psa_data$age_std[j])) + beta[1]*psa_data$vol_std[j]
	psa_data$log_psa[j] <- rnorm(1, mean=lin_pred, sd=sigma_res)}
#summary(psa_data$log_psa)



#get ordered subject variable
pt_data<-pt_data[order(pt_data$obs_eta),]
pt_data$subj<-c(1:n)
psa_data$subj<-rep(0,n_obs_psa)
for(i in 1:n){psa_data$subj[psa_data$id==pt_data$id[i]]<-pt_data$subj[i]}
bx_sim$subj<-rep(0,N)
for(i in 1:n){bx_sim$subj[bx_sim$id==pt_data$id[i]]<-pt_data$subj[i]}









