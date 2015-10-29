
### Rebecca Yates Coley rycoley@gmail.com
### Code for "Bayesian Joint Hierarchical Model for Prediction of Latent Health States with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in data for individuals, MCMC results of trajectories, and make individual-level plots, as seen in Figure 3 


### WORKFLOW: Load libraries, load data characteristics, define functions, load data for indviduals, load posteriors for individuals, make individual-level predicitons, plot individual data and trajectories


### LOAD LIBRARIES
list.of.packages <- c("splines","scales","lme4")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=T)

library("lme4")
library("splines")
library("scales")

### LOAD DATA CHARACTERISTICS, PARAMETERS
load("simulation-data/pars-for-data-sim.RData")


### DEFINE FUNCTIONS
expit<-function(x){return(exp(x)/(1+exp(x)))}


### GET DATA
source("IOP-prep-data-for-jags.R")

rc_pred_ids <- vector(length=12) #reclassification for selected patients
for(i in 1:12){rc_pred_ids[i]<-max(data_use$rc[data_use$subj==pred_ids[i]])}

#quick look at data
par(mfrow=c(4,3), mar=c(2,2,1,1))
for(i in 1:12){
	plot(psa_data$log_psa[psa_data$subj==pred_ids[i]]~psa_data$age[psa_data$subj==pred_ids[i]], type="l", ylim=c(-0.5,3.25)) }


#define age range for each patient
plot_min_yrs<-plot_max_yrs<-vector(length=12)
for(i in 1:12){
	plot_min_yrs[i]<-min(c(psa_data$age[psa_data$subj==pred_ids[i]], pt_data$age_dx[pt_data$subj==pred_ids[i]]))
	plot_max_yrs[i]<- max(c(psa_data$age[psa_data$subj==pred_ids[i]], data_use$age[!is.na(data_use$bx_here) & !is.na(data_use$bx_here) & data_use$subj==pred_ids[i]]))}



### GET POSTERIOR RESULTS

#eta
etahat<-read.csv(paste("results/jags-prediction-iop-eta_hat-1.csv",sep=""))
N_eta<-dim(etahat)[2]
etahat<-as.matrix(etahat[,2:N_eta])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-eta_hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:N_eta])}
eta_mat<-etahat[,(pred_ids-n_eta_known)]
(eta_mean<-as.vector(apply(eta_mat,2,mean)))

#psa data

bet<-read.csv(paste("results/jags-prediction-iop-beta-1.csv",sep=""))
bet<-as.matrix(bet[,2])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-beta-",i,".csv",sep=""))
	bet<-c(bet,res[,2])}

b_vec<-read.csv(paste("results/jags-prediction-iop-b_vec-1-reduced.csv",sep=""))
b_int<-as.matrix(b_vec[,(pred_ids)+1])
b_slope<-as.matrix(b_vec[,(n+1+pred_ids)])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-b_vec-",i,"-reduced.csv",sep=""))
	b_int<-rbind(b_int, res[,(pred_ids)+1])
	b_slope<-rbind(b_slope,res[,(n+1+pred_ids)])}


#reclassification 
gam<-read.csv(paste("results/jags-prediction-iop-gamma_RC-1.csv",sep=""))
gam<-as.matrix(gam[,2:(d_V_RC+2)])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-gamma_RC-",i,".csv",sep=""))
	gam<-rbind(gam,res[,2:(d_V_RC+2)])}


B<-dim(gam)[1]


### MAKE INDIVIDUAL-LEVEL PREDICTIONS

## PSA
psa_ages<-seq(50,90,1)
(npp<-length(psa_ages)) #n psa pred
psa_ages_std<- (psa_ages-psa_age_mean)/psa_age_sd

psa_post_mat<-array(dim=c(B,npp,12))
psa_mean_mat<-matrix(nrow=npp,ncol=12)
for(i in 1:12){
	for(b in 1:B){
		psa_post_mat[b,,i] <- b_int[b,i] + b_slope[b,i]*psa_ages_std + bet[b]*pt_data$vol_std[pt_data$subj==pred_ids[i]]}
	psa_mean_mat[,i] <- apply(psa_post_mat[,,i],2,mean) }


col_psa_grad<-alpha("blue", seq(0, 1, length = 12))[order(c(1:12),decreasing=TRUE)]
psa_grad_up<-psa_grad_low<-array(dim=c(11,2,12) )
quant_up<-c(0.5,seq(0.525,0.975,0.05))
quant_low<-c(seq(0.025,0.475,0.05),0.5)
quant_low<-quant_low[order(quant_low, decreasing=TRUE)]


psa_ages_i<-vector(length=12)

#only show for times after plot_max_yrs
for(i in 1:12){
	psa_ages_i[i]<-psa_ages[psa_ages==min(psa_ages[psa_ages>plot_max_yrs[i]])]
	psa_grad_up[1,1:2,i]<-psa_grad_low[1,1:2,i]<-c(psa_mean_mat[psa_ages==psa_ages_i[i],i], psa_mean_mat[npp,i])
	for(g in 2:11){	
		psa_grad_up[g,1,i]<-quantile(psa_post_mat[,psa_ages==psa_ages_i[i],i],p=quant_up[g])
		psa_grad_up[g,2,i]<-quantile(psa_post_mat[,npp,i],p=quant_up[g])
		psa_grad_low[g,1,i]<-quantile(psa_post_mat[,psa_ages==psa_ages_i[i],i],p=quant_low[g])
		psa_grad_low[g,2,i]<-quantile(psa_post_mat[,npp,i],p=quant_low[g])
	} }	


##RECLASSIFICATION

#may want to use smaller intervals

biopsy_ages<-list()
for(i in 1:12){
	biopsy_ages[[i]]<-seq(max(data_use$age[data_use$subj==pred_ids[i] & !is.na(data_use$bx_here)])+1, pt_data$age_dx[pt_data$subj==pred_ids[i]]+11.5, 1)}


npb<-vector(length=12)
for(i in 1:12){npb[i]<-length(biopsy_ages[[i]])}

biopsy_ages_std<-list()
for(i in 1:12){
	biopsy_ages_std[[i]]<- (biopsy_ages[[i]]-rc_age_mean)/rc_age_sd	}



biopsy_dates<-list()
for(i in 1:12){
	biopsy_dates[[i]]<-seq(max(data_use$date[data_use$subj==pred_ids[i] &  !is.na(data_use$bx_here)]) + 365, pt_data$date_dx[pt_data$subj==pred_ids[i]] + (11.5*365), 365)}

biopsy_dates_ns<-list()
for(i in 1:12){
	biopsy_dates_ns[[i]] <- ns(biopsy_dates[[i]], knots=rc_date_knots, Boundary.knots=rc_date_bknots)}


biopsy_times<-list()
for(i in 1:12){
	biopsy_times[[i]]<-seq(max(data_use$time[data_use$subj==pred_ids[i] &  !is.na(data_use$bx_here)]) + 1, 11, 1)}

biopsy_times_ns<-list()
for(i in 1:12){
	biopsy_times_ns[[i]] <- ns(biopsy_times[[i]], knots=rc_time_knots, Boundary.knots=rc_time_bknots)}

#V_rc_data is intercept, time_ns, date_ns, age_std, eta

#list of matrices
#1 list entry per person
#each matrix has B rows, npb[i] columns 

biopsy_post_mat<-list()
for(i in 1:12){
	biopsy_post_mat[[i]]<-matrix(nrow=B, ncol=npb[i])
	for(j in 1:npb[i]){
		V_RCi<-c(1, biopsy_times_ns[[i]][j,], biopsy_dates_ns[[i]][j,], biopsy_ages_std[[i]][j])
		for(b in 1:B){
			biopsy_post_mat[[i]][b,j] <- as.matrix(gam[b,], nrow=6, ncol=1) %*% as.matrix(c(V_RCi,eta_mat[b,i]), nrow=1, ncol=6) }	}	}



bx.grad.up<-bx.grad.low<-list()

for(i in 1:12){
	bx.grad.up[[i]]<-bx.grad.low[[i]]<-matrix(nrow=11,ncol=npb[i])
	bx.grad.up[[i]][1,]<-bx.grad.low[[i]][1,]<-apply(biopsy_post_mat[[i]],2,median)
	for(g in 2:11){
		bx.grad.up[[i]][g,]<-apply(biopsy_post_mat[[i]],2,quantile,p=quant_up[g])
		bx.grad.low[[i]][g,]<-apply(biopsy_post_mat[[i]],2,quantile,p=quant_low[g])
	}	}
	
bx.traj.mean.pred<-list()
for(i in 1:12){bx.traj.mean.pred[[i]]<-apply(biopsy_post_mat[[i]],2,mean)}


col.bx.grad<-alpha("green3", seq(0, 1, length = 12))[order(c(1:12),decreasing=TRUE)]


bx_data.mats<-list()
for(i in 1:12){
	bx.x <- c(pt_data$age_dx[pt_data$subj==pred_ids[i]], bx_data$age[bx_data$subj==pred_ids[i]])

	if(rc_pred_ids[i]==0){bx_data.mats[[i]]<-as.matrix(cbind(bx.x, rep(-0.9, length(bx.x)))) }
	if(rc_pred_ids[i]==1){bx_data.mats[[i]]<-as.matrix(cbind(bx.x, c(rep(-0.9, (length(bx.x)-1)),3.4))) } }



rc_data.mats<-list()
for(i in 1:12){
	rc.x <- c(pt_data$age_dx[pt_data$subj==pred_ids[i]], rc_data$age[rc_data$subj==pred_ids[i]])
	if(rc_pred_ids[i]==0){
		rc_data.mats[[i]]<-as.matrix(cbind(rc.x, rep(-0.9, length(rc.x))))}
	if(rc_pred_ids[i]==1){
		rc_data.mats[[i]]<-as.matrix(cbind(rc.x, c(rep(-0.9, (length(rc.x)-1)), 3.4)))} }



put.on.scale<-function(value,y.min, y.max){
	return((y.max-y.min)*value + y.min)}



pdf("plots/figure3-individual-data-and-predictions.pdf", height=10, width=8)

par(mfrow=c(4,3), mar=c(2.5,2,3,2.5))
for(i in 1:12){
	
	index<-order(eta_mean,decreasing=TRUE)[i]

	plot(psa_data$log_psa[psa_data$subj==pred_ids[index]]~psa_data$age[psa_data$subj==pred_ids[index]], col="green", bg="blue", ylim=range(-1, 3.5), xlim=c( (plot_min_yrs[index]-0.25) , (plot_max_yrs[index] + 5) ), xlab="", ylab="", xaxt="n", yaxt="n", pch=21, cex=1.5)
	
	Axis(side=1,labels=TRUE, cex.axis=1.25)
	Axis(side=2, at=c(log(1), log(5), log(10)), labels=c("1", "5", "10"), cex.axis=1)
	
	mtext(paste("P(Aggressive PCa)=",round(mean(eta_mat[,index])*100,1),"%",sep=""), 3)


##psa
	lines(psa_mean_mat[,index]~psa_ages, lwd=2, col="blue")
	polygon(x=c(rep(psa_ages_i[index],2), rep(psa_ages[npp],2)), y= c(psa_grad_up[2,1,index], psa_grad_low[2,1,index], psa_grad_low[2,2,index], psa_grad_up[2,2,index]), col=col_psa_grad[1], border=col_psa_grad[1] )
	for(g in 2:10){
		polygon(x=c(rep(psa_ages_i[index],2), rep(psa_ages[npp],2)), y= c(psa_grad_up[(g+1),1,index], psa_grad_up[g,1,index], psa_grad_up[g,2,index], psa_grad_up[(g+1),2,index]), col=col_psa_grad[g], border=NA )
		polygon(x=c(rep(psa_ages_i[index],2), rep(psa_ages[npp],2)), y= c(psa_grad_low[(g+1),1,index], psa_grad_low[g,1,index], psa_grad_low[g,2,index], psa_grad_low[(g+1),2,index]), col=col_psa_grad[g], border=NA )	
	}
	
	points(psa_data$log_psa[psa_data$subj==pred_ids[index]]~psa_data$age[psa_data$subj==pred_ids[index]], col="green", pch=21, cex=1.5, bg="blue")
	
	
##biopsies
if(rc_pred_ids[index]==0){
	
	polygon(x=c(biopsy_ages[[index]], rev(biopsy_ages[[index]])), y=put.on.scale(value=expit(c(bx.grad.up[[index]][2,], rev(bx.grad.low[[index]][2,]))), y.min=-1, y.max=3.5), col=col.bx.grad[1], border=col.bx.grad[1])
for(g in 2:10){
	polygon(x=c(biopsy_ages[[index]], rev(biopsy_ages[[index]])), y=put.on.scale(value=expit(c(bx.grad.up[[index]][(g+1),], rev(bx.grad.up[[index]][g,]))), y.min=-1, y.max=3.5), col=col.bx.grad[g], border=NA )
	polygon(x=c(biopsy_ages[[index]], rev(biopsy_ages[[index]])), y=put.on.scale(value=expit(c(bx.grad.low[[index]][(g+1),], rev(bx.grad.low[[index]][g,]))), y.min=-1, y.max=3.5), col=col.bx.grad[g], border=NA )
	}
	
	}
	

	bx.x <- c(pt_data$age_dx[pt_data$subj==pred_ids[index]], bx_data$age[bx_data$subj==pred_ids[index]])

	if(rc_pred_ids[index]==0){
		points(rep(-0.9, length(bx.x))~bx.x, pch=6, col="blue", cex=2)}
	if(rc_pred_ids[index]==1){
		points( c(rep(-0.9, (length(bx.x)-1)), 3.4) ~bx.x, pch=6, col="blue", cex=2)	}
	
	rc.x <- c(pt_data$age_dx[pt_data$subj==pred_ids[index]], rc_data$age[rc_data$subj==pred_ids[index]])

	if(rc_pred_ids[index]==0){
		points(rep(-0.9, length(rc.x))~rc.x, pch=25, col="blue", bg="green", cex=2)}
	if(rc_pred_ids[index]==1){
		points(c(rep(-0.9, (length(rc.x)-1)), 3.4)~rc.x, pch=25, col="blue", bg="green", cex=2) }
	
	Axis(side=4, at=c(put.on.scale(c(0.25,0.5,0.75), y.min=-1, y.max=3.5)), labels=c("25%","50%","75%"), cex.axis=1)


	
	}

dev.off()


##Plot labels can be added in Powerpoint/Keynote or Adobe Illustrator

