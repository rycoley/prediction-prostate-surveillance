
### Rebecca Yates Coley rycoley@gmail.com
### Code for "Bayesian Joint Hierarchical Model for Prediction of Latent Health States with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in data for individuals, MCMC results of trajectories, and make individual-level plots, as seen in Figure 3 


### WORKFLOW: Load libraries, define functions, load data for indviduals, load posteriors for individuals, make individual-level predicitons, plot individual data and trajectories


### LOAD LIBRARIES
library("splines")
library("scales")


### DEFINE FUNCTIONS
get.ns.basis<-function(obs.data,knots){
#	knots<-quantile(obs.data,p=c(0.25,0.5,0.75))
	od.k1<- obs.data-knots[1]
	od.k1[od.k1<0]<-0
	od.k2<- obs.data-knots[2]
	od.k2[od.k2<0]<-0
	od.k3<- obs.data-knots[3]
	od.k3[od.k3<0]<-0
	return(as.vector((od.k1^3 - od.k3^3)/(knots[3]-knots[1]) - (od.k2^3 - od.k3^3)/(knots[3]-knots[2])))}

expit<-function(x){return(exp(x)/(1+exp(x)))}


### GET DATA
source("IOP-prep-data-for-jags.R")


pred.ids<-c(250:261)

rc.pred.ids <- vector(length=12) #reclassification for selected patients
for(i in 1:12){rc.pred.ids[i]<-max(data.use$rc[data.use$subj==pred.ids[i]])}

#quick look at data
par(mfrow=c(4,3), mar=c(2,2,1,1))
for(i in 1:12){
	plot(psa.data$log.psa[psa.data$subj==pred.ids[i]]~psa.data$age[psa.data$subj==pred.ids[i]], type="l", ylim=c(-0.5,3.25)) }


#define age range for each patient
plot.min.yrs<-plot.max.yrs<-vector(length=12)
for(i in 1:12){
	plot.min.yrs[i]<-min(c(psa.data$age[psa.data$subj==pred.ids[i]], pt.data$age.dx[pt.data$subj==pred.ids[i]]))
	plot.max.yrs[i]<- max(c(psa.data$age[psa.data$subj==pred.ids[i]], data.use$age[!is.na(data.use$bx.here) & data.use$bx.here==1 & data.use$subj==pred.ids[i]]))}


### GET POSTERIOR RESULTS

#eta
etahat<-read.csv(paste("results/jags-prediction-iop-eta.hat-1.csv",sep=""))
N.eta<-dim(etahat)[2]
etahat<-as.matrix(etahat[,2:N.eta])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-eta.hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:N.eta])}
eta.mat<-etahat[,(pred.ids-n_eta_known)]
(eta.mean<-as.vector(apply(eta.mat,2,mean)))

#psa data

bet<-read.csv(paste("results/jags-prediction-iop-beta-1.csv",sep=""))
bet<-as.matrix(bet[,2])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-beta-",i,".csv",sep=""))
	bet<-c(bet,res[,2])}

b.vec<-read.csv(paste("results/jags-prediction-iop-b.vec-1.csv",sep=""))
b.int<-as.matrix(b.vec[,(pred.ids)+1])
b.slope<-as.matrix(b.vec[,(n+1+pred.ids)])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-b.vec-",i,".csv",sep=""))
	b.int<-rbind(b.int, res[,(pred.ids)+1])
	b.slope<-rbind(b.slope,res[,(n+1+pred.ids)])}


#reclassification 
gam<-read.csv(paste("results/jags-prediction-iop-gamma.RC-1.csv",sep=""))
gam<-as.matrix(gam[,2:(d.V.RC+2)])
for(i in 2:5){
	res<-read.csv(paste("results/jags-prediction-iop-gamma.RC-",i,".csv",sep=""))
	gam<-rbind(gam,res[,2:(d.V.RC+2)])}


B<-dim(gam)[1]


### MAKE INDIVIDUAL-LEVEL PREDICTIONS

## PSA
psa.ages<-seq(50,90,1)
(npp<-length(psa.ages)) #n psa pred
psa.ages.std<- (psa.ages-mean(psa.data$age))/sd(psa.data$age)

psa.post.mat<-array(dim=c(B,npp,12))
psa.mean.mat<-matrix(nrow=npp,ncol=12)
for(i in 1:12){
	for(b in 1:B){
		psa.post.mat[b,,i] <- b.int[b,i] + b.slope[b,i]*psa.ages.std + bet[b]*pt.data$std.vol[pt.data$subj==pred.ids[i]]}
	psa.mean.mat[,i] <- apply(psa.post.mat[,,i],2,mean) }


col.psa.grad<-alpha("blue", seq(0, 1, length = 12))[order(c(1:12),decreasing=TRUE)]
psa.grad.up<-psa.grad.low<-array(dim=c(11,2,12) )
quant.up<-c(0.5,seq(0.525,0.975,0.05))
quant.low<-c(seq(0.025,0.475,0.05),0.5)
quant.low<-quant.low[order(quant.low, decreasing=TRUE)]


psa.ages.i<-vector(length=12)

#only show for times after plot.max.yrs
for(i in 1:12){
	psa.ages.i[i]<-psa.ages[psa.ages==min(psa.ages[psa.ages>plot.max.yrs[i]])]
	psa.grad.up[1,1:2,i]<-psa.grad.low[1,1:2,i]<-c(psa.mean.mat[psa.ages==psa.ages.i[i],i], psa.mean.mat[npp,i])
	for(g in 2:11){	
		psa.grad.up[g,1,i]<-quantile(psa.post.mat[,psa.ages==psa.ages.i[i],i],p=quant.up[g])
		psa.grad.up[g,2,i]<-quantile(psa.post.mat[,npp,i],p=quant.up[g])
		psa.grad.low[g,1,i]<-quantile(psa.post.mat[,psa.ages==psa.ages.i[i],i],p=quant.low[g])
		psa.grad.low[g,2,i]<-quantile(psa.post.mat[,npp,i],p=quant.low[g])
	} }	

##save data for shiny app
psa.data.mats<-list()
for(i in 1:12){
	psa.data.mats[[i]]<-as.matrix(cbind(psa.data$age[psa.data$subj==pred.ids[i]], psa.data$log.psa[psa.data$subj==pred.ids[i]]))}


##RECLASSIFICATION

#may want to use smaller intervals

biopsy.ages<-list()
for(i in 1:12){
	biopsy.ages[[i]]<-seq(max(data.use$age[data.use$subj==pred.ids[i] & !is.na(data.use$bx.here)])+1, pt.data$age.dx[pt.data$subj==pred.ids[i]]+11, 1)}

npb<-vector(length=12)
for(i in 1:12){npb[i]<-length(biopsy.ages[[i]])}

biopsy.ages.std<-biopsy.ages
for(i in 1:12){
	for(j in 1:length(biopsy.ages[[i]])){
		biopsy.ages.std[[i]][j]<- data.use$age.std[ abs(data.use$age-biopsy.ages[[i]][j])==min(abs(data.use$age-biopsy.ages[[i]][j])) ][1]	}}


##time between each year is 0.233929
std.yr<-0.1538462


biopsy.sec.times<-list()
for(i in 1:12){
	biopsy.sec.times[[i]]<-seq(max(data.use$sec.time[data.use$subj==pred.ids[i] &  !is.na(data.use$bx.here)])+std.yr, by=std.yr, length.out=npb[i])}


biopsy.bx.times<-list()
for(i in 1:12){
	biopsy.bx.times[[i]]<-seq(max(data.use$time[data.use$subj==pred.ids[i]  & !is.na(data.use$bx.here)])+1, by=1, length.out=npb[i])}
#biopsy.bx.times

#V.RC.data is intercept, age.std, time, time.ns, sec.time.std

#list of matrices
#1 list entry per person
#each matrix has B rows, npb[i] columns 


biopsy.post.mat<-list()
for(i in 1:12){
	biopsy.post.mat[[i]]<-matrix(nrow=B, ncol=npb[i])
	for(j in 1:npb[i]){
		V.RCi<-c(1, biopsy.ages.std[[i]][j], biopsy.bx.times[[i]][j], data.use$time.ns[data.use$time==biopsy.bx.times[[i]][j]][1], biopsy.sec.times[[i]][j])
		for(b in 1:B){
			biopsy.post.mat[[i]][b,j] <- as.matrix(gam[b,], nrow=6, ncol=1) %*% as.matrix(c(V.RCi,eta.mat[b,i]), nrow=1, ncol=6) 
			}
			}
			}


bx.grad.up<-bx.grad.low<-list()

for(i in 1:12){
	bx.grad.up[[i]]<-bx.grad.low[[i]]<-matrix(nrow=11,ncol=npb[i])
	bx.grad.up[[i]][1,]<-bx.grad.low[[i]][1,]<-apply(biopsy.post.mat[[i]],2,median)
	for(g in 2:11){
		bx.grad.up[[i]][g,]<-apply(biopsy.post.mat[[i]],2,quantile,p=quant.up[g])
		bx.grad.low[[i]][g,]<-apply(biopsy.post.mat[[i]],2,quantile,p=quant.low[g])
	}	}
	
bx.traj.mean.pred<-list()
for(i in 1:12){bx.traj.mean.pred[[i]]<-apply(biopsy.post.mat[[i]],2,mean)}


col.bx.grad<-alpha("green3", seq(0, 1, length = 12))[order(c(1:12),decreasing=TRUE)]



bx.data.mats<-list()
for(i in 1:12){
	bx.x <- c(pt.data$age.dx[pt.data$subj==pred.ids[i]], bx.data$age[bx.data$subj==pred.ids[i]])

	if(rc.pred.ids[i]==0){bx.data.mats[[i]]<-as.matrix(cbind(bx.x, rep(-0.9, length(bx.x)))) }
	if(rc.pred.ids[i]==1){bx.data.mats[[i]]<-as.matrix(cbind(bx.x, c(rep(-0.9, (length(bx.x)-1)),3.4))) } }



rc.data.mats<-list()
for(i in 1:12){
	rc.x <- c(pt.data$age.dx[pt.data$subj==pred.ids[i]], rc.data$age[rc.data$subj==pred.ids[i]])
	if(rc.pred.ids[i]==0){
		rc.data.mats[[i]]<-as.matrix(cbind(rc.x, rep(-0.9, length(rc.x))))}
	if(rc.pred.ids[i]==1){
		rc.data.mats[[i]]<-as.matrix(cbind(rc.x, c(rep(-0.9, (length(rc.x)-1)), 3.4)))} }


put.on.scale<-function(value,y.min, y.max){
	return((y.max-y.min)*value + y.min)}



pdf("plots/figure3-individual-data-and-predictions.pdf", height=10, width=8)

par(mfrow=c(4,3), mar=c(2.5,2,3,2.5))
for(i in 1:12){
	
	index<-order(eta.mean,decreasing=TRUE)[i]

	plot(psa.data$log.psa[psa.data$subj==pred.ids[index]]~psa.data$age[psa.data$subj==pred.ids[index]], col="green", bg="blue", ylim=range(-1, 3.5), xlim=c( (plot.min.yrs[index]-0.25) , (plot.max.yrs[index] + 5) ), xlab="", ylab="", xaxt="n", yaxt="n", pch=21, cex=1.5)
	
	Axis(side=1,labels=TRUE, cex.axis=1.25)
	Axis(side=2, at=c(log(1), log(5), log(10)), labels=c("1", "5", "10"), cex.axis=1)
	
	mtext(paste("P(Aggressive PCa)=",round(mean(eta.mat[,index])*100,1),"%",sep=""), 3)



##psa
	lines(psa.mean.mat[,index]~psa.ages, lwd=2, col="blue")
	polygon(x=c(rep(psa.ages.i[index],2), rep(psa.ages[npp],2)), y= c(psa.grad.up[2,1,index], psa.grad.low[2,1,index], psa.grad.low[2,2,index], psa.grad.up[2,2,index]), col=col.psa.grad[1], border=col.psa.grad[1] )
	for(g in 2:10){
		polygon(x=c(rep(psa.ages.i[index],2), rep(psa.ages[npp],2)), y= c(psa.grad.up[(g+1),1,index], psa.grad.up[g,1,index], psa.grad.up[g,2,index], psa.grad.up[(g+1),2,index]), col=col.psa.grad[g], border=NA )
		polygon(x=c(rep(psa.ages.i[index],2), rep(psa.ages[npp],2)), y= c(psa.grad.low[(g+1),1,index], psa.grad.low[g,1,index], psa.grad.low[g,2,index], psa.grad.low[(g+1),2,index]), col=col.psa.grad[g], border=NA )	
	}
	
	points(psa.data$log.psa[psa.data$subj==pred.ids[index]]~psa.data$age[psa.data$subj==pred.ids[index]], col="green", pch=21, cex=1.5, bg="blue")
	
	
##biopsies

if(rc.pred.ids[index]==0){
	
	polygon(x=c(biopsy.ages[[index]], rev(biopsy.ages[[index]])), y=put.on.scale(value=expit(c(bx.grad.up[[index]][2,], rev(bx.grad.low[[index]][2,]))), y.min=-1, y.max=3.5), col=col.bx.grad[1], border=col.bx.grad[1])
for(g in 2:10){
	polygon(x=c(biopsy.ages[[index]], rev(biopsy.ages[[index]])), y=put.on.scale(value=expit(c(bx.grad.up[[index]][(g+1),], rev(bx.grad.up[[index]][g,]))), y.min=-1, y.max=3.5), col=col.bx.grad[g], border=NA )
	polygon(x=c(biopsy.ages[[index]], rev(biopsy.ages[[index]])), y=put.on.scale(value=expit(c(bx.grad.low[[index]][(g+1),], rev(bx.grad.low[[index]][g,]))), y.min=-1, y.max=3.5), col=col.bx.grad[g], border=NA )
	}
	
	}
	

	bx.x <- c(pt.data$age.dx[pt.data$subj==pred.ids[index]], bx.data$age[bx.data$subj==pred.ids[index]])

	if(rc.pred.ids[index]==0){
		points(rep(-0.9, length(bx.x))~bx.x, pch=6, col="blue", cex=2)}
	if(rc.pred.ids[index]==1){
		points( c(rep(-0.9, (length(bx.x)-1)), 3.4) ~bx.x, pch=6, col="blue", cex=2)	}
	
	rc.x <- c(pt.data$age.dx[pt.data$subj==pred.ids[index]], rc.data$age[rc.data$subj==pred.ids[index]])

	if(rc.pred.ids[index]==0){
		points(rep(-0.9, length(rc.x))~rc.x, pch=25, col="blue", bg="green", cex=2)}
	if(rc.pred.ids[index]==1){
		points(c(rep(-0.9, (length(rc.x)-1)), 3.4)~rc.x, pch=25, col="blue", bg="green", cex=2) }
	
	Axis(side=4, at=c(put.on.scale(c(0.25,0.5,0.75), y.min=-1, y.max=3.5)), labels=c("25%","50%","75%"), cex.axis=1)


	
	}


dev.off()

