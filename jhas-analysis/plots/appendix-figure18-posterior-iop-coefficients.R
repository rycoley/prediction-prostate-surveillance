### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
### This code is used to pull in results and plot posterior distributions of IOP components (Figure 7 in Coley et al.)
##These results are not the same as those given in the primary paper for the JHAS cohort analysis. That data is not publically available. These results are from a single simulated dataset and may not reflect paper conclusions.


### DEFINE PRIORS
#possible log-OR on risk of surgery

mu.a<-c(-0.7,0,0.7)
mu.b<-c(-0.7,0,0.7)



### PULL IN RESULTS

# BIOPSY- MAIN EFFECT

res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-nu_BX-1.csv",sep=""))
#dim(res)
nu_bx<-res[,dim(res)[2]] #just need to save main effect of eta on P(BX)
for(i in c(2:5)){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-nu_BX-",i,".csv",sep=""))
	nu_bx<-c(nu_bx,res[,dim(res)[2]])} 
nu_bx_main<-matrix(nrow=length(nu_bx), ncol=10)  #matrix for saving posteriors from non-informative and informative priors
nu_bx_main[,1]<-nu_bx


for(j in 1:3){for(k in 1:3){
res<-read.csv(paste0("jhas-analysis/IOP/results/jags-prediction-iop-sa-",mu.a[j],"-",mu.b[k],"-nu_BX-1.csv"))
nu_bx<-res[,dim(res)[2]]
for(i in c(2:5)){
	res<-read.csv(paste0("jhas-analysis/IOP/results/jags-prediction-iop-sa-",mu.a[j],"-",mu.b[k],"-nu_BX-",i,".csv"))
	nu_bx<-c(nu_bx,res[,dim(res)[2]])}
#i<-5
nu_bx_main[,(3*(j-1)+1+k)]<-nu_bx	
}}




# SURGERY- MAIN EFFECT and INTERACTION

res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-omega_SURG-1.csv",sep=""))
#dim(res)
omega_surg<-as.matrix(res[,(dim(res)[2]-1):dim(res)[2]]) #save main effect of eta and its interaction with prev_RC=1 on P(SURG) 
for(i in c(2:5)){
	res<-read.csv(paste("jhas-analysis/IOP/results/jags-prediction-iop-omega_SURG-",i,".csv",sep=""))
	omega_surg<-rbind(omega_surg, res[,(dim(res)[2]-1):dim(res)[2]])} 
omega_surg_main<-omega_surg_int<-matrix(nrow=dim(omega_surg)[1], ncol=10)  
omega_surg_main[,1]<-omega_surg[,1]
omega_surg_int[,1]<-omega_surg[,2]

for(j in 1:3){for(k in 1:3){
res<-read.csv(paste0("jhas-analysis/IOP/results/jags-prediction-iop-sa-",mu.a[j],"-",mu.b[k],"-omega_SURG-1.csv"))
omega_surg<-as.matrix(res[,(dim(res)[2]-1):dim(res)[2]])
for(i in c(2:5)){
	res<-read.csv(paste0("jhas-analysis/IOP/results/jags-prediction-iop-sa-",mu.a[j],"-",mu.b[k],"-omega_SURG-",i,".csv"))
	omega_surg<-rbind(omega_surg,res[,(dim(res)[2]-1):dim(res)[2]])}
#i<-5
omega_surg_main[,(3*(j-1)+1+k)]<-omega_surg[,1]	
omega_surg_int[,(3*(j-1)+1+k)]<-omega_surg[,2]	
}}




pdf("jhas-analysis/plots/appendix-figure18-posterior-iop-coefficients.pdf", width=18, height=8)

layout(matrix(c(1,2,3,rep(4,3)), nrow=2, ncol=3, byrow=TRUE), width=rep(6,3), height=c(7,1))

par(mar=c(5,3,1,1))

plot(density(nu_bx_main[,1]), lwd=3, ylab="", xlab="", main="", ylim=c(0,3.2), cex.axis=1.5)
for(j in 1:3){for(k in 1:3){
		lines(density(nu_bx_main[,(3*(j-1)+1+k)]), lwd=2, lty=c("dotted", "dotdash", "dashed")[j], col=c("blue", "red", "green3")[k])}}
mtext("log-OR of True State on Biopsy",1,line=3.25, cex=1.5)
#mtext("Density",2,line=2.5, cex=1.25)
abline(v=0, lty="twodash", cex=1.5)


plot(density(omega_surg_main[,1]), lwd=3, ylab="", xlab="", main="", ylim=c(0,1.25), cex.axis=1.5)
for(j in 1:3){for(k in 1:3){
		lines(density(omega_surg_main[,(3*(j-1)+1+k)]), lwd=2, lty=c("dotted", "dotdash", "dashed")[j], col=c("blue", "red", "green3")[k])}}
mtext("log-OR of True State on Surgery",1,line=3.25, cex=1.5)
abline(v=0, lty="twodash", cex=1.5)

plot(density(omega_surg_int[,1]), lwd=3, ylab="", xlab="", main="", ylim=c(0,0.85), cex.axis=1.5)
for(j in 1:3){for(k in 1:3){
		lines(density(omega_surg_int[,(3*(j-1)+1+k)]), lwd=2, lty=c("dotted", "dotdash", "dashed")[j], col=c("blue", "red", "green3")[k])}}
mtext("log-OR of True State x Reclassification on Surgery",1,line=3.25, cex=1.5)
abline(v=0, lty="twodash", cex=1.5)

par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
legend(x=0.2, y=0.7, legend="Flat Prior", lwd=3, lty="solid", bty="n", cex=2)
legend(x=0.35, y=0.9, title="Prior Mean OR for Main Effect on Surgery", legend=seq(-1,1,1), lty=c("dotted","dotdash","dashed"), lwd=rep(2,3), bty="n", horiz=T, cex=2)
legend(x=0.65, y=0.9, title="Prior Mean OR for Interaction on Surgery", legend=seq(-1,1,1), col=c("blue", "red", "green3"), lwd=rep(2,3), bty="n",horiz=T, cex=2)

dev.off()
