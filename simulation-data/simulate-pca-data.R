### This code simulates data in order to demonstrate the joint modeling approach outlined in Coley et al (2015) as well as the importance sampling algorithm outlined in Fisher et al (2015).
### Data is generated using posterior estimates from fitting the joint model to data from the Johns Hopkins Active Surveillance cohort.

### LOAD PACKAGES
library(MASS)
library(splines)

### SET SEED
set.seed(1)


### DEFINE FUNCTIONS
expit<-function(x){return(exp(x)/(1+exp(x)))}

#function to get natural spline basis. (Just an alternate definition. See Ch 11 of Wakefield (2013))
get.ns.basis<-function(obs.data,knots){
	od.k1<- obs.data-knots[1]
	od.k1[od.k1<0]<-0
	od.k2<- obs.data-knots[2]
	od.k2[od.k2<0]<-0
	od.k3<- obs.data-knots[3]
	od.k3[od.k3<0]<-0
	return(as.vector((od.k1^3 - od.k3^3)/(knots[3]-knots[1]) - (od.k2^3 - od.k3^3)/(knots[3]-knots[2])))}


### DEFINE PARAMETER VALUES (similar to posterior estimates)
p_eta <- 0.22

mu_int <- c(1.36, 1.61)
mu_slope <- c(0.26,0.51)
mu.mat <- as.matrix(rbind(mu_int, mu_slope))

Sigma <- matrix(c(0.55^2, 0.04, 0.04, 0.4^2), nrow=2, ncol=2)
sigma_res <- 0.3

beta <- c(0.31)

nu.bx<-c(-1, 0.2, -0.8, 2.3, -1.5, 2.3, -2.7, 0.1, 0.9, 1.2, 2.5, -2.5, -0.5)

gam.rc<-c(-2.5, 0.6, -0.4, 0.1, 0.3, 1.7)

omega.surg<-c(-5.4, -0.4, -1.6, 1.8, 1, 6.5, 2.7, 0.8, -2, -0.8, -0.2, 1.1, 0.6, 2.5)

#from data, for design matrices for biopsy data
mean.age.bx<-69.4
sd.age.bx<-6.5
knots.age.bx<- c(-0.5, 0.1, 0.7)

knots.time.bx<- c(1.3, 3.2, 5.8)

mean.sec.time.bx<-4.5
sd.sec.time.bx<-4.1


### SIMULATE DATA
n <- 1000
id <- c(1:1000)

ages.dx <- rnorm(n, mean=65.5, sd=5.5) #from data
sec.time.dx <- rnorm(n, mean=1.6, sd=4.3)  #secular time, in relation to 2005


pt.data<-as.data.frame(cbind(id,ages.dx, sec.time.dx))
names(pt.data) <- c("id","age.dx","sec.time.dx")

#latent class
pt.data$eta.true <- eta.true <- rbinom(n,1,p_eta)
table(pt.data$eta.true)



### all biopsy data

times<-seq(1,10,1)

ids<-c(rep(1,10))
for(i in 2:n){ids<-c(ids,rep(i,10))}

bx.sim<-as.data.frame(cbind(ids, rep(times,n)))
names(bx.sim)<-c("id","time")

(N<-dim(bx.sim)[1])

bx.sim$eta<-rep(0,N)
for(i in 1:n){
	bx.sim$eta[bx.sim$id==i]<-pt.data$eta.true[pt.data$id==i]}

bx.sim$age<-bx.sim$sec.time<-rep(0,N)
for(i in 1:n){
	bx.sim$age[bx.sim$id==i]<-pt.data$age.dx[i] + bx.sim$time[bx.sim$id==1] + 0.5
	bx.sim$sec.time[bx.sim$id==i]<-pt.data$sec.time.dx[i] + bx.sim$time[bx.sim$id==1] + 0.5}

bx.sim$age.std<-(bx.sim$age-mean.age.bx)/sd.age.bx
bx.sim$age.ns<-get.ns.basis(bx.sim$age.std, knots.age.bx)

bx.sim$sec.time.std<-(bx.sim$sec.time-mean.sec.time.bx)/sd.sec.time.bx

bx.sim$time.ns<-get.ns.basis(bx.sim$time,knots.time.bx)

time.ns.bx.mat<-ns(bx.sim$time, knots=c(2,4,6))
sec.time.ns.bx.mat<-ns(bx.sim$sec.time.std, knots=c(-0.5,0.3,0.9))
sec.time.ns.rrp.mat<-ns(bx.sim$sec.time.std, knots=c(-0.2,0.7))

bx.sim$rm<-rep(0,N)
#bx.sim$rm[bx.sim$sec.time>9.5]<-1
#bx.sim<-bx.sim[bx.sim$rm==0,]
#(N<-dim(bx.sim)[1])


##biopsies
bx.sim$bx.here<-rep(0,N)
bx.sim$num.prev.bx<-rep(1,N)

bx.sub<-bx.sim[bx.sim$time==1,]
(n_bx<-dim(bx.sub)[1])
U.BX<-as.matrix(cbind( rep(1,n_bx), bx.sub$age.std, bx.sub$age.ns, time.ns.bx.mat[bx.sim$time==1,], bx.sub$num.prev.bx, sec.time.ns.bx.mat[bx.sim$time==1,], bx.sub$eta  ))
summary(as.vector(expit(U.BX%*%nu.bx)))

bx.sim$bx.here[bx.sim$time==1]<-rbinom(n,1,as.vector(expit(U.BX%*%nu.bx)))
#table(bx.sim$bx.here[bx.sim$time==1])

for(j in 2:10){
	for(i in 1:n){
		bx.sim$num.prev.bx[bx.sim$id==i & bx.sim$bx.time==j]<-sum(bx.sim$bx.here[bx.sim$id==i & bx.sim$time<j]) + 1}

	bx.sub<-bx.sim[bx.sim$time==j,]
	(n_bx<-dim(bx.sub)[1])
	U.BX<-as.matrix(cbind(rep(1,n_bx), bx.sub$age.std, bx.sub$age.ns, time.ns.bx.mat[bx.sim$time==j,], bx.sub$num.prev.bx, sec.time.ns.bx.mat[bx.sim$time==j,], bx.sub$eta ))
	
	bx.sim$bx.here[bx.sim$time==j]<-rbinom(n,1,as.vector(expit(U.BX%*%nu.bx)))}
#table(bx.sim$bx.here)	

#reclassifications
bx.sim$rc<-bx.sim$prev.G7<-rep(0,N)
rc.sub<-bx.sim[bx.sim$bx.here==1,]
(n_rc<-dim(rc.sub)[1])
V.RC<-as.matrix(cbind(rep(1,n_rc), rc.sub$age.std, rc.sub$time, rc.sub$time.ns, rc.sub$sec.time.std, rc.sub$eta))

bx.sim$rc[bx.sim$bx.here==1]<-rbinom(n_rc,1,as.vector(expit(V.RC%*%gam.rc)))

for(i in 1:n){
	if(sum(bx.sim$rc[bx.sim$id==i]==1)>0){
		rc.time<-min(bx.sim$time[bx.sim$rc==1 & bx.sim$id==i])
		bx.sim$rc[bx.sim$id==i & bx.sim$time>rc.time]<-0
		bx.sim$bx.here[bx.sim$id==i & bx.sim$time>rc.time]<-0
		bx.sim$num.prev.bx[bx.sim$id==i & bx.sim$time>rc.time]<-(bx.sim$num.prev.bx[bx.sim$id==i & bx.sim$time==1] + 1)
		bx.sim$prev.G7[bx.sim$id==i & bx.sim$time>=rc.time]<-1
		bx.sim$rm[bx.sim$id==i & bx.sim$time>(rc.time+2)]<-1}}


bx.sim$rrp<-rep(0,N) #RRP is surgery (radical retropubic prostatectomy)
bx.sim$num.prev.bx.rrp <- bx.sim$num.prev.bx + bx.sim$bx.here

W.SURG<-as.matrix(cbind(rep(1,N), bx.sim$age.std, bx.sim$age.ns, time.ns.bx.mat, sec.time.ns.rrp.mat, bx.sim$num.prev.bx.rrp, bx.sim$prev.G7, bx.sim$eta, (bx.sim$prev.G7*bx.sim$eta) ))

bx.sim$rrp<-rbinom(N,1,as.vector(expit(W.SURG%*%omega.surg)))

#messes up design matrices to delete columns earlier
bx.sim<-bx.sim[bx.sim$rm==0,]
(N<-dim(bx.sim)[1])


pt.data$rc<-pt.data$rrp<-rep(0,n)

for(i in 1:n){
	if(sum(bx.sim$rrp[bx.sim$id==i])>0){
		rrp.time<-min(bx.sim$time[bx.sim$id==i & bx.sim$rrp==1])
		bx.sim$rm[bx.sim$id==i & bx.sim$time>rrp.time]<-1	
		pt.data$rrp[pt.data$id==i]<-1}	}


bx.sim<-bx.sim[bx.sim$rm==0,]
(N<-dim(bx.sim)[1])

for(i in 1:n){
	pt.data$rc[i]<-sum(bx.sim$rc[bx.sim$id==pt.data$id[i]])}
table(pt.data$rc) 

pt.data$obs.eta<-rep(NA,n)
pt.data$obs.eta[pt.data$rrp==1]<-pt.data$eta.true[pt.data$rrp==1]


write.csv(pt.data,"pt-data-sim.csv")
write.csv(bx.sim,"bx-data-sim.csv")



##psa data

psa.time<-seq(-1, max(bx.sim$time[bx.sim$id==1]),0.5)
psa.id<-rep(1, length(psa.time))

for(i in 2:n){
	psa.add<-seq(-1, max(bx.sim$time[bx.sim$id==i]), 0.5)
	psa.time<-c(psa.time,psa.add)
	psa.id<-c(psa.id, rep(i, length(psa.add)))}	
	
psa.data<-as.data.frame(cbind(psa.id, psa.time))
names(psa.data)<-c("id","psa.time")
(n_obs_psa<-dim(psa.data)[1])

psa.data$psa.time<-psa.data$psa.time + runif(n_obs_psa, min=-0.25, max=0.25)
psa.data$age<-vector(length=n_obs_psa)
for(j in 1:n_obs_psa){
	psa.data$age[j] <- psa.data$psa.time[j] + pt.data$age[pt.data$id==psa.data$id[j]]}
	
psa.data$age.std<-(psa.data$age-mean(psa.data$age))/sd(psa.data$age)

pt.data$std.vol<-rnorm(n,0,1)
psa.data$std.vol<-vector(length=n_obs_psa)
for(i in 1:n){
	psa.data$std.vol[psa.data$id==i] <- pt.data$std.vol[i]}

b.vec <- matrix(nrow=n, ncol=2)
for(i in 1:n){
	b.vec[i,] <- mvrnorm(n=1, mu=mu.mat[,(pt.data$eta.true[pt.data$id==i]+1)], Sigma=Sigma)}

write.csv(b.vec,"b-vec-true.csv")



psa.data$log.psa <- vector(length=n_obs_psa)
for(j in 1:n_obs_psa){
	lin.pred <- NULL
	lin.pred <- sum(b.vec[pt.data$id==psa.data$id[j],] * c(1, psa.data$age.std[j])) + beta[1]*psa.data$std.vol[j]
	psa.data$log.psa[j] <- rnorm(1, mean=lin.pred, sd=sigma_res)}
summary(psa.data$log.psa)


write.csv(psa.data,"psa-data-sim.csv")
write.csv(pt.data,"pt-data-sim.csv")




#get ordered subject variable

pt.data<-pt.data[order(pt.data$obs.eta),]
pt.data$subj<-c(1:n)
psa.data$subj<-rep(0,n_obs_psa)
for(i in 1:n){psa.data$subj[psa.data$id==pt.data$id[i]]<-pt.data$subj[i]}
bx.sim$subj<-rep(0,N)
for(i in 1:n){bx.sim$subj[bx.sim$id==pt.data$id[i]]<-pt.data$subj[i]}

write.csv(psa.data,"psa-data-sim.csv")
write.csv(pt.data,"pt-data-sim.csv")
write.csv(bx.sim,"bx-data-sim.csv")









