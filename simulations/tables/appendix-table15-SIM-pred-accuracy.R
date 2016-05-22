### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
#Get predictive accuracy for Table 15 in the appendix


load('simulations/UNADJ/results/sim-UNADJ-res-eta-known.RData')
mean(auc.res)
quantile(auc.res, p=c(0.025, 0.975))
mean(mse.res)
quantile(mse.res, p=c(0.025, 0.975))


load('simulations/UNADJ/results/sim-UNADJ-res-eta.RData')
mean(auc.res)
quantile(auc.res, p=c(0.025, 0.975))
mean(mse.res)
quantile(mse.res, p=c(0.025, 0.975))



load('simulations/IOP-BX/results/sim-IOP-BX-res-eta-known.RData')
mean(auc.res)
quantile(auc.res, p=c(0.025, 0.975))
mean(mse.res)
quantile(mse.res, p=c(0.025, 0.975))


load('simulations/IOP-BX/results/sim-IOP-BX-res-eta.RData')
mean(auc.res)
quantile(auc.res, p=c(0.025, 0.975))
mean(mse.res)
quantile(mse.res, p=c(0.025, 0.975))


load('simulations/IOP-SURG/results/sim-IOP-SURG-res-eta-known.RData')
mean(auc.res)
quantile(auc.res, p=c(0.025, 0.975))
mean(mse.res)
quantile(mse.res, p=c(0.025, 0.975))


load('simulations/IOP-SURG/results/sim-IOP-SURG-res-eta.RData')
mean(auc.res)
quantile(auc.res, p=c(0.025, 0.975))
mean(mse.res)
quantile(mse.res, p=c(0.025, 0.975))



load('simulations/IOP/results/sim-IOP-res-eta-known.RData')
mean(auc.res)
quantile(auc.res, p=c(0.025, 0.975))
mean(mse.res)
quantile(mse.res, p=c(0.025, 0.975))


load('simulations/IOP/results/sim-IOP-res-eta.RData')
mean(auc.res)
quantile(auc.res, p=c(0.025, 0.975))
mean(mse.res)
quantile(mse.res, p=c(0.025, 0.975))





load('simulations/LOGISTIC/results/sim-LOGISTIC-res-eta.RData')

mean(auc.ek.out)
quantile(auc.ek.out, p=c(0.025, 0.975))

mean(mse.ek.out)
quantile(mse.ek.out, p=c(0.025, 0.975))

mean(auc.res)
quantile(auc.res, p=c(0.025, 0.975))

mean(mse.res)
quantile(mse.res, p=c(0.025, 0.975))
