### Rebecca Yates Coley rycoley@gmail.com
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"

##Get results for Tables 10-14 in appendix

load("simulations/UNADJ/results/sim-UNADJ-results.RData")



round(mean(p_eta.est),4)
mean(p_eta.cov)

cbind(round(apply(mu_int.est,2,mean),4), apply(mu_int.cov,2,mean))
cbind(round(apply(mu_slope.est,2,mean),4), apply(mu_slope.cov,2,mean))

round(mean(sigma_int.est),4)
mean(sigma_int.cov)

round(mean(sigma_slope.est),4)
mean(sigma_slope.cov)

round(mean(cov_int_slope.est),4)
mean(cov_int_slope.cov)

round(mean(beta.est),4)
mean(beta.cov)

round(mean(sigma_res.est),4)
mean(sigma_res.cov)


cbind(round(apply(gamma_RC.est,2,mean),4), apply(gamma_RC.cov,2,mean))




load("simulations/IOP-BX/results/sim-IOP-BX-results.RData")


round(mean(p_eta.est),4)
mean(p_eta.cov)

cbind(round(apply(mu_int.est,2,mean),4), apply(mu_int.cov,2,mean))
cbind(round(apply(mu_slope.est,2,mean),4), apply(mu_slope.cov,2,mean))

round(mean(sigma_int.est),4)
mean(sigma_int.cov)

round(mean(sigma_slope.est),4)
mean(sigma_slope.cov)

round(mean(cov_int_slope.est),4)
mean(cov_int_slope.cov)

round(mean(beta.est),4)
mean(beta.cov)

round(mean(sigma_res.est),4)
mean(sigma_res.cov)


cbind(round(apply(nu_BX.est,2,mean),4), apply(gamma_RC.cov,2,mean))

cbind(round(apply(gamma_RC.est,2,mean),4), apply(gamma_RC.cov,2,mean))




load("simulations/IOP-SURG/results/sim-IOP-SURG-results.RData")


round(mean(p_eta.est),4)
mean(p_eta.cov)

cbind(round(apply(mu_int.est,2,mean),4), apply(mu_int.cov,2,mean))
cbind(round(apply(mu_slope.est,2,mean),4), apply(mu_slope.cov,2,mean))

round(mean(sigma_int.est),4)
mean(sigma_int.cov)

round(mean(sigma_slope.est),4)
mean(sigma_slope.cov)

round(mean(cov_int_slope.est),4)
mean(cov_int_slope.cov)

round(mean(beta.est),4)
mean(beta.cov)

round(mean(sigma_res.est),4)
mean(sigma_res.cov)


cbind(round(apply(gamma_RC.est,2,mean),4), apply(gamma_RC.cov,2,mean))

cbind(round(apply(omega_SURG.est,2,mean),4), apply(gamma_RC.cov,2,mean))



load("simulations/IOP/results/sim-iop-results.RData")
#ls()

round(mean(p_eta.est),4)
mean(p_eta.cov)

cbind(round(apply(mu_int.est,2,mean),4), apply(mu_int.cov,2,mean))
cbind(round(apply(mu_slope.est,2,mean),4), apply(mu_slope.cov,2,mean))

round(mean(sigma_int.est),4)
mean(sigma_int.cov)

round(mean(sigma_slope.est),4)
mean(sigma_slope.cov)

round(mean(cov_int_slope.est),4)
mean(cov_int_slope.cov)

round(mean(beta.est),4)
mean(beta.cov)

round(mean(sigma_res.est),4)
mean(sigma_res.cov)


cbind(round(apply(nu_BX.est,2,mean),4), apply(gamma_RC.cov,2,mean))

cbind(round(apply(gamma_RC.est,2,mean),4), apply(gamma_RC.cov,2,mean))

cbind(round(apply(omega_SURG.est,2,mean),4), apply(gamma_RC.cov,2,mean))

