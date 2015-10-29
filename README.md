# prediction-prostate-surveillance
Joint modeling of latent prostate cancer state and observed clinical outcomes, for active surveillance of low risk prostate cancer


October 29, 2015- This is code accompanying Coley et al. (2015) 
"Bayesian Joint Hierarchical Model for Prediction of Latent Health States with Application to Active Surveillance of Prostate Cancer"


—DATA—
Data from Johns Hopkins Active Surveillance cohort is not publicly available, so simulated data is provided in the folder simulation-data. Data was simulated using posterior estimates from fitting the model on the JH AS data. Explanation of variables in codebook.md in simulation-data folder.


—R Scripts for model estimation—
There are R scripts for a combination of two settings: (1) modeling assumptions and (2) code aim:

(1) Modeling assumptions

(a) IOP- Model is fit allowing for an Informative Observation Process (IOP), as described in the paper and

(b) UNADJ- Model assumes biopsy and surgery information is missing as random, so it is unadjusted for the possibility of IOP.


(2) Code Aim

(a) General model fitting, individual-level plots, and calibration plots for observed outcomes- This code demonstrates how the model is fit in JAGS. Some observations are deleted in order to demonstrate individual-level predictions (Figure 3). Model output is also used to make calibration plots for observed outcomes— biopsies performed, reclassification, and surgery (Figure 5). 

(b) CV- The purpose of this code is to run cross-validation on the observed latent state (eta). At each replication of the code, observed true cancer state is masked for one patient. Then, the predicted state (posterior mean of eta) is compared to the true state in PLOTS- predicted-vs-observed-eta.R


There are also R scripts for PLOTS and RESULTS. See explanations below.

—Batch Scripts—
R scripts were written to be run on a SGE computing cluster, in order to run multiple chains in parallel or, in the case of cross-validation, multiple folds. Batch scripts are given in a batch-scripts folder to get users an idea of the memory needed, task id setting, etc, but they will need to be considerably modified for an individual users account and platform. 



-PLOTS-
Figures from paper recreated with simulated data. R scripts entitled PLOTS-… create these figures and save them in the plots folder.


-RESULTS-
RESULTS-prediction-model-jags-inf-obs.R summarizes posterior results and creates trace plots and posterior densities for each parameter (saved in plots/results-check).

RESULTS-reduce-size.R takes posterior samples of patient-level random effects and posterior fitted probabilities for observed outcomes and rounds to three significant digits. (The R/JAGS script samples and saves values up to 7 sig digs, but these files can be cumbersome to upload on github. This code makes them more manageable. Only rounded versions of these results are in this repo, but the code supplied will produce the larger files.)