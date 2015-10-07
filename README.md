# prediction-prostate-surveillance
Joint modeling of latent prostate cancer state and observed clinical outcomes, for active surveillance of low risk prostate cancer


October 6, 2015- This is code accompanying Coley et al. (2015) 
"Bayesian Joint Hierarchical Model for Prediction of Latent Health States with Application to Active Surveillance of Prostate Cancer"


—DATA—
Data from Johns Hopkins Active Surveillance cohort is not publicly available, so simulated data is provided in the folder simulation-data. Data was simulated using posterior estimates from fitting the model on the JH AS data. Explanation of variables in readme.md in this folder.

Note: Since the data generating model is the same as the estimation model, plots of model fit will be overly optimistic. In particular, the CV-AUC for latent state is considerably higher.

—R Scripts—
There are R scripts for a combination of two settings: (1) modeling assumptions and (2) code aim:

(1) Modeling assumptions

(a) IOP- Model is fit allowing for an Informative Observation Process (IOP), as described in the paper and
(b) UNADJ- Model assumes biopsy and surgery information is missing as random, so it is unadjusted for the possibility of IOP.

(2) Code Aim

(a) General model fitting, individual-level plots, and calibration plots for observed outcomes- This code demonstrates how the model is fit in JAGS. Some observations are deleted in order to demonstrate individual-level predictions (Figure 3). Model output is also used to make calibration plots for observed outcomes— biopsies performed, reclassification, and surgery (Figure 5). 

(b) CV- The purpose of this code is to run cross-validation on the observed latent state (eta). At each replication of the code, observed true cancer state is masked for one patient. Then, the predicted state (posterior mean of eta) is compared to the true state in PLOTS- predicted-vs-observed-eta.R

—Batch Scripts—
R scripts were written to be run on a SGE computing cluster, in order to run multiple chains in parallel or, in the case of cross-validation, multiple folds. Batch scripts are given in a batch-scripts folder to get users an idea of the memory needed, task id setting, etc, but they will need to be considerably modified for an individual users account and platform. 

—Results—
.csv files with posterior samples are not posted to github due to size. Email rycoley@gmail.com to request.

-CODEBOOK-
Explains all variables in simulated datasets
