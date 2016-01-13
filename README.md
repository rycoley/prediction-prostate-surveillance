# prediction-prostate-surveillance
Joint modeling of latent prostate cancer state and observed clinical outcomes, for active surveillance of low risk prostate cancer

January 12, 2016- This is code accompanying Coley et al. (2016) 
“A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"


This README contains descriptions of:
Data
R scripts for model estimation
Batch scripts
R scripts for recreating plots
Results and R scripts to summarize results
**Workflow**


—DATA—
Data from Johns Hopkins Active Surveillance cohort is not publicly available, so simulated data is provided in the folder simulation-data. Data was simulated using posterior estimates from fitting the model on the JH AS data. Explanation of variables can be found in codebook.md in the simulation-data folder.


—R Scripts for model estimation—
There are R scripts for a combination of two settings: (1) modeling assumptions and (2) code aim:

(1) Modeling assumptions

(a) IOP- Model is fit allowing for an Informative Observation Process (IOP), as described in the paper and

(b) UNADJ- Model assumes biopsy and surgery information is missing as random, so it is unadjusted for the possibility of IOP.


(2) Code Aim

(a) General model fitting, individual-level plots, and calibration plots for observed outcomes- This code demonstrates how the model is fit in JAGS. Observations after 5 years are deleted for a dozen patients in order to demonstrate individual-level predictions (Figure 6). Model output is also used to make calibration plots for observed outcomes— biopsies performed, reclassification, and surgery (Figure 5). 

(b) CV- The purpose of this code is to run cross-validation on the observed latent state (eta). At each replication of the code, observed true cancer state is masked for one patient. Then, the predicted state (posterior mean of eta) is compared to the true state in PLOTS- predicted-vs-observed-eta.R

(c) Sensitivity Analysis- Model estimation is performed with various informative priors on IOP components. Model output is used to make posterior density plots in Figure 7.


There are also R scripts for PLOTS and RESULTS. See explanations below.

—Batch Scripts—
R scripts were written to be run on a SGE computing cluster, in order to run multiple chains in parallel or, in the case of cross-validation, multiple folds. Batch scripts are given in a batch-scripts folder to get users an idea of the memory needed, task id setting, etc, but they will need to be considerably modified for an individual users account and platform. 



-PLOTS-
Figures from paper recreated with simulated data. R scripts entitled PLOTS-… create these figures and save them in the plots folder.


-RESULTS-
RESULTS-prediction-model-jags-inf-obs.R summarizes posterior results and creates trace plots and posterior densities for each parameter (saved in plots/results-check).

RESULTS-reduce-size.R takes posterior samples of patient-level random effects and posterior fitted probabilities for observed outcomes and rounds to three significant digits. (The R/JAGS script samples and saves values up to 7 sig digs, but these files can be cumbersome to upload on github. This code makes them more manageable. Only rounded versions of these results are in this repo, but the code supplied will produce the larger files.)



-WORKFLOW- 
To reproduce the analysis found in Coley et al (2016) using simulated data, do the following:

I. Reproduce unadjusted (non-IOP) analysis

1. Run  UNADJ-call-jags.R with starting seed 1-5 (as in batch script UNADJ-jags-batch-script.s). This R script calls UNADJ-prep-data-for-jags.R in order to tidy and shape simulated data for estimation and defines settings for posterior sampling, including the JAGS model defined in UNADJ-jags-model.txt, before performing model estimation and saving posterior samples in the “results” folder.

2. Run UNADJ-CV-call-jags.R with starting seeds 1-203 (all patients with true state observed) in order to obtain out-of-sample predictions. (Saved to “results folder.)


II. Reproduce analysis that adjusts for a possible informative observation process (IOP). 

1. Run  IOP-call-jags.R with starting seed 1-5 (as in batch script IOP-jags-batch-script.s). This R script calls IOP-prep-data-for-jags.R in order to tidy and shape simulated data for estimation and defines settings for posterior sampling, including the JAGS model defined in IOP-jags-model.txt, before performing model estimation and saving posterior samples in the “results” folder.

2. Run IOP-CV-call-jags.R with starting seeds 1-203 (all patients with true state observed) in order to obtain out-of-sample predictions. (Saved to “results folder.)




III. Reproduce summary of posterior predictions of true cancer state.

1. PLOTS-fig4-descriptive-predictions-eta.R reproduces stacked histogram and scatter plot in Figure 4 of paper. (Note: Paper gives results for true clinical data. This script will create an analogous plot for simulated data.)

IV. Reproduce performance and goodness-of-fit metrics.

1. PLOTS-fig4-predicted-vs-observed-eta.R reproduces the sensitivity and specificity calculations, AUC curves, and calibration plot for predictions of eta shown in Figure 4.  Estimates and intervals for sensitivity, specificity, and AUC are also reported in the manuscript text. (Note: Paper gives results for true clinical data. This script will create an analogous plot for simulated data.)

2. PLOTS-fig-5-predicted-vs-observed-outcomes.R reproduces calibration plots for observables given in Figure 5. (Note: Paper gives results for true clinical data. This script will create an analogous plot for simulated data.)



V. Reproduce predictions for individual patients.

1. PLOTS-fig6-individual-data-and-predictions plots the data and predictions for a dozen patients who had post-year 5 observations deleted from the analysis datasets. The result is the same as Figure 6 in the manuscript.




VI. Reproduce sensitivity analysis

1. Run IOP-SENSITIVITY-ANALYSIS-call-jags.R to fit the IOP model with prior means for IOP effects (as described in paper). All combinations of COEF_MAIN = -0.7, 0, 0.7 and COEF_INT = -0.7, 0, 0.7 are called via batch scripts. On a SGE cluster, for e.g., -v COEF_MAIN=-0.7,COEF_INT=-0.7. (Saved to “results folder.)

2. Posterior densities shown in Figure 7 can be reproduced by PLOTS-fig7-posterior-iop-coefficients.R (Note: Paper gives results for true clinical data. This script will create an analogous plot for simulated data.)

 
