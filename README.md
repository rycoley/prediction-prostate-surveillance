# prediction-prostate-surveillance
Joint modeling of latent prostate cancer state and observed clinical outcomes, for active surveillance of low risk prostate cancer

May 20, 2016- This is code accompanying ANONYMOUS et al. (2016) 
“A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"


This README contains descriptions of:
Data
R scripts for model estimation (data analysis and simulations)
Batch scripts
R scripts for recreating plots
Results and R scripts to summarize results



——JHAS Data Analysis——

—DATA—
Data from Johns Hopkins Active Surveillance cohort is not publicly available, so simulated data is provided in the folder simulation-data. Data was simulated using posterior estimates from fitting the model on the JH AS data (though all spline knots may not be equal to those given in Appendix Tables 5-7). Explanation of variables can be found in codebook.md in the simulation-data folder.


—R Scripts for model estimation—
There are R scripts for a combination of two settings: (1) modeling assumptions and (2) code aim:

(1) Modeling assumptions

(a) IOP- Model is fit allowing for a Informative Observation Process (IOP) components for biopsy and surgery, as described in the paper; IOP-BX and IOP-SURG only allow for one of the IOP components (biopsy and surgery, respectively)

(b) UNADJ- Model assumes biopsy and surgery information is missing as random, so it is unadjusted for the possibility of IOP.


(2) Code Aim

(a) General model fitting, individual-level plots, and calibration plots for observed outcomes- This code demonstrates how the model is fit in JAGS. Observations after 5 years are deleted for a dozen patients in order to demonstrate individual-level predictions (Figure 22 in Appendix). Model output is also used to make calibration plots for observed outcomes— biopsies performed, reclassification, and surgery (e.g., Figure 6). 

(b) CV- The purpose of this code is to run cross-validation on the observed latent state (eta). At each replication of the code, observed true cancer state is masked for one patient. Then, the predicted state (posterior mean of eta) is compared to the true state to calculate predictive accuracy and visualize calibration.

(c) Robustness of IOP Coefficients- Model estimation is performed with various informative priors on IOP components. Model output is used to make posterior density plots in Appendix Figure 18.


There are also R scripts for PLOTS and RESULTS. See explanations below.

—Batch Scripts—
R scripts were written to be run on a SGE computing cluster, in order to run multiple chains in parallel or, in the case of cross-validation, multiple folds. Batch scripts are given in a batch-scripts folder to get users an idea of the memory needed, task id setting, etc, but they will need to be considerably modified for an individual users account and platform. 


-PLOTS-
Figures from paper recreated with simulated data. R scripts are labeled by the figure they reproduce.


-TABLES-
Tables in the appendix include analysis results for posterior summaries and predictive accuracy. A text file in the folder jhas-analysis/tables directs the users to R scripts to see how those results were obtained.



Notes:

*Results for this simulated dataset are not the same as those given in the primary paper for the JHAS cohort analysis. These results may not reflect paper conclusions. 

*jhas-analysis/IOP/RESULTS-reduce-size.R takes posterior samples of patient-level random effects and posterior fitted probabilities for observed outcomes and rounds to three significant digits. (The R/JAGS script samples and saves values up to 7 sig digs, but these files can be cumbersome to upload on github. This code makes them more manageable. Only rounded versions of these results are in this repo, but the code supplied will produce the larger files.)

*We also use a logistic regression model for comparison. This can be found in simulations/LOGISTIC.





——Simulation studies——


The content of the simulation-studies folder is similar to that of the jhas-analysis folder. Datasets are generated using the parameter values in simulations/simulation-data for four versions of the proposed model (IOP, IOP-SURG, IOP-BX, UNADJ). For each, of 200 simulated datasets, (1) posterior samples are saved, (2) predictions of the true cancer state are saved for patients without surgery, and (3) and importance sampling algorithm (IS) is used to obtain predictions for patients with surgery. R scripts and batch scripts are provided for all of these tasks. Logistic regression (LOGISTIC) is also used for comparison. 

simulations/plots contains R scripts to reproduce figures in the primary paper and appendix. simulations/tables contains R scripts to obtain information in appendix tables.




