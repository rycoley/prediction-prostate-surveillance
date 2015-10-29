# Codebook

### Abbreviations and Acronyms
* PSA: prostate-specific antigen
* PCa: prostate cancer
* BX: biopsy
* RC: reclassification
* SURG: surgical removal of prostate (radical retropubic prostatectomy)
* PT: patient
* DX: diagnosis

### PT data
1 record per patient
* id: pt-specific unique identifier, assigned for data generation
* age_dx: age at PCa diagnosis
* date_dx: calendar date of diagnosis, in numeric form using as.Date() default in R
* eta_true: latent state, either 0 or 1, corresponding to indolent (Gleason <=6) or aggressive (Gleason >=7) cancer (respectively); generated for all patients for data simulation but only observed on a subset with surgery (rap)
* obs_eta: eta_true if surg=1, NA otherwise
* surg: does this subject ever have surgery?
* rc: does this subject ever reclassify, i.e. have a biopsy with a Gleason score 7 or higher
* vol_std: prostate volume, mean and sd standardized
* subj: pt-specific unique identifier, assigned after data simulation to correspond to data sorting. data sorted based on obs_eta as follows: eta observed and =0, eta observed and =1, eta unobserved. (sorting by obs_eta makes posterior sampling in JAGS much easier)

### PSA data
1 record per person per PSA observation
* psa: outcome measured in continuous time. Measurements may occur before initial PCa diagnosis. 
* log_psa: log-transformed PSA used for regression
* age: age at PSA measurement
* age_std: mean and std dev-standardized age at psa (standardized within this dataset); used in X
* vol_std: prostate volume, standardized within psa dataset; used in Z
* subj: pt-specific unique identifier, see definition in pt_data

### BX Data
1 record per annual interval where pt eligible for biopsy or (after rc) surgery.
Patients eligible for BX until RC (as per active surveillance protocol); it is only possible to reclassify once.
Patients eligible for surgery prior to RC or up to 2 years after RC
Patients censored after surgery, 2 years post-RC, or 10 years
No record for initial diagnostic bx 


* eta: true state (for all patients, used to generate data)
* bx_here: biopsy occurs in this interval; NA after a pt reclassifies
* rc: grade reclassification occurs at this annual biopsy
* surg: surgery performed in this annual interval
* subj: pt-specific unique identifier, see definition in pt_data

* time: years since dx for interval
* age: age during this interval (mid-year)
* date: calendar date during this interval (mid-year)

* bx_time_ns, rc_time_ns, surg_time_ns: natural spline representations of time used in covariance matrices U,V,W (respectively) and calculated with ns(time, knots=bx_time_knots, Boundary.knots=bx_time_bknots), etc. 
* bx_date_ns, rc_date_ns, surg_date_ns: natural spline representations of calendar date used in covariance matrices U,V,W (respectively) and calculated with ns(date, knots=bx_date_knots, Boundary.knots=bx_date_bknots), etc. 
* bx_age_ns, surg_age_ns: natural spline representations of age used in covariance matrices U and W (respectively) and calculated with ns(age, knots=bx_age_knots, Boundary.knots=bx_age_bknots), etc. 
* rc_age_std: standardized age for biopsy results, used in covariance matrix V, and calculated with scale(age, rc_age_mean, rc_age_sd)

* num_prev_bx: number of prior bx at beginning of interval (used to predict P(BX) in this interval); everyone has one prior bx (diagnostic) at time=1 
* bs_num_prev_bx_ns: natural spline representation of num_prev_bx for covariate matrix U and calculated with ns(num_prev_bx, knots=bx_num_prev_bx_knots, Boundary.knots=bx_num_prev_bx_bknots)
* surg_num_prev_bx_std: standardized number of previous biopsies at the end of the interval for covariate matrix W and calculated with scale((num_prev_bx + bx_here), surg_num_prev_bx_mean, surg_num_prev_bx_sd)

* prev_G7: bx grade of Gleason 7 or higher (i.e. reclassification) in this interval or previous intervals (included in covariate matrix W)

* rm: indicator for removing a record from dataset (used in data generation); not used in analysis


