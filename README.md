# ZIBGLMM: Zero-Inflated Bivariate Generalized Linear Mixed Model for Meta-Analysis with Double-Zero-Event Studies

## Data
### Real-world data from CDSR
- The `CDSR_data.csv` contain the real-world meta-analyses used in the study.
  - `CDSR.id` Unique identifier for the Cochrane Database of Systematic Reviews (CDSR) entry.
  - `ma.id` Unique identifier for the meta-analysis within the CDSR entry.
  - `study.name` Name or identifier of the individual study included in the meta-analysis.
  - `ID` Unique meta-analysis identifier. (Unique CDSR.id and ma.id combinations.)
  - `counter` A sequential number assigned to each study within a meta-analysis. 
  - `n1` Sample size of the treatment group in the study.
  - `r1` Number of events or successes in the treatment group.
  - `n2` Sample size of the control group in the study.
  - `r2` Number of events or successes in the control group.
  - `total_count` Total number of participants in the study (n1 + n2).
  - `dzs_count` Number of double-zero studies in the meta-analysis (studies where r1 = r2 = 0).
  - `dz_ratio` Ratio of double-zero studies to total studies in the meta-analysis.
### Simulated meta-analyses
- The `simulation_data.csv` contain the simulated meta-analyses used in the study.
  - `metaanalysisid` Unique identifier for each meta-analysis in the dataset.
  - `site` Identifier for the study site or individual study within the meta-analysis.
  - `n1`, `n2`, `y1`, `y2` Sample size and number of events of the treatment/control group in the study.
  - `ZI_study` Zero-inflated study indicator. A binary variable indicating whether the study is generated as a zero-inflated study.
  - `zi_rate` Zero-inflation rate used to generate the meta-analysis. 
  - `n_study` Total number of studies within the meta-analysis. 
  - `b0`, `b1` fixed effects and represent the average risks in the control and treatment groups in the logit scale
  - `yy` Used for SAS data loading. Please safely ignore it. 

## Run the models

### Frequentist BGLMM and ZIBGLMM
- The code for frequentist BGLMM and ZIBGLMM is contained in `ZIBGLMM/code/SAS/` folder
- You can utilize SAS OnDemand for Academics [https://welcome.oda.sas.com/] to run SAS codes.
- Make sure to upload the `simulation_data.csv` and `CDSR_data.csv` to the server along with the codes, make sure to update the path to the files, `datafile="/home/u63419104/simulations_data.csv" `.
- Please make sure to update the path for the output in your SAS code, e.g. `OUTFILE='/home/u63419104/zibglmm_2.csv'`.

### Bayesian BGLMM and ZIBGLMM
- The code for Bayesian BGLMM and ZIBGLMM is contained in `ZIBGLMM/code/bayesian/` folder
- You can set the following arguments 
  - `start_id <- 1`
  - `end_id <- 1`
  - `num_iter = 100000`
  - `thinning = 100`
  - `burnin = 100`
  - `results_folder = "results_bayesian_bglmm"`

### Bayesian ZIBGLMM stan version 
- The code for Bayesian BGLMM and ZIBGLMM implemented in Rstan is contained in `ZIBGLMM/code/rstan_version/` folder
- You can change line 43 to change the model to `zibayesian_bglmm.stan` or `bayesian_bglmm.stan` to run ZIBGLMM or BGLMM by compiling different model files. 
- Note that our original model utilizes an Inverse-Wishart distribution for the priors and adapts the scale matrix dynamically during sampling, which can be complex to replicate in Stan. Consequently, the RStan implementation employs an LKJ correlation prior to the covariance matrix instead of the Inverse-Wishart distribution. The users can explore different prior configurations more conveniently within the RStan framework.

