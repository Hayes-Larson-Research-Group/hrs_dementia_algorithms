# hrs_dementia_algorithms
This repository holds code for the the project titled: "Classification Matters: Divergent Estimates of Dementia Risk Factors in the Health and Retirement Study" 

* **0_paths_jl.R** : Sets the working path references across the project. This script is to help code reviewers to easily run subsequent scripts without manually editing file locations. 

* **function_jl.R** : Reusable utility functions.

* **1_data_construction.R** : Constructing analysis-ready datasets for each dementia classification algorithm. Baseline covariates were obtained from RAND longitudinal file and Tracker file and merged into dementia algorithm datasets. 

* **2_table1.Rmd** : For each algorithm, generate both the unweighted and weighted Table 1. Then, restrict the sample to individuals not classified as having dementia in 2010 by each algorithm, merge in the end-of-follow-up event, and generate the unweighted and weighted Table 1 for these samples as well.

* **3_km_plot.Rmd** : Plot Kaplan-Meier curves for both unweighted and weighted 70+ samples, stratified by algorithms.

* **4_km_plot_race_and_ethnicity.Rmd** : Plot Kaplan-Meier curves for both unweighted and weighted 70+ samples, stratified by race_ethnicity and algorithms 

* **5_incident_rate_Rmd** : Use functions to calculate the weighted and unweighted person-years (PY) and number of cases for each age group, then compute the weighted and unweighted age-specific incidence rates (IR), overall IR, and age-standardized IR along with their confidence intervals (CIs), stratified by algorithms.

* **6_race_and_ethnicity_ir.Rmd** : : Use functions to calculate the weighted and unweighted person-years (PY) and number of cases for each age group, then compute the weighted and unweighted age-specific incidence rates (IR), overall IR, and age-standardized IR along with their confidence intervals (CIs), stratified by algorithms, race_ethnicity.

* **7_main_cox_model_overall_stratified.Rmd** : Fit both age-adjusted and unadjusted models for each algorithm in the 70+ sample using both weighted and unweighted data, then fit models stratified by race/ethnicity.

* **8_main_cox_model_overall_stratified_apoe_sample.Rmd** : In the 70+ sample, this script first restricts to participants with APOE information, then fits age-adjusted and unadjusted Cox models for each algorithm (Langa, Wu, Expert, Hudomiet) using both weighted and unweighted data, and finally fits Cox models stratified by race/ethnicity.

* **9_cox_model_using_apoe_variable_overall_stratified.Rmd** : In the 70+ APOE subsample, fits Cox models with APOE for each algorithm (Langa, Wu, Expert, Hudomiet), using both weighted and unweighted data. Each algorithm gets two models—an APOE + age model and a fully adjusted model—and the same analysis is then repeated stratified by race/ethnicity.

* **10_forest_plot_hr_jl.Rmd** : Builds HR forest plots from the weighted Cox results. Combines the main risk-factor HRs with the APOE-ε4 HR into a single figure, produces a sensitivity-analysis version restricted to the APOE subsample, and produces a race/ethnicity-stratified forest plot faceted by race.

* **11_forest_plot_ir_jl.Rmd** : Builds the age-specific incidence-rate figures from the IR output files. Produces an overall IR-by-age plot comparing the four algorithms, and a race/ethnicity-stratified version faceted by race.

* **12_main_all_adjust_apoe_overall.Rmd** : analysis on the 70+ APOE subsample. Re-fits the main risk-factor Cox models for each algorithm (Langa, Wu, Expert, Hudomiet), weighted and unweighted, but with APOE added as a covariate to every model. Then reshapes the weighted results into long format for forest-plotting.

* **13_reverter_plot.Rmd** : Uses TraMineR sequence-index plots to visualize cognitive-state trajectories over 2010–2020 for each algorithm (Langa-Weir, Hudomiet, Expert, Wu). Plots both 2-category (dementia / no dementia) and 3-category (Normal / CIND / Dementia) trajectories, counts "reverters" (people who move from dementia back to a less-impaired state), and produces final saved plots restricted to the dementia subset who reverted.

* **14_data_construction_censor_2016.R** : data construction for computing dementia incidence rates for each algorithm (Langa-Weir, Hudomiet, Expert, Wu) on the 2016-censored sens2016 sample.

* **15_censor_2016_main_cox_model_overall_stratified.Rmd** : Sensitivity analysis that re-runs the full main Cox pipeline on the 2016-censored dataset. Fits age-adjusted and unadjusted models for each algorithm (Langa, Wu, Expert, Hudomiet), weighted and unweighted, then fits models stratified by race/ethnicity—structurally identical to script #7 but on the sens2016 sample.

* **16_sensitivity_2016_apoe_cox_model.Rmd** : 2016-censoring sensitivity analysis that runs three stages in one file on the sens2016 sample: (1) the main risk-factor Cox models for each algorithm (weighted/unweighted), (2) the APOE-exposure Cox models on the APOE subsample, and (3) reshaping the weighted results into long format and binding the APOE-ε4 HR.

* **17_ir_censor_2016.Rmd** : Computes dementia incidence rates for each algorithm (Langa-Weir, Hudomiet, Expert, Wu) on the 2016-censored sens2016 sample. Produces age-specific, crude, and age-adjusted IRs (standardized to the 2000 US census population), both unweighted (closed-form CIs) and survey-weighted (bootstrap CIs), and writes them out as wide tables by algorithm.

* **18_flowchart.R** : Generates the data for the cohort flow chart, and produce the output needed to build the flow chart.
