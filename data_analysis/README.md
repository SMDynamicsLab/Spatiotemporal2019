Scripts to perform data wrangling, analysis, and visualization.
Auxiliar functions and ready-to-use data are in another folder for clarity

# data_wrangling.m
Takes the raw data of two experiments and generates the mean time series for each condition across subjects, experiments, and size of perturbations (Small and Large) between experiments.

- input: subject_simetric_raw and subject_simetric_xl_raw
- output: subject_simetric-02-Aug-2022.mat, subject_simetric_xl-02-Aug-2022.mat, average_simetric-02-Aug-2022.mat, average_simetric_xl-02-Aug-2022.mat, and average_between_experiments_02-Aug-2022.mat

# statistical_significance.m 
Compute the p-values from generated null distributions of ASYM time series and DIFF time series.

- input: average_between_experiments_02-Aug-2022.mat, subject_simetric-02-Aug-2022.mat, and subject_simetric_xl-02-Aug-2022.mat
- output: clase1_sim-1.mat, clase1_sim-2.mat, clase2_sim-1.mat, and clase2_sim-2.mat

# figure_manager.m
Generates all the figures published in the paper.

- input: average_by_experiment_02-Aug-2022.mat, subject_simetric-02-Aug-2022.mat, subject_simetric_xl-02-Aug-2022.mat, clase1_sim-1.mat, clase1_sim-2.mat, clase2_sim-1.mat, and clase2_sim-2.mat
- output: 8 figures

subject_simetric_raw.mat and subject_simetric_xl_raw.mat contains the data collected in a structure named Subject.

# data folder
Contains the output of data_wrangling.m and statistical_significance.m to run figure_manager.m without the processing of the raw data.

# aux_functions folder
- **f1_variables_trials.m**: calculates variables by trial
- **f2_detects_outliers_trials.m**: detects outliers trials.
- **f3_average_within_subjects.m**: calculates average mean by subject, without outliers trials.
- **f4_detects_subjects.m**: filters subject by number of valid trials.
- **f5_average_between_subjects.m** and **f5_average_between_subjects_xl.m**: average time series between subject of the same class.
- **fdr_bh.m**: performs False Discovery Rate (FDR) correction.
- **columnlegend_RL.m**: fix a problem with the two-column legend and the longitude of the text.
- **ds2nfu.m**: convert data space units into normalized figure units.
