% Data wrangler
% This scripts wrangler the raw data from the experiments published in López, S.L., Laje, R. 
% Spatiotemporal perturbations in paced finger tapping suggest a common 
% mechanism for the processing of time errors. Sci Rep 9, 17814 (2019). 
% https://doi.org/10.1038/s41598-019-54133-x

%% Clean start
clear all
close all

%% Load raw data from both experiments

for i=1:2 
    
    if i==1
        load('subject_simetric_raw.mat')
    else 
        load('subject_simetric_xl_raw.mat')
    end

%% Processing secuence
% 1. Calculates variables by trial
% 2. Detects outliers trials 
% 3. Calculates average mean by subject, without outliers trials
% 4. Filters subject by number of valid trials
% 5. Average time series between subject of the same class

%% Parameters----------------------------------------------
pre_baseline_bips = 7;	    % N° de bips used for the prebaseline
pos_baseline_bips = 5;      % N° de bips used for the posbaseline
out_limit = 145;		    % Maximum asyn permited

%% 1. Calculates variables by trial
[subject]=f1_variables_trials(subject,pre_baseline_bips,pos_baseline_bips,out_limit);

%% 2. Detects outliers trials 
[subject]=f2_detects_outliers_trials(subject);

%% 3. Calculates average mean by subject, without outliers trials
[subject]=f3_average_within_subjects(subject);

%% 4. Filters subject by number of valid trials
n_trial=4;  % Minimum number of valid trials required
[subject,sub_out]=f4_detects_subjects(subject,n_trial);

sub_out_disp=[sub_out];
disp(['ID from subject with few trials = ' num2str(sub_out_disp)])

% Remove one of the authors from the experimental data
if i==1
    sub_au=6;
    subject(6).out=1;
    sub_out=sort([sub_out_disp sub_au]);
else
    sub_au=[1];
    sujeto(1).out=1;
    sub_out=sort([sub_out_disp sub_au]);
end

%% Parameters of analysis
parameters.pre_bs_bips = pre_baseline_bips;
parameters.pos_bs_bips = pos_baseline_bips;
parameters.out_limit = out_limit;

%% Saves 

if i==1
save(['subject_simetric-' date '.mat'],'subject','parameters')
else
save(['subject_simetric_xl-' date '.mat'],'subject','parameters')
end


%% 5. Average time series between subject of the same perturbation

if i==1
    clase_cuatro=[13:15];
    [clase_pert] = f5_average_between_subjects(subject,sub_out,clase_cuatro,parameters);
else
    clase_dos=[13:15];
    [clase_pert] = f5_average_between_subjects_xl(subject,sub_out,clase_dos,parameters);
end

end

%% 6 Average between experiments
exp1=load('average_simetric-02-Aug-2022.mat');
exp2=load('average_simetric_xl-02-Aug-2022.mat');

%% ----------------------- POSTBASELINE -----------------------------------
%% Edito las matrices originales
small=exp1.clase_pert(4).Cond_Mat_Pos;
small(any(any(isnan(small),3),2),:,:) = [];
thirty=exp1.clase_pert(2).Cond_Mat_Pos;
thirty(any(any(isnan(thirty),3),2),:,:) = [];
forty=exp1.clase_pert(3).Cond_Mat_Pos;
forty(any(any(isnan(forty),3),2),:,:) = [];
fifty=exp2.clase_pert(2).Cond_Mat_Pos;
fifty(any(any(isnan(fifty),3),2),:,:) = [];

%% Creates the matrix for each group
clase_exp(1).Cond_Mat_Pos = [small; thirty];
clase_exp(1).origen = [ones(1,size(small,1))*15 ones(1,size(thirty,1))*30]';
clase_exp(2).Cond_Mat_Pos = [forty; fifty];
clase_exp(2).origen = [ones(1,size(forty,1))*45 ones(1,size(fifty,1))*50]';

for exp=1:2
    for c=1:size(clase_exp(exp).Cond_Mat_Pos,3)
        clase_exp(exp).condition(c).n=size(clase_exp(exp).Cond_Mat_Pos,1);
        clase_exp(exp).condition(c).serie_prom_pos(1,:)=mean(clase_exp(exp).Cond_Mat_Pos(:,:,c));
        clase_exp(exp).condition(c).serie_prom_pos(2,:)=[-10:14];
        clase_exp(exp).condition(c).serie_prom_pos_std(1,:)=std(clase_exp(exp).Cond_Mat_Pos(:,:,c));
        clase_exp(exp).condition(c).serie_prom_pos_std(2,:)=[-10:14];
        clase_exp(exp).condition(c).serie_prom_pos_ee(1,:)=clase_exp(exp).condition(c).serie_prom_pos_std(1,:)/sqrt(clase_exp(exp).condition(c).n);
        clase_exp(exp).condition(c).serie_prom_pos_ee(2,:)=[-10:14];
    end
end

%% ----------------------- PreBASELINE -----------------------------------
%% Edito las matrices originales
small=exp1.clase_pert(4).Cond_Mat;
small(any(any(isnan(small),3),2),:,:) = [];
thirty=exp1.clase_pert(2).Cond_Mat;
thirty(any(any(isnan(thirty),3),2),:,:) = [];
forty=exp1.clase_pert(3).Cond_Mat;
forty(any(any(isnan(forty),3),2),:,:) = [];
fifty=exp2.clase_pert(2).Cond_Mat;
fifty(any(any(isnan(fifty),3),2),:,:) = [];

%% Armo y las matrices
clase_exp(1).Cond_Mat = [small; thirty];
clase_exp(1).origen = [ones(1,size(small,1))*15 ones(1,size(thirty,1))*30]';
clase_exp(2).Cond_Mat = [forty; fifty];
clase_exp(2).origen = [ones(1,size(forty,1))*45 ones(1,size(fifty,1))*50]';

for exp=1:2
    for c=1:size(clase_exp(exp).Cond_Mat,3)
        clase_exp(exp).condition(c).n=size(clase_exp(exp).Cond_Mat,1);
        clase_exp(exp).condition(c).serie_prom(1,:)=mean(clase_exp(exp).Cond_Mat(:,:,c));
        clase_exp(exp).condition(c).serie_prom(2,:)=[-10:14];
        clase_exp(exp).condition(c).serie_prom_std(1,:)=std(clase_exp(exp).Cond_Mat(:,:,c));
        clase_exp(exp).condition(c).serie_prom_std(2,:)=[-10:14];
        clase_exp(exp).condition(c).serie_prom_ee(1,:)=clase_exp(exp).condition(c).serie_prom_std(1,:)/sqrt(clase_exp(exp).condition(c).n);
        clase_exp(exp).condition(c).serie_prom_ee(2,:)=[-10:14];
    end
end

save(['average_by_experiment_' date],'clase_exp')
