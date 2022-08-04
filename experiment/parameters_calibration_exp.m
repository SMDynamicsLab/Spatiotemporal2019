function [temp_sizes,esc_sizes,fiteo_up,fiteo_down] = parameters_calibration_exp(cat_M,subject_number,subject,out_limit,mech_sizes)
%% Calibration Data Fit
% Calibration with two linear functions in the extreme positions

s=subject_number;
n_con=length(mech_sizes)+1; % # number of conditions
n_cal=length(subject(s).cal); % # number of calibration trials
pre_baseline_bips=10; % # number of bips included for the pre baseline

ASYN=nan(n_cal,n_con); % Matrix to fill. row: trial, col: ccondition
%1)-M, 2)+M y 3)preBaseline.

for t=1:length(subject(s).cal)
    pert_bip=subject(s).cal(t).mech_bip; % #bip perturbado
    pert_indx=find(subject(s).cal(t).asyn(2,:)==pert_bip); %perturbated bip index
    % pre baseline computed as the mean value of the asynchronies previous to the perturbation bip 
    pre_baseline=mean(subject(s).cal(t).asyn(1,pert_indx-pre_baseline_bips:pert_indx));
    c=find(subject(s).cal(t).mech_size==mech_sizes); % experimental condition
    % When all the trial asynchronies are lower than the outlier limit (asyn-baseline>out_limit)
    % fill the matrix with the value of asyn-bs (-M y +M) and bs (M=0)
            if isempty(find(abs(subject(s).cal(t).asyn(1,pert_indx-10:pert_indx)-pre_baseline)>=out_limit))==1

        ASYN(t,c)=subject(s).cal(t).asyn(1,pert_indx+1)-pre_baseline;
        ASYN(t,n_con)=pre_baseline;
        ASYN(t,n_con+1)=t;
    end
end

asyn=nanmean(ASYN); % mean asynchronies vector

% Linear fit
fitType='poly1'; % Fit Function

x_down=[ang2alt_servo(min(mech_sizes)) 0]';
y_down=[asyn(1) 0]';

x_up=[0 ang2alt_servo(max(mech_sizes))]';
y_up=[0 asyn(2)]';

[fiteo_up,gof_up] = fit(x_up,y_up,fitType);
[fiteo_down,gof_down] = fit(x_down,y_down,fitType);

% limit asynchrony
if abs(fiteo_up.p1)>=abs(fiteo_down.p1) %USUALLY!! If m_up > m_down
    asyn_max=fiteo_down.p1*x_down(1); % asynchrony -M
    % asynchrony category
    for k=length(cat_M):-1:1
        if asyn_max>=cat_M(k)
            asyn_max=cat_M(k);
            break
        end
    end
    
    esc_max=asyn_max/fiteo_down.p1;
    
    asyn_sim=-asyn_max; % symmetric asynchrony
    esc_sim=asyn_sim/fiteo_up.p1;
    
else % When m_up < m_down!
    asyn_max=fiteo_up.p1*x_up(2); % asynchrony +M  
    
    % asynchrony category
    for k=length(cat_M):-1:1
        if abs(asyn_max)>=cat_M(k)
            asyn_max=-cat_M(k);
            break
        end
    end
    
    esc_max=asyn_max/fiteo_up.p1;
    
    asyn_sim=-asyn_max; % symmetric asynchrony
    esc_sim=asyn_sim/fiteo_down.p1;
end

    if abs(asyn_max)==15
        temp_sizes=round([asyn_max asyn_sim]);
        esc_sizes=[esc_max esc_sim];
    elseif abs(asyn_max)<15
        temp_sizes='error';
    else
        temp_sizes=[-15 15];
        if abs(fiteo_up.p1)>=abs(fiteo_down.p1)
        esc_sizes=[asyn_max/fiteo_down.p1 asyn_sim/fiteo_up.p1];
        else
            esc_sizes=[asyn_max/fiteo_up.p1 asyn_sim/fiteo_down.p1];
        end
    end
   subject(s).exp_temp_sizes_final=temp_sizes;
subject(s).exp_mech_sizes_final=esc_sizes; 

end