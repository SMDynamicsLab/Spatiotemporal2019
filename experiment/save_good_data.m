function [ o,trial ] = save_good_data(subject_number,data,subject_folder,exp,bad,step,j,o,stim,resp,asyn,error,t_acc,x_acc,y_acc,z_acc,t_pr,pr,mech_size,mech_bip,trial,B_laps,temp_size)

trial(j).order = o;
trial(j).trial = j;

for h=1:length(B_laps);
    if ismember(j,[B_laps(h)-9:B_laps(h)])==1
        trial(j).block=h;
    else
    end
end

trial(j).mech_size = mech_size;
trial(j).temp_size = temp_size;
trial(j).mech_bip = mech_bip;
trial(j).stim = stim;
trial(j).resp = resp;
trial(j).asyn = asyn;
if isempty(x_acc)==0
trial(j).x_acc(1,:) = x_acc;
trial(j).x_acc(2,:) = t_acc;
trial(j).y_acc(1,:) = y_acc;
trial(j).y_acc(2,:) = t_acc;
trial(j).z_acc(1,:) = z_acc;
trial(j).z_acc(2,:) = t_acc;
trial(j).pr(1,:) = pr;
trial(j).pr(2,:) = t_pr;
end




if step==1
    filename=['temporal_' exp '_ent.mat'];
    step_crudo='_ent';
elseif step==2
    if strcmp('alturas',exp)==1
    filename=['temporal_' exp '_exp.mat'];
    step_crudo='_exp';
    else
    filename=['temporal_' exp '_cal.mat'];
    step_crudo='_cal';
    end
elseif step==3
    filename=['temporal_' exp '_exp.mat'];
    step_crudo='_exp';
end

save(filename,'trial','bad')

% Nombre del archivo .dat
filename_dat = [ num2str(subject_number) '-' exp step_crudo '-trial_' num2str(trial(j).order) '-' num2str(trial(j).trial) '.dat'];
% Ruta completa del archivo .dat
crudo = fullfile(subject_folder,filename_dat);   
        
% Guarda cada trial en un .dat en la carpeta del subject
fid = fopen(crudo, 'wt');                 
fprintf(fid, '%s\n', data{:});           
fclose(fid);

o=o+1;
end