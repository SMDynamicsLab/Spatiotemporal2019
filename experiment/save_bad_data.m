function [ o,r,bad ] = save_bad_data(subject_number,data,subject_folder,exp,trial,step,r,o,stim,resp,asyn,error,t_acc,x_acc,y_acc,z_acc,t_pr,pr,mech_size,mech_bip,bad,B_laps,temp_size)

bad(r).order = o;
bad(r).trial = r;

for h=1:length(B_laps);
    if ismember(r,[B_laps(h)-9:B_laps(h)])==1
        bad(r).block=h;
    else
    end
end

bad(r).mech_size = mech_size;
bad(r).temp_size = temp_size;
bad(r).mech_bip = mech_bip;
bad(r).stim = stim;
bad(r).resp = resp;
bad(r).asyn = asyn;
bad(r).type = error;
if isempty(x_acc)==0
bad(r).x_acc(1,:) = x_acc;
bad(r).x_acc(2,:) = t_acc;
bad(r).y_acc(1,:) = y_acc;
bad(r).y_acc(2,:) = t_acc;
bad(r).z_acc(1,:) = z_acc;
bad(r).z_acc(2,:) = t_acc;
bad(r).pr(1,:) = pr;
bad(r).pr(2,:) = t_pr;
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
filename_dat = [ num2str(subject_number) '-' exp step_crudo '-bad_' num2str(bad(r).order) '-' num2str(bad(r).trial) '.dat'];
% Ruta completa del archivo .dat
crudo = fullfile(subject_folder,filename_dat);   
        
% Guarda cada trial en un .dat en la carpeta del subject
fid = fopen(crudo, 'wt');                 
fprintf(fid, '%s\n', data{:});           
fclose(fid);

o=o+1;
r=r+1;

end

