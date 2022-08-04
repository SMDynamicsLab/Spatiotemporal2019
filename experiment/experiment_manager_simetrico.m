%% Experiment Manager
% Works with 2 arduinos: master_tone_ingka.ino y slave_servo.ino
% Test different fits of the temporal asynchronies for each step

clear all;
delete(instrfind);

%% Experiment Parameters

% Number of bips for trial.
N_stim_ent = 25;
N_stim_cal = 25;
N_stim_exp = 35;

% Repetitions by condition.
n=8;
% n=1;
n_calibracion=8;
% n_calibracion=1;
n_entrenamiento=1;

% Calibration Parameters
out_limit=150;

% Mechanical perturbations
mech_sizes = [33 95]; % size of perturbations
servo_ini_pos=57;
mech_bip_range_ent = [15 18];    % range of bips
mech_bip_range_cal = [15 18];    % range of bips
mech_bip_range_exp = [15 20];    % range of bips
% Step down categories
cat_M=[15 30 45];

cond_mech = max(size(mech_sizes));

%% 1. Register subject data
dir=pwd;
exp='simetrico'; % sufix of the .mat file
filename=['subjects_' exp '.mat'];
[ subject,subject_number ] = f_DatosSujeto(dir, exp);
% Crea la carpeta del subject
subject_folder=fullfile(pwd,'data',num2str(subject_number));
mkdir(subject_folder)

%% Entrenamiento
tent = fopen('temporal_simetrico_ent','w'); % Abre archivos temporales
N_stim = N_stim_ent;
mech_bip_range = mech_bip_range_ent;
temp_sizes=[-30 30];
step=1;
conditions='full';
[duration,trial,bad,mech_sizes,temp_sizes] = loop_central_simetrico(subject_number,subject_folder,conditions,exp,temp_sizes,step,N_stim,n,n_entrenamiento,n_calibracion,servo_ini_pos,cond_mech,mech_sizes,mech_bip_range);
% Save data
subject(subject_number).ent=trial;
subject(subject_number).ent_mech_sizes=mech_sizes;
subject(subject_number).ent_temp_sizes=temp_sizes;
subject(subject_number).ent_duration=duration;
clear trial

if  isempty(bad) == 0
    subject(subject_number).ent_bad=bad;
    clear bad
else
end

% disp('Fin del entrenamiento. Presione una tecla para comenzar con el experimento');
disp('End of the training. Press a key to begin the experiment');
disp(' ');
pause()
fclose(tent);

%% Save data
save(filename,'subject')

%% Calibration
tcal = fopen('temporal_simetrico_cal','w');
step=2;
N_stim = N_stim_cal;
mech_bip_range = mech_bip_range_cal;
temp_sizes=0;
mech_sizes=mech_sizes(mech_sizes~=57);
conditions='mech';
[duration,trial,bad,mech_sizes,temp_sizes] = loop_central_simetrico(subject_number,subject_folder,conditions,exp,temp_sizes,step,N_stim,n,n_entrenamiento,n_calibracion,servo_ini_pos,cond_mech,mech_sizes,mech_bip_range);
% Save data
subject(subject_number).cal=trial;
subject(subject_number).cal_mech_sizes=mech_sizes;
subject(subject_number).cal_temp_sizes=temp_sizes;
subject(subject_number).cal_duration=duration;
clear trial

if  isempty(bad) == 0
    subject(subject_number).cal_bad=bad;
    clear bad
else
end

% disp('Fin de la primera parte.! Presione una tecla para continuar');
disp('End of the first part. Press a key to continue');
disp(' ');
pause()
fclose(tcal);

%% Experiment
texp = fopen('temporal_simetrico_exp','w');
step=3;
N_stim = N_stim_exp;
mech_bip_range = mech_bip_range_exp;
[temp_sizes,esc_sizes,fiteo_up,fiteo_down] = parameters_calibration_exp(cat_M,subject_number,subject,out_limit,mech_sizes);
subject(subject_number).cal_fiteo_up=fiteo_up;
subject(subject_number).cal_fiteo_down=fiteo_down;
altura=esc_sizes;
[ angle ] = hei2ang_servo( height );
mech_sizes=angle;
conditions='full';
[duration,trial,bad,mech_sizes,temp_sizes] = loop_central_simetrico(subject_number,subject_folder,conditions,exp,temp_sizes,step,N_stim,n,n_entrenamiento,n_calibracion,servo_ini_pos,cond_mech,mech_sizes,mech_bip_range);
% Save data
subject(subject_number).exp=trial;
subject(subject_number).exp_mech_sizes=mech_sizes;
subject(subject_number).exp_temp_sizes=temp_sizes;
subject(subject_number).exp_duration=duration;
clear trial

if  isempty(bad) == 0
    subject(subject_number).exp_bad=bad;
    clear bad
else
end

% disp('Fin del experimento.�Muchas gracias por participar!');
disp('End of experiment. Thank you!');
disp(' ');

%% Save data in the .mat file
% Clasiffies trials conditions
[ subject ] = f_assign_condition( subject );
save(filename,'subject')
fclose(texp);
%% Saves the subject data in his own folder
% mkdir(fullfile('data',num2str(subject_number)))
movefile('temporal_simetrico_ent.mat',fullfile('data',num2str(subject_number),[num2str(subject_number) '_data_simetrico_ent.mat']));
movefile('temporal_simetrico_cal.mat',fullfile('data',num2str(subject_number),[num2str(subject_number) '_data_simetrico_cal.mat']));
movefile('temporal_simetrico_exp.mat',fullfile('data',num2str(subject_number),[num2str(subject_number) '_data_simetrico_exp.mat']));