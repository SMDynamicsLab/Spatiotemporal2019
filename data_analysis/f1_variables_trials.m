function [subject]=f1_variables_trials(subject,pre_baseline_bips,pos_baseline_bips,out_limit)

for s=1:length(subject)
    
    if max(subject(s).exp_temp_sizes)<=15
        pre_bips=min(pre_baseline_bips);
        pos_bips=min(pos_baseline_bips);
    else
        pre_bips=max(pre_baseline_bips);
        pos_bips=max(pos_baseline_bips);
    end
    
    for t=1:length(subject(s).exp)
        
        c = subject(s).exp(t).condition;
        
        % indice del primer bip de la serie comun a todos los subjects:
        min_idx = find((subject(s).exp(t).asyn(2,:)-subject(s).exp(t).mech_bip-1)==-10);
        % indice del ultimo bip de la serie comun a todos los subjects:
        max_idx = find((subject(s).exp(t).asyn(2,:)-subject(s).exp(t).mech_bip-1)==14);
        % Serie reducida y relativizada al bip de la perturbacion
        ser_red=[subject(s).exp(t).asyn(1,min_idx:max_idx); subject(s).exp(t).asyn(2,min_idx:max_idx)-subject(s).exp(t).mech_bip-1];
        
        % Ubica a la perturbacion en el bip 0:
        pert_idx = find(ser_red(2,:)==0);
        
        pre_bs=mean(ser_red(1,pert_idx-pre_bips:pert_idx-1));
        pre_std=std(ser_red(1,pert_idx-pre_bips:pert_idx-1));
        pos_bs=mean(ser_red(1,end-pos_bips+1:end));
        pos_std=std(ser_red(1,end-pos_bips+1:end));
        
        new_serie = [ser_red(1,:)-pre_bs; ser_red(2,:)];
        new_serie2 = [ser_red(1,:)-pos_bs; ser_red(2,:)];
        
        subject(s).exp(t).asyn_red = ser_red;
        subject(s).exp(t).asyn_bs = new_serie;
        subject(s).exp(t).asyn_posbs = new_serie2;
                
        subject(s).exp(t).pre_bs = pre_bs;
        subject(s).exp(t).pre_bs_n = pre_bips;
        subject(s).exp(t).pre_bs_std = pre_std;
        subject(s).exp(t).pos_bs = pos_bs;
        subject(s).exp(t).pos_bs_n = pos_bips;
        subject(s).exp(t).pos_bs_std = pos_std;
        
        subject(s).exp(t).out_limit_value = out_limit;
        
        if isempty(find(abs(subject(s).exp(t).asyn_bs(1,:))>=out_limit))==1
            subject(s).exp(t).out_lim=0;
        else
            subject(s).exp(t).out_lim=1;
        end
        
    end
    
end

% save(['subjects_simetrico-' date '.mat'],'subject')

end