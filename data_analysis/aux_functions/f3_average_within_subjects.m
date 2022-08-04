function [subject]=f3_average_within_subjects(subject)

for s=1:length(subject)
    
    clear ASYN
    clear POSASYN
    clear NORM
	clear PRE_BS
	clear POS_BS
	clear PRE_IRI
	clear POS_IRI
    clear meanpreIRI
    clear meanposIRI
    
	cont=ones(1,9);
	
	for t=1:length(subject(s).exp)
		
		c=subject(s).exp(t).condition;
		
		if subject(s).exp(t).out_std==1 || subject(s).exp(t).out_lim==1;
			
			subject(s).exp(t).out=1;
			
			meanpreIRI(t)=nan;
			meanposIRI(t)=nan;
			
			ASYN(cont(c),:,c)=nan(1,length(subject(s).exp(t).asyn_red(1,:)));
            POSASYN(cont(c),:,c)=nan(1,length(subject(s).exp(t).asyn_red(1,:)));
            NORM(cont(c),:,c)=nan(1,length(subject(s).exp(t).asyn_red(1,11:end)));
			PRE_BS(cont(c),c)=nan;
			POS_BS(cont(c),c)=nan;
			PRE_IRI(cont(c),c)=nan;
			POS_IRI(cont(c),c)=nan;
			cont(c) = cont(c) + 1 ;
			
		else
			
			subject(s).exp(t).out=0;
			
			cposIRI=1;cpreIRI=1;
			
			for i=1:length(subject(s).exp(t).resp)-1
				resp_n=subject(s).exp(t).resp(1,i);
				resp_n1=subject(s).exp(t).resp(1,i+1);
				if i<subject(s).exp(t).mech_bip-1
					preIRI(cpreIRI)=resp_n1-resp_n;
					cpreIRI=cpreIRI+1;
				elseif i>=subject(s).exp(t).mech_bip-1+5 %el posIRI se calcula luego del "overshoot"
					posIRI(cposIRI)=resp_n1-resp_n;
					cposIRI=cposIRI+1;
				end
			end
			
			meanpreIRI(t)=mean(preIRI);
			meanposIRI(t)=mean(posIRI);
			
			ASYN(cont(c),:,c)=subject(s).exp(t).asyn_bs(1,:);
            POSASYN(cont(c),:,c)=subject(s).exp(t).asyn_posbs(1,:);

            NORM(cont(c),:,c)=subject(s).exp(t).asyn_red(1,11:end)/(sqrt(sum(subject(s).exp(t).asyn_red(1,11:end).^2)));
			PRE_BS(cont(c),c)=subject(s).exp(t).pre_bs;
			POS_BS(cont(c),c)=subject(s).exp(t).pos_bs;
			PRE_IRI(cont(c),c)=meanpreIRI(t);
			POS_IRI(cont(c),c)=meanposIRI(t);
			cont(c) = cont(c) + 1 ;
			
		end
		
	end
	
	for c=1:9
		
		% Trials incluidos en la serie
		subject(s).condition(c).n=sum(~isnan(ASYN(:,11,c)));
		% Serie promedio
		subject(s).condition(c).serie_prom(1,:)=nanmean(ASYN(:,:,c));
		subject(s).condition(c).serie_prom(2,:)=subject(s).exp(t).asyn_bs(2,:);
        % Serie promedio POST
		subject(s).condition(c).serie_prom_pos(1,:)=nanmean(POSASYN(:,:,c));
		subject(s).condition(c).serie_prom_pos(2,:)=subject(s).exp(t).asyn_bs(2,:);
		% Desvio de la serie promedio
		subject(s).condition(c).serie_prom_pos_std(1,:)=nanstd(POSASYN(:,:,c));
		subject(s).condition(c).serie_prom_pos_std(2,:)=subject(s).exp(t).asyn_bs(2,:);	
        % Serie promedio normalizada
		subject(s).condition(c).serie_prom_norm(1,:)=nanmean(NORM(:,:,c));
		subject(s).condition(c).serie_prom_norm(2,:)=[0:14];
		% Desvio de la serie promedio
		subject(s).condition(c).serie_prom_std(1,:)=nanstd(ASYN(:,:,c));
		subject(s).condition(c).serie_prom_std(2,:)=subject(s).exp(t).asyn_bs(2,:);
        % Desvio de la serie promedio normalizada
		subject(s).condition(c).serie_prom_std_norm(1,:)=nanstd(NORM(:,:,c));
		subject(s).condition(c).serie_prom_std_norm(2,:)=[0:14];
		% Pre baseline promedio de la condition (deber�a ser cercana a cero);
		subject(s).condition(c).mean_prebaseline=nanmean(PRE_BS(:,c));
		% Desvio standard del pre_bs promedio de la condition
		subject(s).condition(c).std_prebaseline=nanstd(PRE_BS(:,c));
		% Post baseline promedio de la condition
		subject(s).condition(c).mean_posbaseline=nanmean(POS_BS(:,c));
		% Desvio standard del post_bs promedio de la condition
		subject(s).condition(c).std_posbaseline=nanstd(POS_BS(:,c));
		% preIRI promedio de la condition
		subject(s).condition(c).mean_preIRI=nanmean(PRE_IRI(:,c));
		% Desvio standard del post IRI de la condition
		subject(s).condition(c).std_preIRI=nanstd(PRE_IRI(:,c));
		% posIRI promedio de la condition
		subject(s).condition(c).mean_posIRI=nanmean(POS_IRI(:,c));
		% Desvio standard del postIRI promedio de la condition
		subject(s).condition(c).std_posIRI=nanstd(POS_IRI(:,c));
		
	end
	
	subject(s).NMA=nanmean(PRE_BS(:,:));
	
end

% save(['subjects_simetrico-' date '.mat'],'subject')
end