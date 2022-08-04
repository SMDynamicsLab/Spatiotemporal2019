function [clase_pert] = f5_average_between_subjects(subject,sub_out,clase_cuatro,parameters)
s=length(subject);
c=length(subject(s).condition);

max_trial=[subject(s).exp(1).asyn_red(2,1):subject(s).exp(1).asyn_red(2,end)]; % coordenadas del trial

% Matriz a completar con las series promedios de cada subject
% fila=subject,columna=#bip,matriz=#condition.
clase_pert(1).Cond_Mat=nan(s,length(max_trial),c);
clase_pert(2).Cond_Mat=nan(s,length(max_trial),c);
clase_pert(3).Cond_Mat=nan(s,length(max_trial),c);
clase_pert(4).Cond_Mat=nan(s,length(max_trial),c);
clase_pert(1).Cond_Mat_Pos=nan(s,length(max_trial),c);
clase_pert(2).Cond_Mat_Pos=nan(s,length(max_trial),c);
clase_pert(3).Cond_Mat_Pos=nan(s,length(max_trial),c);
clase_pert(4).Cond_Mat_Pos=nan(s,length(max_trial),c);

for s=1:length(subject)
    sub_class=0;
    if ismember(s,sub_out)==1 % No incluye los subjects excluidos
    else
        if ismember(max(subject(s).exp_temp_sizes),[clase_cuatro 30 45])==0
        else
            if ismember(max(subject(s).exp_temp_sizes),[clase_cuatro])==1
                p=4;
                clase_pert(p).sizes=clase_cuatro;
                for c=1:length(subject(s).condition);
                    clase_pert(p).Cond_Mat(s,:,c)=subject(s).condition(c).serie_prom(1,:);
                    clase_pert(p).Cond_Mat_Pos(s,:,c)=subject(s).condition(c).serie_prom_pos(1,:);
                end
                if max(subject(s).exp_temp_sizes)==15
                    p=1;
                    clase_pert(p).sizes=15;
                    for c=1:length(subject(s).condition);
                        clase_pert(p).Cond_Mat(s,:,c)=subject(s).condition(c).serie_prom(1,:);
                        clase_pert(p).Cond_Mat_Pos(s,:,c)=subject(s).condition(c).serie_prom_pos(1,:);
                    end
                end
            elseif max(subject(s).exp_temp_sizes)==30
                p=2;
                clase_pert(p).sizes=30;
            elseif max(subject(s).exp_temp_sizes)==45;
                p=3;
                clase_pert(p).sizes=45;
            end
            for c=1:length(subject(s).condition);
                clase_pert(p).Cond_Mat(s,:,c)=subject(s).condition(c).serie_prom(1,:);
                clase_pert(p).Cond_Mat_Pos(s,:,c)=subject(s).condition(c).serie_prom_pos(1,:);
            end
        end
    end
end

for p=1:4
    for c=1:length(subject(s).condition);
        
        % Serie promedio de la condition
        clase_pert(p).condition(c).serie_prom(1,:)=nanmean(clase_pert(p).Cond_Mat(:,:,c));
        clase_pert(p).condition(c).serie_prom(2,:)=max_trial;
        % Desv�o estandar de la serie promedio de la condition
        clase_pert(p).condition(c).serie_prom_std(1,:)=nanstd(clase_pert(p).Cond_Mat(:,:,c));
        clase_pert(p).condition(c).serie_prom_std(2,:)=max_trial;
        % cantidad de subjects incluidos en el promedio
        clase_pert(p).condition(c).n=max(cumsum(~isnan(clase_pert(p).Cond_Mat(:,21,c))));
        % Numero de subjects excluidos en el promedio
        clase_pert(p).condition(c).subject_out=find(isnan(clase_pert(p).Cond_Mat(:,21,1)));
        % Error estandar de la serie promedio de la condition
        clase_pert(p).condition(c).serie_prom_ee(1,:)=clase_pert(p).condition(c).serie_prom_std(1,:)/sqrt(clase_pert(p).condition(c).n);
        clase_pert(p).condition(c).serie_prom_ee(2,:)=max_trial;
        
        % POSTBASELINE
        % Serie promedio de la condition
        clase_pert(p).condition(c).serie_prom_pos(1,:)=nanmean(clase_pert(p).Cond_Mat_Pos(:,:,c));
        clase_pert(p).condition(c).serie_prom_pos(2,:)=max_trial;
        % Desv�o estandar de la serie promedio de la condition
        clase_pert(p).condition(c).serie_prom_pos_std(1,:)=nanstd(clase_pert(p).Cond_Mat_Pos(:,:,c));
        clase_pert(p).condition(c).serie_prom_pos_std(2,:)=max_trial;
        % cantidad de subjects incluidos en el promedio
        clase_pert(p).condition(c).n_pos=max(cumsum(~isnan(clase_pert(p).Cond_Mat_Pos(:,21,c))));
        % Error estandar de la serie promedio de la condition
        clase_pert(p).condition(c).serie_prom_pos_ee(1,:)=clase_pert(p).condition(c).serie_prom_pos_std(1,:)/sqrt(clase_pert(p).condition(c).n);
        clase_pert(p).condition(c).serie_prom_pos_ee(2,:)=max_trial;
        % Numero de subjects excluidos en el promedio
        clase_pert(p).condition(c).subject_out_pos=find(isnan(clase_pert(p).Cond_Mat_Pos(:,21,1)));
        
    end
    
end
% % %% Promedia los NMA entre subjects
% % 
% % clase_pert(1).NMA_mat=nan(3,length(subject),length(subject(s).condition));
% % clase_pert(2).NMA_mat=nan(3,length(subject),length(subject(s).condition));
% % clase_pert(3).NMA_mat=nan(3,length(subject),length(subject(s).condition));
% % clase_pert(4).NMA_mat=nan(3,length(subject),length(subject(s).condition));
% % 
% % for s=1:length(subject)
% %     
% %     if ismember(s,sub_out)==1 % No incluye los subjects excluidos
% %     else
% %         if ismember(max(subject(s).exp_temp_sizes),[clase_cuatro 30 45])==0
% %         else
% %             if ismember(max(subject(s).exp_temp_sizes),[clase_cuatro])==1
% %                 p=4;
% %                 clase_pert(p).sizes=clase_cuatro;
% %                 for c=1:length(subject(s).condition)
% %                     clase_pert(p).NMA_mat(1,s,c)=subject(s).condition(c).mean_prebaseline;
% %                     clase_pert(p).NMA_mat(2,s,c)=subject(s).condition(c).mean_posbaseline;
% %                     clase_pert(p).NMA_mat(3,s,c)=subject(s).condition(c).mean_posbaseline-subject(s).condition(c).mean_prebaseline;
% %                 end
% %                 if max(subject(s).exp_temp_sizes)==15
% %                     p=1;
% %                     clase_pert(p).sizes=15;
% %                     for c=1:length(subject(s).condition)
% %                         clase_pert(p).NMA_mat(1,s,c)=subject(s).condition(c).mean_prebaseline;
% %                         clase_pert(p).NMA_mat(2,s,c)=subject(s).condition(c).mean_posbaseline;
% %                         clase_pert(p).NMA_mat(3,s,c)=subject(s).condition(c).mean_posbaseline-subject(s).condition(c).mean_prebaseline;
% %                     end
% %                 end
% %             elseif max(subject(s).exp_temp_sizes)==30
% %                 p=2;
% %                 clase_pert(p).sizes=30;
% %             elseif max(subject(s).exp_temp_sizes)==45;
% %                 p=3;
% %                 clase_pert(p).sizes=45;
% %             elseif max(subject(s).exp_temp_sizes)>=min(clase_cuatro) && max(subject(s).exp_temp_sizes)<=max(clase_cuatro)
% %                 p=3;
% %                 clase_pert(p).sizes=clase_cuatro;
% %             end
% %             
% %             for c=1:length(subject(s).condition)
% %                 clase_pert(p).NMA_mat(1,s,c)=subject(s).condition(c).mean_prebaseline;
% %                 clase_pert(p).NMA_mat(2,s,c)=subject(s).condition(c).mean_posbaseline;
% %                 clase_pert(p).NMA_mat(3,s,c)=subject(s).condition(c).mean_posbaseline-subject(s).condition(c).mean_prebaseline;
% %             end
% %         end
% %     end
% %     
% % end
% % 
% % for p=1:4
% %     
% %     for c=1:length(subject(s).condition);
% %         clase_pert(p).condition(c).preNMA_prom=nanmean(clase_pert(p).NMA_mat(1,:,c));
% %         clase_pert(p).condition(c).preNMA_std=nanstd(clase_pert(p).NMA_mat(1,:,c));
% %         clase_pert(p).condition(c).posNMA_prom=nanmean(clase_pert(p).NMA_mat(2,:,c));
% %         clase_pert(p).condition(c).posNMA_std=nanstd(clase_pert(p).NMA_mat(2,:,c));
% %         clase_pert(p).condition(c).deltaNMA_prom=nanmean(clase_pert(p).NMA_mat(3,:,c));
% %         clase_pert(p).condition(c).deltaNMA_std=nanstd(clase_pert(p).NMA_mat(3,:,c));
% %         
% %     end
% %     
% % end


%% Grafica las series temporales promedio por condition

for p=1:4
    
    figure()
    
    for i=1:9
        s=i;
        [lineprops] = fg_Colores_presentacion(s);
        mseb(clase_pert(p).condition(i).serie_prom(2,:),clase_pert(p).condition(i).serie_prom(1,:),clase_pert(p).condition(i).serie_prom_ee(1,:),lineprops,0.1);
        hold all
        xlim([-10 15])      
        l=legend('-M-T','-M','-M+T','-T','iso','+T','+M-T','+M','+M+T');
        set(l,'Location','Southwest')
        thand = get(gca,'title');
        set(thand,'string',{'Mean Time Series (+ee)',['SC=' num2str(clase_pert(p).sizes) 'ms']},'fontsize',10);
        
        set(gca,'FontSize',8)
        xlhand = get(gca,'xlabel');
        set(xlhand,'string','step {\itn} relative to perturbation','fontsize',10);
        ylhand = get(gca,'ylabel');
        set(ylhand,'string','e_n [ms]','fontsize',10);

    end
end

save(['average_simetric-' date '.mat'],'clase_pert','parameters')

end