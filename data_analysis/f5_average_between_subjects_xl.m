function [clase_pert] = f5_average_between_subjects_xl(subject,sub_out,clase_dos,parameters)

s=length(subject);
c=length(subject(s).condition);

max_trial=[subject(s).exp(1).asyn_red(2,1):subject(s).exp(1).asyn_red(2,end)]; % coordenadas del trial

% Matriz a completar con las series promedios de cada sujeto
% fila=sujeto,columna=#bip,matriz=#condition.
clase_pert(1).Cond_Mat=nan(s,length(max_trial),c);
clase_pert(2).Cond_Mat=nan(s,length(max_trial),c);
clase_pert(1).Cond_Mat_Pos=nan(s,length(max_trial),c);
clase_pert(2).Cond_Mat_Pos=nan(s,length(max_trial),c);

for s=1:length(subject)
    pert=abs(subject(s).cal_asyn_max);
    if subject(s).out==0 % No incluye los sujetos outliers
        if ismember(pert,clase_dos)==1
            p=2;
            clase_pert(p).sizes=clase_dos;
            for c=1:length(subject(s).condition);
                clase_pert(p).Cond_Mat_Pos(s,:,c)=subject(s).condition(c).serie_prom_pos(1,:);
                clase_pert(p).Cond_Mat(s,:,c)=subject(s).condition(c).serie_prom(1,:);
            end
        end
        if pert==15
            p=1;
            clase_pert(p).sizes=15;
            for c=1:length(subject(s).condition);
                clase_pert(p).Cond_Mat_Pos(s,:,c)=subject(s).condition(c).serie_prom_pos(1,:);
                clase_pert(p).Cond_Mat(s,:,c)=subject(s).condition(c).serie_prom(1,:);
            end
        end
    end
end

for p=1:2
    for c=1:length(subject(s).condition);
        
        % Serie promedio de la condition
        clase_pert(p).condition(c).serie_prom(1,:)=nanmean(clase_pert(p).Cond_Mat(:,:,c));
        clase_pert(p).condition(c).serie_prom(2,:)=max_trial;
        % Desv�o estandar de la serie promedio de la condition
        clase_pert(p).condition(c).serie_prom_std(1,:)=nanstd(clase_pert(p).Cond_Mat(:,:,c));
        clase_pert(p).condition(c).serie_prom_std(2,:)=max_trial;
        % cantidad de sujetos incluidos en el promedio
        clase_pert(p).condition(c).n=max(cumsum(~isnan(clase_pert(p).Cond_Mat(:,21,c))));
        % Numero de sujetos excluidos en el promedio
        clase_pert(p).condition(c).sujeto_out=find(isnan(clase_pert(p).Cond_Mat(:,21,1)));
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
        % cantidad de sujetos incluidos en el promedio
        clase_pert(p).condition(c).n_pos=max(cumsum(~isnan(clase_pert(p).Cond_Mat_Pos(:,21,c))));
        % Error estandar de la serie promedio de la condition
        clase_pert(p).condition(c).serie_prom_pos_ee(1,:)=clase_pert(p).condition(c).serie_prom_pos_std(1,:)/sqrt(clase_pert(p).condition(c).n);
        clase_pert(p).condition(c).serie_prom_pos_ee(2,:)=max_trial;
        % Numero de sujetos excluidos en el promedio
        clase_pert(p).condition(c).sujeto_out_pos=find(isnan(clase_pert(p).Cond_Mat_Pos(:,21,1)));
        
    end
    
end

%% Grafica las series temporales promedio por condition

for p=1:2
    
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

save(['average_simetric_xl-' date '.mat'],'clase_pert','parameters')

end