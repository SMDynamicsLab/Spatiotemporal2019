%% Figures Manager
% This scripts generates the figures published in López, S.L., Laje, R. 
% Spatiotemporal perturbations in paced finger tapping suggest a common 
% mechanism for the processing of time errors. Sci Rep 9, 17814 (2019). 
% https://doi.org/10.1038/s41598-019-54133-x

%% Clean start
clear all
close all

%% General figures settings
fsize=7;
ftype='Helvetica';
fsize_leg=6;
fsize_letra=7;
msize = 2;
lwidth=1;
comparaciones=5;
alpha=0.05;
grey=[0.8267 0.8267 0.8267];
colores=ametrine(5);

%% Data loading

load('average_by_experiment_02-Aug-2022.mat')

filename_in_1_subj = 'subject_simetric-02-Aug-2022.mat';
aux_1 = load(filename_in_1_subj);
filename_in_2_subj = 'subject_simetric_xl-02-Aug-2022.mat';
aux_2 = load(filename_in_2_subj);

%% Figure 2.
% Resynchronization time series for every condition averaged across subjects.

% Settings
fig_size_cm = [14 5];
comparaciones=5;
alpha=0.05;
colores=ametrine(9);
letras={'a)','b)'};
titles={'Small perturbations','Large perturbations'};
lgd_str={'-S-T','-S','-S+T','-T','iso','+S','+S-T','+S','+S+T'};

figure(2);
clf(2);

for pan=1:2
    
    exp=pan;
    y_lim_inf = -70;
    y_lim_sup = 70;
    y_lim_inf_tick = -60;
    y_lim_sup_tick = 60;
    y_int = 20;
    
    subplot(1,2,pan)
    plot([-2 14],[0 0],'k-','handlevisibility','off');
    hold on;
    plot([0 0],[y_lim_inf y_lim_sup],'k-','handlevisibility','off');
   
    for i=1:9
        lineprops.col={colores(i,:)};
        lineprops.width=lwidth;
        lineprops.style='-';
        ptm(i) = mseb(clase_exp(exp).condition(i).serie_prom(2,:),clase_exp(exp).condition(i).serie_prom(1,:),clase_exp(exp).condition(i).serie_prom_ee(1,:),lineprops,1);
        hold all;
        
    end
    
    % Axes
    xlabel('Beep {\itn} relative to perturbation','Fontsize',fsize,'FontName',ftype)
    xlim([-2 14])
    set(gca,'XTick',[0:5:14])
    set(gca,'FontSize',fsize);
    
    ylim([y_lim_inf y_lim_sup])
    ylabel( 'e_n (ms)','Fontsize',fsize,'FontName',ftype)
    set(gca,'YTick',[y_lim_inf_tick:y_int:y_lim_sup_tick])
    set(gca,'FontSize',fsize);
    
    % Title
    title(titles(pan),'fontweight','bold','Fontsize',fsize,'FontName',ftype)
      
    if pan==1
        %Legend
        lh = columnlegend_RL(3,lgd_str,'padding',-0.2,'fontsize',fsize_leg,'location','Southeast');
        set(lh,'fontsize',5);
        legend('boxoff');
        
    end    
    box off;
end

% Letters
l_pos=[0.07 0.98 0 0; 0.51 0.98 0 0 ];

for i=1:length(letras)
    l=annotation('textbox', l_pos(i,:), 'string',letras(i));
    l_g=get(l);
    set(l,'FontSize',fsize_letra,'FontWeight','bold');
end

% Print
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

% print('-dpng','-r600','figure_all.png')

%% Figure 3. 
% The response to large temporal perturbations is asymmetric.

fig_size_cm = [14 10];
patrones=[4 4 2 2];
invertidas=[6 6 8 8];
clases=[1 2 1 2];
simms=[1 1 2 2];
letras={'a)','b)','c)','d)'};

figure(3);
clf(3);

for pan=1:4
    
    pat=patrones(pan);
    inv=invertidas(pan);
%     sim=simms(pan);
    clase=clases(pan);
    contS=1;
    
    %% Carga los datos del bootstrap para el FDR
    if pan==1 || pan==2
        sim = 1;
    else
        sim = 2;
    end
    
    filename=['clase',num2str(clase),'_sim-',num2str(sim),'.mat'];
    load(filename);
      
%     % FDR
%         p_test=p_val_emp(sim,2:comparaciones+1);
%         p_test(p_test==0)=0.001;
% 
%     [h,crit_p,adj_ci_cvrg,adj_p]=fdr_bh(p_test,0.05,'pdep','yes');
%     sig_bips=find(adj_p<alpha);
    
    % FDR
        if pan==2
        p_test=p_val_emp(1,2:comparaciones+1);
        p_test(p_test==0)=0.001;
    else
       p_test=p_val_emp(1,2:comparaciones+1);
    end
    [h,crit_p,adj_ci_cvrg,adj_p]=fdr_bh(p_test,0.05,'pdep','yes');
    sig_bips=find(adj_p<alpha);

    
    %% Grafico
    
    % Colores
    if pan<=2
        color=colores(1:2,:);
    elseif pan>2
        color=colores(3:4,:);
    end
        
    subplot(2,2,pan)

if pan==3 || pan==4
    y_lim_inf = -25;
    y_lim_sup = 25;
    y_lim_inf_tick = -20;
    y_lim_sup_tick = 20;
    y_int = 10;
elseif pan==1 || pan==2
    y_lim_inf = -60;
    y_lim_sup = 60;
    y_lim_inf_tick = -60;
    y_lim_sup_tick = 60;
    y_int = 20;
end
    %Rectangulo
    if pan==2
    ptr = rectangle('Position',[sig_bips(1)-0.5 y_lim_inf+1 sig_bips(end)-sig_bips(1)+1 y_lim_sup*2],'FaceColor',grey,'EdgeColor',grey);
    end
    hold on;
    plot([-2 14],[0 0],'k-');
    plot([0 0],[y_lim_inf y_lim_sup],'k-');

    lineprops(1).col={color(1,:)};
    lineprops(1).width=lwidth;
    lineprops(1).style='-';
    % Serie patron
    x_p=clase_exp(clase).condition(pat).serie_prom_pos(2,:);
    y_p=clase_exp(clase).condition(pat).serie_prom_pos(1,:);
    e_p=clase_exp(clase).condition(pat).serie_prom_pos_ee(1,:);
    ptm1 = mseb(x_p,y_p,e_p,lineprops(1),1);
    % Serie invertida
    x_i = clase_exp(clase).condition(inv).serie_prom_pos(2,:);
    y_i = clase_exp(clase).condition(inv).serie_prom_pos(1,:);
    e_i = clase_exp(clase).condition(inv).serie_prom_pos_ee(1,:);
    lineprops(2).col={color(2,:)};
    lineprops(2).width=lwidth;
    lineprops(2).style='-';
    ptm2 = mseb(x_i,y_i,e_i,lineprops(2),1);
    % Resta
    x_r = clase_exp(clase).condition(inv).serie_prom_pos(2,:);
    y_r = y_p+y_i;
    e_r = sqrt((e_p).^2+(e_i).^2);
    lineprops(3).col={[0 0 0]};
    lineprops(3).width=lwidth;
    lineprops(3).style='--';
    ptm3 = mseb(x_i,y_r,e_i,lineprops(3),1);
    
if pan==1
    y_lab = ylabel( {'{\bf{Temporal perturbations}}','','e_n (ms)'},'Fontsize',fsize,'FontName',ftype);
    set(y_lab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
elseif pan==3
   y_lab = ylabel( {'{\bf{Spatial perturbations}}','','e_n (ms)'},'Fontsize',fsize,'FontName',ftype);
   set(y_lab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
end

% Perturbation type
if pan==1
title( 'Small perturbations','fontweight','bold','Fontsize',fsize,'FontName',ftype)
elseif pan==2
title( 'Large perturbations','fontweight','bold','Fontsize',fsize,'FontName',ftype)
end

% X axes
xlim([-2 14])
set(gca,'XTick',[0:5:14])
set(gca,'FontSize',fsize);
if pan>2
xlabel( 'Beep {\itn} relative to perturbation','Fontsize',fsize,'FontName',ftype)
end

% Y axes
ylim([y_lim_inf y_lim_sup])
set(gca,'YTick',[y_lim_inf_tick:y_int:y_lim_sup_tick])
set(gca,'FontSize',fsize);

    %Asteriscos
    for bip=1:length(adj_p)
        if adj_p(1,bip)<alpha
            plot(bip,y_lim_inf+5,'*k','MarkerSize',msize)
        end
    end

% Legend
if pan==1 || pan==2
    lh = legend([ptm1.mainLine ptm2.mainLine ptm3.mainLine],'-T','+T','ASYM');
elseif pan==3 || pan==4
    lh = legend([ptm1.mainLine ptm2.mainLine ptm3.mainLine],'-S','+S','ASYM');
end
    set(lh,'Location','Southeast');
    leg_pos = get(lh,'position');
    aux = ds2nfu([0 y_lim_inf+0.5*y_int 0 0]);
    leg_pos(2) = aux(2);
    set(lh,'position',leg_pos);
    legend('boxoff')

    box off;
end

% Letters
l_pos=[0.05 0.95 0 0; 0.49 0.95 0 0 ;0.05 0.48 0 0;0.49 0.48 0 0];
for i=1:length(letras)
    
    l=annotation('textbox', l_pos(i,:), 'string',letras(i));
    l_g=get(l);
    set(l,'FontSize',fsize_letra,'FontWeight','bold');
    
end

% Save
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

% print('-dpng','-r600','figure_simples.png');

%% Figure 4. 
% Large combined perturbations produce asymmetric responses.

patrones=[1 1 3 3];
invertidas=[9 9 7 7];
clases=[1 2 1 2];
letras={'a)','b)','c)','d)'};
simms=[4 4 3 3];
titles={'Small perturbations','Large perturbations'};

figure(4);
clf(4);

for pan=1:4
    
    pat=patrones(pan);
    inv=invertidas(pan);
    clase=clases(pan);
    sim=simms(pan);
    contS=1;    
    
    %% Carga los datos del bootstrap para el FDR
    if pan==1 || pan==2
        sim = 4;
    else
        sim = 3;
    end
    
    filename=['clase',num2str(clase),'_sim-',num2str(sim) '.mat'];
    load(filename);
    
    % FDR
    [h,crit_p,adj_ci_cvrg,adj_p]=fdr_bh(p_val_emp(1,2:comparaciones+1),0.05,'pdep','yes');
    sig_bips=find(adj_p<alpha);

    %% Grafico
    
    % Colores
    if pan<=2
        color=colores(1:2,:);
    elseif pan>2
        color=colores(3:4,:);
    end
    
    subplot(2,2,pan)

    % Eje Y
    y_int=20;
    y_lim_inf_tick = -60;
    y_lim_sup_tick = 60;
    if pan==1 || pan==2
        y_lim_inf = -70;
        y_lim_sup = 70;
    elseif pan==3 || pan==4
        y_lim_inf = -60;
        y_lim_sup = 60;
    end
    
    %Rectangulo
    if pan==2
        rectangle('Position',[sig_bips(1)-0.5 y_lim_inf+1 sig_bips(end)-sig_bips(1)+1 y_lim_sup*2],'FaceColor',grey,'EdgeColor',grey)
    end
    hold on;
    plot([-2 14],[0 0],'k-');
    plot([0 0],[y_lim_inf y_lim_sup],'k-');
       
    lineprops(1).col={color(1,:)};
    lineprops(1).width=lwidth;
    lineprops(1).style='-';
    % Serie patron
    x_p=clase_exp(clase).condition(pat).serie_prom_pos(2,:);
    y_p=clase_exp(clase).condition(pat).serie_prom_pos(1,:);
    e_p=clase_exp(clase).condition(pat).serie_prom_pos_ee(1,:);
    ptm1 = mseb(x_p,y_p,e_p,lineprops(1),1);
    hold on
    % Serie invertida
    x_i = clase_exp(clase).condition(inv).serie_prom_pos(2,:);
    y_i = clase_exp(clase).condition(inv).serie_prom_pos(1,:);
    e_i = clase_exp(clase).condition(inv).serie_prom_pos_ee(1,:);

    lineprops(2).col={color(2,:)};
    lineprops(2).width=lwidth;
    lineprops(2).style='-';
    ptm2 = mseb(x_i,y_i,e_i,lineprops(2),1);
    % Resta
    x_r = clase_exp(clase).condition(inv).serie_prom_pos(2,:);
    y_r = y_p+y_i;
    e_r = sqrt((e_p).^2+(e_i).^2);

    lineprops(3).col={[0 0 0]};
    lineprops(3).width=lwidth;
    lineprops(3).style='--';
    ptm3 = mseb(x_i,y_r,e_i,lineprops(3),1);
    
    if pan>2
        xlabel( 'Beep {\itn} relative to perturbation','Fontsize',fsize,'FontName',ftype)
    end
    
    if pan==1
        y_lab = ylabel( {'{\bf{Analogous perturbations}}','','e_n (ms)'},'Fontsize',fsize,'FontName',ftype);
        set(y_lab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    elseif pan==3
        y_lab = ylabel( {'{\bf{Opposite perturbations}}','','e_n (ms)'},'Fontsize',fsize,'FontName',ftype);
        set(y_lab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    end
    
    % Title
    if pan==1 || pan==2
    title(titles(pan),'fontweight','bold','Fontsize',fsize,'FontName',ftype)
    end
    
    % Eje X
    xlim([-2 14])
    set(gca,'XTick',[0:5:14])
    set(gca,'FontSize',fsize);
    
    ylim([y_lim_inf y_lim_sup])
    set(gca,'YTick',[y_lim_inf_tick:y_int:y_lim_sup_tick])
    set(gca,'FontSize',fsize);
    
    %Asteriscos
    for bip=1:length(adj_p)
        if adj_p(1,bip)<alpha
            plot(bip,y_lim_inf+8,'*k','MarkerSize',msize)
        end
    end
    
    % Legend
    if pan==1 || pan==2
        lh = legend([ptm1.mainLine ptm2.mainLine ptm3.mainLine],'-S-T','+S+T','ASYM');
    elseif pan==3 || pan==4
        lh = legend([ptm1.mainLine ptm2.mainLine ptm3.mainLine],'-S+T','+S-T','ASYM');
    end
    set(lh,'Location','Southeast')
    leg_pos = get(lh,'position');
    aux = ds2nfu([0 y_lim_inf+0.5*y_int 0 0]);
    leg_pos(2) = aux(2);
    set(lh,'position',leg_pos);
    legend('boxoff')
    
   
    box off;
    
end

% Letters
l_pos=[0.05 0.95 0 0; 0.49 0.95 0 0 ;0.05 0.48 0 0;0.49 0.48 0 0];
for i=1:length(letras)
    
    l=annotation('textbox', l_pos(i,:), 'string',letras(i));
    l_g=get(l);
    set(l,'FontSize',fsize_letra,'FontWeight','bold');
    
end

% Save
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

% print('-dpng','-r600','figure_combined.png');

%% Figure 5.
% Simple perturbations are additive.

fig_size_cm = [14 20];

combinadas=[1 1 3 3 7 7 9 9];
mecanicas=[2 2 2 2 8 8 8 8];
temporales=[4 4 6 6 4 4 6 6];
clases=[1 2 1 2 1 2 1 2];
letras={'a)','b)','c)','d)','e)','f)','g)','h)'};

figure(5)
clf(5)

for pan=1:8
    
    comb=combinadas(pan);
    mec=mecanicas(pan);
    temp=temporales(pan);
    clase=clases(pan);
    contS=1;
    
    %Colores
    if pan<3
    color=colores(1,:);
    elseif pan==3 || pan==4
    color=colores(2,:);
    elseif pan==5 || pan==6
        color=colores(3,:);
    elseif pan==7 || pan==8
         color=colores(4,:);
    end
    subplot(4,2,pan)
    
    % Eje Y
    y_int = 20;
    if pan==1 || pan==2
        y_lim_inf = -30;
        y_lim_sup = 70;
        y_lim_inf_tick = -20;
        y_lim_sup_tick = 60;
    elseif pan==3 || pan==4
        y_lim_inf = -60;
        y_lim_sup = 40;
        y_lim_inf_tick = -60;
        y_lim_sup_tick = 40;
    elseif pan==5 || pan==6
        y_lim_inf = -30;
        y_lim_sup = 60;
        y_lim_inf_tick = -20;
        y_lim_sup_tick = 60;
    elseif pan==7 || pan==8
        y_lim_inf = -75;
        y_lim_sup = 30;
        y_lim_inf_tick = -60;
        y_lim_sup_tick = 20;
    end

    plot([-2 14],[0 0],'k-');
    hold on;
    plot([0 0],[y_lim_inf y_lim_sup],'k-');
    
    %Suma experimental
    lineprops(1).col={color};
    lineprops(1).width=lwidth;
    lineprops(1).style='-';
    x_exp=clase_exp(clase).condition(comb).serie_prom_pos(2,:);
    y_exp=clase_exp(clase).condition(comb).serie_prom_pos(1,:);
    e_exp=clase_exp(clase).condition(comb).serie_prom_pos_ee(1,:);
    ptm1 = mseb(x_exp,y_exp,e_exp,lineprops(1),1);
%     hold on
    % Suma aritmetica
    x_sum=clase_exp(clase).condition(comb).serie_prom_pos(2,:);
    y_sum=clase_exp(clase).condition(mec).serie_prom_pos(1,:)+clase_exp(clase).condition(temp).serie_prom_pos(1,:);
    e_sum=sqrt(clase_exp(clase).condition(mec).serie_prom_pos_ee(1,:).^2+clase_exp(clase).condition(temp).serie_prom_pos_ee(1,:).^2);
    lineprops(2)=lineprops(1);
    lineprops(2).style=':';
    ptm2 = mseb(x_sum,y_sum,e_sum,lineprops(2),1);
    % Resta
%     [lineprops] = fg_Colores_presentacion(20);
    lineprops(3)=lineprops(1);
    lineprops(3).style='--';
    lineprops(3).col={[0 0 0]};
    x_resta=clase_exp(clase).condition(comb).serie_prom_pos(2,:);
    y_resta=y_exp-y_sum;
    e_resta=sqrt((e_exp).^2+(e_sum).^2);
    ptm3 = mseb(x_resta,y_resta,e_resta,lineprops(3),1);
    
    if pan>6
        xlabel( 'Beep {\itn} relative to perturbation','Fontsize',fsize,'FontName',ftype)
    end
    ylabel('e_n (ms)','Fontsize',fsize,'FontName',ftype);
    
    % Tipo de perturbacion
    if pan==1
        title( 'Small perturbations','fontweight','bold','Fontsize',fsize,'FontName',ftype)
    elseif pan==2
        title( 'Large perturbations','fontweight','bold','Fontsize',fsize,'FontName',ftype)
    end

    
    % Eje X
    xlim([-2 14])
    set(gca,'XTick',[0:5:14])
    set(gca,'FontSize',fsize);
    
    % Eje Y
    ylim([y_lim_inf y_lim_sup])
    set(gca,'YTick',[y_lim_inf_tick:y_int:y_lim_sup_tick])
    set(gca,'FontSize',fsize);
   
    
%     %Asteriscos
%         for bip=1:length(adj_p)
%             if adj_p(1,bip)<alpha
%                 plot(bip,y_lim_inf+5,'*r','MarkerSize',msize)
%             end
%         end
    
    % Legend
    if pan==1
        lh = legend([ptm1.mainLine ptm2.mainLine ptm3.mainLine],'-S-T','(-S)+(-T)','DIFF');
        set(lh,'Location','Northeast');
        legend('boxoff');
    elseif pan==3
        lh = legend([ptm1.mainLine ptm2.mainLine ptm3.mainLine],'-S+T','(-S)+(+T)','DIFF');
        set(lh,'Location','Southeast');
        legend('boxoff');
    elseif pan==5
        lh = legend([ptm1.mainLine ptm2.mainLine ptm3.mainLine],'+S-T','(+S)+(-T)','DIFF');
        set(lh,'Location','Northeast');
        legend('boxoff');
    elseif pan==7
        lh = legend([ptm1.mainLine ptm2.mainLine ptm3.mainLine],'+S+T','(+S)+(+T)','DIFF');
        set(lh,'Location','Southeast');
        legend('boxoff');
    end
    
    box off;
end

% Letras
l_pos=[0.06 0.943 0 0; 0.5 0.943 0 0;...
       0.06 0.723 0 0; 0.5 0.723 0 0;...
       0.06 0.505 0 0; 0.5 0.505 0 0;...
       0.06 0.285 0 0; 0.5 0.285 0 0];
for i=1:length(letras)
    l=annotation('textbox', l_pos(i,:), 'string',letras(i));
    l_g=get(l);
     set(l,'FontSize',fsize_letra,'FontWeight','bold');
end

%% Uncomment to save
% set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
% set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);
% print('-dpng','-r600','figure_additivity.png');

%% Figure 1.
% Resynchronization time series for every condition averaged across subjects.

period_pre = 500;
perturb_sizes = [15 30 45 50];
period_post_delta = [-1 0 1 -1 0 1 -1 0 1];
cond_names = {'-S-T','-S','-S+T','-T','iso','+T','+S-T','+S','+S+T'};
cond_types = {'T','S','opuestas','analogas'};
nbr_pert_sizes = length(perturb_sizes);
nbr_conds = length(cond_names);
nbr_cond_types = length(cond_types);
ceroM = 57;

aux_1.subject(1) = [];	% subject=1 tiene una pert_size especial, queda afuera

nbr_subjects_1 = length(aux_1.subject);	% pert_size 15,30,45
nbr_subjects_2 = length(aux_2.subject);	% pert_size 50
nbr_subjects = nbr_subjects_1 + nbr_subjects_2;
stims_series = aux_1.subject(1).exp(1).asyn_posbs(2,:);
nbr_stims_series = length(stims_series);
perturb_bip_series = find(stims_series==0);
isi_series = period_pre*ones(nbr_pert_sizes,nbr_conds,nbr_stims_series);
isi_series_cell = cell(nbr_pert_sizes,nbr_conds);
delta_isi_series = zeros(nbr_pert_sizes,nbr_conds,nbr_stims_series);
isi_series(:,:,1:perturb_bip_series-1) = period_pre;
delta_isi_series(:,:,1:perturb_bip_series-1) = 0;

% initialize structure
data_series_subj = [];
for pert = 1:nbr_pert_sizes
	for cond = 1:nbr_conds
		data_series_subj(pert,cond).subject = [];
		data_series_subj(pert,cond).list_of_subjects = [];
		delta_isi_series(pert,cond,perturb_bip_series) = perturb_sizes(pert)*period_post_delta(cond);
 		isi_series(pert,cond,perturb_bip_series:end) = period_pre + delta_isi_series(pert,cond,perturb_bip_series);
		isi_series_cell{pert,cond} = squeeze(isi_series(pert,cond,:));
		for subj = 1:nbr_subjects
			data_series_subj(pert,cond).subject(subj).asyns = [];
		end
	end
end
% isi_series = period_pre + delta_isi_series;
list_of_subjects = cell(nbr_pert_sizes,nbr_conds);

% load data: pert_size=[15 30 45]
for subj = 1:nbr_subjects_1
	if aux_1.subject(subj).out==0		% subject no outlier
		temp_sizes = [aux_1.subject(subj).exp(:).temp_size];
		pert_size = max(temp_sizes);
		if abs(pert_size)==45
			pert = 3;
		elseif abs(pert_size)==30
			pert = 2;
		elseif abs(pert_size)==15	% pert_size=9 es para el subject 1
			pert = 1;
		end
		nbr_trials = length(aux_1.subject(subj).exp);
		for trial = 1:nbr_trials
			if aux_1.subject(subj).exp(trial).out==0		% trial no outlier
				cond = aux_1.subject(subj).exp(trial).condition;
				asyns = aux_1.subject(subj).exp(trial).asyn_posbs(1,:);
				if isempty(data_series_subj(pert,cond).subject(subj).asyns)
					data_series_subj(pert,cond).subject(subj).asyns = asyns;
				else
					data_series_subj(pert,cond).subject(subj).asyns = cat(1,data_series_subj(pert,cond).subject(subj).asyns,asyns);
				end
				list_of_subjects{pert,cond} = [list_of_subjects{pert,cond} subj];
			end
		end
	end
end
% load data: pert_size=[50]
for subj = 1:nbr_subjects_2
	subj_tot = nbr_subjects_1 + subj;
	if aux_2.subject(subj).out==0		% subject no outlier
		pert = 4;
		nbr_trials = length(aux_2.subject(subj).exp);
		for trial = 1:nbr_trials
			if aux_2.subject(subj).exp(trial).out==0		% trial no outlier
				cond = aux_2.subject(subj).exp(trial).condition;
				asyns = aux_2.subject(subj).exp(trial).asyn_posbs(1,:);
				if isempty(data_series_subj(pert,cond).subject(subj_tot).asyns)
					data_series_subj(pert,cond).subject(subj_tot).asyns = asyns;
				else
					data_series_subj(pert,cond).subject(subj_tot).asyns = cat(1,data_series_subj(pert,cond).subject(subj_tot).asyns,asyns);
				end
				list_of_subjects{pert,cond} = [list_of_subjects{pert,cond} subj_tot];
			end
		end
	end
end

% wrapping up
for pert = 1:nbr_pert_sizes
	for cond = 1:nbr_conds
		list_of_subjects{pert,cond} = unique(list_of_subjects{pert,cond});
	end
end
data_series_subj_ave = nan(nbr_pert_sizes,nbr_conds,nbr_subjects,nbr_stims_series);
data_series_alltrials = cell(nbr_pert_sizes,nbr_conds);
for pert = 1:nbr_pert_sizes
	for cond = 1:nbr_conds
		aux = [];
		for subj = 1:nbr_subjects
			if find(list_of_subjects{pert,cond}==subj)
				data_series_subj_ave(pert,cond,subj,:) = mean(data_series_subj(pert,cond).subject(subj).asyns,1);
				aux = cat(1,aux,data_series_subj(pert,cond).subject(subj).asyns);
			end
		end
		data_series_alltrials{pert,cond} = aux;
	end
end
data_series_subj_ave_pertcond_ave = squeeze(nanmean(data_series_subj_ave,3));


% embedding trajectories
data_trajec_2D = repmat(data_series_subj_ave_pertcond_ave,[1 1 1 2]);
data_trajec_alltrials_2D = cellfun(@(x) repmat(x,[1 1 2]),data_series_alltrials,'UniformOutput',false);
data_trajec_subj_ave_2D = repmat(data_series_subj_ave,[1 1 1 1 2]);
% suma y diferencia con el paso anterior: ((e_n + e_{n-1})/2; (e_n - e_{n-1})/2)
aux_1 = 0.5*(data_trajec_2D(:,:,:,2) + circshift(data_trajec_2D(:,:,:,2),1,3));
aux_2 = 0.5*(data_trajec_2D(:,:,:,2) - circshift(data_trajec_2D(:,:,:,2),1,3));
data_trajec_2D(:,:,:,1) = aux_1;
data_trajec_2D(:,:,:,2) = aux_2;
aux_1 = 0.5*(data_trajec_subj_ave_2D(:,:,:,:,2) + circshift(data_trajec_subj_ave_2D(:,:,:,:,2),1,4));
aux_2 = 0.5*(data_trajec_subj_ave_2D(:,:,:,:,2) - circshift(data_trajec_subj_ave_2D(:,:,:,:,2),1,4));
data_trajec_subj_ave_2D(:,:,:,:,1) = aux_1;
data_trajec_subj_ave_2D(:,:,:,:,2) = aux_2;
data_trajec_alltrials_2D = cellfun(@(x) cat(3,0.5*(x(:,:,1)+circshift(x(:,:,1),1,2)),0.5*(x(:,:,2)-circshift(x(:,:,2),1,2))),data_trajec_alltrials_2D,'UniformOutput',false);
xlabel_embed = '(e_n + e_{n-1})/2 (ms)';
ylabel_embed = '(e_n - e_{n-1})/2 (ms)';
xlabel_embed_norm = '(e_n + e_{n-1})/2 (nondim)';
ylabel_embed_norm = '(e_n - e_{n-1})/2 (nondim)';
data_trajec_2D(:,:,1,:) = [];
data_trajec_alltrials_2D = cellfun(@(x) x(:,2:end,:),data_trajec_alltrials_2D,'UniformOutput',false);
data_trajec_subj_ave_2D(:,:,:,1,:) = [];
stims_trajec = stims_series(2:end);
nbr_stims_trajec = length(stims_trajec);
perturb_bip_trajec = find(stims_trajec==0);

%% figure: step change time series

cond = 4;
pert_size = 4;
subj = 1;
trial = 7;

subj_nbr = list_of_subjects{pert_size,cond}(subj);
datos = data_series_subj(pert_size,cond).subject(subj_nbr).asyns(trial,:);
isi = 500*ones(size(datos));
deltaT = -50;
isi(perturb_bip_series:end) = isi(perturb_bip_series:end) + deltaT;
time_axis = [1:nbr_stims_series] - perturb_bip_series;

fsize = 7.5;
fsize_legend = 6;
msize = 8;
lwidth = 1;
color_map = {'c','m','b','k'};
linestyle_map = {'--','-'};
marker_map = {'d','.'};
markersize_map = [10 50];
linewidth_map = [5 1];
fig_size_cm = [10 5];


new_cmap = ametrine(9);

figure(1);
clf(1);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);


subplot(3,1,[1 2]);
plot(time_axis([1 end]),[0 0],'k-');
hold on;
plot([0 0],[-25 75],'k-');
% plot(time_axis,datos,'.-','color',new_cmap(1,:),'linewidth',lwidth,'markersize',msize);
plot(time_axis,datos,'k.-','linewidth',lwidth,'markersize',msize);
set(gca,'xtick',[-10:5:15],'ytick',[-25:25:75],'xticklabel',[]);
set(gca,'fontsize',fsize);
xlim([-10 14]);
ylim([-25 75]);
ylabel('e_n (ms)');
posf = ds2nfu(gca,[-1.8 42.5 1.3 0]);
pta = annotation('textarrow',[posf(1) posf(1)+posf(3)],[posf(2) posf(2)+posf(4)],'string',{'forced';'error'});
set(pta,'headstyle','plain','headwidth',3,'headlength',5,'linewidth',1.3,'fontsize',fsize);
ptt = text(-13,75,'a)');
set(ptt,'fontsize',1.25*fsize,'fontweight','bold');

subplot(3,1,3);
plot([0 0],[440 510],'k-');
hold on;
% plot(time_axis,isi,'.-','color',new_cmap(1,:),'linewidth',lwidth,'markersize',msize);
plot(time_axis,isi,'k.-','linewidth',lwidth,'markersize',msize);
set(gca,'xtick',[-10:5:15],'ytick',[450:50:500]);
set(gca,'fontsize',fsize);
xlim([-10 14]);
ylim([440 510]);
xlabel('Beep n relative to perturbation');
ylabel({'Stimulus';'period (ms)'});


%%% Uncomment to save figure
% print_filename = 'figure_stepchange_timeseries';
% print('-dpdf',[print_filename '.pdf']);

%% Figure 6.
% Experimental reconstruction of the phase space by means of an embedding.

data_series = data_series_subj_ave_pertcond_ave;
selected_conds = [4 6];

data_series_T = data_series(:,selected_conds,:);
data_series_50_noS = permute(cat(1,data_series(4,[4 7 1],:),data_series(4,[6 3 9],:)),[2 1 3]);

data_trajec_2D_T = data_trajec_2D(:,[4 6],:,:);
data_trajec_2D_M = data_trajec_2D(:,[2 8],:,:);
data_trajec_2D_opuestas = data_trajec_2D(:,[7 3],:,:);
data_trajec_2D_analogas = data_trajec_2D(:,[1 9],:,:);

data_trajec_2D_50 = [];
data_trajec_2D_50(1,:,:,:) = data_trajec_2D_T(4,:,:,:);
data_trajec_2D_50(2,:,:,:) = data_trajec_2D_M(4,:,:,:);
data_trajec_2D_50(3,:,:,:) = data_trajec_2D_opuestas(4,:,:,:);
data_trajec_2D_50(4,:,:,:) = data_trajec_2D_analogas(4,:,:,:);

data_trajec_2D_50_noS = data_trajec_2D_50([1 3 4],:,:,:);

lims_x_T_series = [-2 14];
lims_y_T_series = [-60 60];
lims_x_T_embed = [-60 60];
lims_y_T_embed = [-20 25];
lims_x_50_series = [-2 14];
lims_y_50_series = [-60 60];
lims_x_50_embed = [-50 50];
lims_y_50_embed = [-30 30];
init_bip = perturb_bip_trajec;
end_bip = length(stims_trajec);

fsize = 7;
fsize_legend = 6;
msize = 8;
lwidth = 1;
color_map = {'c','m','b','k'};
linestyle_map = {'--','-'};
marker_map = {'d','.'};
markersize_map = [10 50];
linewidth_map = [5 1];
% fig_size_cm = [20 8];
fig_size_cm = [14 5];%[12 5];

new_cmap = ametrine(4);

figure(6);
clf(6);

set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

%%% Panel A
subplot(1,2,1);
ptp = nan(nbr_pert_sizes,2);
set(gca,'fontsize',fsize);
plot(lims_x_T_series,[0 0],'k');
hold on;
plot([0 0],lims_y_T_series,'k');
for pert_size = 1:nbr_pert_sizes
	for cond = 1:2
		ptp(pert_size,cond) = plot(stims_series,squeeze(data_series_T(pert_size,cond,:)),'color',new_cmap(pert_size,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
	end
end
xlim(lims_x_T_series);
ylim(lims_y_T_series);
set(gca,'xtick',[0:5:15],'ytick',[lims_y_T_series(1):20:lims_y_T_series(end)]);
xlabel('Beep n relative to perturbation');
ylabel('e_n (ms)');
% data_series_T_legend = {'-T 15','-T 30','-T 45','-T 50','+T 15','+T 30','+T 45','+T 50'};
data_series_T_legend = {'15','30','45','50'};
% legend(reshape(ptp,1,[]),data_series_T_legend,'location','southeast','fontsize',fsize_legend);
legend(ptp(:,2),data_series_T_legend,'location','southeast','fontsize',fsize_legend);
legend boxoff;
ptt = text(7,50,{'continuous: +T';'dashed: -T'});
set(ptt,'fontsize',fsize_legend);
ptt = text(-5.7,60,'a)');
set(ptt,'fontsize',fsize,'fontweight','bold');
box off;

%%% Panel B
subplot(1,2,2);
ptp = nan(nbr_pert_sizes,2);
set(gca,'fontsize',fsize);
plot([0 0],lims_y_T_embed,'k');
hold on;
plot(lims_x_T_embed,[0 0],'k');
for pert_size = 1:nbr_pert_sizes
	for cond = 1:2
		ptp(pert_size,cond) = plot(squeeze(data_trajec_2D_T(pert_size,cond,init_bip:end_bip,1)),squeeze(data_trajec_2D_T(pert_size,cond,init_bip:end_bip,2)),'color',new_cmap(pert_size,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
		plot(squeeze(data_trajec_2D_T(pert_size,cond,init_bip,1)),squeeze(data_trajec_2D_T(pert_size,cond,init_bip,2)),'.','color',new_cmap(pert_size,:),'markersize',msize);
	end
end
xlim(lims_x_T_embed);
ylim(lims_y_T_embed);
set(gca,'xtick',[lims_x_T_embed(1):20:lims_x_T_embed(end)],'ytick',[lims_y_T_embed(1):10:lims_y_T_embed(end)]);
xlabel(xlabel_embed);
ylabel(ylabel_embed);
ptt = text(-90,25,'b)');
set(ptt,'fontsize',fsize,'fontweight','bold');

%%% Uncomment to save figure
figure(6);
% print_filename = 'figure_embed_T';
% print('-dpdf',[print_filename '.pdf']);
% print('-dpng','-r600',[print_filename '.png']);
% print('-dtiff','-r600',[print_filename '.tiff']);
% print('-depsc','-r600',[print_filename '.eps']);

%% Figure 7
% Response to combined perturbations of ±50 ms and experimental reconstruction
% of phase space. (a) Averaged time series of traditional simple temporal 
% perturbations (+T and ?T), combined opposite perturbations (?S+T and +S?T),
% and combined analogous perturbations (+S+T and ?S?T). 

fsize_legend = 5;
new_cmap = ametrine(3);

figure(7);
clf(7);

set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

subplot(1,2,1);
ptp = nan(nbr_cond_types-1,2);
set(gca,'fontsize',fsize);
plot(lims_x_T_series,[0 0],'k');
hold on;
plot([0 0],lims_y_50_series,'k');
for cond_type = 1:nbr_cond_types-1
	for cond = 1:2
		ptp(cond_type,cond) = plot(stims_series,squeeze(data_series_50_noS(cond_type,cond,:)),'color',new_cmap(cond_type,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
	end
end
xlim(lims_x_50_series);
ylim(lims_y_50_series);
set(gca,'xtick',[0:5:15],'ytick',[lims_y_50_series(1):20:lims_y_50_series(end)]);
xlabel('Beep n relative to perturbation');
ylabel('e_n (ms)');
data_series_50_noS_legend = {'-T','-T+S','-T-S','+T','+T-S','+T+S'};
legend(reshape(ptp,1,[]),data_series_50_noS_legend,'location','southeast','fontsize',fsize_legend);
legend boxoff;
ptt = text(-5.7,60,'a)');
set(ptt,'fontsize',fsize,'fontweight','bold');
box off;

%%% Panel A
subplot(1,2,2);
ptp = nan(nbr_cond_types-1,2);
set(gca,'fontsize',fsize);
plot([0 0],lims_y_50_embed,'k');
hold on;
plot(lims_x_50_embed,[0 0],'k');
for cond_type = 1:nbr_cond_types-1
	for cond = 1:2
		ptp(cond_type,cond) = plot(squeeze(data_trajec_2D_50_noS(cond_type,cond,init_bip:end_bip,1)),squeeze(data_trajec_2D_50_noS(cond_type,cond,init_bip:end_bip,2)),'color',new_cmap(cond_type,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
		plot(squeeze(data_trajec_2D_50_noS(cond_type,cond,init_bip,1)),squeeze(data_trajec_2D_50_noS(cond_type,cond,init_bip,2)),'.','color',new_cmap(cond_type,:),'markersize',msize);
	end
end
xlim(lims_x_50_embed);
ylim(lims_y_50_embed);
set(gca,'xtick',[lims_x_50_embed(1):20:lims_x_50_embed(end)]-10,'ytick',[lims_y_50_embed(1):10:lims_y_50_embed(end)]);
xlabel(xlabel_embed);
ylabel(ylabel_embed);
ptt = text(-75,30,'b)');
set(ptt,'fontsize',fsize,'fontweight','bold');

%%% Uncomment to save figure
figure(7);
% print_filename = 'figure_embed_ST';
% print('-dpdf',[print_filename '.pdf']);
% print('-dpng','-r600',[print_filename '.png']);
% print('-dtiff','-r600',[print_filename '.tiff']);

%% Figure 8
% Access to previously unexplored system states.

theta = 0:0.01:2*pi;

selected_perts = [3 4];%[1 2 3 4];
norm_pert = selected_perts(end);

data_trajec_subj_ave_2D_T = data_trajec_subj_ave_2D(:,[4 6],:,:,:);
data_trajec_subj_ave_2D_M = data_trajec_subj_ave_2D(:,[2 8],:,:,:);
data_trajec_subj_ave_2D_opuestas = data_trajec_subj_ave_2D(:,[7 3],:,:,:);
data_trajec_subj_ave_2D_analogas = data_trajec_subj_ave_2D(:,[1 9],:,:,:);

% normalizar todos los tamaÃ±os de perturbacion de una misma condition
max_T_norm = repmat(max(abs(data_trajec_2D_T(selected_perts,:,:,:)),[],3),[1 1 nbr_stims_trajec 1]);
max_M_norm = repmat(max(abs(data_trajec_2D_M(selected_perts,:,:,:)),[],3),[1 1 nbr_stims_trajec 1]);
max_op_norm = repmat(max(abs(data_trajec_2D_opuestas(selected_perts,:,:,:)),[],3),[1 1 nbr_stims_trajec 1]);
max_an_norm = repmat(max(abs(data_trajec_2D_analogas(selected_perts,:,:,:)),[],3),[1 1 nbr_stims_trajec 1]);

% con esto llevo todas al tamaÃ±o que corresponde a la pert que elegi como normalizacion
max_T_up = repmat(max_T_norm(find(selected_perts==norm_pert),:,:,:),[length(selected_perts),1,1,1]);
max_M_up = repmat(max_M_norm(find(selected_perts==norm_pert),:,:,:),[length(selected_perts),1,1,1]);
max_op_up = repmat(max_op_norm(find(selected_perts==norm_pert),:,:,:),[length(selected_perts),1,1,1]);
max_an_up = repmat(max_an_norm(find(selected_perts==norm_pert),:,:,:),[length(selected_perts),1,1,1]);

max_M_up = max_M_up./max_T_up;
max_op_up = max_op_up./max_T_up;
max_an_up = max_an_up./max_T_up;
max_T_up = max_T_up./max_T_up;

data_trajec_2D_T_scaleup = data_trajec_2D_T(selected_perts,:,:,:).*max_T_up./max_T_norm;
data_trajec_2D_T_scaleup_ave = squeeze(mean(data_trajec_2D_T_scaleup,1));
data_trajec_2D_M_scaleup = data_trajec_2D_M(selected_perts,:,:,:).*max_M_up./max_M_norm;
data_trajec_2D_M_scaleup_ave = squeeze(mean(data_trajec_2D_M_scaleup,1));
data_trajec_2D_op_scaleup = data_trajec_2D_opuestas(selected_perts,:,:,:).*max_op_up./max_op_norm;
data_trajec_2D_op_scaleup_ave = squeeze(mean(data_trajec_2D_op_scaleup,1));
data_trajec_2D_an_scaleup = data_trajec_2D_analogas(selected_perts,:,:,:).*max_an_up./max_an_norm;
data_trajec_2D_an_scaleup_ave = squeeze(mean(data_trajec_2D_an_scaleup,1));

max_T_up_subj = permute(repmat(max_T_up,[1 1 1 1 nbr_subjects]),[1 2 5 3 4]);
max_T_norm_subj = permute(repmat(max_T_norm,[1 1 1 1 nbr_subjects]),[1 2 5 3 4]);
data_trajec_subj_ave_2D_T_scaleup = data_trajec_subj_ave_2D_T(selected_perts,:,:,:,:).*max_T_up_subj./max_T_norm_subj;
% concateno todos los subjects de todas las pert_sizes
data_trajec_subj_ave_2D_T_scaleup_all = squeeze(data_trajec_subj_ave_2D_T_scaleup(1,:,:,:,:));
for pert = 2:length(selected_perts)
	data_trajec_subj_ave_2D_T_scaleup_all = cat(2,data_trajec_subj_ave_2D_T_scaleup_all,squeeze(data_trajec_subj_ave_2D_T_scaleup(pert,:,:,:,:)));
end
data_trajec_subj_ave_2D_T_scaleup_all_init = data_trajec_subj_ave_2D_T_scaleup_all(:,:,perturb_bip_trajec,:);
data_trajec_subj_ave_2D_T_scaleup_all_init_center = squeeze(nanmean(data_trajec_subj_ave_2D_T_scaleup_all_init,2));
data_trajec_subj_ave_2D_T_scaleup_all_init_std = squeeze(nanstd(data_trajec_subj_ave_2D_T_scaleup_all_init,0,2));
data_trajec_subj_ave_2D_T_scaleup_all_init_radius = data_trajec_subj_ave_2D_T_scaleup_all_init_std;
data_trajec_subj_ave_2D_T_scaleup_all_init_xcircle = nan(2,length(theta));
data_trajec_subj_ave_2D_T_scaleup_all_init_ycircle = nan(2,length(theta));
aux_x = arrayfun(@(x,y) x*cos(theta) + y,data_trajec_subj_ave_2D_T_scaleup_all_init_radius(:,1),data_trajec_subj_ave_2D_T_scaleup_all_init_center(:,1),'UniformOutput',false);
aux_y = arrayfun(@(x,y) x*sin(theta) + y,data_trajec_subj_ave_2D_T_scaleup_all_init_radius(:,2),data_trajec_subj_ave_2D_T_scaleup_all_init_center(:,2),'UniformOutput',false);
for pert_sign = 1:2
	data_trajec_subj_ave_2D_T_scaleup_all_init_xcircle(pert_sign,:) = aux_x{pert_sign};
	data_trajec_subj_ave_2D_T_scaleup_all_init_ycircle(pert_sign,:) = aux_y{pert_sign};
end

max_M_up_subj = permute(repmat(max_M_up,[1 1 1 1 nbr_subjects]),[1 2 5 3 4]);
max_M_norm_subj = permute(repmat(max_M_norm,[1 1 1 1 nbr_subjects]),[1 2 5 3 4]);
data_trajec_subj_ave_2D_M_scaleup = data_trajec_subj_ave_2D_M(selected_perts,:,:,:,:).*max_M_up_subj./max_M_norm_subj;
% concateno todos los subjects de todas las pert_sizes
data_trajec_subj_ave_2D_M_scaleup_all = squeeze(data_trajec_subj_ave_2D_M_scaleup(1,:,:,:,:));
for pert = 2:length(selected_perts)
	data_trajec_subj_ave_2D_M_scaleup_all = cat(2,data_trajec_subj_ave_2D_M_scaleup_all,squeeze(data_trajec_subj_ave_2D_M_scaleup(pert,:,:,:,:)));
end
data_trajec_subj_ave_2D_M_scaleup_all_init = data_trajec_subj_ave_2D_M_scaleup_all(:,:,perturb_bip_trajec,:);
data_trajec_subj_ave_2D_M_scaleup_all_init_center = squeeze(nanmean(data_trajec_subj_ave_2D_M_scaleup_all_init,2));
data_trajec_subj_ave_2D_M_scaleup_all_init_std = squeeze(nanstd(data_trajec_subj_ave_2D_M_scaleup_all_init,0,2));
data_trajec_subj_ave_2D_M_scaleup_all_init_err = data_trajec_subj_ave_2D_M_scaleup_all_init_std/sqrt((length(selected_perts)*nbr_subjects));
data_trajec_subj_ave_2D_M_scaleup_all_init_radius = data_trajec_subj_ave_2D_M_scaleup_all_init_std;
data_trajec_subj_ave_2D_M_scaleup_all_init_xcircle = nan(2,length(theta));
data_trajec_subj_ave_2D_M_scaleup_all_init_ycircle = nan(2,length(theta));
aux_x = arrayfun(@(x,y) x*cos(theta) + y,data_trajec_subj_ave_2D_M_scaleup_all_init_radius(:,1),data_trajec_subj_ave_2D_M_scaleup_all_init_center(:,1),'UniformOutput',false);
aux_y = arrayfun(@(x,y) x*sin(theta) + y,data_trajec_subj_ave_2D_M_scaleup_all_init_radius(:,2),data_trajec_subj_ave_2D_M_scaleup_all_init_center(:,2),'UniformOutput',false);
for pert_sign = 1:2
	data_trajec_subj_ave_2D_M_scaleup_all_init_xcircle(pert_sign,:) = aux_x{pert_sign};
	data_trajec_subj_ave_2D_M_scaleup_all_init_ycircle(pert_sign,:) = aux_y{pert_sign};
end

max_op_up_subj = permute(repmat(max_op_up,[1 1 1 1 nbr_subjects]),[1 2 5 3 4]);
max_op_norm_subj = permute(repmat(max_op_norm,[1 1 1 1 nbr_subjects]),[1 2 5 3 4]);
data_trajec_subj_ave_2D_opuestas_scaleup = data_trajec_subj_ave_2D_opuestas(selected_perts,:,:,:,:).*max_op_up_subj./max_op_norm_subj;
% concateno todos los subjects de todas las pert_sizes
data_trajec_subj_ave_2D_opuestas_scaleup_all = squeeze(data_trajec_subj_ave_2D_opuestas_scaleup(1,:,:,:,:));
for pert = 2:length(selected_perts)
	data_trajec_subj_ave_2D_opuestas_scaleup_all = cat(2,data_trajec_subj_ave_2D_opuestas_scaleup_all,squeeze(data_trajec_subj_ave_2D_opuestas_scaleup(pert,:,:,:,:)));
end
data_trajec_subj_ave_2D_opuestas_scaleup_all_init = data_trajec_subj_ave_2D_opuestas_scaleup_all(:,:,perturb_bip_trajec,:);
data_trajec_subj_ave_2D_opuestas_scaleup_all_init_center = squeeze(nanmean(data_trajec_subj_ave_2D_opuestas_scaleup_all_init,2));
data_trajec_subj_ave_2D_opuestas_scaleup_all_init_std = squeeze(nanstd(data_trajec_subj_ave_2D_opuestas_scaleup_all_init,0,2));
data_trajec_subj_ave_2D_opuestas_scaleup_all_init_err = data_trajec_subj_ave_2D_opuestas_scaleup_all_init_std/sqrt((length(selected_perts)*nbr_subjects));
data_trajec_subj_ave_2D_opuestas_scaleup_all_init_radius = data_trajec_subj_ave_2D_opuestas_scaleup_all_init_std;
data_trajec_subj_ave_2D_opuestas_scaleup_all_init_xcircle = nan(2,length(theta));
data_trajec_subj_ave_2D_opuestas_scaleup_all_init_ycircle = nan(2,length(theta));
aux_x = arrayfun(@(x,y) x*cos(theta) + y,data_trajec_subj_ave_2D_opuestas_scaleup_all_init_radius(:,1),data_trajec_subj_ave_2D_opuestas_scaleup_all_init_center(:,1),'UniformOutput',false);
aux_y = arrayfun(@(x,y) x*sin(theta) + y,data_trajec_subj_ave_2D_opuestas_scaleup_all_init_radius(:,2),data_trajec_subj_ave_2D_opuestas_scaleup_all_init_center(:,2),'UniformOutput',false);
for pert_sign = 1:2
	data_trajec_subj_ave_2D_opuestas_scaleup_all_init_xcircle(pert_sign,:) = aux_x{pert_sign};
	data_trajec_subj_ave_2D_opuestas_scaleup_all_init_ycircle(pert_sign,:) = aux_y{pert_sign};
end

max_an_up_subj = permute(repmat(max_an_up,[1 1 1 1 nbr_subjects]),[1 2 5 3 4]);
max_an_norm_subj = permute(repmat(max_an_norm,[1 1 1 1 nbr_subjects]),[1 2 5 3 4]);
data_trajec_subj_ave_2D_analogas_scaleup = data_trajec_subj_ave_2D_analogas(selected_perts,:,:,:,:).*max_an_up_subj./max_an_norm_subj;
% concateno todos los subjects de todas las pert_sizes
data_trajec_subj_ave_2D_analogas_scaleup_all = squeeze(data_trajec_subj_ave_2D_analogas_scaleup(1,:,:,:,:));
for pert = 2:length(selected_perts)
	data_trajec_subj_ave_2D_analogas_scaleup_all = cat(2,data_trajec_subj_ave_2D_analogas_scaleup_all,squeeze(data_trajec_subj_ave_2D_analogas_scaleup(pert,:,:,:,:)));
end
data_trajec_subj_ave_2D_analogas_scaleup_all_init = data_trajec_subj_ave_2D_analogas_scaleup_all(:,:,perturb_bip_trajec,:);
data_trajec_subj_ave_2D_analogas_scaleup_all_init_center = squeeze(nanmean(data_trajec_subj_ave_2D_analogas_scaleup_all_init,2));
data_trajec_subj_ave_2D_analogas_scaleup_all_init_std = squeeze(nanstd(data_trajec_subj_ave_2D_analogas_scaleup_all_init,0,2));
data_trajec_subj_ave_2D_analogas_scaleup_all_init_err = data_trajec_subj_ave_2D_analogas_scaleup_all_init_std/sqrt((length(selected_perts)*nbr_subjects));
data_trajec_subj_ave_2D_analogas_scaleup_all_init_radius = data_trajec_subj_ave_2D_analogas_scaleup_all_init_std;
data_trajec_subj_ave_2D_analogas_scaleup_all_init_xcircle = nan(2,length(theta));
data_trajec_subj_ave_2D_analogas_scaleup_all_init_ycircle = nan(2,length(theta));
aux_x = arrayfun(@(x,y) x*cos(theta) + y,data_trajec_subj_ave_2D_analogas_scaleup_all_init_radius(:,1),data_trajec_subj_ave_2D_analogas_scaleup_all_init_center(:,1),'UniformOutput',false);
aux_y = arrayfun(@(x,y) x*sin(theta) + y,data_trajec_subj_ave_2D_analogas_scaleup_all_init_radius(:,2),data_trajec_subj_ave_2D_analogas_scaleup_all_init_center(:,2),'UniformOutput',false);
for pert_sign = 1:2
	data_trajec_subj_ave_2D_analogas_scaleup_all_init_xcircle(pert_sign,:) = aux_x{pert_sign};
	data_trajec_subj_ave_2D_analogas_scaleup_all_init_ycircle(pert_sign,:) = aux_y{pert_sign};
end

lightcyan = [0.75 1 1];
lightblue = [0.75 0.85 1];
lightblack = [0.75 0.75 0.75];
lightmagenta = [1 0.65 1];

fsize = 14;
fsize_legend = 11;
msize = 15;
lwidth = 2;
new_cmap = ametrine(4);
new_cmap_area = brighten(ametrine(4),0.8);
new_cmap_area(end,3) = 0.75;
linestyle_map = {'--','-'};
marker_map = {'d','.'};
markersize_map = [10 50];
linewidth_map = [5 1];
fig_size_cm = [30 12];
lims_x = [-2 1.5];
lims_y = [-2 2];

figure(8);
clf(8);

set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);


%%% Panel A
subplot(1,2,2);
ptp = nan(nbr_pert_sizes,2);
set(gca,'fontsize',fsize);
plot([0 0],lims_y,'k-','linewidth',1);
hold on;
plot(lims_x,[0 0],'k-','linewidth',1);

cond = 1;
ptf = fill(squeeze(data_trajec_subj_ave_2D_T_scaleup_all_init_xcircle(cond,:)),squeeze(data_trajec_subj_ave_2D_T_scaleup_all_init_ycircle(cond,:)),new_cmap_area(1,:));
set(ptf,'facecolor',new_cmap_area(1,:),'edgecolor',new_cmap_area(1,:));
cond = 2;
ptf = fill(squeeze(data_trajec_subj_ave_2D_T_scaleup_all_init_xcircle(cond,:)),squeeze(data_trajec_subj_ave_2D_T_scaleup_all_init_ycircle(cond,:)),new_cmap_area(1,:));
set(ptf,'facecolor',new_cmap_area(1,:),'edgecolor',new_cmap_area(1,:));
cond = 2;
ptf = fill(squeeze(data_trajec_subj_ave_2D_opuestas_scaleup_all_init_xcircle(cond,:)),squeeze(data_trajec_subj_ave_2D_opuestas_scaleup_all_init_ycircle(cond,:)),new_cmap_area(3,:));
set(ptf,'facecolor',new_cmap_area(3,:),'edgecolor',new_cmap_area(3,:));
cond = 1;
ptf = fill(squeeze(data_trajec_subj_ave_2D_opuestas_scaleup_all_init_xcircle(cond,:)),squeeze(data_trajec_subj_ave_2D_opuestas_scaleup_all_init_ycircle(cond,:)),new_cmap_area(3,:));
set(ptf,'facecolor',new_cmap_area(3,:),'edgecolor',new_cmap_area(3,:));
cond = 1;
ptf = fill(squeeze(data_trajec_subj_ave_2D_analogas_scaleup_all_init_xcircle(cond,:)),squeeze(data_trajec_subj_ave_2D_analogas_scaleup_all_init_ycircle(cond,:)),new_cmap_area(4,:));
set(ptf,'facecolor',new_cmap_area(4,:),'edgecolor',new_cmap_area(4,:));
cond = 2;
ptf = fill(squeeze(data_trajec_subj_ave_2D_analogas_scaleup_all_init_xcircle(cond,:)),squeeze(data_trajec_subj_ave_2D_analogas_scaleup_all_init_ycircle(cond,:)),new_cmap_area(4,:));
set(ptf,'facecolor',new_cmap_area(4,:),'edgecolor',new_cmap_area(4,:));
cond = 1;
ptf = fill(squeeze(data_trajec_subj_ave_2D_M_scaleup_all_init_xcircle(cond,:)),squeeze(data_trajec_subj_ave_2D_M_scaleup_all_init_ycircle(cond,:)),new_cmap_area(2,:));
set(ptf,'facecolor',new_cmap_area(2,:),'edgecolor',new_cmap_area(2,:));
cond = 2;
ptf = fill(squeeze(data_trajec_subj_ave_2D_M_scaleup_all_init_xcircle(cond,:)),squeeze(data_trajec_subj_ave_2D_M_scaleup_all_init_ycircle(cond,:)),new_cmap_area(2,:));
set(ptf,'facecolor',new_cmap_area(2,:),'edgecolor',new_cmap_area(2,:));
for cond = 1:2
	ptp_T(cond) = plot(squeeze(data_trajec_2D_T_scaleup_ave(cond,init_bip:end_bip,1)),squeeze(data_trajec_2D_T_scaleup_ave(cond,init_bip:end_bip,2)),'color',new_cmap(1,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
	plot(squeeze(data_trajec_2D_T_scaleup_ave(cond,init_bip,1)),squeeze(data_trajec_2D_T_scaleup_ave(cond,init_bip,2)),'.','color',new_cmap(1,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
	ptp_M(cond) = plot(squeeze(data_trajec_2D_M_scaleup_ave(cond,init_bip:end_bip,1)),squeeze(data_trajec_2D_M_scaleup_ave(cond,init_bip:end_bip,2)),'color',new_cmap(2,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
	plot(squeeze(data_trajec_2D_M_scaleup_ave(cond,init_bip,1)),squeeze(data_trajec_2D_M_scaleup_ave(cond,init_bip,2)),'.','color',new_cmap(2,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
	ptp_op(cond) = plot(squeeze(data_trajec_2D_op_scaleup_ave(cond,init_bip:end_bip,1)),squeeze(data_trajec_2D_op_scaleup_ave(cond,init_bip:end_bip,2)),'color',new_cmap(3,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
	plot(squeeze(data_trajec_2D_op_scaleup_ave(cond,init_bip,1)),squeeze(data_trajec_2D_op_scaleup_ave(cond,init_bip,2)),'.','color',new_cmap(3,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
	ptp_an(cond) = plot(squeeze(data_trajec_2D_an_scaleup_ave(cond,init_bip:end_bip,1)),squeeze(data_trajec_2D_an_scaleup_ave(cond,init_bip:end_bip,2)),'color',new_cmap(4,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
	plot(squeeze(data_trajec_2D_an_scaleup_ave(cond,init_bip,1)),squeeze(data_trajec_2D_an_scaleup_ave(cond,init_bip,2)),'.','color',new_cmap(4,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
end
xlim(lims_x);
ylim(lims_y);
set(gca,'xtick',[lims_x(1):1:lims_x(end)],'ytick',[lims_y(1):1:lims_y(end)]);
xlabel(xlabel_embed_norm);
ylabel(ylabel_embed_norm);
data_legend = {'-T','-S','-T+S','-T-S','+T','+S','+T-S','+T+S'};
legend(reshape([ptp_T; ptp_M; ptp_op; ptp_an],1,[]),data_legend,'location','northwest','fontsize',fsize_legend);
legend boxoff;
ptt = text(-2.6,2,'b)');
set(ptt,'fontsize',fsize,'fontweight','bold');
set(gca,'linewidth',1);

%%% Panel B
subplot(1,2,1);
ptp = nan(nbr_pert_sizes,2);
set(gca,'fontsize',fsize,'linewidth',1);
plot([0 0],lims_y,'k-','linewidth',1);
hold on;
plot(lims_x,[0 0],'k-','linewidth',1);
cond = 1;
ptf = fill(squeeze(data_trajec_subj_ave_2D_T_scaleup_all_init_xcircle(cond,:)),squeeze(data_trajec_subj_ave_2D_T_scaleup_all_init_ycircle(cond,:)),new_cmap_area(1,:));
set(ptf,'facecolor',new_cmap_area(1,:),'edgecolor',new_cmap_area(1,:));
cond = 2;
ptf = fill(squeeze(data_trajec_subj_ave_2D_T_scaleup_all_init_xcircle(cond,:)),squeeze(data_trajec_subj_ave_2D_T_scaleup_all_init_ycircle(cond,:)),new_cmap_area(1,:));
set(ptf,'facecolor',new_cmap_area(1,:),'edgecolor',new_cmap_area(1,:));
for cond = 1:2
	ptp_T(cond) = plot(squeeze(data_trajec_2D_T_scaleup_ave(cond,init_bip:end_bip,1)),squeeze(data_trajec_2D_T_scaleup_ave(cond,init_bip:end_bip,2)),'color',new_cmap(1,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
	plot(squeeze(data_trajec_2D_T_scaleup_ave(cond,init_bip,1)),squeeze(data_trajec_2D_T_scaleup_ave(cond,init_bip,2)),'.','color',new_cmap(1,:),'linestyle',linestyle_map{cond},'linewidth',lwidth,'markersize',msize);
end
xlim(lims_x);
ylim(lims_y);
set(gca,'xtick',[lims_x(1):1:lims_x(end)],'ytick',[lims_y(1):1:lims_y(end)]);
xlabel(xlabel_embed_norm);
ylabel(ylabel_embed_norm);
data_legend = {'-T','+T'};
legend(reshape(ptp_T,1,[]),data_legend,'location','northwest','fontsize',fsize_legend);
legend boxoff;
ptt = text(-2.6,2,'a)');
set(ptt,'fontsize',fsize,'fontweight','bold');
set(gca,'linewidth',1);

%%% Uncomment to save figure
% print_filename = 'figure_areas';
% print('-dpdf',[print_filename '.pdf']);
% print('-dpng','-r600',[print_filename '.png']);