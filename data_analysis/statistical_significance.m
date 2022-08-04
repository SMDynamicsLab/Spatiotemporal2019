%% Statistical significance
% Script to compute the p-values for the null distributions of ASYM time
% series and DIFF time series as published in López, S.L., Laje, R. 
% Spatiotemporal perturbations in paced finger tapping suggest a common 
% mechanism for the processing of time errors. Sci Rep 9, 17814 (2019). 
% https://doi.org/10.1038/s41598-019-54133-x

%% Clean Start
clear all
close all

load('average_between_experiments_02-Aug-2022.mat')

n_perm=1000;
cont=1;
contS=1;
fontsize=16;
bip_comp=5;
clase= input('Small (15+30): 1; Large (45+50): 2? ');
sim = input('Simetria a analizar: 1)Temporal 2)Spatial 3)Opposite 4)Analogous ');
n_perm=input('Permutation number? (default: 1000) ');
if isempty(n_perm)
    n_perm = 1000;
end
alfa=input('alpha? (default: 0.05) ');
if isempty(alfa)
    alfa = 0.05;
end
comparaciones=[4 6;2 8;3 7;1 9];
titulos = {'Temporal','Spatial','Opposite','Analogous'};

if clase==1
    g1=load('subject_simetric-02-Aug-2022.mat');
    g2=load('subject_simetric-02-Aug-2022.mat');
    cond_g1=[13:15];
    cond_g2=[30];
elseif clase==2
    g1=load('subject_simetric-02-Aug-2022.mat');
    g2=load('subject_simetric_xl-02-Aug-2022.mat');
    cond_g1=[45];
    cond_g2=[50];
end
    
%% Generates the DIFF time serie for each subject within each group of interest

for s=1:size(clase_exp(clase).Cond_Mat_Pos,1)    
        c1=comparaciones(sim,1);
        c2=comparaciones(sim,2);
        resta(contS,:)=clase_exp(clase).Cond_Mat_Pos(s,11:end,c1)+clase_exp(clase).Cond_Mat_Pos(s,11:end,c2);
        contS=contS+1;
end

for bip=1:size(resta,2)
    [h_real,p_real,ci_real,stats_real] = ttest(resta(:,bip));
    t_values_real(bip)=stats_real.tstat;
end

t_max_real=t_values_real(find(abs(t_values_real)==max(abs(t_values_real))));

%% Samples from each group

%Group1
for s=1:length(g1.subject)    
    if g1.subject(s).out==0 && ismember(max(g1.subject(s).exp_temp_sizes),cond_g1)
%         disp(s)
        for t=1:length(g1.subject(s).exp)
            c=g1.subject(s).exp(t).condition;
            if g1.subject(s).exp(t).out==0
                if c==c1 
                    Mat(cont,:)=g1.subject(s).exp(t).asyn_posbs(1,11:end);
                    cont=cont+1;
                elseif c==c2
                    Mat(cont,:)=-g1.subject(s).exp(t).asyn_posbs(1,11:end);
                    cont=cont+1;
                end
            end
        end
    end
end

%Group2
for s=1:length(g2.subject)    
    if g2.subject(s).out==0 && ismember(max(g2.subject(s).exp_temp_sizes),cond_g2)
%         disp(s)
        for t=1:length(g2.subject(s).exp)
            c=g2.subject(s).exp(t).condition;
            if g2.subject(s).exp(t).out==0
                if c==c1 
                    Mat(cont,:)=g2.subject(s).exp(t).asyn_posbs(1,11:end);
                    cont=cont+1;
                elseif c==c2
                    Mat(cont,:)=-g2.subject(s).exp(t).asyn_posbs(1,11:end);
                    cont=cont+1;
                end
            end
        end
    end
end

%% Form the subjects

n_subj=size(clase_exp(clase).Cond_Mat_Pos,1);
n_trials=size(Mat,1);
serie_diff_surr=[];

idx=round((n_trials-1).*rand(16,n_subj,n_perm) + 1);
for i=1:16
    for j=1:n_subj
        for k=1:1000
            surr(i,j,k,:)=Mat(idx(i,j,k),:);
        end
    end
end

% idx=round((n_trials-1).*rand(16,8,n_perm) + 1);
% for i=1:16
%     for j=1:8
%         for k=1:1000
%             surr(i,j,k,:)=Mat(idx(i,j,k),:);
%         end
%     end
% end

minust_surr=squeeze(mean(surr(1:8,:,:,:),1));
plust_surr=squeeze(mean(surr(9:16,:,:,:),1));
% serie_diff=minust+plust;
serie_diff_surr=minust_surr-plust_surr;
serie_diff_surr_media=squeeze(mean(serie_diff_surr,1));

% minust_idx=round((n_trials-1).*rand(8,1) + 1);
% plust_idx=round((n_trials-1).*rand(8,1) + 1);

%% SUBTRACTION VS DIST of 1000 DIFFERENCES---------------------------------------
figure(1)
clf(1)
plot([0:14],serie_diff_surr_media,'b')
hold on
plot([0:14],mean(resta),'r')

%% Calculate the empirical p_values
serie_resta_prom=mean(resta);

for bip=1:size(serie_diff_surr_media,2)
    asyn_diff(bip)=serie_resta_prom(bip);
    n_emp(bip)=size(serie_diff_surr_media,1);
%     p_val_emp(bip)=size(find(serie_diff_surr_media(:,1)>=asyn_diff(bip)),1)/(n_emp(bip)+1);
pos_size=size(find(serie_diff_surr_media(:,1)>=abs(asyn_diff(bip))),1);
neg_size=size(find(serie_diff_surr_media(:,1)<=-abs(asyn_diff(bip))),1);  
p_val_emp(bip)=(pos_size + neg_size)/(n_emp(bip)+1);
    
end

[h,crit_p,adj_ci_cvrg,adj_p]=fdr_bh(p_val_emp(1,2:bip_comp+1),alfa,'pdep','yes');
sig_bip=adj_p<alfa;

%%% Uncomment to save
% save(['clase' num2str(clase) '_sim-' num2str(sim) '-' num2str(date)],'p_val_emp')

%% DISTRIBUTION OF DIFFERENCES------------------------------------------------
for bip=1:size(serie_diff_surr_media,2)
    q_sup(bip) = quantile(serie_diff_surr_media(:,bip),1-alfa/2);
    q_inf(bip) = quantile(serie_diff_surr_media(:,bip),alfa/2);
end

figure(2)
clf(2)
plot([0:14],q_sup,':r','Linewidth',3)
hold on
plot([0:14],mean(resta),'b','Linewidth',3)
plot([0:14],q_inf,':r','Linewidth',3)

for bip=1:length(adj_p)
   if adj_p(1,bip)<alfa
       plot(bip,0,'*g')
   end
end

legend(['Percentiles ' num2str(1-alfa/2) ' y ' num2str(alfa/2)],'average time serie')
xlabel('step n relative to perturbation','FontSize',fontsize)
ylabel('e_{n} [ms]','FontSize',fontsize)
if clase==1
    title(['Average - Small Pert.: Simmetry ' titulos(sim)],'FontSize',fontsize)
elseif clase==2
    title(['Average - Large Pert.: Simmetry ' titulos(sim)],'FontSize',fontsize)
end

%% DISTRIBUTION TMAX
for j=1:size(serie_diff_surr,2)
    
    for bip=1:size(serie_diff_surr,3)
        
        [h,p,ci,stats] = ttest(serie_diff_surr(:,j,bip));
        t_values(j,bip)=stats.tstat;
        
    end
    
    t_max(j)=t_values(j,(find(abs(t_values(j,:))==max(abs(t_values(j,:))))));

end

q_sup = quantile(t_max,1-alfa/2);
q_inf = quantile(t_max,alfa/2);