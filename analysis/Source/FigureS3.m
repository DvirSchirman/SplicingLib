clear
close all
clc

addpath('../Data/');
addpath(genpath('Functions/'));

Figures_str = '../Figures/FigureS3/';
if ~exist(Figures_str)
    mkdir(Figures_str)
end

load('SplicingLib_db.mat')

get_GC =@(x) (sum(x=='C')+sum(x=='G'))./length(x);
get_Y = @(x) (sum(x=='T')+sum(x=='C'))/length(x);
get_T = @(x) (sum(x=='T'))/length(x);
get_C = @(x) (sum(x=='C'))/length(x);

inds = find(cellfun(@isempty, splicing_lib_tbl.SS5_inds));    
splicing_lib_tbl(inds,:)=[];
inds = find(cellfun(@isempty, splicing_lib_tbl.SS3_inds));    
splicing_lib_tbl(inds,:)=[];

splicing_lib_tbl(cell2mat(cellfun(@(x) ~isempty(strfind(x,'ss control')),splicing_lib_tbl.type,'un',0)),:)=[];
splicing_lib_tbl(cell2mat(cellfun(@(x) ~isempty(strfind(x,'double')),splicing_lib_tbl.type,'un',0)),:)=[];
splicing_lib_tbl(cellfind(splicing_lib_tbl.type,'synthetic control'),:)=[];
splicing_lib_tbl(cellfind(splicing_lib_tbl.type,'synthetic alternative background control'),:)=[];
splicing_lib_tbl(cellfind(splicing_lib_tbl.type,'synthetic hairpin - control'),:)=[];


splicing_lib_tbl.GC=cellfun(@(x) get_GC(x),splicing_lib_tbl.seq);
splicing_lib_tbl.intron_seq=cellfun(@(x,y,z) x(y(1)+42:z(end)+42),splicing_lib_tbl.seq,splicing_lib_tbl.SS5_inds,splicing_lib_tbl.SS3_inds,'un',0);
splicing_lib_tbl.intron_GC=cellfun(@(x) get_GC(x(1:end-22)), splicing_lib_tbl.intron_seq);
splicing_lib_tbl.polyY=cellfun(@(x) get_Y(x(end-19:end)), splicing_lib_tbl.intron_seq);
splicing_lib_tbl.polyT=cellfun(@(x) get_T(x(end-19:end)), splicing_lib_tbl.intron_seq);
splicing_lib_tbl.polyC=cellfun(@(x) get_C(x(end-19:end)), splicing_lib_tbl.intron_seq);

inds = cellfind(splicing_lib_tbl.SS_5,'GTATGT');
inds = intersect(inds, cellfind(splicing_lib_tbl.branch_seq,'TACTAAC'));
inds = intersect(inds, [cellfind(splicing_lib_tbl.SS_3,'TAG'); cellfind(splicing_lib_tbl.SS_3,'CAG')]);

splicing_lib_tbl = splicing_lib_tbl(inds,:);


%% Figure S3A

bins=[0.2,0.3,0.4,0.45,0.5,0.55,0.6,1];
bins_sp_cell={};
bins_str={};
for i=2:length(bins)
    inds=find(splicing_lib_tbl.polyY>=bins(i-1) & splicing_lib_tbl.polyY<bins(i) & splicing_lib_tbl.polyC>0.3*mean(bins(i-1:i)));
    bins_sp_cell=[bins_sp_cell, {splicing_lib_tbl.splicing_eff_median(inds)}];
    if i==2
        bins_str{i-1}=sprintf('0-%.2f',bins(i));
    elseif i==length(bins)
        bins_str{i-1}=sprintf('%.2f-1',bins(i-1));
    else
        bins_str{i-1}=sprintf('%.2f-%.2f',bins(i-1),bins(i));
    end
end
figure('units','centimeters','outerposition',[2 2 32 18])
bins_sp_cell_nozero=cellfun(@(x) x(x>0),bins_sp_cell,'un',0);
nozero_frac=cellfun(@(x,y) length(x)/length(y),bins_sp_cell_nozero,bins_sp_cell);

CT = cbrewer('qual','Paired',10);
myviolinplot_mean(bins_sp_cell,repmat(CT(1,:),length(bins_sp_cell_nozero),1))
ax = gca;
ax.XTickLabel =bins_str;
ax.FontSize=18;
ylabel('Splicing efficiency','FontSize',20)
xlabel('Y-enrichment','FontSize',20)


mean_vec=cellfun(@nanmean,bins_sp_cell_nozero);
if mod(length(bins),2) %%odd
    bins_mean=[reshape(bins(1:end-1),2,(length(bins)-1)/2), reshape(bins(2:end),2,(length(bins)-1)/2)];
else %%even
    bins_mean=[reshape(bins,2,length(bins)/2), reshape(bins(2:end-1),2,(length(bins)-2)/2)];
end
[~,inds]=sort(bins_mean(1,:));
bins_mean=bins_mean(:,inds);
bins_mean=mean(bins_mean);
[r,p]=nancorr(bins_mean',mean_vec');


export_fig(sprintf('%sA - polyY_violin.png',Figures_str),'-png','-r100','-transparent');


%% Figure S3B

bins=[0,0.25,0.35,0.45,1];
bins=[0 linspace(0.3,0.45,7) 1];
bins_sp_cell={};
bins_str={};
for i=2:length(bins)
    inds=find(splicing_lib_tbl.intron_GC>=bins(i-1) & splicing_lib_tbl.intron_GC<bins(i));
    bins_sp_cell=[bins_sp_cell, {splicing_lib_tbl.splicing_eff_median(inds)}];
    if i==2
        bins_str{i-1}=sprintf('0-%.2f',bins(i));
    elseif i==length(bins)
        bins_str{i-1}=sprintf('%.2f-1',bins(i-1));
    else
        bins_str{i-1}=sprintf('%.2f-%.2f',bins(i-1),bins(i));
    end
end
figure('units','centimeters','outerposition',[2 2 32 18])
bins_sp_cell_nozero=cellfun(@(x) x(x>0),bins_sp_cell,'un',0);
nozero_frac=cellfun(@(x,y) length(x)/length(y),bins_sp_cell_nozero,bins_sp_cell);

mean_vec=cellfun(@nanmean,bins_sp_cell);
if mod(length(bins),2) %%odd
    bins_mean=[reshape(bins(1:end-1),2,(length(bins)-1)/2), reshape(bins(2:end),2,(length(bins)-1)/2)];
else %%even
    bins_mean=[reshape(bins,2,length(bins)/2), reshape(bins(2:end-1),2,(length(bins)-2)/2)];
end
[~,inds]=sort(bins_mean(1,:));
bins_mean=bins_mean(:,inds);
bins_mean=mean(bins_mean);
[r,p]=nancorr(bins_mean',mean_vec','Pearson');

CT = cbrewer('qual','Paired',10);
myviolinplot_mean(bins_sp_cell,repmat(CT(1,:),length(bins_sp_cell_nozero),1))
ax = gca;
ax.XTickLabel =bins_str;
ax.FontSize=16;
ylabel(sprintf('Splicing efficiency'),'FontSize',30)
xlabel('GC content','FontSize',30)

export_fig(sprintf('%sB - GC_violin.png',Figures_str),'-png','-r100','-transparent');


%% Figure S3C

bins = [0,73,88,105,121,150];
bins_sp_cell={};
bins_str={};
for i=2:length(bins)
    inds=find(cell2mat(splicing_lib_tbl.intron_len)>=bins(i-1) & cell2mat(splicing_lib_tbl.intron_len)<bins(i));
    bins_sp_cell=[bins_sp_cell, {splicing_lib_tbl.splicing_eff_median(inds)}];
    if i==2
        bins_str{i-1}=sprintf('0-%d',bins(i));
    elseif i==length(bins)
        bins_str{i-1}=sprintf('%d-148',bins(i-1));
    else
        bins_str{i-1}=sprintf('%d-%d',bins(i-1),bins(i));
    end
    
end
figure('units','centimeters','outerposition',[2 2 32 18])
bins_sp_cell_nozero=cellfun(@(x) x(x>0),bins_sp_cell,'un',0);
nozero_frac=cellfun(@(x,y) length(x)/length(y),bins_sp_cell_nozero,bins_sp_cell);

p_distribution_mat = get_T_p_val_mat(bins_sp_cell_nozero);
p_distribution_mat_h = p_distribution_mat < 0.05/(size(p_distribution_mat,1)*(size(p_distribution_mat,1)-1)/2);

CT = cbrewer('qual','Paired',10);
myviolinplot_mean(bins_sp_cell,repmat(CT(1,:),length(bins_sp_cell),1))
ax = gca;
ax.XTickLabel = bins_str;
ax = gca;
ax.XAxis.Limits=[0.5 length(bins_sp_cell)+0.5];
ax.XTick=1:length(bins_sp_cell);
ax.XTickLabel =bins_str;
ax.FontSize = 25;
ylabel(sprintf('Splicing efficiency'),'FontSize',35)
xlabel('Intron length','FontSize',35)
% ylim([0 100])

export_fig(sprintf('%sC - intron_length.png',Figures_str),'-png','-r100','-transparent');



%% Figure S3D
load('SplicingLib_db.mat')
tmp_splicing_lib_tbl = splicing_lib_tbl([cellfind(splicing_lib_tbl.type,'synthetic'); cellfind(splicing_lib_tbl.type,'synthetic mutated')],:);
consensus_inds=cellfind(tmp_splicing_lib_tbl.SS_5,'GTATGT');
consensus_inds = intersect(consensus_inds,cellfind(tmp_splicing_lib_tbl.branch_seq,'TACTAAC'));
tmp_tbl = tmp_splicing_lib_tbl(consensus_inds,:);

inds_AG_3 = [cellfind(tmp_tbl.SS_3,'AAG'); cellfind(tmp_tbl.SS_3,'CAG'); cellfind(tmp_tbl.SS_3,'TAG')];

inds_AAG = cellfind(tmp_tbl.SS_3,'AAG');
inds_CAG = cellfind(tmp_tbl.SS_3,'CAG');
inds_TAG = cellfind(tmp_tbl.SS_3,'TAG');
figure('units','centimeters','outerposition',[2 2 24 13.5])
myviolinplot_mean({tmp_tbl.cryptic_splicing_eff_log(inds_AAG),tmp_tbl.cryptic_splicing_eff_log(inds_CAG),tmp_tbl.cryptic_splicing_eff_log(inds_TAG),tmp_tbl.cryptic_splicing_eff_log(setdiff(1:size(tmp_tbl,1),inds_AG_3))})
ax = gca;
ylabel('Cryptic splicing efficiency')
ax.XTickLabel = {'AAG','CAG','TAG','Other'};
ax.FontSize = 20;

export_fig(sprintf('%sD - cryptic_intended_3SS_synthetic.png',Figures_str),'-png','-r100','-transparent');
%% Figure S3D

% inds=find(cellfun(@length,splicing_lib_tbl.branch_pos_from_3)==1);
% tbl=splicing_lib_tbl(inds,:);
% CT = cbrewer('qual','Paired',10);
% 
% bins=[0,20,30,40,50,100];
% bins_sp_cell={};
% bins_str={};
% for i=2:length(bins)
%     inds=find(cell2mat(tbl.branch_pos_from_3)>=bins(i-1) & cell2mat(tbl.branch_pos_from_3)<bins(i));
%     bins_sp_cell=[bins_sp_cell, {tbl.splicing_eff_median(inds)}];
%     if i==2
%         bins_str{i-1}=sprintf('0-%d',bins(i));
%     elseif i==length(bins)
%         bins_str{i-1}=sprintf('%d-100',bins(i-1));
%     else
%         bins_str{i-1}=sprintf('%d-%d',bins(i-1),bins(i));
%     end
% end
% figure('units','centimeters','outerposition',[2 2 32 18])
% bins_sp_cell_nozero=cellfun(@(x) x(x>0),bins_sp_cell,'un',0);
% 
% p_distribution_mat = get_T_p_val_mat(bins_sp_cell_nozero);
% p_distribution_mat_h = p_distribution_mat < 0.05/(size(p_distribution_mat,1)*(size(p_distribution_mat,1)-1)/2);
% 
% myviolinplot_mean(bins_sp_cell,repmat(CT(1,:),length(bins_sp_cell),1))
% ax = gca;
% ax.XTick=1:length(bins_sp_cell);
% ax.XTickLabel =bins_str;
% ax.XAxis.Limits=[0.5 length(bins_sp_cell)+0.5];
% ax.FontSize=25;
% ylabel(sprintf('Splicing efficiency'),'FontSize',30)
% 
% export_fig(sprintf('%sD - branch_to_3SS.png',Figures_str),'-png','-r100','-transparent');
% 

