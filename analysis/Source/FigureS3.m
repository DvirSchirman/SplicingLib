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

alt_inds = find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic alternative background')),splicing_lib_tbl.type,'un',0)));
syn_inds = find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic')),splicing_lib_tbl.type,'un',0)) & ...
    strcmp(splicing_lib_tbl.SS_5,'GTATGT') & strcmp(splicing_lib_tbl.branch_seq,'TACTAAC') & ...
    (strcmp(splicing_lib_tbl.SS_3,'CAG') | strcmp(splicing_lib_tbl.SS_3,'TAG')));

syn_backgrounds_tbl=splicing_lib_tbl([syn_inds;alt_inds],:);

alt_bg_SE_tbl=table();

n=1;
for i=1:size(syn_backgrounds_tbl,1)
    if strcmp(syn_backgrounds_tbl.SS_3{i},'CAG')
        id3=1;
    elseif strcmp(syn_backgrounds_tbl.SS_3{i},'TAG')
        id3=2;
    else
        error('unsupported 3'' site');
    end
    
    id_vec=[syn_backgrounds_tbl.intron_len{i},syn_backgrounds_tbl.branch_pos_from_3{i},syn_backgrounds_tbl.W_tail_len{i} ,id3];
    syn_backgrounds_tbl.id_vec{i}=num2str(id_vec);
    
    if i>1
        ind = cellfind(alt_bg_SE_tbl.id_vec,syn_backgrounds_tbl.id_vec{i});
        if isempty(ind)
            alt_bg_SE_tbl.id_vec{n}=syn_backgrounds_tbl.id_vec{i};
            alt_bg_SE_tbl.intron_len{n}=syn_backgrounds_tbl.intron_len{i};
            alt_bg_SE_tbl.branch_pos_from_3{n}=syn_backgrounds_tbl.branch_pos_from_3{i};
            alt_bg_SE_tbl.W_tail_len{n}=syn_backgrounds_tbl.W_tail_len{i};
            alt_bg_SE_tbl.id3(n)=id3;
            n=n+1;
        end
    else
        alt_bg_SE_tbl.id_vec{n}=syn_backgrounds_tbl.id_vec{i};
        alt_bg_SE_tbl.intron_len{n}=syn_backgrounds_tbl.intron_len{i};
        alt_bg_SE_tbl.branch_pos_from_3{n}=syn_backgrounds_tbl.branch_pos_from_3{i};
        alt_bg_SE_tbl.W_tail_len{n}=syn_backgrounds_tbl.W_tail_len{i};
        alt_bg_SE_tbl.id3(n)=id3;
        n=n+1;
    end   
end

alt_bg_dSE_tbl=alt_bg_SE_tbl;
alt_bg_SE_tbl{:,6:15}=zeros(size(alt_bg_SE_tbl,1),10);
for i=1:size(syn_backgrounds_tbl,1)
    ind = cellfind(alt_bg_SE_tbl.id_vec,syn_backgrounds_tbl.id_vec{i});
    bg = syn_backgrounds_tbl.background{i};
    alt_bg_SE_tbl{ind,5+bg}=syn_backgrounds_tbl.splicing_eff_median(i);
end
alt_bg_dSE_tbl{:,6:14}=alt_bg_SE_tbl{:,7:15}-alt_bg_SE_tbl{:,6};
bg_names={'MUD1','UBC9','SNC1','rand1','rand2','rand3','rand4','rand5','rand6','rand7'};
alt_bg_SE_tbl.Properties.VariableNames(6:end)=bg_names;
alt_bg_dSE_tbl.Properties.VariableNames(6:end)=bg_names(2:end);

%%
figure('units','centimeters','outerposition',[5 5 32 18])
violin_cell={};
for i=2:length(bg_names)
    violin_cell{1,i-1}=alt_bg_dSE_tbl.(bg_names{i});
end
myviolinplot_mean(violin_cell)
ax = gca;
ax.XTickLabel = bg_names(2:end);
% xlabel('Bacground sequence')
% ylabel('\Delta Splicing efficiency')
ax.FontSize = 20;


[~,p] = cellfun(@ttest, violin_cell);

p_mat = get_Tpaired_p_val_mat(violin_cell);

export_fig(sprintf('%sA - Alternative_backgrounds.png',Figures_str),'-png','-r100','-transparent');

%% Figure S3B

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


export_fig(sprintf('%sB - polyY_violin.png',Figures_str),'-png','-r100','-transparent');


%% Figure S3C

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

export_fig(sprintf('%sC - GC_violin.png',Figures_str),'-png','-r100','-transparent');


%% Figure S3D

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

export_fig(sprintf('%sD - intron_length.png',Figures_str),'-png','-r100','-transparent');


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

