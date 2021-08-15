clear
close all
clc

addpath('../Data/');
addpath(genpath('Functions/'));

Figures_str = '../Figures/FigureS1/';
if ~exist(Figures_str)
    mkdir(Figures_str)
end

load('SplicingLib_db.mat')
load('double_sp_eff_tbl.mat')

seed=20200310;
rng(seed);

get_GC =@(x) (sum(x=='C')+sum(x=='G'))./length(x);
get_Y = @(x) (sum(x=='T')+sum(x=='C'))/length(x);
get_T = @(x) (sum(x=='T'))/length(x);
get_C = @(x) (sum(x=='C'))/length(x);

inds = find(cellfun(@isempty, splicing_lib_tbl.SS5_inds));    
splicing_lib_tbl(inds,:)=[];
inds = find(cellfun(@isempty, splicing_lib_tbl.SS3_inds));    
splicing_lib_tbl(inds,:)=[];

ctrl_inds = [cellfind(splicing_lib_tbl.type,'synthetic control'); cellfind(splicing_lib_tbl.type,'synthetic alternative background control'); ...
    cellfind(splicing_lib_tbl.type,'synthetic hairpin - control'); cellfind(splicing_lib_tbl.type,'endogenous ss control'); ...
    cellfind(splicing_lib_tbl.type,'orthologous ss control')];
ctrl_tbl = splicing_lib_tbl(ctrl_inds,:);

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

% load('ddG_tbl.mat')
% [~,inds1,inds2]=intersect(splicing_lib_tbl.BC,ddG_tbl.BC);
% splicing_lib_tbl.ddG=nan*ones(size(splicing_lib_tbl,1),1);
% splicing_lib_tbl.ddG(inds1)=ddG_tbl.ddG(inds2);

%% Figure S1A
figure
histogram(splicing_lib_tbl.splicing_eff_median(splicing_lib_tbl.splicing_eff_median>0),0:0.05:1)
xlabel('Splicing efficiency');
ylabel('Number of variants');
ax = gca;
ax.FontSize = 18;

export_fig(sprintf('%sA - SE_histogram.png',Figures_str),'-png','-r100','-transparent');

%% Figure S1B
inds = find(splicing_lib_tbl.BC_num==4);
multiple_BC_tbl=splicing_lib_tbl(inds,:);
designed_seq=cellfun(@(x) x(43:end-30),multiple_BC_tbl.seq,'un',0);
designed_seq_u=unique(designed_seq);
multiple_BC_SE_mat=nan*ones(length(designed_seq_u),8);
multiple_BC_cryptic_SE_mat=nan*ones(length(designed_seq_u),8);
n=1;
for i=1:length(designed_seq_u)
    if ~mod(i,100)
        disp(i)
    end
    inds=cellfind(designed_seq,designed_seq_u{i});
    if length(inds)<2
        continue
    end
    SE=multiple_BC_tbl.splicing_eff_median(inds);
    multiple_BC_SE_mat(n,1:length(inds))=SE';
    cryptic_SE=multiple_BC_tbl.cryptic_splicing_eff_log(inds);
    multiple_BC_cryptic_SE_mat(n,1:length(inds))=cryptic_SE';
    n=n+1;
end
multiple_BC_SE_mat(n:end,:)=[];
multiple_BC_cryptic_SE_mat(n:end,:)=[];

multiple_BC_noZero_mat=multiple_BC_SE_mat(nansum(multiple_BC_SE_mat,2)>0,:);
multiple_BC_var=nanvar(multiple_BC_SE_mat,[],2);
multiple_BC_noZero_var=nanvar(multiple_BC_noZero_mat,[],2);

multiple_BC_cryptic_noZero_mat=multiple_BC_cryptic_SE_mat(nansum(multiple_BC_cryptic_SE_mat,2)>0,:);
multiple_BC_cryptic_var=nanvar(multiple_BC_cryptic_SE_mat,[],2);
multiple_BC_cryptic_noZero_var=nanvar(multiple_BC_cryptic_noZero_mat,[],2);

multiple_BC_total_SE_mat = multiple_BC_SE_mat + multiple_BC_cryptic_SE_mat;
multiple_BC_total_noZero_mat=multiple_BC_total_SE_mat(nansum(multiple_BC_total_SE_mat,2)>0,:);
multiple_BC_total_var=nanvar(multiple_BC_total_SE_mat,[],2);
multiple_BC_total_noZero_var=nanvar(multiple_BC_total_noZero_mat,[],2);

multiple_BC_noZero_vec = multiple_BC_noZero_mat(isfinite(multiple_BC_noZero_mat));
multiple_BC_ctrl_noZero = reshape(randsample(multiple_BC_noZero_vec,floor(length(multiple_BC_noZero_vec)/4)*4),floor(length(multiple_BC_noZero_vec)/4),4);

multiple_BC_ctrl_noZero_var=nanvar(multiple_BC_ctrl_noZero,[],2);

multiple_BC_total_noZero_vec = multiple_BC_total_noZero_mat(isfinite(multiple_BC_total_noZero_mat));
multiple_BC_ctrl_total_noZero = reshape(randsample(multiple_BC_total_noZero_vec,floor(length(multiple_BC_total_noZero_vec)/4)*4),floor(length(multiple_BC_total_noZero_vec)/4),4);

multiple_BC_ctrl_total_noZero_var=nanvar(multiple_BC_ctrl_total_noZero,[],2);

load('multiple_BC_ctrl_stats.mat')
%%
figure
histogram(randsample(multiple_BC_total_var_mean,1e4),'Normalization','probability')
hold on
stem(mean(multiple_BC_total_noZero_var),0.08,'^','linewidth',1.5)
xlim([0 0.2])
h=legend('10^4 randomized sets','multiple barcode sets');
h.Location='southwest';
h.FontSize=12;
ax = gca;
ax.FontSize=18;
xlabel('Mean variance','fontsize',20)
ylabel('Frequency','fontsize',20)

export_fig(sprintf('%sB - multiple_BC.png',Figures_str),'-png','-r100','-transparent');

%% Figure S1C

mycorrplot(multiple_BC_total_noZero_mat(:,1:4))
r = nancorr(multiple_BC_total_noZero_mat(:,1:4));
r = nancorr(multiple_BC_total_noZero_mat(:,1:4),'Spearman');
r(r==1)=nan;
nanmean(r(:))

export_fig(sprintf('%sC - multiple_BC_total_corrplot.png',Figures_str),'-png','-r100','-transparent');


%% Figure S1D

unspliced_RNA = randn(5e3,1);
spliced_RNA = randn(5e3,1);

SE = 10.^spliced_RNA./(10.^spliced_RNA+10.^unspliced_RNA);

figure
scatter(SE,log10(10.^spliced_RNA+10.^unspliced_RNA))
xlabel('Splicing efficiency','fontsize',16)
ylabel('Total RNA abundance','fontsize',16)
ax = gca;
ax.FontSize = 20;

[r,p]=nancorr(SE,log10(10.^spliced_RNA+10.^unspliced_RNA))


export_fig(sprintf('%sD - SE_expression_scatter_randomized.png',Figures_str),'-png','-r100','-transparent');

