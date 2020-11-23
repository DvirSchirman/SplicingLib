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


%% Figure S1A
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

export_fig(sprintf('%sA - multiple_BC.png',Figures_str),'-png','-r100','-transparent');

%% Figure S1B

mycorrplot(multiple_BC_total_noZero_mat(:,1:4))
r = nancorr(multiple_BC_total_noZero_mat(:,1:4));
r(r==1)=nan;
nanmean(r(:))

export_fig(sprintf('%sB - multiple_BC_total_corrplot.png',Figures_str),'-png','-r100','-transparent');


%% Figure S1C

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


export_fig(sprintf('%sB - SE_expression_scatter_randomized.png',Figures_str),'-png','-r100','-transparent');

