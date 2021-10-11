clear
close all
clc

addpath('../Data/');
addpath(genpath('Functions/'));

Figures_str = '../Figures/Figure3/';
if ~exist(Figures_str)
    mkdir(Figures_str)
end

load('SplicingLib_db.mat')

get_GC =@(x) (sum(x=='C')+sum(x=='G'))./length(x);
get_Y = @(x) (sum(x=='T')+sum(x=='C'))/length(x);
get_T = @(x) (sum(x=='T'))/length(x);

%% remove negative control variants
splicing_lib_tbl(cell2mat(cellfun(@(x) ~isempty(strfind(x,'ss control')),splicing_lib_tbl.type,'un',0)),:)=[];
splicing_lib_tbl(cell2mat(cellfun(@(x) ~isempty(strfind(x,'double')),splicing_lib_tbl.type,'un',0)),:)=[];
splicing_lib_tbl(cellfind(splicing_lib_tbl.type,'synthetic control'),:)=[];
splicing_lib_tbl(cellfind(splicing_lib_tbl.type,'synthetic alternative background control'),:)=[];
splicing_lib_tbl(cellfind(splicing_lib_tbl.type,'synthetic hairpin - control'),:)=[];

%% Figure 3A
figure
inds = find(splicing_lib_tbl.splicing_eff_median>0 | splicing_lib_tbl.cryptic_splicing_eff_log>0);
scatter(splicing_lib_tbl.splicing_eff_median+splicing_lib_tbl.cryptic_splicing_eff_log,splicing_lib_tbl.expression_log)
ax = gca;
ax.FontSize = 20;
splicing_lib_tbl.expression_log(~isfinite(splicing_lib_tbl.expression_log))=nan;
[r,p] = nancorr(splicing_lib_tbl.splicing_eff_median+splicing_lib_tbl.cryptic_splicing_eff_log,splicing_lib_tbl.expression_log);
display(sprintf('total splicing efficiency VS expression r=%.2f',r))

export_fig(sprintf('%sA - SE_Expression_Scatter.png',Figures_str),'-png','-r100','-transparent');

%% Figure 3C

inds=find(cellfun(@isempty,splicing_lib_tbl.intron_len));
splicing_lib_tbl(inds,:)=[];
inds=find(cellfun(@length,splicing_lib_tbl.branch_pos_from_3)~=1);
splicing_lib_tbl(inds,:)=[];
inds=find(cellfun(@length,splicing_lib_tbl.W_tail_len)~=1);
splicing_lib_tbl(inds,:)=[];


splicing_lib_tbl.intron_len = cell2mat(splicing_lib_tbl.intron_len);
splicing_lib_tbl.branch_pos_from_3 = cell2mat(splicing_lib_tbl.branch_pos_from_3);
splicing_lib_tbl.W_tail_len = cell2mat(splicing_lib_tbl.W_tail_len);

syn_inds=find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic')),splicing_lib_tbl.type,'un',0)));
syn_inds=[syn_inds;  find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic control')),splicing_lib_tbl.type,'un',0)))];
syn_tbl=splicing_lib_tbl(syn_inds,:);
consensus_inds = 1:size(syn_tbl,1);
consensus_inds=cellfind(syn_tbl.SS_5,'GTATGT');
consensus_inds = intersect(consensus_inds,cellfind(syn_tbl.branch_seq,'TACTAAC'));
consensus_inds = intersect(consensus_inds,[cellfind(syn_tbl.SS_3,'TAG'); cellfind(syn_tbl.SS_3,'CAG')] );

inds = setdiff(1:size(syn_tbl,1),consensus_inds);
ctrl_inds = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'control')),syn_tbl.type,'un',0)));
ctrl_inds = [ctrl_inds; cellfind(syn_tbl.orth_organism,'agos')];
inds = setdiff(inds,ctrl_inds);

figure('units','centimeters','outerposition',[2 2 18 16])
CT=cbrewer('qual','Set1',9);
CT = [CT(5,:); CT(1,:); CT(9,:)];
myviolinplot_mean({syn_tbl.splicing_eff_median(consensus_inds),-1,syn_tbl.splicing_eff_median(inds)},CT)
ylabel(sprintf('Splicing\nefficiency'),'FontSize',30)
ax = gca;
% ax.FontSize = 30;
xlim([0.5 3.5])
ylim([0 1])
ax.XTickLabel = {'Consensus sites','','Non-consensus sites'};
ax.FontSize = 20;
% ax.XTickLabelRotation = -45;
% xlabel('5'' Splice site','FontSize',30)

export_fig(sprintf('%sB - consensus_sites_distribution_combinatorial.png',Figures_str),'-png','-r100','-transparent');

%% Figure 3D


SS_5=['GTATGT';'GTACGT';'GTATGA';'GTAAGT';'GTATGC';'GCATGT'];
SS_5_freq=[198; 24; 10; 9; 9 ; 4];
SS_5_prob=SS_5_freq/sum(SS_5_freq);
SS_5_prob=SS_5_prob([4,2,3,5,1]);

ss5=unique(syn_tbl.SS_5);
% ss5(2)=[];
% ss5(1)=[];

consensus_inds = 1:size(syn_tbl,1);
% consensus_inds=cellfind(syn_tbl.SS_5,'GTATGT');
consensus_inds = intersect(consensus_inds,cellfind(syn_tbl.branch_seq,'TACTAAC'));
consensus_inds = intersect(consensus_inds,[cellfind(syn_tbl.SS_3,'TAG'); cellfind(syn_tbl.SS_3,'CAG')] );
syn_tbl = syn_tbl(consensus_inds,:);

consensus_inds = cellfind(syn_tbl.SS_5,'GTATGT');
Consensus_tbl = syn_tbl(consensus_inds,:);
NonConsensus_tbl = syn_tbl(setdiff(1:size(syn_tbl,1),consensus_inds),:);

for i = 1:size(NonConsensus_tbl)
    tmp_inds = find(Consensus_tbl.intron_len==NonConsensus_tbl.intron_len(i));
    tmp_inds = intersect(tmp_inds,find(Consensus_tbl.branch_pos_from_3==NonConsensus_tbl.branch_pos_from_3(i)));
    tmp_inds = intersect(tmp_inds,find(Consensus_tbl.W_tail_len==NonConsensus_tbl.W_tail_len(i)));
    tmp_inds = intersect(tmp_inds,cellfind(Consensus_tbl.SS_3,NonConsensus_tbl.SS_3{i}));
    
    ref_sp_eff = mean(Consensus_tbl.splicing_eff_median(tmp_inds));
    
    NonConsensus_tbl.sp_eff_diff(i) = NonConsensus_tbl.splicing_eff_median(i)-ref_sp_eff;
end


figure('units','centimeters','outerposition',[2 2 32 16])
ss5(cellfind(ss5,'GTATGT'))=[];
splicing_eff_cell={};
for i=1:length(ss5)
    splicing_eff_cell=[splicing_eff_cell, {NonConsensus_tbl.sp_eff_diff(find(strcmp(NonConsensus_tbl.SS_5,ss5{i})))}];
end
p_distribution_mat = get_T_p_val_mat(splicing_eff_cell);
p_distribution_mat_h = p_distribution_mat < 0.05/(size(p_distribution_mat,1)*(size(p_distribution_mat,1)-1)/2);
myviolinplot_mean(splicing_eff_cell);
ax = gca;
ax.FontSize = 25;
ylabel(sprintf('\\DeltaSplicing\nefficiency'),'FontSize',30)
xlim([0.5 length(ss5)+0.5])
ax = gca;
ax.XTick=1:length(splicing_eff_cell);
ss5=cellfun(@(x) strrep(x,'T','U'),ss5,'un',0);
ax.XTickLabel = ss5;
ax.FontSize = 30;
ax.XTickLabelRotation = -45;
xlabel('5'' Splice site','FontSize',30)
xlim([0.5 length(ss5)+0.5])
ylim([-1 1])

export_fig(sprintf('%sD1 - SS_5_dSE.png',Figures_str),'-png','-r100','-transparent');

figure
% SS_5_prob=[0; 0; SS_5_prob];
CT=cbrewer('qual','Set1',4);
CT_temp = cbrewer('qual','Set1',9);
CT = [CT; CT_temp(5,:)];
h=pie(SS_5_prob/sum(SS_5_prob),repmat({''},length(SS_5_prob),1));%,ss5(2:end));
for i=1:2:length(h)
    h(i).FaceColor=CT((i+1)/2,:);
end
for i = 2:2:length(h)
    h(i).FontSize = 22;
end

export_fig(sprintf('%sD2 - SS_5_abundance.png',Figures_str),'-png','-r100','-transparent');

%% Figure 3E



syn_tbl=splicing_lib_tbl(syn_inds,:);

branch=['NNCTAAC';'NNCTAAT';'NNTTAAC';'TACTAAC'];
branch_freq=[253-226,1,4,226];
branch_prob=branch_freq/sum(branch_freq(1:3));
branch=unique(syn_tbl.branch_seq);
% branch(1)=[];
% branch(1)=[];

tmp_branch_inds=find(cell2mat(cellfun(@(y) (strcmp(y,'NNCTAAC')),splicing_lib_tbl.branch_seq,'un',0)));
tmp_branch_tbl=splicing_lib_tbl(tmp_branch_inds,:);
tmp_branch_NN=cellfun(@(x,y) x(y(1)+40:y(1)+41),tmp_branch_tbl.seq,tmp_branch_tbl.branch_inds,'un',0);
TA_inds=find(strcmp(tmp_branch_NN,'TA'));

splicing_lib_tbl.branch_seq(tmp_branch_inds(TA_inds))=repmat({'TACTAAC'},length(TA_inds),1);

syn_inds=find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic')),splicing_lib_tbl.type,'un',0)));
syn_inds=[syn_inds;  find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic control')),splicing_lib_tbl.type,'un',0)))];
syn_tbl=splicing_lib_tbl(syn_inds,:);

consensus_inds=cellfind(syn_tbl.SS_5,'GTATGT');
% consensus_inds = intersect(consensus_inds,cellfind(syn_tbl.branch_seq,'TACTAAC'));
consensus_inds = intersect(consensus_inds,[cellfind(syn_tbl.SS_3,'TAG'); cellfind(syn_tbl.SS_3,'CAG')] );
syn_tbl = syn_tbl(consensus_inds,:);

consensus_inds = cellfind(syn_tbl.branch_seq,'TACTAAC');
Consensus_tbl = syn_tbl(consensus_inds,:);
NonConsensus_tbl = syn_tbl(setdiff(1:size(syn_tbl,1),consensus_inds),:);

for i = 1:size(NonConsensus_tbl)
    tmp_inds = find(Consensus_tbl.intron_len==NonConsensus_tbl.intron_len(i));
    tmp_inds = intersect(tmp_inds,find(Consensus_tbl.branch_pos_from_3==NonConsensus_tbl.branch_pos_from_3(i)));
    tmp_inds = intersect(tmp_inds,find(Consensus_tbl.W_tail_len==NonConsensus_tbl.W_tail_len(i)));
    tmp_inds = intersect(tmp_inds,cellfind(Consensus_tbl.SS_3,NonConsensus_tbl.SS_3{i}));
    
    ref_sp_eff = mean(Consensus_tbl.splicing_eff_median(tmp_inds));
    
    NonConsensus_tbl.sp_eff_diff(i) = NonConsensus_tbl.splicing_eff_median(i)-ref_sp_eff;
end

figure('units','centimeters','outerposition',[2 2 32 16])
branch(cellfind(branch,'TACTAAC'))=[];
splicing_eff_cell={};
for i=1:length(branch)
    splicing_eff_cell=[splicing_eff_cell, {NonConsensus_tbl.sp_eff_diff(find(strcmp(NonConsensus_tbl.branch_seq,branch{i})))}];
end
p_distribution_mat = get_T_p_val_mat(splicing_eff_cell);
p_distribution_mat_h = p_distribution_mat < 0.05/(size(p_distribution_mat,1)*(size(p_distribution_mat,1)-1)/2);
myviolinplot_mean(splicing_eff_cell);
ax = gca;
tmp_xlim=ax.XLim;
xlim([1 tmp_xlim(2)])
ax.XTick = 1:length(branch);
ax.FontSize = 30;
ylabel(sprintf('\\DeltaSplicing\nefficiency'),'FontSize',30)
xlim([0.5 length(branch)+0.5])
ax = gca;
ax.XTick=1:length(splicing_eff_cell);
branch=cellfun(@(x) strrep(x,'T','U'),branch,'un',0);
ax.XTickLabel =branch;
ax.XTickLabelRotation = -45;
ax.FontSize = 30;
xlabel('Branch site','fontsize',30)
xlim([0.5 length(branch)+0.5])
ylim([-1 1])

export_fig(sprintf('%sE1 - branch_dSE.png',Figures_str),'-png','-r100','-transparent');

figure
CT=cbrewer('qual','Set1',3);
CT = [CT; CT_temp(5,:)];
h=pie(branch_prob/sum(branch_prob),repmat({''},length(branch_prob),1));%,branch(2:end));
for i=1:2:length(h)
    h(i).FaceColor=CT((i+1)/2,:);
end
for i = 2:2:length(h)
    h(i).FontSize = 22;
end

export_fig(sprintf('%sE2 - branch_abundance.png',Figures_str),'-png','-r100','-transparent');

%% Figure 3F

syn_tbl=splicing_lib_tbl(syn_inds,:);

ss3=unique(syn_tbl.SS_3);
% ss3=ss3([2,4,1,3,5]);
% ss3(1)=[];
% ss3(1)=[];

SS_3=['TAG';'CAG';'AAG'];
SS_3_freq=[137;116;7];
SS_3_prob=SS_3_freq/sum(SS_3_freq);
SS_3_prob=SS_3_prob([3,2,1]);

consensus_inds=cellfind(syn_tbl.SS_5,'GTATGT');
consensus_inds = intersect(consensus_inds,cellfind(syn_tbl.branch_seq,'TACTAAC'));
% consensus_inds = intersect(consensus_inds,[cellfind(syn_tbl.SS_3,'TAG'); cellfind(syn_tbl.SS_3,'CAG')] );
syn_tbl = syn_tbl(consensus_inds,:);

consensus_inds = [cellfind(syn_tbl.SS_3,'TAG'); cellfind(syn_tbl.SS_3,'CAG')];
Consensus_tbl = syn_tbl(consensus_inds,:);
NonConsensus_tbl = syn_tbl(setdiff(1:size(syn_tbl,1),consensus_inds),:);

for i = 1:size(syn_tbl)
    switch syn_tbl.SS_3{i}
        case 'TAG'
            tmp_inds = find(Consensus_tbl.intron_len==syn_tbl.intron_len(i));
            tmp_inds = intersect(tmp_inds,find(Consensus_tbl.branch_pos_from_3==syn_tbl.branch_pos_from_3(i)));
            tmp_inds = intersect(tmp_inds,find(Consensus_tbl.W_tail_len==syn_tbl.W_tail_len(i)));
            tmp_inds = intersect(tmp_inds,cellfind(Consensus_tbl.SS_3,'CAG'));
        case 'CAG'
            tmp_inds = find(Consensus_tbl.intron_len==syn_tbl.intron_len(i));
            tmp_inds = intersect(tmp_inds,find(Consensus_tbl.branch_pos_from_3==syn_tbl.branch_pos_from_3(i)));
            tmp_inds = intersect(tmp_inds,find(Consensus_tbl.W_tail_len==syn_tbl.W_tail_len(i)));
            tmp_inds = intersect(tmp_inds,cellfind(Consensus_tbl.SS_3,'TAG'));
        otherwise        
            tmp_inds = find(Consensus_tbl.intron_len==syn_tbl.intron_len(i));
            tmp_inds = intersect(tmp_inds,find(Consensus_tbl.branch_pos_from_3==syn_tbl.branch_pos_from_3(i)));
            tmp_inds = intersect(tmp_inds,find(Consensus_tbl.W_tail_len==syn_tbl.W_tail_len(i)));
        %     tmp_inds = intersect(tmp_inds,cellfind(Consensus_tbl.SS_3,NonConsensus_tbl.SS_3{i}));
    end
    
    ref_sp_eff = mean(Consensus_tbl.splicing_eff_median(tmp_inds));
    
    syn_tbl.sp_eff_diff(i) = syn_tbl.splicing_eff_median(i)-ref_sp_eff;
end

figure('units','centimeters','outerposition',[2 2 32 16])
splicing_eff_cell={};
for i=1:length(ss3)
    splicing_eff_cell=[splicing_eff_cell, {syn_tbl.sp_eff_diff(find(strcmp(syn_tbl.SS_3,ss3{i})))}];
end
p_distribution_mat = get_T_p_val_mat(splicing_eff_cell);
p_distribution_mat_h = p_distribution_mat < 0.05/(size(p_distribution_mat,1)*(size(p_distribution_mat,1)-1)/2);
myviolinplot_mean(splicing_eff_cell);
ax = gca;
ax.FontSize = 30;
ylabel(sprintf('\\DeltaSplicing\nefficiency'),'FontSize',30)
xlim([0.5 length(ss3)+0.5])
ax = gca;
ax.XTick=1:length(splicing_eff_cell);
ss3=cellfun(@(x) strrep(x,'T','U'),ss3,'un',0);
ax.XTickLabel =ss3;
ax.XTickLabelRotation = -45;
ax.FontSize = 30;
xlabel('3'' splice site','fontsize',30)
xlim([0.5 length(ss3)+0.5])
ylim([-1 1])

export_fig(sprintf('%sF1 - SS_3_dSE.png',Figures_str),'-png','-r100','-transparent');

figure

CT=cbrewer('qual','Set1',3);
h=pie(SS_3_prob/sum(SS_3_prob),repmat({''},length(SS_3_prob),1));%,ss3(2:end));
for i=1:2:length(h)
    h(i).FaceColor=CT((i+1)/2,:);
end
for i = 2:2:length(h)
    h(i).FontSize = 22;
end

export_fig(sprintf('%sF2 - SS_3_abundance.png',Figures_str),'-png','-r100','-transparent');

%% Figure 3G
load('SplicingLib_db.mat')

tmp_splicing_lib_tbl = splicing_lib_tbl([cellfind(splicing_lib_tbl.type,'endogenous'); cellfind(splicing_lib_tbl.type,'orthologous')],:);
consensus_inds=cellfind(tmp_splicing_lib_tbl.SS_5,'GTATGT');
consensus_inds = intersect(consensus_inds,cellfind(tmp_splicing_lib_tbl.branch_seq,'TACTAAC'));
tmp_tbl = tmp_splicing_lib_tbl(consensus_inds,:);

inds_AG_3 = [cellfind(tmp_tbl.SS_3,'AAG'); cellfind(tmp_tbl.SS_3,'CAG'); cellfind(tmp_tbl.SS_3,'TAG')];

inds_AAG = cellfind(tmp_tbl.SS_3,'AAG');
inds_CAG = cellfind(tmp_tbl.SS_3,'CAG');
inds_TAG = cellfind(tmp_tbl.SS_3,'TAG');

figure('units','centimeters','outerposition',[2 2 24 13.5])
myviolinplot_mean({tmp_tbl.splicing_eff_median(inds_AAG),tmp_tbl.splicing_eff_median(inds_CAG),tmp_tbl.splicing_eff_median(inds_TAG),tmp_tbl.splicing_eff_median(setdiff(1:size(tmp_tbl,1),inds_AG_3))})
ax = gca;
ylabel('Cryptic splicing efficiency')
ax.XTickLabel = {'AAG','CAG','TAG','Other'};
ax.FontSize = 20;

export_fig(sprintf('%sG - intended_3SS_endogenous.png',Figures_str),'-png','-r100','-transparent');

