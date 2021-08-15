clear
close all
clc

addpath('../Data/');
addpath(genpath('Functions/'));

Figures_str = '../Figures/FigureS2/';
if ~exist(Figures_str)
    mkdir(Figures_str)
end

load('SplicingLib_db.mat')

%% Figure S2A
inds=find(cellfun(@isempty,splicing_lib_tbl.intron_len));
splicing_lib_tbl(inds,:)=[];
inds=find(cellfun(@length,splicing_lib_tbl.branch_pos_from_3)~=1);
splicing_lib_tbl(inds,:)=[];
inds=find(cellfun(@length,splicing_lib_tbl.W_tail_len)~=1);
splicing_lib_tbl(inds,:)=[];
inds = cellfind(splicing_lib_tbl.type,'R loops - intron');
splicing_lib_tbl(inds,:)=[];

splicing_lib_tbl.intron_len = cell2mat(splicing_lib_tbl.intron_len);
splicing_lib_tbl.branch_pos_from_3 = cell2mat(splicing_lib_tbl.branch_pos_from_3);
splicing_lib_tbl.W_tail_len = cell2mat(splicing_lib_tbl.W_tail_len);

%%%% short introns %%%%

inds = find(splicing_lib_tbl.intron_len<105);
splicing_lib_tbl = splicing_lib_tbl(inds,:);

SS_5={'GTATGT';'GTACGT';'GTATGA';'GTAAGT';'GTATGC'};
SS_5_freq=[198; 24; 10; 9; 9 ; 4];
SS_5_prob=SS_5_freq/sum(SS_5_freq);
SS_5_prob=SS_5_prob([4,2,3,5,1]);

syn_inds=find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic')),splicing_lib_tbl.type,'un',0)));
syn_inds=[syn_inds;  find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic control')),splicing_lib_tbl.type,'un',0)))];
syn_tbl=splicing_lib_tbl(syn_inds,:);

ss5=SS_5;

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
    a=1;
end

ss5(cellfind(ss5,'GTATGT'))=[];
sp_eff_cell={};
for i=1:length(ss5)
    sp_eff_cell=[sp_eff_cell, {NonConsensus_tbl.sp_eff_diff(find(strcmp(NonConsensus_tbl.SS_5,ss5{i})))}];
end
splicing_eff_cell_ss5_short=sp_eff_cell;

%%%% long introns %%%%

load('SplicingLib_db.mat')
inds=find(cellfun(@isempty,splicing_lib_tbl.intron_len));
splicing_lib_tbl(inds,:)=[];
inds=find(cellfun(@length,splicing_lib_tbl.branch_pos_from_3)~=1);
splicing_lib_tbl(inds,:)=[];
inds=find(cellfun(@length,splicing_lib_tbl.W_tail_len)~=1);
splicing_lib_tbl(inds,:)=[];
inds = cellfind(splicing_lib_tbl.type,'R loops - intron');
splicing_lib_tbl(inds,:)=[];

splicing_lib_tbl.intron_len = cell2mat(splicing_lib_tbl.intron_len);
splicing_lib_tbl.branch_pos_from_3 = cell2mat(splicing_lib_tbl.branch_pos_from_3);
splicing_lib_tbl.W_tail_len = cell2mat(splicing_lib_tbl.W_tail_len);

inds = find(splicing_lib_tbl.intron_len>105);
splicing_lib_tbl = splicing_lib_tbl(inds,:);

SS_5={'GTATGT';'GTACGT';'GTATGA';'GTAAGT';'GTATGC'};
SS_5_freq=[198; 24; 10; 9; 9 ; 4];
SS_5_prob=SS_5_freq/sum(SS_5_freq);
SS_5_prob=SS_5_prob([4,2,3,5,1]);

syn_inds=find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic')),splicing_lib_tbl.type,'un',0)));
syn_inds=[syn_inds;  find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic control')),splicing_lib_tbl.type,'un',0)))];
syn_tbl=splicing_lib_tbl(syn_inds,:);

ss5=SS_5;

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
    a=1;
end

ss5(cellfind(ss5,'GTATGT'))=[];
sp_eff_cell={};
for i=1:length(ss5)
    sp_eff_cell=[sp_eff_cell, {NonConsensus_tbl.sp_eff_diff(find(strcmp(NonConsensus_tbl.SS_5,ss5{i})))}];
end
splicing_eff_cell_ss5_long=sp_eff_cell;

figure('units','centimeters','outerposition',[2 2 32 16])
null_dist = randn(10,1)-100;
splicing_eff_cell = {};
for i=1:length(splicing_eff_cell_ss5_short)
    splicing_eff_cell = [splicing_eff_cell, splicing_eff_cell_ss5_short(i), splicing_eff_cell_ss5_long(i), null_dist];
end
CT = cbrewer('qual','Paired',10);
CT = CT([5,6,1:4,9,10],:);
CT2 = zeros(length(splicing_eff_cell),3);
CT2(1:3:end,:)=CT(1:2:end,:);
CT2(2:3:end,:)=CT(2:2:end,:);
myviolinplot_mean(splicing_eff_cell,CT2);
ax = gca;
ax.FontSize = 25;
ylabel(sprintf('\\Delta Splicing\nefficiency'),'FontSize',30)
xlim([0.5 length(ss5)+0.5])
ax = gca;
ax.XTick=1:length(splicing_eff_cell);
ss5=cellfun(@(x) strrep(x,'T','U'),ss5,'un',0);
ax.XTickLabel = repmat({''},length(splicing_eff_cell),1);
ax.XTickLabel(1:3:end) = ss5;
ax.FontSize = 30;
ax.XTickLabelRotation = -45;
xlabel('5'' Splice site','FontSize',30)
xlim([0.5 length(splicing_eff_cell)+0.5])
ylim([-1 1])
for i=1:length(splicing_eff_cell_ss5_short)
    [~, p_ss5(i)] = ttest2(splicing_eff_cell_ss5_short{i},splicing_eff_cell_ss5_long{i});
end

export_fig(sprintf('%sA - ss5_diff_length.png',Figures_str),'-png','-r100','-transparent'); 

%% Figure S2B
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

export_fig(sprintf('%sB - cryptic_intended_3SS_synthetic.png',Figures_str),'-png','-r100','-transparent');

%% Figure S2C

inds = cellfind(splicing_lib_tbl.type,'synthetic');
inds = [inds; cellfind(splicing_lib_tbl.type,'synthetic mutated')];
splicing_lib_tbl=splicing_lib_tbl(inds,:);
inds = find(cellfun(@isempty,splicing_lib_tbl.branch_seq));
splicing_lib_tbl(inds,:)=[];
consensus_inds=cellfind(splicing_lib_tbl.SS_5,'GTATGT');
consensus_inds = intersect(consensus_inds,cellfind(splicing_lib_tbl.branch_seq,'TACTAAC'));
consensus_inds = intersect(consensus_inds,[cellfind(splicing_lib_tbl.SS_3,'TAG'); cellfind(splicing_lib_tbl.SS_3,'CAG')] );

consensus_spliced_mean = nanmean(splicing_lib_tbl.splicing_eff_median(consensus_inds));

Alphabet='ACGT';

consensus_ss5='GTATGT';
non_consensus_mat_ss5=zeros(4,length(consensus_ss5));
non_consensus_mat_ss5_p=nan*ones(4,length(consensus_ss5));
for i=1:length(consensus_ss5)
    consensus_inds = cellfind(splicing_lib_tbl.branch_seq,'TACTAAC');
    consensus_inds = intersect(consensus_inds,[cellfind(splicing_lib_tbl.SS_3,'TAG'); cellfind(splicing_lib_tbl.SS_3,'CAG')] );
    tmp_inds=setdiff(1:length(consensus_ss5),i);
    consensus_inds = intersect(consensus_inds,find(cellfun(@(x) strcmp(x(tmp_inds),consensus_ss5(tmp_inds)),splicing_lib_tbl.SS_5)));
    loci_vec=cellfun(@(x) x(i),splicing_lib_tbl.SS_5);
    for j=1:4
        if Alphabet(j)==consensus_ss5(i)
            continue;
        end
        inds = find(loci_vec==Alphabet(j));
        inds = intersect(inds,consensus_inds);
        if isempty(inds)
            non_consensus_mat_ss5(j,i)=nan;
        else
            non_consensus_mat_ss5(j,i)=nanmean(splicing_lib_tbl.splicing_eff_median(inds));
            [~, non_consensus_mat_ss5_p(j,i)] = ttest(splicing_lib_tbl.splicing_eff_median(inds));
        end
    end
end
       
non_consensus_mat_ss5_h=non_consensus_mat_ss5_p<0.025/sum(sum(isfinite(non_consensus_mat_ss5_p)));
for i=1:length(consensus_ss5)
    non_consensus_mat_ss5(find(Alphabet==consensus_ss5(i)),i)=consensus_spliced_mean;
end

figure
imagesc(non_consensus_mat_ss5);
ax = gca;
ax.CLim=[0 0.5];
colormap(jet(256))
colorbar;

export_fig(sprintf('%sC1 - ss5_mutations_imagesc.png',Figures_str),'-png','-r100','-transparent');

consensus_branch='TACTAAC';
non_consensus_mat_branch=zeros(4,length(consensus_branch));
non_consensus_mat_branch_p=nan*ones(4,length(consensus_branch));
for i=1:length(consensus_branch)
    consensus_inds=cellfind(splicing_lib_tbl.SS_5,'GTATGT');
    consensus_inds = intersect(consensus_inds,[cellfind(splicing_lib_tbl.SS_3,'TAG'); cellfind(splicing_lib_tbl.SS_3,'CAG')] );
    tmp_inds=setdiff(1:length(consensus_branch),i);
    consensus_inds = intersect(consensus_inds,find(cellfun(@(x) strcmp(x(tmp_inds),consensus_branch(tmp_inds)),splicing_lib_tbl.branch_seq)));
    loci_vec=cellfun(@(x) x(i),splicing_lib_tbl.branch_seq);
    for j=1:4
        if Alphabet(j)==consensus_branch(i)
            continue;
        end
        inds = find(loci_vec==Alphabet(j));
        inds = intersect(inds,consensus_inds);
        if isempty(inds)
            non_consensus_mat_branch(j,i)=nan;
        else
            non_consensus_mat_branch(j,i)=nanmean(splicing_lib_tbl.splicing_eff_median(inds));
            [~, non_consensus_mat_branch_p(j,i)] = ttest(splicing_lib_tbl.splicing_eff_median(inds));
        end
    end
end

non_consensus_mat_branch_h=non_consensus_mat_branch_p<0.025/sum(sum(isfinite(non_consensus_mat_branch_p)));
for i=1:length(consensus_branch)
    non_consensus_mat_branch(find(Alphabet==consensus_branch(i)),i)=consensus_spliced_mean;
end
figure
h=imagesc(non_consensus_mat_branch);
ax = gca;
ax.CLim=[0 0.5];
colormap(jet(256))
colorbar;

export_fig(sprintf('%sC2 - branch_mutations_imagesc.png',Figures_str),'-png','-r100','-transparent');


consensus_ss3_1='TAG';
consensus_ss3_2='CAG';
non_consensus_mat_ss3=zeros(4,length(consensus_ss3_1));
non_consensus_mat_ss3_p=nan*ones(4,length(consensus_ss3_1));
for i=1:length(consensus_ss3_1)
    consensus_inds = cellfind(splicing_lib_tbl.SS_5,'GTATGT');
    consensus_inds = intersect(consensus_inds,cellfind(splicing_lib_tbl.branch_seq,'TACTAAC') );
    loci_vec=cellfun(@(x) x(i),splicing_lib_tbl.SS_3);
    for j=1:4
        if i>1
            if Alphabet(j)==consensus_ss3_1(i)
                continue;
            end
            tmp_inds=setdiff(1:length(consensus_ss3_1),i);
            consensus_inds = intersect(consensus_inds,find(cellfun(@(x) strcmp(x(tmp_inds),consensus_ss3_1(tmp_inds))+strcmp(x(tmp_inds),consensus_ss3_2(tmp_inds)),splicing_lib_tbl.SS_3)));
        else
            tmp_inds=setdiff(1:length(consensus_ss3_1),i);
            consensus_inds = intersect(consensus_inds,find(cellfun(@(x) strcmp(x(tmp_inds),consensus_ss3_1(tmp_inds)),splicing_lib_tbl.SS_3)));
        end
        inds = find(loci_vec==Alphabet(j));
        inds = intersect(inds,consensus_inds);
        if isempty(inds)
            non_consensus_mat_ss3(j,i)=nan;
        else
            non_consensus_mat_ss3(j,i)=nanmean(splicing_lib_tbl.splicing_eff_median(inds));
            [~, non_consensus_mat_ss3_p(j,i)] = ttest(splicing_lib_tbl.splicing_eff_median(inds));
        end
    end
end

non_consensus_mat_ss3_h=non_consensus_mat_ss3_p<0.025/sum(sum(isfinite(non_consensus_mat_ss3_p)));
non_consensus_mat_ss3(non_consensus_mat_ss3==0)=consensus_spliced_mean;
non_consensus_mat_ss3(isnan(non_consensus_mat_ss3))=0;
figure
imagesc(non_consensus_mat_ss3);
ax = gca;
ax.CLim=[0 0.5];
colormap(jet(256))
colorbar;
    
export_fig(sprintf('%sC3 - ss3_mutations_imagesc.png',Figures_str),'-png','-r100','-transparent');

%% Figure S2D
load('SplicingLib_db.mat')

load('site_dG_GC_30.mat')

[~, inds1, inds2] = intersect(splicing_lib_tbl.BC, sliding_dG_GC_tbl.BC);
splicing_lib_tbl.ss5_dG=[];
splicing_lib_tbl.ss3_dG=[];
splicing_lib_tbl.branch_dG=[];
splicing_lib_tbl.ss5_dG(inds1) = sliding_dG_GC_tbl.ss5_dG(inds2);
splicing_lib_tbl.ss3_dG(inds1) = sliding_dG_GC_tbl.ss3_dG(inds2);
splicing_lib_tbl.branch_dG(inds1) = sliding_dG_GC_tbl.branch_dG(inds2);

get_GC =@(x) (sum(x=='C')+sum(x=='G'))./length(x);
get_Y = @(x) (sum(x=='T')+sum(x=='C'))/length(x);
get_T = @(x) (sum(x=='T'))/length(x);

neg_ctrl_inds = cell2mat(cellfun(@(x) ~isempty(strfind(x,'ss control')),splicing_lib_tbl.type,'un',0));
neg_ctrl_inds = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'ss control')),splicing_lib_tbl.type,'un',0)));
neg_ctrl_inds = [neg_ctrl_inds; find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'double')),splicing_lib_tbl.type,'un',0)))];
neg_ctrl_inds = [neg_ctrl_inds; cellfind(splicing_lib_tbl.type,'synthetic control')];
neg_ctrl_inds = [neg_ctrl_inds; cellfind(splicing_lib_tbl.type,'synthetic alternative background control')];
neg_ctrl_inds = [neg_ctrl_inds; cellfind(splicing_lib_tbl.type,'synthetic hairpin - control')];

neg_ctrl_tbl = splicing_lib_tbl(neg_ctrl_inds,:);
splicing_lib_tbl(neg_ctrl_inds,:)=[];

syn_mut_inds = find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic hairpin' )),splicing_lib_tbl.type,'un',0)));
% syn_mut_inds = [syn_mut_inds; find(cell2mat(cellfun(@(y) (strcmp(y,'endogenous - mutated sites' )),splicing_lib_tbl.type,'un',0)))];
% syn_mut_inds = [syn_mut_inds; find(cell2mat(cellfun(@(y) (strcmp(y,'endogenous pombe - mutated sites' )),splicing_lib_tbl.type,'un',0)))];

syn_mut_tbl = splicing_lib_tbl(syn_mut_inds,:);
% stem_len_max=cellfun(@(x,y,z) nanmax([x,y,z]),syn_mut_tbl.ss5_stem_len,syn_mut_tbl.branch_stem_len,syn_mut_tbl.ss3_stem_len);
% inds=find(stem_len_max>14);
% syn_mut_tbl=syn_mut_tbl(inds,:);

syn_inds=find(cell2mat(cellfun(@(y) (strcmp(y,'synthetic')),splicing_lib_tbl.type,'un',0)));
consensus_inds=cellfind(splicing_lib_tbl.SS_5,'GTATGT');
consensus_inds = intersect(consensus_inds,cellfind(splicing_lib_tbl.branch_seq,'TACTAAC'));
consensus_inds = intersect(consensus_inds,[cellfind(splicing_lib_tbl.SS_3,'TAG'); cellfind(splicing_lib_tbl.SS_3,'CAG')] );
syn_consesus_inds = intersect(syn_inds,consensus_inds);

consensus_tbl = splicing_lib_tbl(syn_consesus_inds,:);

hairpin_SE_tbl = table();
n=1;
for i=1:size(consensus_tbl,1)
    if strcmp(consensus_tbl.SS_3{i},'CAG')
        id3=1;
    elseif strcmp(consensus_tbl.SS_3{i},'TAG')
        id3=2;
    else
        error('unsupported 3'' site');
    end
    
    id_vec=[consensus_tbl.intron_len{i},consensus_tbl.branch_pos_from_3{i},consensus_tbl.W_tail_len{i} ,id3];
    consensus_tbl.id_vec{i}=num2str(id_vec);
end

consensus_SE_CAG = consensus_tbl.splicing_eff_median(cellfind(consensus_tbl.id_vec, '121   30    8    1'));
consensus_SE_TAG = consensus_tbl.splicing_eff_median(cellfind(consensus_tbl.id_vec, '121   30    8    2'));

for i=1:size(syn_mut_tbl,1)
    if strcmp(syn_mut_tbl.SS_3{i},'CAG')
        syn_mut_tbl.dSE(i) = (syn_mut_tbl.splicing_eff_median(i) - consensus_SE_CAG)/consensus_SE_CAG;
    elseif strcmp(syn_mut_tbl.SS_3{i},'TAG')
        syn_mut_tbl.dSE(i) = (syn_mut_tbl.splicing_eff_median(i) - consensus_SE_TAG)/consensus_SE_TAG;
    else
        error('unsupported 3'' site');
    end
end

states = [  1,0,0;...
            0,1,0;...
            0,0,1;...
            1,1,0;...
            1,0,1;...
            0,1,1;...
            1,1,1];

violin_cell={};
states_str=[];
for i=1:7
    violin_cell{i}=syn_mut_tbl.dSE(sum(cell2mat(syn_mut_tbl.site_hairpin)==states(i,:),2)==3);
    states_str=[states_str;mat2str(states(i,:))];
end

figure('units','centimeters','outerposition',[5 5 32 18])
myviolinplot_mean(violin_cell)
ax = gca;
ylabel('\Delta Splicing efficiency')
% ylim([0 1])
% ax.YTick=-1:0.2:0;
% ax.YTickLabel={'0';'0.2';'0.4';'0.6';'0.8';'1'};
ax.FontSize=25;
% ax.XTickLabel=[{'[0 0 0]'}, states_str];
% xlabel('Mutated sites ([5'', branch, 3''])','FontSize',12)
% ylabel(sprintf('Non-Zero\nSplicing efficiency'),'FontSize',12)
xlim([0.5 length(violin_cell)+0.5]);

p = get_T_p_val_mat(violin_cell);

export_fig(sprintf('%sD - hairpin_sites.png',Figures_str),'-png','-r100','-transparent');





