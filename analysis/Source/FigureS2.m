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






